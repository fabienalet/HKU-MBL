#ifndef OBS_H
#define OBS_H

#ifdef USE_MKL
#include <mkl.h>
#include <mkl_cblas.h>
#include <mkl_lapack.h>
// Some work-around in case MKL is not available
#else
typedef struct __MKL_Complex16 {
  double real;
  double imag;
} MKL_Complex16;
typedef int MKL_INT;
typedef MKL_Complex16 lapack_complex_double;
#include <cblas.h>
#include <lapacke.h>
#endif

class observable {
  typedef std::vector<unsigned short int> Conf;

 private:
  basis *basis_pointer;

 public:
  observable(basis *this_basis_pointer, int this_number_threads);
  observable() {}
  ~observable() {}

  // useful to know for openMP computation
  int number_threads;
  void init();

  // Entanglement entropy
  double entang_entropy(double q);
  void compute_entanglement_spectrum(PetscScalar *state);
  std::vector<double> entanglement_spectrum;

  // Local magnetization and imbalance
  std::vector<double> sz_local;
  void compute_local_magnetization(PetscScalar *state);
  double product_state_imbalance(int, int, int);
  // copied from the basis
  std::vector<double> sz_for_basis_state;
};

observable::observable(basis *this_basis_pointer, int this_number_threads) {
  number_threads = this_number_threads;
  basis_pointer = this_basis_pointer;

  // copy from basis
  sz_for_basis_state.resize(basis_pointer->NUMBER_OF_STATES, 0);
  sz_for_basis_state = basis_pointer->sz_for_basis_state;
}

// For the state provided in parameters, this function computes
// <state|S_i^z|state> for all sites i, and store it in sz_local
// state is entirely owned by this processor
void observable::compute_local_magnetization(PetscScalar *state) {
  int L = basis_pointer->L;
  int LA = basis_pointer->LA;
  int LB = basis_pointer->LB;
  // local magnetization
  sz_local.resize(L, 0.);
  for (int i = 0; i < L; i++) sz_local[i] = 0.;

  std::vector<unsigned short int> config(L, 0);

  int i = 0;  // configuration index
  for (int nsa = 0; nsa < basis_pointer->valid_sectors; ++nsa) {
    for (int ca = 0; ca < basis_pointer->Confs_in_A[nsa].size(); ++ca) {
      for (int r = 0; r < LA; ++r) {
        config[r] = basis_pointer->Confs_in_A[nsa][ca][r];
      }
      for (int cb = 0; cb < basis_pointer->Confs_in_B[nsa].size(); ++cb) {
        for (int r = 0; r < LB; ++r) {
          config[r + LA] = basis_pointer->Confs_in_B[nsa][cb][r];
        }

        for (int d = 0; d < L; ++d) {
#ifdef PETSC_USE_COMPLEX
          sz_local[d] += sz_for_basis_state[config[d]] *
                         PetscRealPart(state[i] * PetscConjComplex(state[i]));
#else
          sz_local[d] += sz_for_basis_state[config[d]] * state[i] * state[i];
#endif
        }

        i++;
      }
    }
  }
}

// Return the observable <S_i^z(t)S_i^z(0)> assuming initial state is the
// product state index0, corresponding to the configuration (nsa0, ca0, cb0) in
// the computational basis. CAUTION : compute_local_magnetization has to be
// called before.

double observable::product_state_imbalance(int nsa0, int ca0, int cb0) {
  double imbalance = 0;
  int L = basis_pointer->L;
  int LA = basis_pointer->LA;
  imbalance = 0;
  std::vector<int> config(L, 0);
  for (int d = 0; d < L; ++d) {
    if (d < LA) {
      config[d] = basis_pointer->Confs_in_A[nsa0][ca0][d];
    } else {
      config[d] = basis_pointer->Confs_in_B[nsa0][cb0][d - LA];
    }
    imbalance += sz_for_basis_state[config[d]] * sz_local[d];
  }
  imbalance /= (L / 4.);
  return imbalance;
}

// Computes Renyi entanglement entropy S_q = 1/(1-q) log \sum_i lambda_i^q
// with lambda_i eigenvalues of the reduced density matrix
// For q = 1, us von Neumann entropy S_1 = - sum_i \lambda_i log(\lambda_i)
// CAUTION: compute_entanglement_spectrum has to be called before !

double observable::entang_entropy(double q) {
  double entropy = 0.;
  if (q == 1) {
    for (int i = 0; i < entanglement_spectrum.size(); i++) {
      double ai = entanglement_spectrum[i];
      if (ai != 0) {
        entropy += -ai * log(ai);
      }
    }
  } else {
    for (int i = 0; i < entanglement_spectrum.size(); i++) {
      double ai = entanglement_spectrum[i];
      entropy += pow(ai, q);
    }
    entropy = log(entropy) / (1. - q);
  }
  return entropy;
}

/*
 Computes the entanglement spectrum of the state provided as a parameter
 state is entirely owned by this processor

 The reduced density matrix is block-diagonal
A singular value decomposition (SVD) is performed with Lapack, sector by sector
Indeed the basis has been constructed such that basis states in the same
symmetry sector are consecutive state = [ sector 0 | sector 1 | ... | last
sector ]

 Uses openMP, and if available, MKL to speed up the SVD

 CAUTION : state will unfortunately be completely messed up by the call to
lapack, can't use it again!!
*/
void observable::compute_entanglement_spectrum(PetscScalar *state) {
  entanglement_spectrum.resize(0);

// We use the maximal number of threads available as we u
#ifdef USE_MKL
  mkl_set_num_threads(number_threads);
#else
  omp_set_num_threads(number_threads);
#endif

  for (int nsa = 0; nsa < basis_pointer->valid_sectors;
       ++nsa) {  // loop over all sectors
    int sizeA = basis_pointer->Confs_in_A[nsa].size();
    int sizeB = basis_pointer->Confs_in_B[nsa].size();
    int sectorsize = sizeA * sizeB;
    int start = basis_pointer->starting_conf[nsa];

    if ((sectorsize > 0)) {  // only for non-void sectors
      if (sectorsize ==
          1) {  // don't run the full machinery if the sector size is 1
        double ame = PetscAbsScalar(state[start]);
        entanglement_spectrum.push_back(ame * ame);
      } else {  // sectorsize>1, do the SVD

        /*
        Here we do the SVD using lapack
        The call to lapack are different whether the state is complex or
        real, and whether we use MKL or not (the later being checked by
        USE_MKL
        This part is pretty ugly but should be efficient
        */
        int minsize = min(sizeA, sizeB);
        // Will contain the singular values for this sector
        std::vector<double> local_svd_spectrum(minsize, 0.);

        {
          MKL_INT m = sizeB;
          MKL_INT n = sizeA;
          MKL_INT lda = m;
          MKL_INT ldu = m;
          MKL_INT ldvt = n;
          MKL_INT info, lwork;
          MKL_INT iwork[8 * minsize];
          lwork = -1;
#ifdef PETSC_USE_COMPLEX
#ifdef USE_MKL
          MKL_Complex16 wkopt;
          MKL_Complex16 *work;
          MKL_Complex16 u[ldu * m], vt[ldvt * n];
          double rwork[5 * m * m + 7 * m];
          zgesdd("N", &m, &n, ((MKL_Complex16 *)&(state[start])), &lda,
                 &local_svd_spectrum[0], u, &ldu, vt, &ldvt, &wkopt, &lwork,
                 rwork, iwork, &info);
#else
          double __complex__ wkopt;
          double __complex__ *work;
          double __complex__ u[ldu * m], vt[ldvt * n];
          double rwork[5 * m * m + 7 * m];
          zgesdd_("N", &m, &n, ((double __complex__ *)&(state[start])), &lda,
                  &local_svd_spectrum[0], u, &ldu, vt, &ldvt, &wkopt, &lwork,
                  rwork, iwork, &info);
#endif
          lwork = (MKL_INT)wkopt.real;
          work = (MKL_Complex16 *)malloc(lwork * sizeof(MKL_Complex16));
#ifdef USE_MKL
          zgesdd("N", &m, &n, ((MKL_Complex16 *)&(state[start])), &lda,
                 &local_svd_spectrum[0], u, &ldu, vt, &ldvt, work, &lwork,
                 rwork, iwork, &info);
#else
          zgesdd_("N", &m, &n, ((double __complex__ *)&(state[start])), &lda,
                  &local_svd_spectrum[0], u, &ldu, vt, &ldvt, work, &lwork,
                  rwork, iwork, &info);
#endif
#else
          double wkopt;
          double *work;
          double u[ldu * m], vt[ldvt * n];
#ifdef USE_MKL
          dgesdd("N", &m, &n, &state[start], &lda, &local_svd_spectrum[0], u,
                 &ldu, vt, &ldvt, &wkopt, &lwork, iwork, &info);
#else
          dgesdd_("N", &m, &n, &state[start], &lda, &local_svd_spectrum[0], u,
                  &ldu, vt, &ldvt, &wkopt, &lwork, iwork, &info);
#endif
          lwork = (MKL_INT)wkopt;
          work = (double *)malloc(lwork * sizeof(double));
#ifdef USE_MKL
          dgesdd("N", &m, &n, &state[start], &lda, &local_svd_spectrum[0], u,
                 &ldu, vt, &ldvt, work, &lwork, iwork, &info);
#else
          dgesdd_("N", &m, &n, &state[start], &lda, &local_svd_spectrum[0], u,
                  &ldu, vt, &ldvt, work, &lwork, iwork, &info);
#endif
#endif
          free((void *)work);
        }

        // now stores the entanglement spectrum by squating the singular values
        double s;
        for (int rr = 0; rr < local_svd_spectrum.size(); ++rr) {
          s = local_svd_spectrum[rr];
          entanglement_spectrum.push_back(s * s);
        }

      }  // sectorsize>1
    }    // sectorsize>0
  }      // loop over sectors

/*
  // It could be cleaner to sort the entanglement spectrum ...
  std::sort(entanglement_spectrum.begin(), entanglement_spectrum.end());

  // and to check that it sums to 1
  double sum = 0.;
  for (int pp = 0; pp < entanglement_spectrum.size(); ++pp) {
    sum += entanglement_spectrum[pp];
  }
  */

// Restores the number of threads to 1, just in case this interferes with future
// calls to PETSc or the linear solver...
#ifdef USE_MKL
  mkl_set_num_threads(1);
#else
  omp_set_num_threads(1);
#endif
}

#endif
