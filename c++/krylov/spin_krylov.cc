static char help[] =
    "Demo code of Krylov time evolution for XXZ chain with random fields \n"
    "For Les Houches summer school 2019 \n"
    "(C) 2019, Fabien Alet and Nicolas Macé\n";

#define PETSC_DESIRE_COMPLEX
#define PETSC_USE_COMPLEX 1
#include <omp.h>
#include <complex>

#include <slepcmfn.h>

#include <slepceps.h>
#include <algorithm>
#include <iostream>
#include <map>
#include <random>

using namespace std;

#include "lib/krylov_parameters.h"

#include "lib/Spin_basis.h"

#include "lib/SpinOneHalfXXZ_disorder.h"

#include "lib/Spin_observable.h"

int main(int argc, char **argv) {
  cout.precision(20);
  /************* Init Petsc and Slepc *********/
  SlepcInitialize(&argc, &argv, "krylov.options", help);

  /************* Init parallel work *********************************/
  // For parallelization on node (openMP)
#ifdef USE_MKL
  int ENV_NUM_THREADS = mkl_get_max_threads();
  mkl_set_num_threads(ENV_NUM_THREADS);
#else
  int ENV_NUM_THREADS = omp_get_max_threads();
#endif
  omp_set_num_threads(
      ENV_NUM_THREADS); 

  // For parallelization between nodes (MPI)
  int myrank, mpisize;
  MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
  MPI_Comm_size(MPI_COMM_WORLD, &mpisize);

  /******************* Init basis and observables ***************************/
  // spin 1/2 = 2 states per site
  int number_of_states = 2;

  // Sz=0 sector by default
  double Sz = 0;
  PetscOptionsGetReal(NULL, NULL, "-Sz", &Sz, NULL);
  // chain size L=8 by default
  int L = 8;
  PetscOptionsGetInt(NULL, NULL, "-L", &L, NULL);

  // init basis
  basis mybasis(L, Sz, myrank, number_of_states);
  int nconf = mybasis.total_number_of_confs;
  if (myrank == 0) {
    std::cout << "# L = " << L << " number of states=" << nconf << std::endl;
  }
  // init observable
  observable myobservable(&mybasis, ENV_NUM_THREADS);

  // Petsc data structure for the matrix
  PetscInt Istart, Iend;
  Mat H;
  /************************** Hamiltonian **********************************/
  {
    // distribute evenly the matrix lines between processes
    std::vector<int> local_block_sizes(mpisize, nconf / mpisize);
    for (size_t i = 0; i < nconf % mpisize; i++) local_block_sizes[i]++;
    Istart = 0;
    for (size_t i = 0; i < myrank; i++) Istart += local_block_sizes[i];
    Iend = Istart + local_block_sizes[myrank];

    // init Hamiltonian
    Hamiltonian myHamiltonian(&mybasis, &H);
    myHamiltonian.init_matrix_sizes(Istart, Iend);
    MatSetUp(H);
    myHamiltonian.create_matrix(Istart, Iend);
  }

  //  Matrix assembly
  MatSetOption(H, MAT_SYMMETRIC, PETSC_TRUE);
  MatSetOption(H, MAT_SYMMETRY_ETERNAL, PETSC_TRUE);
  MatAssemblyBegin(H, MAT_FINAL_ASSEMBLY);
  MatAssemblyEnd(H, MAT_FINAL_ASSEMBLY);
  if (myrank == 0) {
    std::cout << "#Hamiltonian matrix assembly done." << std::endl;
  }

  // MatView(H, PETSC_VIEWER_STDOUT_WORLD);

  /******************* Now get extremal energies ***************************/
  // SLEPc data structure for the extremal eigenvalue problem
  EPS eps;
  PetscScalar Eminc, Emaxc;
  bool do_ev = 0;
  if (do_ev) {  // just for checking
    EPSCreate(PETSC_COMM_WORLD, &eps);
    EPSSetOperators(eps, H, NULL);
    EPSSetProblemType(eps, EPS_HEP);
    EPSSetDimensions(eps, 1, PETSC_DECIDE, PETSC_DECIDE);

    // get the ground-state energy (in this sector)
    EPSSetWhichEigenpairs(eps, EPS_SMALLEST_REAL);
    EPSSolve(eps);
    EPSGetEigenvalue(eps, 0, &Eminc, NULL);
    double Emin = PetscRealPart(Eminc);
    if (0 == myrank) std::cout << "Emin of H = " << Emin << "\n";

    // get the maximal energy (in this sector)
    EPSSetWhichEigenpairs(eps, EPS_LARGEST_REAL);
    EPSSolve(eps);
    EPSGetEigenvalue(eps, 0, &Emaxc, NULL);
    double Emax = PetscRealPart(Emaxc);
    if (0 == myrank) std::cout << "Emax of H = " << Emax << "\n";
  }
  /************************* Krylov computation *************************/
  Parameters myparameters(myrank);
  /**** Initialize Krylov from Slepc *****/
  MatScale(H, -1.0 * PETSC_i);
  MFN mfn;
  FN fct;
  MFNCreate(PETSC_COMM_WORLD, &mfn);
  MFNSetOperator(mfn, H);

  MFNGetFN(mfn, &fct);
  FNSetType(fct, FNEXP);
  MFNSetTolerances(mfn, 1e-12, PETSC_DEFAULT);
  MFNSetFromOptions(mfn);

  int ineel = mybasis.ineel;

  /*****  Get and initialize vectors *****/
  Vec Psi_t, res;
  MatCreateVecs(H, PETSC_NULL, &Psi_t);
  MatCreateVecs(H, PETSC_NULL, &res);

  /********* Initial state ***************/
  // list of initial states (each of these consists of one basis vector)
  std::vector<int> init_states;
  if (myparameters.product_state_start) {
    // If no number of sample is specified, we start from just 1
    int seed3 = 74310;
    PetscOptionsGetInt(NULL, NULL, "-seed_inistates", &seed3, NULL);

    int Nsamp = myparameters.num_product_states;
    std::random_device device;
    std::mt19937 generator(device());
    generator.seed(seed3);
    std::uniform_int_distribution<> distr(0, nconf - 1);
    for (int i = 0; i < Nsamp; i++) {
      init_states.push_back(distr(generator));
    }
    if (myrank == 0) {
      std::cout << "# " << Nsamp
                << " random product states initialized with seed_inistates = "
                << seed3 << std::endl;
    }
  } else if (myparameters.cdw_start) {
    if (myrank == 0) {
      std::cout << "# Starting from Neel state = " << ineel << endl;
    }
    init_states.push_back(mybasis.ineel);
  } else {
    if (myrank == 0) {
      std::cout << "# Nothing specified. Starting from Neel state = " << ineel
                << endl;
    }
    init_states.push_back(mybasis.ineel);
  }

  if (myrank == 0) {
    std::cout << "# Using initial states ";
    for (auto i0 : init_states) std::cout << i0 << " ";
    std::cout << std::endl;
  }

  /*********** loop over initial states **********/
  for (auto i0 : init_states) {
    // parameters of the initial configuration
    int nsa0, ca0, cb0;
    mybasis.get_conf_coefficients(i0, nsa0, ca0, cb0);

    // Initialise Psi(0) and Res(0)
    VecSet(Psi_t, 0.0);
    VecSetValue(Psi_t, i0, 1.0, INSERT_VALUES);
    VecAssemblyBegin(Psi_t);
    VecAssemblyEnd(Psi_t);
    VecAssemblyBegin(res);
    VecAssemblyEnd(res);
    // just to make sure: compute the norm
    PetscReal norm;
    VecNormalize(Psi_t, &norm);
    VecCopy(Psi_t, res);

    if (myrank == 0) {
      std::cout << "IMBALANCE ENTANGLEMENT_ENTROPY_S1 #STATE " << i0 << endl;
    }
    /****** Time loop ******/
    int t_index;
    double dt_measure =
        (myparameters.TEEmax - myparameters.TEEmin) / myparameters.nmeasures;
    double time_next_measure = myparameters.TEEmin;
    int each_measurement = myparameters.num_times / myparameters.nmeasures;
    for (t_index = 0; t_index <= myparameters.num_times; ++t_index) {
      double t = myparameters.time_points[t_index];
      if (myparameters.delta_t_points[t_index] != 0) {
        FNSetScale(fct, myparameters.delta_t_points[t_index], 1.0);
        // result of time evolution written in res
        MFNSolve(mfn, Psi_t, res);
      }

      //  if (myrank == 0)
      //    std::cout << "... Solved time t=" << t << std::flush << std::endl;

      /************** Measurements ************/
      if ((t_index % each_measurement) == 0) {
        // Will repatriate the distributed vector res into a local vector
        Vec res_local;
        // Different strategies depending on MPI-distributed or not
        // only 1 mpi proc
        if (mpisize == 1) {
          VecCreateSeq(PETSC_COMM_SELF, nconf, &res_local);
          VecCopy(res, res_local);
        }
        // more than 1 proc
        else {
          VecScatter ctx;
          VecScatterCreateToZero(res, &ctx, &res_local);
          VecScatterBegin(ctx, res, res_local, INSERT_VALUES, SCATTER_FORWARD);
          VecScatterEnd(ctx, res, res_local, INSERT_VALUES, SCATTER_FORWARD);
          VecScatterDestroy(&ctx);
        }

        // only the processor 0 will do the measurement job ...
        if (myrank == 0) {
          PetscScalar *state;
          // provide the whole res_local vector to processor 0
          VecGetArray(res_local, &state);
          if (myparameters.measure_local || myparameters.measure_imbalance) {
            myobservable.compute_local_magnetization(state);
            if (myparameters.measure_imbalance) {
              double Imb = myobservable.product_state_imbalance(nsa0, ca0, cb0);
              std::cout << "IMBALANCE " << i0 << " " << t << " " << Imb
                        << std::endl;
            }
            if (myparameters.measure_local) {
              std::vector<double> Siz(L, 0.);
              Siz = myobservable.sz_local;
              for (int r = 0; r < L; ++r) {
                std::cout << "SZ " << i0 << " " << t << " " << r << " "
                          << Siz[r] << std::endl;
              }
            }
          }
          if (myparameters.measure_entanglement) {
            myobservable.compute_entanglement_spectrum(state);
            double S1 = myobservable.entang_entropy(1);

            std::cout << "ENTANGLEMENT_ENTROPY_S1 " << i0 << " " << t << " "
                      << S1 << endl;
          }
          // reaffect the data to each processor
          VecRestoreArray(res_local, &state);
        }  // end of 0processor
        VecDestroy(&res_local);
      }  // end measurements

      /************** End of  Measurements ***********/
      // put back Res in Psi_t for next time evolution
      VecCopy(res, Psi_t);
    }  // end of t_index
  }    // end of io loop

  SlepcFinalize();
  return 0;
}
