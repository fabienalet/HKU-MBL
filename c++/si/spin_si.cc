static char help[] =
    "Demo code of Shift-Invert ED for XXZ chain with random fields \n"
    "For Les Houches summer school 2019 \n"
    "(C) 2019, Fabien Alet and Nicolas Macé\n";

#include <omp.h>
#include <slepceps.h>
// Don't use boost for this demo code
// #include <boost/dynamic_bitset.hpp>
#include <algorithm>
#include <iostream>
#include <map>
#include <random>

using namespace std;

#include "lib/Spin_basis.h"

#include "lib/SpinOneHalfXXZ_disorder.h"

#include "lib/Spin_observable.h"

int main(int argc, char **argv) {
  cout.precision(20);
  /************* Init Petsc and Slepc *********/
  SlepcInitialize(&argc, &argv, "si.options", help);

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

  /************* Shift-invert computation *********/

  /********** Get target energies *************/
  // Get the targeted energies : targets is a vector in units of 0 <= epsilon <=
  // 1 (where 0=Emin and 1=Emax)
  std::vector<double> targets;
  PetscBool targets_set;
  char *targets_c_string = new char[1000];
  // example pf command line option : -epsilon "0.5 0.52"
  PetscOptionsGetString(NULL, NULL, "-epsilon", targets_c_string, 1000,
                        &targets_set);
  if (targets_set) {
    std::string targets_string(targets_c_string);
    std::stringstream tsstr;
    tsstr.str(targets_string);
    double tgt;
    while (tsstr >> tgt) {
      targets.push_back(tgt);
    }
    delete[] targets_c_string;
  } else {
    // if nothing specified, use infinite temperature
    PetscScalar E_infinite_temperature;
    MatGetTrace(H, &E_infinite_temperature);
    E_infinite_temperature /= nconf;
    if (0 == myrank) {
      std::cout << "#Targeting E(beta = 0) = " << E_infinite_temperature
                << " (epsilon = "
                << (E_infinite_temperature - Eminc) / (Emaxc - Eminc) << ")\n";
    }
    targets.push_back((E_infinite_temperature - Eminc) / (Emaxc - Eminc));
  }

  for (double &renorm_target : targets) {  // loop over all targets wanted
    // target is in units of H
    double target = renorm_target * (Emaxc - Eminc) + Eminc;
    if (0 == myrank) {
      std::cout << "#Processing target " << target
                << " [ epsilon=" << renorm_target << " ] ... ";
    }
    /****** Do the shift-invert computation for this target ********/
    // SLEPc structure for the shift-invert computation
    EPS eps2;
    EPSCreate(PETSC_COMM_WORLD, &eps2);
    EPSSetOperators(eps2, H, NULL);
    EPSSetProblemType(eps2, EPS_HEP);
    ST st;
    EPSGetST(eps2, &st);

    EPSSetWhichEigenpairs(eps2, EPS_TARGET_REAL);
    // This will allow to choose the number of targeted states, as well as the
    // linear solver, linear solver options etc in the file si.options
    EPSSetFromOptions(eps2);
    EPSSetTarget(eps2, target);
    EPSSolve(eps2);
    if (0 == myrank) std::cout << "#Solved done. \n";

    /********* Now do measurements *******************/

    // Initialize a vector in the same processor structure as H to store an
    // eigenvector
    Vec xr;
    MatCreateVecs(H, PETSC_NULL, &xr);

    // Check if entanglement entropy is wanted
    PetscBool measure_entanglement = PETSC_TRUE;
    PetscOptionsGetBool(NULL, NULL, "-measure_entanglement",
                        &measure_entanglement, NULL);

    // Get the number of converged states and loop over them
    PetscInt nconv = 0;
    EPSGetConverged(eps2, &nconv);

    std::vector<double> energies;
    std::vector<double> rgap;
    for (int i = 0; i < nconv; i++) {
      PetscScalar Er;
      // retrieve the i-th eigenvalue and eigenstate
      EPSGetEigenpair(eps2, i, &Er, PETSC_NULL, xr, PETSC_NULL);
      energies.push_back(Er);
      if (myrank == 0) {
        std::cout << "E_" << i << " = " << Er << std::endl;
      }
      /********** Do entanglement entropy ***********/
      if (measure_entanglement) {
        // if we want entanglement, get back the eigenstate locally on processor
        // 0 (Vec_local) and compute its entanglement_spectrum
        Vec Vec_local;

        if (mpisize ==
            1) {  // if serial, a simple Sequential Vector and copy will do
          VecCreateSeq(PETSC_COMM_SELF, nconf, &Vec_local);
          VecCopy(xr, Vec_local);
        } else {  // otherwise we use a Petsc scatter objet, and scatter the
                  // eigenstate (stored on several processors) into the local
                  // vector
          VecScatter ctx;
          VecScatterCreateToZero(xr, &ctx, &Vec_local);
          VecScatterBegin(ctx, xr, Vec_local, INSERT_VALUES, SCATTER_FORWARD);
          VecScatterEnd(ctx, xr, Vec_local, INSERT_VALUES, SCATTER_FORWARD);
          VecScatterDestroy(&ctx);
        }

        if (myrank == 0) {  // Now do entanglement computation on node 0
          // for some reason we need to reallocate the SeqVector to a simple
          // array (called state)
          PetscScalar *state;
          VecGetArray(Vec_local, &state);
          myobservable.compute_entanglement_spectrum(state);
          double S1 = myobservable.entang_entropy(1);
          std::cout << "Entanglement entropy =  " << S1 << "\n";
          VecRestoreArray(Vec_local, &state);
        }
        VecDestroy(&Vec_local);
      }
    }  // loop over converged states

    /********************** Compute gap ratios *******************/
    std::sort(energies.begin(), energies.end());
    for (int r = 1; r < (nconv - 1); ++r) {
      double e1 = energies[r];
      double g0 = e1 - energies[r - 1];
      double g1 = energies[r + 1] - e1;
      if (g0 > g1) {
        rgap.push_back(g1 / g0);
      } else {
        rgap.push_back(g0 / g1);
      }
    }

    if ((myrank == 0) && (rgap.size() > 0)) {
      std::cout << "*** Gap ratios\n";
      double sum = 0.;
      for (int i = 0; i < rgap.size(); ++i) {
        std::cout << rgap[i] << "\n";
        sum += rgap[i];
      }
      std::cout << "< r > = " << sum / rgap.size() << "\n";
    }
  }  // loop over targets

  // End gracefully
  SlepcFinalize();
  return 0;
}
