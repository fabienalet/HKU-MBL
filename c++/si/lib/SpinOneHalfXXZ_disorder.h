#ifndef HAMILTONIAN_H
#define HAMILTONIAN_H

class Hamiltonian {
 private:
  basis *basis_pointer;
  Mat *pointer_to_H;

 public:
  //  Hamiltonian() { get_parameters(); };
  Hamiltonian(basis *this_basis_pointer, Mat *matrix_pointer);
  ~Hamiltonian() {}

  void get_parameters();

  void init_matrix_sizes(PetscInt Istart, PetscInt Iend);
  void create_matrix(PetscInt Istart, PetscInt Iend);

  // This will be duplicated from the basis
  int L;
  int LA;
  int LB;
  int myrank;

  // Random-field XXZ spin chain parameters
  double Delta;
  double field_default;
  double disorder;
  double J;
  std::vector<double> field;
  std::vector<double> coupling;
  PetscBool pbc;
  PetscInt seed;
};

Hamiltonian::Hamiltonian(basis *this_basis_pointer, Mat *matrix_pointer) {
  basis_pointer = this_basis_pointer;
  pointer_to_H = matrix_pointer;

  // copy information from the basis for simplicity
  L = basis_pointer->L;
  LA = basis_pointer->LA;
  LB = basis_pointer->LB;
  myrank = basis_pointer->myrank;

  get_parameters();
}

void Hamiltonian::get_parameters() {
  J = 1.0;
  disorder = 0.;
  Delta = 1.;
  field_default = 0.;
  pbc = PETSC_FALSE;
  seed = 74310;

  PetscOptionsGetReal(NULL, NULL, "-J", &J, NULL);
  PetscOptionsGetReal(NULL, NULL, "-Delta", &Delta, NULL);
  PetscOptionsGetReal(NULL, NULL, "-field_default", &field_default, NULL);
  PetscOptionsGetBool(NULL, NULL, "-pbc", &pbc, NULL);
  PetscOptionsGetReal(NULL, NULL, "-disorder", &disorder, NULL);
  PetscOptionsGetInt(NULL, NULL, "-seed", &seed, NULL);  // CHKERRQ(ierr);

  coupling.resize(L, J);
  if (!(pbc)) {
    coupling[L - 1] = 0.;
  }

  field.resize(L, field_default);

  if (disorder > 0) {
    std::random_device device;
    std::mt19937 generator(device());
    generator.seed(seed);
    std::uniform_real_distribution<double> box(-disorder, disorder);
    for (int i = 0; i < L; i++) {
      field[i] += box(generator);
    }
  }

  if (myrank == 0) {
    std::cout << "# field= { ";
    for (int i = 0; i < L; i++) {
      std::cout << field[i] << " ";
    }
    std::cout << " }" << endl;
  }
}

// Istart and Iend are the two indices delimiting the indices managed by the
// current processor. For each configuration index managed by the current
// processor (Istart <= i <= Iend), and for each non-zero matrix element H_ij,
// this function computes whether j is also managed by the current processor (
// Istart <= j <= Iend) or not For optimal memory management, PETSc prefers to
// know the number of elements which are stored on the same processor (d_nnz) or
// on another one (o_nzz)

void Hamiltonian::init_matrix_sizes(PetscInt Istart, PetscInt Iend) {
  // the number of configurations (or matrix lines) handled by this processor
  size_t local_size = Iend - Istart;
  // for each configuration, the number of matrix-elements which are on this
  // processor too
  std::vector<PetscInt> d_nnz(local_size, 0);
  // on another one
  std::vector<PetscInt> o_nnz(local_size, 0);

  int i = 0;
  for (int nsa = 0; nsa < basis_pointer->valid_sectors;
       ++nsa) {  // loop over all sectors
    for (int ca = 0; ca < basis_pointer->Confs_in_A[nsa].size();
         ++ca) {  // loop over all A configurations in this sector
      // total configuration
      std::vector<unsigned short int> config(L, 0);
      // copy the configuration A in the total configuration
      std::vector<unsigned short int> confA =
          basis_pointer->Confs_in_A[nsa][ca];
      for (int r = 0; r < LA; ++r) {
        config[r] = confA[r];
      }
      for (int cb = 0; cb < basis_pointer->Confs_in_B[nsa].size();
           ++cb) {  // loop over all B configurations in this sector
        if ((i >= Istart) && (i < Iend)) {
          // construct B-part only if in the good range
          // copy the configuration B in the total configuration
          std::vector<unsigned short int> confB =
              basis_pointer->Confs_in_B[nsa][cb];
          for (int r = 0; r < LB; ++r) {
            config[r + LA] = confB[r];
          }

          // will store the configuration after the eventual spin-flip
          std::vector<unsigned short int> newconfig = config;
          std::vector<unsigned short int> newconfA = confA;
          std::vector<unsigned short int> newconfB = confB;
          int j;
          // value of diagonal element H_ii
          double diag = 0.;
          for (int r = 0; r < L; ++r) {
            // First consider the exchange part
            if (coupling[r]) {
              int rb = r;
              int rb2 = (r + 1) % L;
              //
              if (config[rb] ==
                  config[rb2]) {  // nearest-neighbor spins identical
                diag += 0.25 * Delta * coupling[r];
              } else {  // nearest-neighbors spins of opposite orientation
                diag -= 0.25 * Delta * coupling[r];
                // flip them in the newconfig
                newconfig[rb] = 1 - config[rb];
                newconfig[rb2] = 1 - config[rb2];

                // compute the sector to which this newconfig belongs, by
                // finding its (nleft,nright)
                int nleft = 0;
                int nright = 0;
                for (int p = 0; p < LA; ++p) {
                  newconfA[p] = newconfig[p];
                  nleft += newconfig[p];
                }
                for (int p = 0; p < LB; ++p) {
                  newconfB[p] = newconfig[p + LA];
                  nright += newconfig[p + LA];
                }
                // the new sector
                int new_nsa =
                    basis_pointer
                        ->particle_sector[std::make_pair(nleft, nright)];

                // find the index of this new configuratino
                j = basis_pointer->starting_conf[new_nsa] +
                    basis_pointer->InverseMapA[new_nsa][newconfA] *
                        basis_pointer->Confs_in_B[new_nsa].size() +
                    basis_pointer->InverseMapB[new_nsa][newconfB];

                if ((j >= Iend) or
                    (j < Istart)) {  // j not on the same processor
                  o_nnz[i - Istart]++;
                } else {  // j  on the same processor
                  d_nnz[i - Istart]++;
                }

                // undo the spin flip
                newconfig[rb] = config[rb];
                newconfig[rb2] = config[rb2];
              }
            }

            // Now the ield part
            if (config[r]) {
              diag += field[r] * 0.5;
            } else {
              diag -= field[r] * 0.5;
            }
          }

          if (diag != 0) {  // Make sure H_ii is non-zero
            d_nnz[i - Istart]++;
          }
        }
        i++;
      }  // loop over cb
    }    // over ca
  }      // over sectors

  // send the d_nnz, o_nnz information to PETSc
  MatCreateAIJ(PETSC_COMM_WORLD, Iend - Istart, PETSC_DECIDE,
               basis_pointer->total_number_of_confs,
               basis_pointer->total_number_of_confs, 0, d_nnz.data(), 0,
               o_nnz.data(), &(*pointer_to_H));
}

// Now this function will compute and allocate matrix elements in PETSc
// Its logic is identical to the previous one

void Hamiltonian::create_matrix(PetscInt Istart, PetscInt Iend) {
  int i = 0;
  for (int nsa = 0; nsa < basis_pointer->valid_sectors; ++nsa) {
    for (int ca = 0; ca < basis_pointer->Confs_in_A[nsa].size(); ++ca) {
      std::vector<unsigned short int> config(L, 0);
      std::vector<unsigned short int> confA =
          basis_pointer->Confs_in_A[nsa][ca];
      for (int r = 0; r < LA; ++r) {
        config[r] = confA[r];
      }
      for (int cb = 0; cb < basis_pointer->Confs_in_B[nsa].size(); ++cb) {
        if ((i >= Istart) && (i < Iend)) {
          std::vector<unsigned short int> confB =
              basis_pointer->Confs_in_B[nsa][cb];
          for (int r = 0; r < LB; ++r) {
            config[r + LA] = confB[r];
          }

          std::vector<unsigned short int> newconfig = config;
          std::vector<unsigned short int> newconfA = confA;
          std::vector<unsigned short int> newconfB = confB;
          int j;
          double diag = 0.;
          for (int r = 0; r < L; ++r) {
            if (coupling[r]) {
              int rb = r;
              int rb2 = (r + 1) % L;
              if (config[rb] == config[rb2]) {
                diag += 0.25 * Delta * coupling[r];
              } else {
                diag -= 0.25 * Delta * coupling[r];
                newconfig[rb] = 1 - config[rb];
                newconfig[rb2] = 1 - config[rb2];

                int nleft = 0;
                int nright = 0;
                for (int p = 0; p < LA; ++p) {
                  newconfA[p] = newconfig[p];
                  nleft += newconfig[p];
                }
                for (int p = 0; p < LB; ++p) {
                  newconfB[p] = newconfig[p + LA];
                  nright += newconfig[p + LA];
                }

                int new_nsa =
                    basis_pointer
                        ->particle_sector[std::make_pair(nleft, nright)];
                j = basis_pointer->starting_conf[new_nsa] +
                    basis_pointer->InverseMapA[new_nsa][newconfA] *
                        basis_pointer->Confs_in_B[new_nsa].size() +
                    basis_pointer->InverseMapB[new_nsa][newconfB];

                // H_ij + = 1/2 J
                MatSetValue(*pointer_to_H, i, j, (PetscScalar)0.5 * coupling[r],
                            ADD_VALUES);

                newconfig[rb] = config[rb];
                newconfig[rb2] = config[rb2];
              }
            }
            // field part
            if (config[r]) {
              diag += field[r] * 0.5;
            } else {
              diag -= field[r] * 0.5;
            }
          }
          if (diag != 0) {  // H_ii + = diag
            MatSetValue(*pointer_to_H, i, i, (PetscScalar)diag, ADD_VALUES);
          }
        }
        i++;
      }  // loop over cb
    }    // over ca
  }      // over nsA sectors
}

#endif
