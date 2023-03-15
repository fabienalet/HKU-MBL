#ifndef BASIS_H
#define BASIS_H
#include <bitset>
class basis {
  typedef std::vector<unsigned short int> Conf;

 private:
 public:
  basis(int L_, double Sz, int myrank_, int num_states_);
  basis(int L_, int LA_, double Sz, int myrank_, int num_states_);
  basis(){};
  ~basis() {}

  short int NUMBER_OF_STATES;
  std::vector<int> numberparticle_for_basis_state;
  std::vector<double> sz_for_basis_state;

  void int_to_occupation(unsigned long int num, int powof2, int L,
                         Conf& occupation);
  bool is_config_valid_inA(int num, int powof2, int L, int Nmax);

  std::vector<unsigned int> number_conf_in_A;
  std::vector<std::vector<Conf> > Confs_in_A;
  std::vector<unsigned int> starting_conf;
  std::vector<std::map<Conf, unsigned long int> > InverseMapA;

  std::vector<unsigned int> number_conf_in_B;
  std::vector<std::vector<Conf> > Confs_in_B;
  std::vector<std::map<Conf, unsigned long int> > InverseMapB;

  double TotalSz;
  int Nparticles;
  std::map<std::pair<int, int>, int> particle_sector;
  int valid_sectors;
  int LA, LB;
  int L;
  int myrank;
  int total_number_of_confs;

  void clean_basis();
  void init();
  unsigned long int index(Conf conf);
  void init_vectors_product_state(Vec& Psi_t, Vec& Psi_t2, int i0);
  void get_conf_coefficients(unsigned long int index, int& nsa, int& ca,
                             int& cb);

  // saving the indices of the Neel configurations
  int ineel;
  int ineel2;
};

basis::basis(int L_, double Sz_, int myrank_, int num_states_) {
  // size of the Hilbert space of 1 spin
  NUMBER_OF_STATES = num_states_;
  L = L_;
  myrank = myrank_;
  TotalSz = Sz_;

  LA = L / 2;
  LB = L - LA;
  init();
}

basis::basis(int L_, int LA_, double Sz_, int myrank_, int num_states_) {
  // size of the Hilbert space of 1 spin
  NUMBER_OF_STATES = num_states_;
  L = L_;
  myrank = myrank_;
  TotalSz = Sz_;

  LA = LA_;
  LB = L - LA;
  init();
}

void basis::clean_basis() {
  number_conf_in_A.clear();
  Confs_in_A.clear();
  starting_conf.clear();
  InverseMapA.clear();
  number_conf_in_B.clear();
  Confs_in_B.clear();
  InverseMapB.clear();

  particle_sector.clear();
}

void basis::init() {
  // NUMBER_OF_STATES = 2 for spin-1/2 model

  // A few simple definitions to navigate from occupation to spin notations
  //
  double Sz_quantum = 0.5 * (NUMBER_OF_STATES - 1);
  Nparticles = TotalSz + L * Sz_quantum;

  numberparticle_for_basis_state.resize(NUMBER_OF_STATES, 0);
  sz_for_basis_state.resize(NUMBER_OF_STATES, 0);
  for (int r = 0; r < NUMBER_OF_STATES; ++r) {
    numberparticle_for_basis_state[r] = r;
    sz_for_basis_state[r] = r - Sz_quantum;
  }

  // Count the number of valid sectors of the form (na,nb) with na+nb=Nparticles
  // , where na (nb) is the number of particles on A (B)
  valid_sectors = 0;
  for (int nleft = 0; nleft <= (LA * NUMBER_OF_STATES); ++nleft) {
    int nright = Nparticles - nleft;
    if ((nright >= 0) && (nright <= (LB * NUMBER_OF_STATES))) {
      particle_sector[std::make_pair(nleft, nright)] = valid_sectors;
      valid_sectors++;
    }
  }

  // now
  number_conf_in_A.resize(valid_sectors);
  number_conf_in_B.resize(valid_sectors);
  InverseMapA.resize(valid_sectors);
  Confs_in_A.resize(valid_sectors);
  InverseMapB.resize(valid_sectors);
  Confs_in_B.resize(valid_sectors);

  //  This vector will contain the # of the first configuration in each valid
  //  sector (which is the sum of all valid configurations in previous sectors)
  starting_conf.resize(valid_sectors + 1);

  // basic trick to find the minimal number of bits needed to encode a basis
  // state on one site for spin 1/2, correct_power_of_two = 1 (2 states), for
  // spin 1, correct_power_of_two = 2 (3 states) etc
  int correct_power_of_two = 1;
  while (pow(2, correct_power_of_two) < NUMBER_OF_STATES) {
    correct_power_of_two++;
  }

  // Now we iterate over all integers from 0 to 2^LA (for spin 1/2) and check if
  // their bit representation is a valid configuration for LA

  for (unsigned long int r = 0; r < pow(pow(2, correct_power_of_two), LA);
       ++r) {
    if (is_config_valid_inA(r, correct_power_of_two, LA,
                            LA * NUMBER_OF_STATES)) {
      // If yes, we find in which sector they belong (that is, we compute na)
      std::vector<unsigned short int> this_conf(LA, 0);
      int_to_occupation(r, correct_power_of_two, LA, this_conf);
      int nleft = 0;
      for (int p = 0; p < LA; ++p) {
        nleft += numberparticle_for_basis_state[this_conf[p]];
      }
      // of course the corresponding number of particles on the right is fixed
      // by the magnetization sector
      int nright = Nparticles - nleft;
      // find the corresponding (nleft,nright) sector
      std::map<std::pair<int, int>, int>::const_iterator it =
          particle_sector.find(std::make_pair(nleft, nright));
      //
      if (it != particle_sector.end()) {  // make sure this sector exists ...
        if ((InverseMapA[it->second]).count(this_conf) ==
            0) {  // make sure we did not already find it for some reason
          // Add this conf to the configuration list in this sector
          Confs_in_A[it->second].push_back(this_conf);
          // Add its index to the Inverse map of this sector
          InverseMapA[it->second][this_conf] = number_conf_in_A[it->second];
          // increase this sector size by 1
          number_conf_in_A[it->second] += 1;
        }
      }
    }
  }

  // Repeat the same strategy for the B subsystem
  for (unsigned long int r = 0; r < pow(pow(2, correct_power_of_two), LB);
       ++r) {
    if (is_config_valid_inA(r, correct_power_of_two, LB,
                            LB * NUMBER_OF_STATES)) {
      std::vector<unsigned short int> this_conf(LB, 0);
      int_to_occupation(r, correct_power_of_two, LB, this_conf);
      int nright = 0;
      for (int p = 0; p < LB; ++p) {
        nright += numberparticle_for_basis_state[this_conf[p]];
      }
      int nleft = Nparticles - nright;

      std::map<std::pair<int, int>, int>::const_iterator it =
          particle_sector.find(std::make_pair(nleft, nright));
      if (it != particle_sector.end()) {
        if ((InverseMapB[it->second]).count(this_conf) == 0) {
          Confs_in_B[it->second].push_back(this_conf);
          InverseMapB[it->second][this_conf] = number_conf_in_B[it->second];
          number_conf_in_B[it->second] += 1;
        }
      }
    }
  }
  // The total number of valid configurations is then given by the # of confs in
  // A times the # in B, summed over all sectors
  total_number_of_confs = 0;
  for (int nsa = 0; nsa < valid_sectors; ++nsa) {
    total_number_of_confs += Confs_in_A[nsa].size() * Confs_in_B[nsa].size();
  }

  // Computing the index of all first configurations in each sector
  starting_conf[0] = 0;
  for (int nsa = 1; nsa < valid_sectors; ++nsa) {
    starting_conf[nsa] =
        starting_conf[nsa - 1] +
        Confs_in_A[nsa - 1].size() * Confs_in_B[nsa - 1].size();
  }

  // detect the Neel states index
  int nleft = 0, nright = 0;
  std::vector<unsigned short int> neelA(LA, 0);
  for (int p = 1; p < LA; p += 2) {
    neelA[p] = NUMBER_OF_STATES - 1;
    nleft += 1;
  }
  std::vector<unsigned short int> neelB(LB, 0);
  /* if LA is even, the Néel state is up on the odd-numbered sites of B,
  otherwise it's the opposite */
  int p0 = LA + 1;
  if (LA % 2 == 1) p0 = LA;
  for (int p = p0; p < L; p += 2) {
    neelB[p - LA] = NUMBER_OF_STATES - 1;
    nright += 1;
  }
  int nsa = particle_sector[std::make_pair(nleft, nright)];
  ineel = starting_conf[nsa] +
          InverseMapA[nsa][neelA] * Confs_in_B[nsa].size() +
          InverseMapB[nsa][neelB];

  // now do the opposite Neel
  nleft = 0;
  nright = 0;
  std::vector<unsigned short int> neel2A(LA, 0);
  for (int p = 0; p < LA; p += 2) {
    neel2A[p] = NUMBER_OF_STATES - 1;
    nleft += 1;
  }
  std::vector<unsigned short int> neel2B(LB, 0);
  /* if LA is even, the Néel state is up on the odd-numbered sites of B,
  otherwise it's the opposite */
  p0 = LA;
  if (LA % 2 == 1) p0 = LA + 1;
  for (int p = p0; p < L; p += 2) {
    neel2B[p - LA] = NUMBER_OF_STATES - 1;
    nright += 1;
  }
  nsa = particle_sector[std::make_pair(nleft, nright)];
  int ineel2 = starting_conf[nsa] +
               InverseMapA[nsa][neel2A] * Confs_in_B[nsa].size() +
               InverseMapB[nsa][neel2B];
}

unsigned long int basis::index(Conf conf) {
  /*
   * Return the index (in the basis) of a given conf
   */
  Conf confA(conf.begin(), conf.begin() + LA);
  Conf confB(conf.begin() + LA, conf.end());
  // count the number of particles in the left and right subsystems
  int nleft = 0, nright = 0;
  for (auto& Sz : confA) nleft += Sz;
  for (auto& Sz : confB) nright += Sz;
  // sector number
  int nsa = particle_sector[std::make_pair(nleft, nright)];
  // index
  unsigned long int iconf = starting_conf[nsa] +
                            InverseMapA[nsa][confA] * Confs_in_B[nsa].size() +
                            InverseMapB[nsa][confB];
  return iconf;
}

// Never used ?
void basis::get_conf_coefficients(unsigned long int index, int& nsa, int& ca,
                                  int& cb) {
  // This function is ugly and non-optimal.  Fortunately we rarely use it

  // using the fact that index = cb + LB*(ca + LA*nsa)
  // int reduced_index = index/LB;
  // cb = index%LB;
  // nsa = reduced_index/LA;
  // ca = reduced_index%LA;
  unsigned long int i = 0;
  for (nsa = 0; nsa < valid_sectors; ++nsa) {
    for (ca = 0; ca < Confs_in_A[nsa].size(); ++ca) {
      for (cb = 0; cb < Confs_in_B[nsa].size(); ++cb) {
        if (i == index) return;
        i++;
      }
    }
  }
}

// Never used ?
void basis::init_vectors_product_state(Vec& Psi_t, Vec& Psi_t2, int i0) {
  // add 1 at position i0 in vector Psi_t
  VecSetValue(Psi_t, i0, 1.0, INSERT_VALUES);

  int i = 0;

  for (int nsa = 0; nsa < valid_sectors;
       ++nsa) {  // int nsb=partner_sector[nsa];
    for (int ca = 0; ca < Confs_in_A[nsa].size(); ++ca) {
      for (int cb = 0; cb < Confs_in_B[nsa].size(); ++cb) {
        if (i == i0) {
          // prefactor = state of the first spin of configuration B(nsa, cb)
          PetscScalar prefactor =
              numberparticle_for_basis_state[Confs_in_B[nsa][cb][0]];
          VecSetValue(Psi_t2, i0, prefactor, INSERT_VALUES);
        }
        i++;
      }
    }
  }
}

// For this demo code, we don't want to use boost::dynamic_bitset to avoid
// installing boost so we use the ugly fixed-size bitset, assuming a maximum of
// 32 bits to encode basis states for a subsystem.
// The commented lines use boost::dynamic_bitset and provide a more flexible
// code

// void basis::int_to_occupation(unsigned long int num, std::vector<unsigned
// short int>& occupation) {
void basis::int_to_occupation(unsigned long int num, int correct_power_of_two,
                              int L,
                              std::vector<unsigned short int>& occupation) {
  std::bitset<32> bb(num);
  // boost::dynamic_bitset<> bb(correct_power_of_two * L, num);
  for (int r = 0; r < L; ++r) {
    occupation[r] = 0;
    for (int x = 0; x < correct_power_of_two; ++x) {
      occupation[r] += pow(2, x) * bb[correct_power_of_two * r + x];
    }
  }
}

bool basis::is_config_valid_inA(int num, int correct_power_of_two, int L,
                                int Nmax) {
  std::bitset<32> bb(num);
  // boost::dynamic_bitset<> bb(correct_power_of_two * L, num);
  int total_n = 0;

  for (int r = 0; r < L; ++r) {
    int nn = 0;
    for (int x = 0; x < correct_power_of_two; ++x) {
      nn += pow(2, x) * bb[correct_power_of_two * r + x];
    }
    total_n += nn;
    if (total_n > Nmax) {
      return 0;
    }
  }
  return 1;
}

#endif
