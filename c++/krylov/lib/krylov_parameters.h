#ifndef PARAM_H
#define PARAM_H

using namespace std;

class Parameters {
 private:
 public:
  Parameters(int myrank_);
  ~Parameters() {}

  void Initialize_timegrid();
  std::vector<double> time_points;
  std::vector<double> delta_t_points;

  int myrank;

  PetscBool measure_entanglement;
  PetscBool measure_local;
  PetscBool measure_imbalance;

  PetscInt num_times;
  PetscReal Tmin;
  PetscReal Tmax;
  PetscBool use_linear_timegrid;
  PetscBool loggrid;
  PetscReal dt;

  PetscBool cdw_start;
  PetscBool product_state_start;
  PetscInt num_product_states;
  PetscReal TEEmin;
  PetscReal TEEmax;
  PetscInt nmeasures;
};

Parameters::Parameters(int myrank_) {
  myrank = myrank_;

  num_times = 100;
  Tmin = 0.;
  Tmax = 1.;
  use_linear_timegrid = PETSC_TRUE;
  loggrid = PETSC_FALSE;
  // i01=-1; i02=-1;

  cdw_start = PETSC_FALSE;
  product_state_start = PETSC_FALSE;
  TEEmin = Tmin;
  TEEmax = Tmax;
  nmeasures = num_times;

  measure_entanglement = PETSC_FALSE;
  measure_imbalance = PETSC_FALSE;
  measure_local = PETSC_FALSE;
  dt = 1.0;

  PetscOptionsGetInt(NULL, NULL, "-num_times", &num_times, NULL);
  PetscOptionsGetBool(NULL, NULL, "-loggrid", &loggrid, NULL);
  PetscOptionsGetReal(NULL, NULL, "-Tmax", &Tmax, NULL);
  PetscOptionsGetReal(NULL, NULL, "-Tmin", &Tmin, NULL);
  if (loggrid) {
    use_linear_timegrid = PETSC_FALSE;
    if (Tmin == 0) {
      Tmin = 1.;
    }
  }
  PetscOptionsGetBool(NULL, NULL, "-product_state_start", &product_state_start,
                      NULL);
  num_product_states = 1;
  PetscOptionsGetInt(NULL, NULL, "-num_product_states", &num_product_states,
                     NULL);
  PetscOptionsGetBool(NULL, NULL, "-cdw_start", &cdw_start, NULL);
  TEEmin = Tmin;
  TEEmax = Tmax;
  nmeasures = num_times;
  PetscOptionsGetInt(NULL, NULL, "-num_measures", &nmeasures, NULL);
  PetscOptionsGetBool(NULL, NULL, "-measure_entanglement",
                      &measure_entanglement, NULL);
  PetscOptionsGetBool(NULL, NULL, "-measure_imbalance", &measure_imbalance,
                      NULL);
  PetscOptionsGetBool(NULL, NULL, "-measure_local", &measure_local, NULL);

  int each_measurement = num_times / nmeasures;

  Initialize_timegrid();
}

void Parameters::Initialize_timegrid() {
  double b = 1;
  for (int kk = 0; kk <= num_times; ++kk) {
    double Tj, Deltat;
    if (false == use_linear_timegrid) {
      b = exp(log(Tmax / Tmin) / static_cast<double>(num_times));
      Tj = Tmin * pow(b, kk);
      Deltat;
      if (kk == 0) {
        Deltat = Tj;
      } else {
        Deltat = Tj * (1. - 1. / b);
      }
    } else {
      if (kk == 0) {
        Deltat = Tmin;
        Tj = Tmin;
      } else {
        Deltat = (Tmax - Tmin) / static_cast<double>(num_times);
        Tj = Tmin + kk * Deltat;
      }
    }
    time_points.push_back(Tj);
    delta_t_points.push_back(Deltat);
  }

  if (myrank == 0) {
    std::cout << "#Time evolution from " << Tmin << " to " << Tmax << ", with "
              << num_times << " points on a ";
    if (use_linear_timegrid) {
      std::cout << "linear";
    } else {
      std::cout << "logarithmic";
    }
    std::cout << " time grid\n";
  }
}

#endif
