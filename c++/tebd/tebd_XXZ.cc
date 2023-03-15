static char help[] =
    "Demo code of TEBD time evolution for XXZ chain with random fields \n"
  "Uses Itensor, and largely inspired from Itensor tutorials \n"
    "For Les Houches summer school 2019 \n"
    "2019, Fabien Alet and Nicolas Macé\n";

#include "itensor/all.h"


using std::vector;
using std::move;
using namespace itensor;


BondGate CreateXXZrandomfield_Gate(const SpinHalf &sites,
			       int i, int j, int N, double hi, double hj, double localtimestep)
{
  // be careful with site term for the boundary spins
  if ( (i==1 || i==N) ) { hi*=2;} 
  if ( (j==1 || j==N) ) { hj*=2;} 
auto hterm = sites.op("Sz",i)*sites.op("Sz",j);
	      hterm += 0.5*sites.op("S+",i)*sites.op("S-",j);
	      hterm += 0.5*sites.op("S-",i)*sites.op("S+",j);
	      hterm -= 0.5*hi*sites.op("Sz",i)*sites.op("Id",j);
	      hterm -= 0.5*hj*sites.op("Id",i)*sites.op("Sz",j);

	      return BondGate(sites,i,j,BondGate::tReal,localtimestep,hterm);
}
	     

int
main()
    {

int N = 20; //number of sites
// Assume N even for Trotter-suzuki below ...

Real tstep = 0.005; //time step (smaller is generally more accurate)
Real ttotal = 100; //total time to evolve
 Real dt_measures=1; // We measure every dt_measures
Real cutoff = 1E-12; //truncation error cutoff when restoring MPS form
 int trotter=2; // order of trotter decomposition (2 and 4 available at the moment)

 double disorder=5.0;
 int seed = 74310;
 std::vector<double> fields(N+1,0.);
 
if (disorder > 0) {
    std::random_device device;
    std::mt19937 generator(device());
    generator.seed(seed);
    std::uniform_real_distribution<double> box(-disorder, disorder);
    for (int i = 1; i <=N; i++) {
      fields[i] += box(generator);
    }
 }


 std::cout << "#### Starting TEBD computation for XXZ random field chain\n";
 std::cout << "# Total evolution time = " << ttotal << "\n";
std::cout << "# TEBD time step = " <<  tstep<< "\n";
std::cout << "# Measures every dt = " << dt_measures << "\n";
std::cout << "# Cutoff for SVD's  = " << cutoff << "\n";
std::cout << "# Trotter Suzuki order  = " << trotter << "\n";

std::cout << "# N = " << N << "\n";
std::cout << "# disorder = " << disorder << " , seed = " << seed << "\n";
std::cout << "# Fields = {";
for (int i = 1; i <=N; i++) { std::cout << fields[i] << " ";} std::cout << "}\n\n";

std::cout << "t Imbalance S_1\n";

auto sites = SpinHalf(N);

//Make initial MPS psi to be in the Neel state
auto state = InitState(sites);
for(auto j : range1(N))
    {
    state.set(j,j%2==1?"Up":"Dn");
    }

auto psi = MPS(state);

auto gates = vector<BondGate>();
	if (trotter==2) {
	  // first to the right
	  for(int b = 1; b <= N-1; ++b)
	    {
	      
	      gates.push_back(CreateXXZrandomfield_Gate(sites,b,b+1,N,fields[b],fields[b+1],tstep/2));
		
	    }
	  // then back
	  for(int b = N-1; b >= 1; --b)
    {
      
		  gates.push_back(CreateXXZrandomfield_Gate(sites,b,b+1,N,fields[b],fields[b+1],tstep/2));
		 
    }

	}

	if (trotter==4) { // Forest-Ruth formula
	  // https://doi.org/10.1016/0167-2789(90)90019-L

	  double theta=1./(2.-pow(2.,1./3.));
	  // exp(A theta dt/2)
	  for(int b = 1; b <= N-1; b+=2)
	    {
	      gates.push_back(CreateXXZrandomfield_Gate(sites,b,b+1,N,fields[b],fields[b+1],0.5*theta*tstep));
	    } 
	  // exp(B theta dt)
	  for(int b = N-2; b >=2; b-=2)
	    {
	      gates.push_back(CreateXXZrandomfield_Gate(sites,b,b+1,N,fields[b],fields[b+1],theta*tstep));
	    }		   
	  // exp(A (1-theta) dt/2)
		for(int b = 1; b <= N-1; b+=2)
	    {
	      gates.push_back(CreateXXZrandomfield_Gate(sites,b,b+1,N,fields[b],fields[b+1],0.5*(1-theta)*tstep));
	    }	   
   // exp(B (1-2theta) dt)
	  for(int b = N-2; b >=2; b-=2)
	    {
	      gates.push_back(CreateXXZrandomfield_Gate(sites,b,b+1,N,fields[b],fields[b+1],(1-2*theta)*tstep));
	    }
	     // exp(A (1-theta) dt/2)
		for(int b = 1; b <= N-1; b+=2)
	    {
	      gates.push_back(CreateXXZrandomfield_Gate(sites,b,b+1,N,fields[b],fields[b+1],0.5*(1-theta)*tstep));
	    } 	   
	   // exp(B theta dt)
	  for(int b = N-2; b >=2; b-=2)
	    {
	      gates.push_back(CreateXXZrandomfield_Gate(sites,b,b+1,N,fields[b],fields[b+1],theta*tstep));
	    }
		// exp(A theta dt/2)
	  for(int b = 1; b <= N-1; b+=2)
	    {
	      gates.push_back(CreateXXZrandomfield_Gate(sites,b,b+1,N,fields[b],fields[b+1],0.5*theta*tstep));
	    }	  
			   }

	
	double t=0;
	while(t<ttotal) {
	  
	  // start with measurements
	  
	  // do imbalance measurement ...
	  auto psip = psi; Real imbalance=0.;
	  
	  for(int b=1; b <= N; ++b)  {
	    psip.position(b); 
	    ITensor ket = psip.A(b);
	    ITensor bra = dag(prime(ket,"Site"));
	    ITensor Szop = sites.op("Sz",b);
	    auto szb = eltC(bra*Szop*ket).real();
	    if (b%2) { imbalance+=szb; } else { imbalance-=szb;}
	    
	  }

	  imbalance*=2./N;
	  
	  
	  // now Entanglement Entropy
	   psip = psi; psip.position(N/2);
	  ITensor wf = psip.A(N/2)*psip.A(N/2+1);
	  auto U = psip.A(N/2);
	  ITensor S,V;
	  auto spectrum = svd(wf,U,S,V);
	  Real SvN = 0.;
	  for(auto p : spectrum.eigs())
	    {
	      if(p > 1E-12) SvN += -p*log(p);
	    }
	  
	  std::cout << t << " " << imbalance << " " << SvN << "\n";
	  
	  
	  // do time-evolution
	  gateTEvol(gates,dt_measures,tstep,psi,{"Cutoff=",cutoff,"ShowPercent",false});
	  // we should try to get back the discarded weight...
	  t+=dt_measures;
	}
	  


    return 0;
    }
