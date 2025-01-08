#include <string>
#include <cmath>

#define  MODE   0                    //0:pp collisions; 1:AA collisions
#define  GLB    0                    //0:turn off Glauber model; 1:turn on Glauber model
#define  MULTI  0                    //0:turn off multiplicity; 1:turn on mean multiplicity; 2:turn on multiplicity with probability distribution
#define  TAUS   0                    //0:turn off start time (default tau0=0.6); 0:turn on start time with qhat and jet pT dependences


namespace Input{

//Initial conditions
#if   MODE == 1
const int A=208;                     //Target
const int B=208;                     //Projectile
const double Cen=0.0;                //Centrality
#endif
const double Ecm=5020.;              //GeV, collisional energy in CoM frame
const std::string name="parton";     //particle name in hard scattering: parton,d,u,s,c,b,t,g,photon,Z0,Wp,Higgs
const int scale=0;                   //0,+1,-1: pT scale to control the pT errer bar in pdf and alphas from LHAPDF
const unsigned int nf=3;             //quark flavors
const unsigned int nloop=1;          //1:LO calculation; 2:NLO calculation

//Monte Carlo calls
const size_t calls=20000;

//Range of integration
const int    pTnum=100;
const double pTmin=40.0;             //GeV, transverse momentum ranges of final observed jet
const double pTmax=1040.;
const double rap3 =2.8;              //rapidity of final triggered jet
const double rap4 =2.8;              //rapidity of final observed jet 
#if   MODE == 1
const double elossmin=0.0;
const double elossmax=100.;
#if   GLB == 1
const double xmin =-15.;             //fm, QGP size, same as Hydro
const double xmax = 15.;
const double ymin =-15.;
const double ymax = 15.;
const double thetamin=0.0;
const double thetamax=2.*M_PI;
#endif
#endif
 
//Energy loss
#if   MODE == 1
const double Rsize=0.4;              //jet cone size
const double alphaSmed=0.2;          //medium scale in BDMPS
const double Qmed=15.;               //medium scale in Multiplicity, used in MULTI=1
#if   MULTI != 2
const int    nEloss=1;               //numbers of energy loss
#elif MULTI == 2
const int    nEloss=10;              //numbers of energy loss
#endif
#if   GLB == 0
const double omegaC=6.5;             //GeV; fixed energy loss in BDMPS; 14.:MULTI=0; (4.5,4.3,4.0,2.6):MULTI=1; (15,12,10):MULTI=2 
#elif GLB == 1
const double qhat0=2.5;              //GeV/fm^2; 6.7:MULTI=0; 1.2:MULTI=1; ???:MULTI=2
#endif
#endif

}


////////////////////////////////////////////////////////////
//Instructions
//MODE=0  and GLB=0:     make pp.exe
//MODE=1  and GLB=0:     make AAFixed.exe
//MODE=1  and GLB=1:     make AA.exe"
//MULTI=0 and nEloss=1:  single parton energy loss
//MULTI=1 and nEloss=1:  single parton with multi mean energy loss
//MULTI=2 and nEloss=10: multi-partons with different energy loss
