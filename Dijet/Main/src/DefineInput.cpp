#include <string>
#include <cmath>

#define MODE  0                      //0:pp collisions; 1:AA collisions
#define GLB   0                      //0:turn off Glauber model; 1:turn on Glauber model
#define MULTI 0                      //0:turn off multiplicity effect; 1:turn on mean multiplicity effect; 2:turn on multiplicity with probability distribution
#define TAUS  0                      //0:turn off start time (default tau0=0.6); 0:turn on start time with qhat and jet pT dependences


namespace Input{

//Initial conditions
#if   MODE == 1
const int A=208;                     //Target
const int B=208;                     //Projectile
const double Cen=0.0;                //Centrality
#endif
const double Ecm=5020.;              //GeV
const std::string name="parton";     //particle name: parton,d,u,s,c,b,t,g,photon,Z0,Wp,Higgs
const int scale=0;                   //0,+1,-1: pT scale to control the pT errer bar in pdf and alphas from LHAPDF
const unsigned int nf=3;             //quark flavor
const unsigned int nloop=1;          //1:LO calculation; 2:NLO calculation

//Monte Carlo calls
const size_t calls=20000;

//Range of integration
const int    pTnum=100;              //unit:GeV
const double pTmin=40.0;
const double pTmax=1040.;
const double rap3 =2.8;              //triggered jet rapidity
const double rap4 =2.8;              //observed  jet rapidity
#if   MODE == 1
const double elossmin=0.0;
const double elossmax=100.;
#if   GLB == 1
const double xmin =-15.;             //unit:fm
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
const double Qmed=10.;               //medium scale in MultiplicityBessel
#if   MULTI != 2 
const int    nEloss=1;               //numbers of energy loss
#elif MULTI == 2
const int    nEloss=10;              //numbers of energy loss
#endif
#if   GLB == 0
const double omegaC=14.;             //unit:GeV; fixed energy loss in BDMPS; 14.:MULTI=0; (4.5,4.3,4.0,2.6):MULTI=1; (15,12,10):MULTI=2 
#elif GLB == 1
const double qhat0 =1.0;             //unit:GeV/fm^2; 6.6:MULTI=0; 1.2:MULTI=1; ???:MULTI=2
#endif
#endif

}


////////////////////////////////////////////////////////////
//Instructions
//MODE=0 needs to "make pp.exe"
//MODE=1 and GLB=0 need to "make AAFixed.exe"
//MODE=1 and GLB=1 need to "make AA.exe"
//MULTI=0 and nEloss=1  mean single parton energy loss
//MULTI=1 and nEloss=1  mean single parton with multi same energy loss
//MULTI=2 and nEloss=15 mean multi-partons with different energy loss
