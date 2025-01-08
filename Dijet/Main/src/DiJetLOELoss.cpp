#include "DiJetLOELoss.h"


DiJetLOELoss::DiJetLOELoss(): EnergyLoss(){
  _glauber = new Glauber();
}

DiJetLOELoss::DiJetLOELoss(const int A, const int B, const double Cen, const double Ecm, Particle *Parton): EnergyLoss(A,B,Cen,Ecm,Parton){
  _glauber = new Glauber(A,B,Cen);
}

DiJetLOELoss::~DiJetLOELoss(){
  delete _glauber;
}


//####################################################################################################
//Thickness function in Glauber model
void DiJetLOELoss::setGlauber(const double x, const double y){
  _glauber->setThicknessTAB(x,y);
}

double DiJetLOELoss::getGlauber(){
  // if(GLB==1) return _glauber->getThicknessTAB()/(2.*M_PI);  //get _TAB0
  if(GLB==1) return _glauber->getThicknessTAB()/(2.*M_PI*_glauber->getThicknessTAB0());  //normalization
  if(GLB==0) return 1.;  //without the geometry effect
}


//####################################################################################################
//Energy loss in BDMPS
void DiJetLOELoss::setEnergyLoss(const std::string &type, const double *range, const size_t dim){
  if(MODE==1){
    int idim;
    if(GLB==0) idim=3;                                     //without Glauber geometric effect
    if(GLB==1) idim=6;                                     //with Glauber geometric effect
    setMul(idim,dim);                                      //initial and final integral dimension in jet energy loss
    setEpsilon(range);                                     //energy loss ε1,ε2,ε3,...,εn
    if(MULTI==1) setMultiplicity(type,range[0]);           //mean multiplicity
    setpT(range[0]);                                       //jet pT with energy loss
    if(GLB==1){
      if(TAUS==0){setOmegaC(range[3],range[4],range[5]);}  //OmegaC with initial time (tau0=0.6)
      if(TAUS==1){
        setTauStart();                                     //start time with qhat and jet pT dependences
        setOmegaCTau(range[3],range[4],range[5]);          //OmegaC with start time dependence
      }
    }
    if(MULTI==0){setDepsilon(type);}                       //primary energy loss with BDMPS D(ε1)
    // if(MULTI==0){setDepsilonGluonMean();}                  //primary energy loss with BDMPS D(ε1)
    // if(MULTI==1){setDepsilon(type);}                       //mean-multi gluon energy loss when setMultiplicity(type,range[0]) function is used before
    if(MULTI==1){setDepsilonGluonMean();}                  //mean-multi gluon energy loss when setMultiplicity(type,range[0]) function is used before
    if(MULTI==2){
      setDepsilon(type);                                   //first  energy loss with BDMPS D(ε1)
      setDepsilonGluon();                                  //second energy loss to final energy loss with BDMPS D(ε2)*D(ε3)*...*D(εn)
    }
  }

}

double DiJetLOELoss::getEnergyLoss(){
  return getDepsilon();
}


//####################################################################################################
//Kinematic parameters in leading order
void DiJetLOELoss::setParameters(const double yTrig3, const double yAsso4){
  return setLOParameters(yTrig3, yAsso4);
}


//####################################################################################################
//Leading order cross section
double DiJetLOELoss::getSigmaLO(){
  return SigmaLO();
}
double DiJetLOELoss::getSigmaLOQuark(){
  return SigmaLOQuark();
}
double DiJetLOELoss::getSigmaLOGluon(){
  return SigmaLOGluon();
}


//####################################################################################################
double DiJetLOELoss::SigmaLOELoss(){
  return getSigmaLO()*getEnergyLoss()*getGlauber();
}

double DiJetLOELoss::SigmaLOELossQuark(){
  return getSigmaLOQuark()*getEnergyLoss()*getGlauber();
}

double DiJetLOELoss::SigmaLOELossGluon(){
  return getSigmaLOGluon()*getEnergyLoss()*getGlauber();
}
