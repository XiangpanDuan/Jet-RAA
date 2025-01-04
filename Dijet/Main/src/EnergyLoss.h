#ifndef ENERGYLOSS_H
#define ENERGYLOSS_H

#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <cmath>
#include <gsl/gsl_sf_bessel.h>
#include "DiJetLO.h"
#include "EnergyLossOmegaC.h"


class EnergyLoss: public DiJetLO{
 
 private:

  EnergyLossOmegaC *_elossOmegaC;
  int    _A,_B;                            //A:Target, B:Projectile
  double _Cen;                             //centrality
  double _Ecm;                             //collisional energy in CoM frame
  double _Rsize;                           //jet cone size
  int    _iMul,_fMul;                      //initial and final integral dimension in jet energy loss
  double _nMul;;                           //number of same jet energy loss
  double _qhat0;                           //BDMPS: 2-8 GeV^2/fm, High Twist or GLV: 0.3-1 GeV^2/fm
  double _alphaSmed;                       //medium scale in BDMPS
  double _Qmed;                            //medium scale in MultiplicityBessel
  double _nDepsilon;                       //total 
  std::vector<double> _epsilon,_Depsilon;  //jet energy loss for multi-partons
  int    _xnum,_ynum,_taunum,_thetanum;
  double _xmin,_xmax;
  double _ymin,_ymax;
  double _omegaC;
  std::vector<std::vector<std::vector<double>>> _omegaCTable;  //omegaC table
  double _tau0;  //start time for jet energy loss
 
  
 public:

  EnergyLoss();
  EnergyLoss(const int A, const int B, const double Cen, const double Ecm, Particle *Parton);
  ~EnergyLoss();

  void setpT(const double pT);
  inline void setqhat0(const double qhat0)   {_qhat0=qhat0;}                    //jet quenching parameter
  inline void setOmegaC(const double omegaC) {_omegaC=omegaC;}                  //energy loss parameter
  inline void setalphaSmed(const double alphaSmed) {_alphaSmed=alphaSmed;}      //medium scale in BDMPS
  inline void setRsize(const double Rsize)   {_Rsize=Rsize;}                    //jet cone size
  inline void setQmed(const double Qmed)     {_Qmed=Qmed;}                      //medium scale in MultiplicityBessel
  inline void setMul(const int iMul, const int fMul) {_iMul=iMul; _fMul=fMul;}  //initial and final jet multiplicity

  //Start time with qhat and jet pT dependences
  void setTauStart();
  //Initial data (x,y,tau,temp) produced by 2+1D VISHydro
  void setTempTable();
  //QAG integration to get omegaC
  void setOmegaCTau(const double x0, const double y0, const double theta);
  //Initial data (x,y,theta,omegaC) produced by 2+1D VISHydro
  void setOmegaCTable();
  //Trilinear interpolation to get omegaC
  bool setOmegaC(const double x, const double y, const double theta);
  inline double getOmegaC() const {return _omegaC;}

  //Jet energy loss when traversing QGP in BDMPS formalism
  void setEpsilon(const double *epsilon);
  inline double getEpsilon(const double i) const {return _epsilon[i];}
  //One parton (quark or gluon) jet energy loss
  void setDepsilon(const std::string &type);
  void setDepsilonGluon();
  void setDepsilonGluonMean();
  //Multi-gluons jet energy loss
  inline double getDepsilon() const {return _nDepsilon;}
  
  //Multiplicity in jet energy loss
  void setMultiplicity(const std::string &type, const double pT);
  void setMultiplicityBessel(const std::string &type, const double pT);
  inline double getMultiplicit() const {return _nMul;}

};

#endif  //ENERGYLOSS_H
