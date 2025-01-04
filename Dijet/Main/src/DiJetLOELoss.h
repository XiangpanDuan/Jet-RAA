#ifndef DIJETLOELOSS_H
#define DIJETLOELOSS_H

#include <iostream>
#include <cmath>
#include "Glauber.h"
#include "EnergyLoss.h"


class DiJetLOELoss: public EnergyLoss{
 
 private:

  Glauber *_glauber;  //glauber model
  

 public:

  DiJetLOELoss();
  DiJetLOELoss(const int A, const int B, const double Cen, const double Ecm, Particle *Parton);
  ~DiJetLOELoss();

  void   setGlauber(const double x, const double y);
  double getGlauber();
  void   setEnergyLoss(const std::string &type, const double *range, const size_t dim);
  double getEnergyLoss();
  void   setParameters(const double yTrig3, const double yAsso4);

  void   setSigmaLO(const double yTrig3, const double yAsso4);
  double getSigmaLO();
  double getSigmaLOQuark();
  double getSigmaLOGluon();

  double SigmaLOELoss();
  double SigmaLOELossQuark();
  double SigmaLOELossGluon();
  
};

#endif  //DIJETLOELOSS_H
