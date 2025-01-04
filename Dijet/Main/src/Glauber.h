#ifndef GLAUBER_H
#define GLAUBER_H

#include <iostream>
#include <vector>
#include <cmath>
#include <gsl/gsl_math.h>
#include <gsl/gsl_integration.h>


class Glauber{

 private:

  int    _A,_B;  //A:Target, B:Projectile
  double _Cen;
  double _rho0_A,_R_A,_a_A,_omega_A;
  double _rho0_B,_R_B,_a_B,_omega_B;
  double _rho_A,_rho_B;
  double _TA,_TB,_TAB,_TAB0;
  int    _numGauss;
  std::vector<double> _zi,_wi;
  

 public:
  
  Glauber();
  Glauber(const int A, const int B, const double Cen);
  ~Glauber(){}
  void setNuclearParameters(const int A, const int B);
  void setGaussianQuadrature();  //Gaussian-Legendre Quadrature
  void setThicknessTAB0(const int A, const int B);
  inline double getThicknessTAB0() const {return _TAB0;}
  
  void setDensityRhoA(const double x,const  double y, const double z);  //Woods-Saxon density distributions
  void setDensityRhoB(const double x,const  double y, const double z);  //Woods-Saxon density distributions
  inline double getDensityRhoA() const {return _rho_A;}
  inline double getDensityRhoB() const {return _rho_B;}
  void setThicknessTA(const double sx, const double sy);   //TA
  void setThicknessTB(const double sx, const double sy);   //TB
  inline double getThicknessTA() const {return _TA;}
  inline double getThicknessTB() const {return _TB;}
  void setThicknessTAB(const double sx, const double sy);  //TAB with certain b=0
  void setThicknessTAB(const double sx, const double sy, const double bx, const double by);  //TAB: Thickness function
  inline double getThicknessTAB() const {return _TAB;}

};

#endif  //GLAUBER_H
