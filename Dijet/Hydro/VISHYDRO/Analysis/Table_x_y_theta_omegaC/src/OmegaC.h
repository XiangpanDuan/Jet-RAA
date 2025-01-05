#ifndef OMEGAC_H
#define OMEGAC_H

#include <iostream>
#include <fstream>
#include <vector>
#include <cmath>
#include <gsl/gsl_integration.h>


class OmegaC{

 private:

  int    _xnum,_ynum,_taunum,_thetanum;  //table size from 2+1D VisHydro
  double _xmin,_xmax;
  double _ymin,_ymax;
  double _xbin,_ybin,_taubin;
  double _x0,_y0,_tau;
  double _x,_y,_theta;
  std::vector<std::vector<std::vector<double>>> _TempTable;
  double _Temp0,_Temp;
  double _omegaC;
  double _qhat0;  //jet quenching parameter (BDMPS: 2-8 GeV^2/fm; High Twist or GLV: 0.3-1 GeV^2/fm)

  //Variables for QAG adaptive integration
  gsl_integration_workspace *_workspace;
  gsl_function _Fun;
  double _rangemin, _rangemax;
  

 public:

  OmegaC();
  ~OmegaC(){}
  void setInitialData();

  //Hydro cut
  bool HydroTaucut (const double tau);
  bool HydroTempcut(const double temp);
  //Input setting
  inline void setqhat0(const double qhat0) {_qhat0=qhat0;}
  inline void setxytheta(const double x0, const double y0, const double theta) {_x0=x0; _y0=y0; _theta=theta;}
  //Temperature
  double TriInterpolation(const double x, const double y, const double tau);
  double Temperature(const double tau);
  //OmegaC
  void setOmegaC(const double tau);
  inline double getOmegaC() const {return _omegaC;}

  //QAG adaptive integration
  void setupQAG(double (*fun)(double, void *), double rangemin, double rangemax, void *para);
  void calculateQAG(double &res, double &err);
  void cleanupQAG();
};

#endif  //OMEGAC_H
