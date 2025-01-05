#ifndef OMEGAC_H
#define OMEGAC_H

#include <iostream>
#include <fstream>
#include <vector>
#include <cmath>
#include <gsl/gsl_integration.h>


class OmegaC{

 private:

  int _xnum,_ynum,_taunum,_thetanum;  //table size from 2+1D VisHydro (_xnum=nx,_ynum=ny,_taunum=ntau,_thetanum=ntheta)
  double _xmin,_xmax;
  double _ymin,_ymax;
  double _xbin,_ybin,_taubin;
  std::vector<std::vector<std::vector<double>>> _Temp;
  double _Temp0;
  double _TempC,_omegaC;
  double _x0,_y0,_tau;
  double _x,_y,_theta;
  double _qhat0;  //BDMPS: 2-8 GeV^2/fm, High Twist or GLV: 0.3-1 GeV^2/fm

  //Variables for QAG adaptive integration
  gsl_integration_workspace *_workspace;
  gsl_function _Fun;
  double _rangemin, _rangemax;  //integration boundary
  

 public:

  OmegaC();
  OmegaC(const double qhat0);
  ~OmegaC(){}
  void setInitialData();

  bool   HydroTaucut (const double tau);
  bool   HydroTempcut(const double temp);
  void   setInputData(const double ix, const double iy, const double itheta);
  double TriInterpolation(const double xfin, const double yfin, const double tauin);
  double Temperature(const double tau);
  void   setOmegaC(const double tau);
  inline double getOmegaC() const {return _omegaC;};

  //QAG adaptive integration
  void setupQAG(double (*fun)(double, void *), double rangemin, double rangemax, void *para);
  void calculateQAG(double &res, double &err);
  void cleanupQAG();
};

#endif  //OMEGAC_H
