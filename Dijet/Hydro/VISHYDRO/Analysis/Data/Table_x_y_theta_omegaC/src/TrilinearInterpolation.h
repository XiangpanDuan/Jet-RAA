#ifndef TRILINEARINTERPOLATION_H
#define TRILINEARINTERPOLATION_H

#include <iostream>
#include <fstream>
#include <sstream>
#include <iomanip>
#include <cmath>
#include <vector>
#include <gsl/gsl_math.h>
#include <gsl/gsl_integration.h>
#include <gsl/gsl_monte.h>
#include <gsl/gsl_monte_vegas.h>
// #define MC_MISER
// #include <gsl/gsl_monte_miser.h>


class TrilinearInterpolation{

 private:

  int _xnum,_ynum,_taunum,_thetanum;
  double _xmin,_xmax;
  double _ymin,_ymax;
  double _xbin,_ybin,_taubin;
  std::vector<std::vector<std::vector<double>>> _Temp;
  double _qhat0,_Temp0,_TempC,_omegaC;
  double _x0,_y0,_tau;
  double _x,_y,_theta;

  //Variables for QAG adaptive integration
  gsl_integration_workspace *_workspace;
  gsl_function _Fun;
  double _rangemin, _rangemax;  //integration boundary
  

 public:

  TrilinearInterpolation();
  TrilinearInterpolation(const double qhat0);
  ~TrilinearInterpolation(){}
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

#endif  //TRILINEARINTERPOLATION_H
