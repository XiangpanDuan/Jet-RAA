#include <iostream>
#include <fstream>
#include <sstream>
#include <iomanip>
#include <cmath>
#include <vector>
#include <gsl/gsl_math.h>
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


  //Variables for Monte Carlo
  const gsl_rng_type *_T;
  gsl_rng *_r;
#ifdef MC_MISER
  gsl_monte_miser_state * _mcStatus;
#else
  gsl_monte_vegas_state * _mcStatus;
#endif
  gsl_monte_function _mcFun;
  size_t _dim, _calls;
  double *_rangemin, *_rangemax; //integration boundary
  

public:

  TrilinearInterpolation();
  TrilinearInterpolation(const double qhat0);
  void setInitialData();

  bool   HydroTaucut(const double tau);
  bool   HydroTempcut(const double temp);
  void   setInputData(const double ix, const double iy, const double itheta);
  double TriInterpolation(const double xfin, const double yfin, const double tauin);
  double Temperature(const double tau);
  void   setOmegaC(const double tau);
  double getOmegaC();


  //Monte Carlo
  void setupMC(double (*pds)(double *, size_t, void *), size_t dim, double *rangemin, double *rangemax, void *para);
  void cleanupMC();
  void setCalls(size_t calls);
  void calculate(double &res, double &err);
};

