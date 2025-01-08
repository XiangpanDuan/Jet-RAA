#ifndef ENERGYLOSSOMEGAC_H
#define ENERGYLOSSOMEGAC_H

#include <iostream>
#include <vector>
#include <cmath>
#include <gsl/gsl_integration.h>


class EnergyLossOmegaC{
  
 private:

  int    _xnum,_ynum,_taunum,_thetanum;  //table size from 2+1D VISHydro
  double _xmin,_xmax;
  double _ymin,_ymax;
  double _xbin,_ybin,_taubin;
  double _x0,_y0,_tau;
  double _x,_y,_theta;
  std::vector<std::vector<std::vector<double>>> _TempTable;
  double _Temp0,_Temp;
  double _omegaC;


 public:

  EnergyLossOmegaC();
  ~EnergyLossOmegaC(){}
  
  //Hydro cut
  bool HydroTaucut(const double tau);
  bool HydroTempcut(const double temp);
  //Temperature table
  inline void setTemp0(const double temp0) {_Temp0=temp0;}
  inline void setTemp(const int ix, const int iy, const int itau, const double temp) {_TempTable[ix][iy][itau]=temp;}
  //Temperature
  inline void setxytheta(const double x0, const double y0, const double theta) {_x0=x0; _y0=y0; _theta=theta;}
  double TriInterpolation(const double x, const double y, const double tau);
  double Temperature(const double tau);
  //OmegaC
  void setOmegaC(const double tau);
  inline double getOmegaC() const {return _omegaC;}
  
  
  //Integration function must keeps in class with static member function
  static double FunctionOmegaC(double tau, void *params){
    EnergyLossOmegaC *param = (EnergyLossOmegaC*) params;
    param->setOmegaC(tau);
    return param->getOmegaC();
  }

  //Calculate omegaC with start time
  double CalculateOmegaC(const double taustart){
    //QAG adaptive integration
    double taumin=taustart;
    double taumax=20.;
    gsl_function Fun;
    Fun.function = &FunctionOmegaC;
    Fun.params = this;
    gsl_integration_workspace *workspace = gsl_integration_workspace_alloc(10000);
    double res=0.0, err=0.0;
    gsl_integration_qag(&Fun, taumin, taumax, 0, 1e-3, 10000, 6, workspace, &res, &err);
    gsl_integration_workspace_free(workspace);
    return res;
  }

};

#endif  //ENERGYLOSSOMEGAC_H
