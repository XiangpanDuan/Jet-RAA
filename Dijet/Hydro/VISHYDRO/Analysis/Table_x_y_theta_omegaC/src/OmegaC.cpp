#include "OmegaC.h"


OmegaC::OmegaC(){
  _xnum=301;_ynum=301;_taunum=201;_thetanum=36;
  _xmin=-15.;_xmax=15.;_ymin=-15.;_ymax=15.;
  _xbin=(_xmax-_xmin)/(_xnum-1);
  _ybin=(_ymax-_ymin)/(_ynum-1);
  _taubin=0.1;
  _qhat0=1.0;
  _Temp0=0.0;
  _Temp=0.0;
  _omegaC=0.0;
  _TempTable.resize(_xnum, std::vector<std::vector<double>>(_ynum, std::vector<double>(_taunum, 0.0)));
  setInitialData();
}

//Read initial data (x,y,tau,temp) produced by 2+1D VisHydro
void OmegaC::setInitialData(){
  std::ifstream InputFile;
  InputFile.open("../../Table_x_y_tau_temp/Output/x_y_tau_temp.dat");
  double x_val,y_val,tau_val,temp_val;
  for(int ix=0; ix<_xnum; ix++){
    for(int iy=0; iy<_ynum; iy++){
      for(int itau=0; itau<_taunum; itau++){
        InputFile >> x_val >> y_val >> tau_val >> temp_val;
        _TempTable[ix][iy][itau]=temp_val;
      }
    }
  }
  InputFile.close();
  int xhalf=(int)(std::floor(_xnum/2));
  int yhalf=(int)(std::floor(_ynum/2));
  _Temp0=_TempTable[xhalf][yhalf][6];  //T0: x=0, y=0, tau0=0.6
  // std::cout << "Temp0 = " << _Temp0 << std::endl;
}


//####################################################################################################
//Hydro cut
//Initial time cut
bool OmegaC::HydroTaucut(const double tau){
  double tau0=0.6;  //initial time in hydro
  return tau>=tau0;
  // if(tau<tau0) return false;
  // if(tau>=tau0) return true;
}

//Freeze-out temperature cut
bool OmegaC::HydroTempcut(const double temp){
  double tempcut=0.154;  //from lattice QCD (PRD 90 (2014) 094503)
  return temp>=tempcut;
  // if(temp<tempcut) return false;
  // if(temp>=tempcut) return true;
}

double OmegaC::TriInterpolation(const double x, const double y, const double tau){
  int    xlow,xhigh,ylow,yhigh,taulow,tauhigh;
  double dx,dy,dtau;
  double C,C0,C1,C00,C01,C10,C11,C000,C001,C010,C011,C100,C101,C110,C111;
  xlow=(int)(std::floor(x/_xbin)); xhigh=xlow+1;
  ylow=(int)(std::floor(y/_ybin)); yhigh=ylow+1;
  taulow=(int)(std::floor(tau/_taubin)); tauhigh=taulow+1;
  if(xlow<(_xmin/_xbin) || xlow>=(_xmax/_xbin)) return 0.0;
  if(ylow<(_ymin/_ybin) || ylow>=(_ymax/_ybin)) return 0.0;
  if(taulow>=_taunum) return 0.0;
  dx=(x/_xbin-xlow)/(xhigh-xlow);
  dy=(y/_ybin-ylow)/(yhigh-ylow);
  dtau=(tau/_taubin-taulow)/(tauhigh-taulow);
  //attention x and y range
  int xhalf=(int)(std::floor(_xnum/2));
  int yhalf=(int)(std::floor(_ynum/2));
  C000=_TempTable[xlow +xhalf][ylow +yhalf][taulow];
  C001=_TempTable[xlow +xhalf][ylow +yhalf][tauhigh];
  C010=_TempTable[xlow +xhalf][yhigh+yhalf][taulow];
  C011=_TempTable[xlow +xhalf][yhigh+yhalf][tauhigh];
  C100=_TempTable[xhigh+xhalf][ylow +yhalf][taulow];
  C101=_TempTable[xhigh+xhalf][ylow +yhalf][tauhigh];
  C110=_TempTable[xhigh+xhalf][yhigh+yhalf][taulow];
  C111=_TempTable[xhigh+xhalf][yhigh+yhalf][tauhigh];
  //interpolate along x axis
  C00=C000*(1.-dx)+C100*dx;
  C01=C001*(1.-dx)+C101*dx;
  C10=C010*(1.-dx)+C110*dx;
  C11=C011*(1.-dx)+C111*dx;
  //interpolate along y axis
  C0=C00*(1.-dy)+C10*dy;
  C1=C01*(1.-dy)+C11*dy;
  //interpolate along tau axis
  C=C0*(1.-dtau)+C1*dtau;
  // double CC=C000*(1.-dx)*(1.-dy)*(1.-dtau)+C001*(1.-dx)*(1.-dy)*dtau+C010*(1.-dx)*dy*(1.-dtau)+C011*(1.-dx)*dy*dtau+
  //         C100*dx*(1.-dy)*(1.-dtau)+C101*dx*(1.-dy)*dtau+C110*dx*dy*(1.-dtau)+C111*dx*dy*dtau;
  // std::cout << "C=" << C << " CC=" << CC << std::endl;
  return C;
}

//Temperature
double OmegaC::Temperature(const double tau){
  _tau=tau;
  _x=_x0+_tau*std::cos(_theta);
  _y=_y0+_tau*std::sin(_theta);
  double Temp=TriInterpolation(_x,_y,_tau);
  // std::cout << "_x=" << _x << " _y=" << _y << " Temp=" << Temp << std::endl;
  return Temp;
}

//OmegaC
void OmegaC::setOmegaC(const double tau){
  _omegaC=0.0;
  _Temp=0.0;
  if(HydroTaucut(tau)==true) _Temp=Temperature(tau);
  if(HydroTaucut(tau)==true && HydroTempcut(_Temp)==true){
    _omegaC=_qhat0/std::pow(_Temp0,3)*std::pow(_Temp,3)*tau;  //actually _qhat0 = 1
  }
}


//####################################################################################################
//QAG adaptive integration
void OmegaC::setupQAG(double (*fun)(double, void *), double rangemin, double rangemax, void *para){
  _Fun.function=fun;
  _Fun.params=para;
  _rangemin=rangemin; _rangemax=rangemax;
  _workspace=gsl_integration_workspace_alloc(10000);
}

void OmegaC::calculateQAG(double &res, double &err){
  gsl_integration_qag(&_Fun, _rangemin, _rangemax, 0, 1e-3, 10000, 6, _workspace, &res, &err);   //QAG adaptive integration
  // gsl_integration_qags(&_Fun, _rangemin, _rangemax, 0, 1e-3, 10000, _workspace, &res, &err);  //QAGS adaptive integration with singularities
}

void OmegaC::cleanupQAG(){
  gsl_integration_workspace_free(_workspace);
}