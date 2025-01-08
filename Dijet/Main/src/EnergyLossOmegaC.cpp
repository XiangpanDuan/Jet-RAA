#include "EnergyLossOmegaC.h"


EnergyLossOmegaC::EnergyLossOmegaC(){
  _xnum=301;_ynum=301;_taunum=201;_thetanum=36;
  _xmin=-15.;_xmax=15.;
  _ymin=-15.;_ymax=15.;
  _xbin=(_xmax-_xmin)/(_xnum-1);
  _ybin=(_ymax-_ymin)/(_ynum-1);
  _taubin=0.1;
  _Temp0=0.0;
  _Temp=0.0;
  _omegaC=0.0;
  _TempTable.resize(_xnum, std::vector<std::vector<double>>(_ynum, std::vector<double>(_taunum, 0.0)));
}


//####################################################################################################
//Hydro cut
//Initial time cut
bool EnergyLossOmegaC::HydroTaucut(const double tau){
  double tau0=0.6;  //initial time in hydro medium
  return tau>=tau0;
}

//Freeze-out temperature cut
bool EnergyLossOmegaC::HydroTempcut(const double temp){
  double tempcut=0.154;  //from lattice QCD [PRD 90 (2014) 094503]
  return temp>=tempcut;
}


//####################################################################################################
//Trilinear interpolation to get temperature
double EnergyLossOmegaC::TriInterpolation(const double x, const double y, const double tau){
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
  //Attention x and y range
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
  //Interpolate along x axis
  C00=C000*(1.-dx)+C100*dx;
  C01=C001*(1.-dx)+C101*dx;
  C10=C010*(1.-dx)+C110*dx;
  C11=C011*(1.-dx)+C111*dx;
  //Interpolate along y axis
  C0=C00*(1.-dy)+C10*dy;
  C1=C01*(1.-dy)+C11*dy;
  //Interpolate along tau axis
  C=C0*(1.-dtau)+C1*dtau;
  // double CC=C000*(1.-dx)*(1.-dy)*(1.-dtau)+C001*(1.-dx)*(1.-dy)*dtau+C010*(1.-dx)*dy*(1.-dtau)+C011*(1.-dx)*dy*dtau+
  //         C100*dx*(1.-dy)*(1.-dtau)+C101*dx*(1.-dy)*dtau+C110*dx*dy*(1.-dtau)+C111*dx*dy*dtau;
  // std::cout << "Test2 " << C000 << " " << C001 << " " << C010 << " " << C011 << " " << C100 << " " << C101 << " " << C110 << " " << C111 << std::endl;
  // std::cout << "Test2 " << C00 << " " << C01 << " " << C10 << " " << C11 << std::endl;
  // std::cout << "Test2 " << C0 << " " << C1 << std::endl;
  // std::cout << "Test2 " << C << std::endl;
  return C;
}

//Temperature
double EnergyLossOmegaC::Temperature(const double tau){
  _tau=tau;
  _x=_x0+_tau*std::cos(_theta);
  _y=_y0+_tau*std::sin(_theta);
  double Temp=TriInterpolation(_x,_y,_tau);
  // std::cout << _x0 << " " << _y0 << " " << _theta << " " << _x << " " << _y << " " << _tau << " " << Temp << std::endl;
  return Temp;
}


//####################################################################################################
//OmegaC
void EnergyLossOmegaC::setOmegaC(const double tau){
  _omegaC=0.0;
  _Temp=0.0;
  if(HydroTaucut(tau)==true) _Temp=Temperature(tau);
  if(HydroTaucut(tau)==true && HydroTempcut(_Temp)==true){
    double qhat0=1.0;
    _omegaC=qhat0/std::pow(_Temp0,3)*std::pow(_Temp,3)*tau;
  }
}
