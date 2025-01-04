#include "Glauber.h"


Glauber::Glauber(){
  _A=208;
  _B=208;
  _Cen=0.0;
  _TAB0=304.152;
  setNuclearParameters(_A,_B);
  setGaussianQuadrature();
}

Glauber::Glauber(const int A, const int B, const double Cen){
  _A=A;
  _B=A;
  _Cen=Cen;
  setThicknessTAB0(_A,_B);
  setNuclearParameters(_A,_B);
  setGaussianQuadrature();
}

void Glauber::setThicknessTAB0(const int A, const int B){
  if(A==208 && B==208) _TAB0=304.152;  //Pb+Pb
}

//In 115 pages of [Ultrarelativistic Heavy-Ion Collisions] and 3-11 chapter of [A short course on Relativistic Heavy Ion Collision]
void Glauber::setNuclearParameters(const int A, const int B){
  //_rho0: the nucleon density in the center of the nucleus, _R: the nuclear radius, _a: the skin depth, _omega: deviations from a spherical shape
  if(A==208) {_rho0_A=0.16;   _R_A=6.624; _a_A=0.549; _omega_A=0.0;}     //Pb _A=207, _rho0_A=0.159241
  if(A==197) {_rho0_A=0.1693; _R_A=6.38;  _a_A=0.535; _omega_A=0.0;}     //Au
  if(A==16 ) {_rho0_A=0.1654; _R_A=2.608; _a_A=0.513; _omega_A=-0.051;}  //O
  //
  if(B==208) {_rho0_B=0.16;   _R_B=6.624; _a_B=0.549; _omega_B=0.0;}     //Pb _B=207, _rho0_B=0.159241
  if(B==197) {_rho0_B=0.1693; _R_B=6.38;  _a_B=0.535; _omega_B=0.0;}     //Au
  if(B==16 ) {_rho0_B=0.1654; _R_B=2.608; _a_B=0.513; _omega_B=-0.051;}  //O
}

//Gaussian-Legendre Quadrature
void Glauber::setGaussianQuadrature(){
  _numGauss=100;  //data size
  _zi.resize(_numGauss);
  _wi.resize(_numGauss);
  gsl_integration_glfixed_table *table=gsl_integration_glfixed_table_alloc(_numGauss);
  double zmin=-20.0, zmax=20.0;  //range of z
  for(int i=0; i<_numGauss; i++){
    gsl_integration_glfixed_point(zmin, zmax, i, &_zi[i], &_wi[i], table);   //Gaussian-Legendre roots and weights
    // std::cout << i+1 << " " << _wi[i] << " " << _zi[i] << std::endl;
  }
  gsl_integration_glfixed_table_free(table);
}


//####################################################################################################
//Woods-Saxon density distributions
void Glauber::setDensityRhoA(const double x, const double y, const double z){
  double r=std::sqrt(x*x+y*y+z*z);
  _rho_A=_rho0_A*((1.+_omega_A*std::pow(r/_R_A,2))/(1.+std::exp((r-_R_A)/_a_A)));
}

void Glauber::setDensityRhoB(const double x, const double y, const double z){
  double r=std::sqrt(x*x+y*y+z*z);
  _rho_B=_rho0_B*((1.+_omega_B*std::pow(r/_R_B,2))/(1.+std::exp((r-_R_B)/_a_B)));
}


//####################################################################################################
//TA
void Glauber::setThicknessTA(const double sx, const double sy){
  // double zmin=-20.,zmax=20.;  //range of z
  _TA=0.0;
  for(int i=0; i<_numGauss; i++){
    setDensityRhoA(sx,sy,_zi[i]);  //updata current _rho_A
    _TA+=_rho_A*_wi[i];            //Legendre-Gaussian Quadrature
    // std::cout << i+1 << "  " << _wi[i] << "  " << _zi[i] << "  " << _rho_A << "  " << _TA << std::endl;
  }
}

//TB
void Glauber::setThicknessTB(const double sx, const double sy){
  // double zmin=-20.,zmax=20.;  //range of z
  _TB=0.0;
  for(int i=0; i<_numGauss; i++){
    setDensityRhoB(sx,sy,_zi[i]);  //updata current _rho_B
    _TB+=_rho_B*_wi[i];            //Legendre-Gaussian Quadrature
  }
}


//####################################################################################################
//TAB:Thickness function
void Glauber::setThicknessTAB(const double sx, const double sy){
  //rotate impact parameter b to x-axis
  // double bx=0.0, by=0.0;
  // if(_A==208 && _B==208){
  //   if(_Cen>=0.00 && _Cen<0.05)  {bx=0.;  by=0.;}   //0-5%
  //   if(_Cen>=0.05 && _Cen<0.10)  {bx=4.;  by=4.;}   //5-10%
  //   if(_Cen>=0.10 && _Cen<0.20)  {bx=6.;  by=6.;}   //10-20%
  //   if(_Cen>=0.20 && _Cen<0.30)  {bx=8.;  by=8.;}   //20-30%
  //   if(_Cen>=0.30 && _Cen<0.50)  {bx=10.; by=10.;}  //30-50%
  //   if(_Cen>=0.50 && _Cen<0.100) {bx=12.; by=12.;}  //50-100%
  // }
  double bx=0.0, by=0.0;
  setThicknessTA(sx+bx/2.,sy+by/2.);
  setThicknessTB(sx-bx/2.,sy-by/2.);
  _TAB=_TA*_TB;
}

void Glauber::setThicknessTAB(const double sx, const double sy, const double bx, const double by){
  setThicknessTA(sx+bx/2.,sy+by/2.);
  setThicknessTB(sx-bx/2.,sy-by/2.);
  _TAB=_TA*_TB;
}
