#include "ElectroWeak.h"


ElectroWeak::ElectroWeak(){
  // //All the values are taken from PDG 2018
  // _mZ=91.1876;    //in GeV
  // _mW=80.398;     //in GeV
  // _mHiggs=125.1;  //in GeV
  // // _GF=1.1663787e-5;  //in GeV^{-2}

  //gmu scheme shall be used by taking the above three parameters as an input
  _sw2=0.22264585;
  _cw2=1-_sw2;
  _alpha=1.0/132.33843228;

  // //gmu scheme shall be used by taking the above three parameters as an input
  // _cw2=_mW/_mZ; _cw2*=_cw2;
  // _sw2=1.0-_cw2;
  // _alpha=_mW*_mW*sqrt(2.0)*_GF*_sw2/M_PI;
}
