#ifndef ELECTROWEAK_H
#define ELECTROWEAK_H

#include <cmath>


class ElectroWeak{
 
 private:
  
  // double _mZ;      //Z boson mass in GeV
  // double _mW;      //W boson mass in GeV
  // double _mHiggs;  //Higgs boson mass in GeV
  // // double _GF;   //G_F in GeV^{-2}

  double _cw2;    //cos(theta_W)^2
  double _sw2;    //sin(theta_W)^2
  double _alpha;  //electromagnetic alpha


 public:

  ElectroWeak();
  ~ElectroWeak(){}

  // inline double mZ()     const {return _mZ;}
  // inline double mW()     const {return _mW;}
  // inline double mHiggs() const {return _mHiggs;}
  inline double cw2()   const {return _cw2;}
  inline double sw2()   const {return _sw2;}
  inline double alpha() const {return _alpha;}
  
};

#endif  //ELECTROWEAK_H
