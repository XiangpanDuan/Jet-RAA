#include "EnergyLoss.h"


EnergyLoss::EnergyLoss(): DiJetLO(){
  _A=208; _B=208; _Cen=0.; _Ecm=5.02e3;
  _Rsize=0.4;
  _iMul=1.;
  _fMul=1.;
  _nMul=1.;
  _qhat0=1.;
  _alphaSmed=0.2;
  _Qmed=0.5;
  _omegaC=0.0;
  _nDepsilon=1.0;
  _epsilon.resize(_fMul,0.0);
  _Depsilon.resize(_fMul,1.0);
  _xnum=301,_ynum=301,_taunum=201,_thetanum=36;
  _xmin=-15.,_xmax=15.;
  _ymin=-15.,_ymax=15.;
  _omegaCTable.resize(_xnum, std::vector<std::vector<double>>(_ynum, std::vector<double>(_thetanum, 0.0)));
  _elossOmegaC = new EnergyLossOmegaC();
}

EnergyLoss::EnergyLoss(const int A, const int B, const double Cen, const double Ecm, Particle *Parton): DiJetLO(Ecm,Parton){
  _A=A; _B=B; _Cen=Cen; _Ecm=Ecm;
  _Rsize=0.4;
  _iMul=1.;
  _fMul=1.;
  _nMul=1.;
  _qhat0=1.;
  _alphaSmed=0.2;
  _Qmed=0.5;
  _omegaC=0.0;
  _nDepsilon=1.0;
  _epsilon.resize(_fMul,0.0);
  _Depsilon.resize(_fMul,1.0);
  _xnum=301,_ynum=301,_taunum=201,_thetanum=36;
  _xmin=-15.,_xmax=15.;
  _ymin=-15.,_ymax=15.;
  _omegaCTable.resize(_xnum, std::vector<std::vector<double>>(_ynum, std::vector<double>(_thetanum, 0.0)));
  _elossOmegaC = new EnergyLossOmegaC();
}

EnergyLoss::~EnergyLoss(){
  delete _elossOmegaC;
}


void EnergyLoss::setpT(const double pT){
  _pT=pT;
  if(MULTI==0) _pT+=_epsilon[_iMul];        //_iMul=_fMul
  if(MULTI==1) _pT+=_nMul*_epsilon[_iMul];  //_iMul=_fMul, _nMul: mean multiplicity
  if(MULTI==2){
    for(int i=_iMul; i<_fMul; i++){
      _pT+=_epsilon[i];
      // std::cout << "_pT=" << _pT << ",  i=" << i << ",  _fMul=" << _fMul << ",  _epsilon[i]=" << _epsilon[i] << std::endl;
    }
  }
  // std::cout << "_pT=" << _pT << ",  _iMul=" << _iMul << ",  _fMul=" << _fMul << ",  _nMul=" << _nMul << std::endl;
}


//####################################################################################################
//Start time with qhat and jet pT dependences
void EnergyLoss::setTauStart(){
  // _tau0=0.6;  //test
  _tau0=0.2*std::pow(_qhat0/25.,-3./7.)*std::pow(_pT,2./7.);  //unit conversion: 1 ℏ𝑐 = 0.1973269804 GeV fm
  if(_tau0<0.6) _tau0=0.6;
  // std::cout << _pT << " " << _tau0 << std::endl;
}

//Initial data (x,y,tau,temp) produced by 2+1D VISHydro
void EnergyLoss::setTempTable(){
  std::string Astr,Bstr;
  if(_A==1  ) Astr="p";
  if(_A==208) Astr="Pb";
  if(_B==1  ) Bstr="p";
  if(_B==208) Bstr="Pb";
  std::stringstream ss;
  ss << "./InputData/EnergyLoss_" << Astr << Bstr << "_" << _Ecm << "GeV_cen" << _Cen << "_xytautemp.dat";
  std::string InputString=ss.str();
  std::ifstream InputFile;
  InputFile.open(InputString);
  std::cout << "Reading data from " << InputString << std::endl;
  double xval,yval,tauval,tempval;
  for(int ix=0; ix<_xnum; ix++){
    for(int iy=0; iy<_ynum; iy++){
      for(int itau=0; itau<_taunum; itau++){
        InputFile >> xval >> yval >> tauval >> tempval;
        _elossOmegaC->setTemp(ix,iy,itau,tempval);
        if(ix==(int)(std::floor(_xnum/2)) && iy==(int)(std::floor(_ynum/2)) && itau==6) _elossOmegaC->setTemp0(tempval);
        // if(xval==0 && yval==0 && tauval==0.6) std::cout << ix << "  " << iy << "  " << itau << std::fixed << std::setprecision(16) << " Temp0=" << tempval << std::endl;  //zero point test
      }
    }
  }
  InputFile.close();
}

//QAG integration to get omegaC
void EnergyLoss::setOmegaCTau(const double x0, const double y0, const double theta){
  _omegaC=0.0;
  _elossOmegaC->setxytheta(x0,y0,theta);
  _omegaC=_elossOmegaC->CalculateOmegaC(_tau0)*_qhat0;  //attention qhat0
  // std::cout << _pT << " " << _tau0 << " " << _omegaC << std::endl;
}


//####################################################################################################
//Initial data (x,y,theta,omegaC) produced by 2+1D VISHydro
void EnergyLoss::setOmegaCTable(){
  std::string Astr,Bstr;
  if(_A==1  ) Astr="p";
  if(_A==208) Astr="Pb";
  if(_B==1  ) Bstr="p";
  if(_B==208) Bstr="Pb";
  std::stringstream ss;
  ss << "./InputData/EnergyLoss_" << Astr << Bstr << "_" << _Ecm << "GeV_cen" << _Cen << "_xythetaomegaC.dat";
  std::string InputString=ss.str();
  std::ifstream InputFile;
  InputFile.open(InputString);
  std::cout << "Reading data from " << InputString << std::endl;
  double xval,yval,thetaval,omegaCval;
  for(int ix=0; ix<_xnum; ix++){
    for(int iy=0; iy<_ynum; iy++){
      for(int itheta=0; itheta<_thetanum; itheta++){
        InputFile >> xval >> yval >> thetaval >> omegaCval;
        _omegaCTable[ix][iy][itheta]=omegaCval*_qhat0;  //attention _qhat0
      }
    }
  }
  InputFile.close();
  // std::cout << _omegaCTable[_xnum/2][_ynum/2][0] << std::endl; //x=0, y=0, theth=0
}

//Trilinear interpolation to get omegaC
bool EnergyLoss::setOmegaC(const double x, const double y, const double theta){
  _omegaC=0.0;
  int    xlow,xhigh,ylow,yhigh,thetalow,thetahigh;
  double xbin,ybin,thetabin;
  double C,C0,C1,C00,C01,C10,C11,C000,C001,C010,C011,C100,C101,C110,C111;
  xbin=(_xmax-_xmin)/(_xnum-1);
  ybin=(_ymax-_ymin)/(_ynum-1);
  thetabin=2.*M_PI/_thetanum;
  xlow=(int)(std::floor(x/xbin)); xhigh=xlow+1;
  ylow=(int)(std::floor(y/ybin)); yhigh=ylow+1;
  thetalow=(int)(std::floor(theta/thetabin)); thetahigh=thetalow+1;
  if(thetalow==(_thetanum-1)) {thetahigh=0;}            //attention location angle = 0, 2PI
  if(thetalow==_thetanum) {thetalow=0;  thetahigh=1;}   //attention location angle = 0, 2PI
  if(thetalow>_thetanum)  return false;
  if(xlow<_xmin/xbin || xlow>=_xmax/xbin) return false;
  if(ylow<_ymin/ybin || ylow>=_ymax/ybin) return false;
  double dx=(x/xbin-xlow)/1.;  //double dx=(x/xbin-xlow)/(xhigh-xlow);
  double dy=(y/ybin-ylow)/1.;  //double dy=(y/ybin-ylow)/(yhigh-ylow);
  double dtheta=(theta/thetabin-thetalow)/1.;  //double dtheta=(theta/thetabin-thetalow)/(thetahigh-thetalow);
  // if(thetahigh==0) dtheta=(theta/thetabin-thetalow)/1.;
  //Attention x and y range
  int xhalf=(int)(std::floor(_xnum/2));
  int yhalf=(int)(std::floor(_ynum/2));
  C000=_omegaCTable[xlow +xhalf][ylow +yhalf][thetalow];  //get real lacation omegaCerature
  C001=_omegaCTable[xlow +xhalf][ylow +yhalf][thetahigh];
  C010=_omegaCTable[xlow +xhalf][yhigh+yhalf][thetalow];
  C011=_omegaCTable[xlow +xhalf][yhigh+yhalf][thetahigh];
  C100=_omegaCTable[xhigh+xhalf][ylow +yhalf][thetalow];
  C101=_omegaCTable[xhigh+xhalf][ylow +yhalf][thetahigh];
  C110=_omegaCTable[xhigh+xhalf][yhigh+yhalf][thetalow];
  C111=_omegaCTable[xhigh+xhalf][yhigh+yhalf][thetahigh];
  //Interpolate along x axis
  C00=C000*(1.-dx)+C100*dx;
  C01=C001*(1.-dx)+C101*dx;
  C10=C010*(1.-dx)+C110*dx;
  C11=C011*(1.-dx)+C111*dx;
  //Interpolate along y axis
  C0=C00*(1.-dy)+C10*dy;
  C1=C01*(1.-dy)+C11*dy;
  //Interpolate along theta axis
  C=C0*(1.-dtheta)+C1*dtheta;
  // C=C000*(1.-dx)*(1.-dy)*(1.-dtheta)+C001*(1.-dx)*(1.-dy)*dtheta+C010*(1.-dx)*dy*(1.-dtheta)+C011*(1.-dx)*dy*dtheta+ \
  //   C100*dx*(1.-dy)*(1.-dtheta)+C101*dx*(1.-dy)*dtheta+C110*dx*dy*(1.-dtheta)+C111*dx*dy*dtheta;
  _omegaC=C;
  // std::cout << _omegaC << std::endl;
  // std::cout << x << " " << y << " " << theta << " " << _omegaC << std::endl;
  // std::cout << thetalow << " " << thetahigh << std::endl;
  // std::cout << theta/thetabin-thetalow << " " << thetahigh-thetalow << " " << dtheta << std::endl;
  return true;
}


//####################################################################################################
//Jet energy loss when traversing QGP in BDMPS formalism
void EnergyLoss::setEpsilon(const double *epsilon){
  _epsilon.resize(_fMul,0.0);
  for(int i=_iMul; i<_fMul; i++){
    _epsilon[i]=epsilon[i];
    // std::cout << _iMul << " " << i << " " << epsilon[i] << std::endl;
  }
}

//One parton (quark or gluon) jet energy loss
void EnergyLoss::setDepsilon(const std::string &type){
  double CR,QGratio;
  if(type=="quark") {CR=_qcd->CF(); QGratio=1.0;}                    //quark
  if(type=="gluon") {CR=_qcd->Nc(); QGratio=_qcd->Nc()/_qcd->CF();}  //gluon
  double alphaS=_alphaSmed;                                          //medium scale, alphaS=0.2 (default)
  double alpha=2.*CR*alphaS/M_PI;
  double omegaC=QGratio*_omegaC;
  _Depsilon.resize(_fMul,1.0);
  _nDepsilon=1.0;
  //First energy loss
  _Depsilon[_iMul]=1./_epsilon[_iMul]*std::sqrt(alpha*alpha*omegaC/(2.*_epsilon[_iMul]))*std::exp(-M_PI*alpha*alpha*omegaC/(2.*_epsilon[_iMul]));
  _nDepsilon*=_Depsilon[_iMul];
  // std::cout << type << " " << alphaS << " " << _Depsilon[_iMul] << " " << _nDepsilon << std::endl;
}

//Multi-gluons jet energy loss
void EnergyLoss::setDepsilonGluon(){
  double CR=_qcd->Nc();                  //gluon
  double QGratio=_qcd->Nc()/_qcd->CF();  //gluon
  double alphaS=_alphaSmed;              //medium scale, alphaS=0.2 (default)
  double alpha=2.*CR*alphaS/M_PI;
  double omegaC=QGratio*_omegaC;
  //Start loop from second energy loss
  for(int i=(_iMul+1); i<_fMul; i++){
    _Depsilon[i]=1./_epsilon[i]*std::sqrt(alpha*alpha*omegaC/(2.*_epsilon[i]))*std::exp(-M_PI*alpha*alpha*omegaC/(2.*_epsilon[i]));
    _nDepsilon*=_Depsilon[i];
    // if(i==7 || i==8) std::cout << "gluon " << _Depsilon[i] << " " << _nDepsilon << std::endl;
  }
}

//Multi-gluons jet mean energy loss
void EnergyLoss::setDepsilonGluonMean(){
  double CR=_qcd->Nc();                  //gluon
  double QGratio=_qcd->Nc()/_qcd->CF();  //gluon
  double alphaS=_alphaSmed;              //medium scale, alphaS=0.2 (default)
  double alpha=2.*CR*alphaS/M_PI;
  double omegaC=QGratio*_omegaC;
  _Depsilon.resize(_fMul,1.0);
  _nDepsilon=1.0;
  //Mean energy loss
  _Depsilon[_iMul]=1./_epsilon[_iMul]*std::sqrt(alpha*alpha*omegaC/(2.*_epsilon[_iMul]))*std::exp(-M_PI*alpha*alpha*omegaC/(2.*_epsilon[_iMul]));
  _nDepsilon*=_Depsilon[_iMul];
}


//####################################################################################################
//Multiplicity in jet energy loss
// //Formula 2.17(P6) in "Quenching of hadron spectra in media" [JHEP 09 (2001) 033]
// void EnergyLoss::setMultiplicity(const std::string &type, const double pT){
//   double QGratio;
//   if(type=="quark") QGratio=1.0;  //quark
//   if(type=="gluon") QGratio=1.6;  //gluon
//   double TR=_qcd->TF();
//   double b0=(11.*_qcd->Nc()-4.*_qcd->Nf()*TR)/(12.*M_PI);
//   double a0=1./4.+5.*_qcd->Nf()/(54.*M_PI*b0);
//   // double varLL=0.016965;  //LLA  from pythia8 simulation
//   double varMLL=0.0403164; //MLLA from pythia8 simulation
//   // //Multiplicity
//   // //with Qmed
//   // double Qmed=0.5; //pT,track cut
//   // double nVal;
//   // nVal=QGratio*std::exp(std::pow(2.*_qcd->Nc()/(M_PI*b0),0.5)*(std::sqrt(std::log(((pT*_Rsize)*(pT*_Rsize))/(_lambdaQCD*_lambdaQCD)))-std::sqrt(std::log((Qmed*Qmed)/(_lambdaQCD*_lambdaQCD)))));
//   // // _nMul=varLL*nVal; //LLA varLL with Qmed
//   // _nMul=varMLL*nVal*std::pow(std::log((Qmed*Qmed)/(_lambdaQCD*_lambdaQCD))/std::log(((pT*_Rsize)*(pT*_Rsize))/(_lambdaQCD*_lambdaQCD)),a0); //MLLA varMLL with Qmed
//   //without Qmed
//   double alphaS=1./(b0*std::log((pT*_Rsize)*(pT*_Rsize)/(_lambdaQCD*_lambdaQCD)));
//   // _nMul=QGratio*varLL*std::exp((1./b0)*std::pow(2.*_qcd->Nc()/(M_PI*alphaS),0.5)); //LLA c1
//   _nMul=QGratio*varMLL*std::exp((1./b0)*std::pow(2.*_qcd->Nc()/(M_PI*alphaS),0.5)+a0*std::log(alphaS)); //MLLA c2

//   //Get initial pT
//   double pTin=pT;
//   double pTout=pT+_nMul*_epsilon[_iMul];
//   while((pTout-pTin)>0.1){
//     pTin=pTout;
//     // //with Qmed
//     // nVal=QGratio*std::exp(std::pow(2.*_qcd->Nc()/(M_PI*b0),0.5)*(std::sqrt(std::log((pTin*_Rsize*pTin*_Rsize)/(_lambdaQCD*_lambdaQCD)))-std::sqrt(std::log((Qmed*Qmed)/(_lambdaQCD*_lambdaQCD)))));
//     // // _nMul=varLL*nVal; //LLA varLL with Qmed
//     // _nMul=varMLL*nVal*std::pow(std::log((Qmed*Qmed)/(_lambdaQCD*_lambdaQCD))/std::log((pTin*_Rsize*pTin*_Rsize)/(_lambdaQCD*_lambdaQCD)),a0); //MLLA varMLL with Qmed
//     //without Qmed
//     alphaS=1./(b0*std::log((pTin*_Rsize)*(pTin*_Rsize)/(_lambdaQCD*_lambdaQCD)));
//     // _nMul=QGratio*varLL*std::exp((1./b0)*std::pow(2.*_qcd->Nc()/(M_PI*alphaS),0.5)); //LLA c1
//     _nMul=QGratio*varMLL*std::exp((1./b0)*std::pow(2.*_qcd->Nc()/(M_PI*alphaS),0.5)+a0*std::log(alphaS)); //MLLA c2
//     pTout=pT+_nMul*_epsilon[_iMul];
//     // std::cout << pT << " " << pTin << " " << pTout << " " << _nMul <<std::endl;
//   }
// }

//Formula 5.51(P124) for LL and 7.34(P174) for LL and MLL in [Basics of Perturbative QCD]
void EnergyLoss::setMultiplicity(const std::string &type, const double pT){
  double varLL=1.0;  //LL
  // double varMLL=0.24;  //MLL
  double QGratio=1.0;
  if(type=="quark") QGratio=1.0;  //quark
  if(type=="gluon") QGratio=1.6;  //gluon
  // double a=11./3.*_qcd->Nc()+2.*_qcd->Nf()/(3.*_qcd->Nc()*_qcd->Nc());  //for MLL
  double b=11./3.*_qcd->Nc()-2./3.*_qcd->Nf();
  double A=std::sqrt(16.*_qcd->Nc()/b);
  // double B=a/b;     //for MLL
  double Qmed=_Qmed;   //medium scale in Multiplicity, Qmed=0.5 (default), Qmed=_lambdaQCD+1e-10
  double Q=pT*_Rsize;  //R dependence
  double ypT=std::log(Q/Qmed);
  double lambdapT=std::log(Qmed/_lambdaQCD);
  double YpT=ypT+lambdapT;
  double xx1=A*std::sqrt(YpT);
  double xx2=A*std::sqrt(lambdapT);
  _nMul=QGratio*varLL*xx1*(gsl_sf_bessel_I1(xx1)*gsl_sf_bessel_K0(xx2)+gsl_sf_bessel_K1(xx1)*gsl_sf_bessel_I0(xx2));  //LL
  // _nMul=QGratio*varMLL*xx1*std::pow(xx2/xx1,B)*(gsl_sf_bessel_Inu(B+1,xx1)*gsl_sf_bessel_Knu(B,xx2)+gsl_sf_bessel_Knu(B+1,xx1)*gsl_sf_bessel_Inu(B,xx2));  //MLL
  // std::cout << pT << " " << _nMul << std::endl;

  //Get initial pT
  double pTin=pT;
  double pTout=pT+_nMul*_epsilon[_iMul];
  while((pTout-pTin)>0.1){
    pTin=pTout;
    Q=pTin*_Rsize;
    ypT=std::log(Q/Qmed);
    lambdapT=std::log(Qmed/_lambdaQCD);
    YpT=ypT+lambdapT;
    xx1=A*std::sqrt(YpT);
    xx2=A*std::sqrt(lambdapT);
    _nMul=QGratio*varLL*xx1*(gsl_sf_bessel_I1(xx1)*gsl_sf_bessel_K0(xx2)+gsl_sf_bessel_K1(xx1)*gsl_sf_bessel_I0(xx2));  //LL
    // _nMul=QGratio*varMLL*xx1*std::pow(xx2/xx1,B)*(gsl_sf_bessel_Inu(B+1,xx1)*gsl_sf_bessel_Knu(B,xx2)+gsl_sf_bessel_Knu(B+1,xx1)*gsl_sf_bessel_Inu(B,xx2));  //MLL
    pTout=pT+_nMul*_epsilon[_iMul];
  }
}
