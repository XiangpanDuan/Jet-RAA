#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <cmath>
#include <gsl/gsl_math.h>
#include <gsl/gsl_sf_bessel.h>


//Mean multiplicity in jet energy losss
//Formula 5.51(P124) for LLA and 7.34(P174) for MLLA in [Basics of Perturbative QCD]
int main()
{
  double Qmed=35.;  //medium scale
  int    mode=0;    //0:LL; 1:MLL
  std::stringstream ss;
  if(mode==0) ss << "../Output/MultiLL_Qmed"  << Qmed << ".dat";
  if(mode==1) ss << "../Output/MultiMLL_Qmed" << Qmed << ".dat";
  std::string OutputString=ss.str();
  std::ofstream OutputFile;
  OutputFile.open(OutputString);
  OutputFile << "# pT Multi_Quark Multi_Gluon" << std::endl;

  int    nf=3.;
  int    Nc=3;
  double LambdaQCD=0.245748;  //keep consistent with the _lambdaQCD value in the Main/src/QCD.cpp
  double Rsize=0.4;
  double a=11./3.*Nc+2.*nf/(3.*Nc*Nc);
  double b=11./3.*Nc-2./3.*nf;
  double A=std::sqrt(16.*Nc/b);
  double B=a/b;
  double ypT,lambdapT,YpT;
  double x1,x2;
  double Multi_Q,Multi_G;
  double varLL=1.0;       //LL
  double varMLL=0.24;     //MLL
  // double varMLL=1.0;   //MLL
  double QGratio_Q=1.00;  //quark
  double QGratio_G=1.62;  //gluon
  double pT,Q;
  double pTmin=50.,pTbin=10.;
  int    pTnum=100;
  for(int i=0; i<pTnum; i++){
    pT=pTmin+i*pTbin;
    Q=pT*Rsize;  //R dependence
    ypT=std::log(Q/Qmed);
    lambdapT=std::log(Qmed/LambdaQCD);
    YpT=ypT+lambdapT;
    x1=A*std::sqrt(YpT);
    x2=A*std::sqrt(lambdapT);
    if(mode==0){
      Multi_Q=QGratio_Q*varLL*x1*(gsl_sf_bessel_I1(x1)*gsl_sf_bessel_K0(x2)+gsl_sf_bessel_K1(x1)*gsl_sf_bessel_I0(x2));  //LL
      Multi_G=QGratio_G*varLL*x1*(gsl_sf_bessel_I1(x1)*gsl_sf_bessel_K0(x2)+gsl_sf_bessel_K1(x1)*gsl_sf_bessel_I0(x2));  //LL
    }
    if(mode==1){
      Multi_Q=QGratio_Q*varMLL*x1*std::pow(x2/x1,B)*(gsl_sf_bessel_Inu(B+1,x1)*gsl_sf_bessel_Knu(B,x2)+gsl_sf_bessel_Knu(B+1,x1)*gsl_sf_bessel_Inu(B,x2));  //MLL
      Multi_G=QGratio_G*varMLL*x1*std::pow(x2/x1,B)*(gsl_sf_bessel_Inu(B+1,x1)*gsl_sf_bessel_Knu(B,x2)+gsl_sf_bessel_Knu(B+1,x1)*gsl_sf_bessel_Inu(B,x2));  //MLL
    }
    OutputFile << pT << " " << Multi_Q << " " << Multi_G << std::endl;
  }
  OutputFile.close();


  // //
  // double TR=0.5;
  // double b0=(11.*Nc-4.*nf*TR)/(12.*M_PI);
  // double a0=1./4.+5.*nf/(54.*M_PI*b0);
  // double varLL=0.016965;    //LLA  from pythia8 simulation
  // double varMLL=0.0403164;  //MLLA from pythia8 simulation
  // for(int i=0; i<pTnum; i++){
  //   pT=pTmin+i*pTbin;
  //   Q=pT*Rsize;
  //   double alphaS=1./(b0*std::log(Q*Q/(LambdaQCD*LambdaQCD)));
  //   if(mode==0){
  //     Multi_Q=QGratio_Q*varLL*std::exp(std::pow(2.*Nc/(M_PI*b0)*std::log(Q*Q/(LambdaQCD*LambdaQCD)),0.5));  //LLA
  //     Multi_G=QGratio_G*varLL*std::exp(std::pow(2.*Nc/(M_PI*b0)*std::log(Q*Q/(LambdaQCD*LambdaQCD)),0.5));  //LLA
  //   }
  //   if(mode==1){
  //     Multi_Q=QGratio_Q*varMLL*std::exp(std::pow(2.*Nc/(M_PI*b0)*std::log(Q*Q/(LambdaQCD*LambdaQCD)),0.5)+a0*std::log(alphaS));  //MLLA
  //     Multi_G=QGratio_G*varMLL*std::exp(std::pow(2.*Nc/(M_PI*b0)*std::log(Q*Q/(LambdaQCD*LambdaQCD)),0.5)+a0*std::log(alphaS));  //MLLA
  //   }
  // }


  return 0;
}