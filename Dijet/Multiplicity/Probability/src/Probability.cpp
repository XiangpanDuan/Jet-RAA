#include <iostream>
#include <fstream>
#include <sstream>
#include <iomanip>
#include <vector>
#include <cmath>
#include <gsl/gsl_integration.h>


const  int pTnum=100;
double pTmin=40.;
double pTmax=1040;
double pT,pTbin,pTbinmin,pTbinmax;
const  int Pnum=10;
const  int num=100;
std::vector<std::vector<double>> P (pTnum, std::vector<double>(Pnum, 0.0));
std::vector<std::vector<double>> Su(pTnum, std::vector<double>(Pnum, 0.0));
std::vector<std::vector<std::vector<double>>> PTable (pTnum, std::vector<std::vector<double>>(Pnum, std::vector<double>(num, 0.0)));
std::vector<std::vector<std::vector<double>>> SuTable(pTnum, std::vector<std::vector<double>>(Pnum, std::vector<double>(num, 0.0)));
double Pval=0.0;
int    ival,jval;
double Rsize=0.4;
double Qmed=19.9;  //25->10
double LambdaQCD=0.245748;  //keep consistent with the _lambdaQCD value in the Main/src/QCD.cpp
double lambda=std::log(Qmed/LambdaQCD);
unsigned int Nc=3;
unsigned int nf=3;
double b=11./3.*Nc-2./3.*nf;
double ymin,ymax,ybin;
double PSum,nMul;

//------------------------------------------------------------
int Factorial(const int n){
    if(n<=1) return 1;        //0!=1, 1!=1
    return n*Factorial(n-1);  //n!
}

//------------------------------------------------------------
double LinearInterpolation(const double yp){
  int ylow=(int)(std::floor(yp/ybin));
  int yhigh=ylow+1;
  if(ylow<0 || yhigh>=num) return 0.0;
  return PTable[ival][jval][ylow]+(PTable[ival][jval][yhigh]-PTable[ival][jval][ylow])*(yp-ylow*ybin)/ybin;
}

//------------------------------------------------------------
double Function(double yp, void *params){
  // (void)(params);
  double gamma0square=4.*Nc/(b*(yp+lambda));
  double Pval=LinearInterpolation(yp);
  int    nFactorial=Factorial(jval+1);
  return (ymax-yp)*gamma0square*Pval*nFactorial;
}

//------------------------------------------------------------
double FunctionBase(double yp, void *params){
  // (void)(params);
  double gamma0square=4.*Nc/(b*(yp+lambda));
  return (ymax-yp)*gamma0square;
}


//------------------------------------------------------------
void CalProbability(){

  std::stringstream ss;
  ss << "../Output/Probability_Qmed" << Qmed << ".dat";
  std::string OutputString=ss.str();
  std::ofstream OutputFile;
  OutputFile.open(OutputString);
  OutputFile << "# pT P1 P2 P3 P4 P5 P6 P7 P8 P9 P10 <Mul>" << std::endl;

  pTbin=(pTmax-pTmin)/pTnum;
  for(int i=0; i<pTnum; i++){
    ival=i;
    pT=pTmin+pTbin/2.+pTbin*i;
    //pTbinmin=pT-pTbin/2.;
    pTbinmax=pT+pTbin/2.;
    pTbinmax=pTbinmax*Rsize;  //attention Rsize
    OutputFile << pT+pTbin/2. << " ";

    PSum=0.0;
    nMul=0.0;
    for(int j=0; j<Pnum; j++){
      jval=j;
      ybin=std::log(pTbinmax/Qmed)/num;

      double res,err;

      //Index
      if(j==0){
        for(int k=0; k<num; k++){
          ymin=0.0;
          ymax=ybin*(k+1);  //update current y value
          //QAG adaptive integration
          gsl_function FunBase;
          FunBase.function = &FunctionBase;
          FunBase.params = nullptr;
          gsl_integration_workspace *workspace = gsl_integration_workspace_alloc(10000);
          res=0.0; err=0.0;
          gsl_integration_qag(&FunBase, ymin, ymax, 0, 1e-3, 10000, 6, workspace, &res, &err);
          gsl_integration_workspace_free(workspace);
          PTable[i][j][k]=exp(-res);
        }
      }
      //Index
      for(int k=0; k<num; k++){
        if(j>0){
          for(int jj=0; jj<j; jj++){
            PTable[i][j][k]+=(double)(jj+1)/j*PTable[i][j-1-jj][k]*SuTable[i][jj][k]/Factorial(jj+1);
          }
        }
        // if(j>0){
        //   if(j==1)  PTable[i][j][k]=1./Factorial(j)*PTable[i][0][k]*(pow(SuTable[i][0][k],j));
        //   if(j==2)  PTable[i][j][k]=1./Factorial(j)*PTable[i][0][k]*(pow(SuTable[i][0][k],j)+SuTable[i][1][k]);
        //   if(j==3)  PTable[i][j][k]=1./Factorial(j)*PTable[i][0][k]*(pow(SuTable[i][0][k],j)\
        //                                                           +3*pow(SuTable[i][0][k],1)*SuTable[i][1][k]+SuTable[i][2][k]);
        //   if(j==4)  PTable[i][j][k]=1./Factorial(j)*PTable[i][0][k]*(pow(SuTable[i][0][k],j)\
        //                                                           +6*pow(SuTable[i][0][k],2)*SuTable[i][1][k]+3*pow(SuTable[i][1][k],2)\
        //                                                           +4*pow(SuTable[i][0][k],1)*SuTable[i][2][k]+SuTable[i][3][k]);
        //   if(j==5)  PTable[i][j][k]=1./Factorial(j)*PTable[i][0][k]*(pow(SuTable[i][0][k],j)\
        //                                                          +10*pow(SuTable[i][0][k],3)*SuTable[i][1][k]\
        //                                                          +10*pow(SuTable[i][0][k],2)*SuTable[i][2][k]+10*SuTable[i][1][k]*SuTable[i][2][k]\
        //                                                          +5 *pow(SuTable[i][0][k],1)*(3*pow(SuTable[i][1][k],2)+SuTable[i][3][k])+SuTable[i][4][k]);
        //   if(j==6)  PTable[i][j][k]=1./Factorial(j)*PTable[i][0][k]*(pow(SuTable[i][0][k],j)\
        //                                                          +15*pow(SuTable[i][0][k],4)*SuTable[i][1][k]+15*pow(SuTable[i][1][k],3)\
        //                                                          +20*pow(SuTable[i][0][k],3)*SuTable[i][2][k]+10*pow(SuTable[i][2][k],2)+15*SuTable[i][1][k]*SuTable[i][3][k]\
        //                                                          +15*pow(SuTable[i][0][k],2)*(3*pow(SuTable[i][1][k],2)+SuTable[i][3][k])\
        //                                                          +6 *pow(SuTable[i][0][k],1)*(10*SuTable[i][1][k]*SuTable[i][2][k]+SuTable[i][4][k])+SuTable[i][5][k]);
        //   if(j==7)  PTable[i][j][k]=1./Factorial(j)*PTable[i][0][k]*(pow(SuTable[i][0][k],j)\
        //                                                          +21*pow(SuTable[i][0][k],5)*SuTable[i][1][k]\
        //                                                          +35*pow(SuTable[i][0][k],4)*SuTable[i][2][k]+105*pow(SuTable[i][1][k],2)*SuTable[i][2][k]+35*SuTable[i][2][k]*SuTable[i][3][k]\
        //                                                          +35*pow(SuTable[i][0][k],3)*(3*pow(SuTable[i][1][k],2)+SuTable[i][3][k])+21*SuTable[i][1][k]*SuTable[i][4][k]\
        //                                                          +21*pow(SuTable[i][0][k],2)*(10*SuTable[i][1][k]*SuTable[i][2][k]+SuTable[i][4][k])\
        //                                                          +7 *pow(SuTable[i][0][k],1)*(15*pow(SuTable[i][1][k],3)+10*pow(SuTable[i][2][k],2)+15*SuTable[i][1][k]*SuTable[i][3][k]+SuTable[i][5][k])+SuTable[i][6][k]);
        //   if(j==8)  PTable[i][j][k]=1./Factorial(j)*PTable[i][0][k]*(pow(SuTable[i][0][k],j)\
        //                                                          +28*pow(SuTable[i][0][k],6)*SuTable[i][1][k]+105*pow(SuTable[i][1][k],4)\
        //                                                          +56*pow(SuTable[i][0][k],5)*SuTable[i][2][k]+210*pow(SuTable[i][1][k],2)*SuTable[i][3][k]+35*pow(SuTable[i][3][k],2)\
        //                                                          +70*pow(SuTable[i][0][k],4)*(3*pow(SuTable[i][1][k],2)+SuTable[i][3][k])+56*SuTable[i][2][k]*SuTable[i][4][k]\
        //                                                          +56*pow(SuTable[i][0][k],3)*(10*SuTable[i][1][k]*SuTable[i][2][k]+SuTable[i][4][k])+28*SuTable[i][1][k]*(10*pow(SuTable[i][2][k],2)+SuTable[i][5][k])\
        //                                                          +28*pow(SuTable[i][0][k],2)*(15*pow(SuTable[i][1][k],3)+10*pow(SuTable[i][2][k],2)+15*SuTable[i][1][k]*SuTable[i][3][k]+SuTable[i][5][k])\
        //                                                          +8 *pow(SuTable[i][0][k],1)*(105*pow(SuTable[i][1][k],2)*SuTable[i][2][k]+35*SuTable[i][2][k]*SuTable[i][3][k]+21*SuTable[i][1][k]*SuTable[i][4][k]+SuTable[i][6][k])+SuTable[i][7][k]);
        //   if(j==9)  PTable[i][j][k]=1./Factorial(j)*PTable[i][0][k]*(pow(SuTable[i][0][k],j)\
        //                                                         +36 *pow(SuTable[i][0][k],7)*SuTable[i][1][k]\
        //                                                         +84 *pow(SuTable[i][0][k],6)*SuTable[i][2][k]+1260*pow(SuTable[i][1][k],3)*SuTable[i][2][k]+280*pow(SuTable[i][2][k],3)\
        //                                                         +126*pow(SuTable[i][0][k],5)*(3*pow(SuTable[i][1][k],2)+SuTable[i][3][k])+378*pow(SuTable[i][1][k],2)*SuTable[i][4][k]+126*SuTable[i][3][k]*SuTable[i][4][k]\
        //                                                         +126*pow(SuTable[i][0][k],4)*(10*SuTable[i][1][k]*SuTable[i][2][k]+SuTable[i][4][k])+84*SuTable[i][2][k]*SuTable[i][5][k]\
        //                                                         +84 *pow(SuTable[i][0][k],3)*(15*pow(SuTable[i][1][k],3)+10*pow(SuTable[i][2][k],2)+15*SuTable[i][1][k]*SuTable[i][3][k]+SuTable[i][5][k])+36*SuTable[i][1][k]*(35*SuTable[i][2][k]*SuTable[i][3][k]+SuTable[i][6][k])\
        //                                                         +36 *pow(SuTable[i][0][k],2)*(105*pow(SuTable[i][1][k],2)*SuTable[i][2][k]+35*SuTable[i][2][k]*SuTable[i][3][k]+21*SuTable[i][1][k]*SuTable[i][4][k]+SuTable[i][6][k])\
        //                                                         +9  *pow(SuTable[i][0][k],1)*(105*pow(SuTable[i][1][k],4)+210*pow(SuTable[i][1][k],2)*SuTable[i][3][k]+35*pow(SuTable[i][3][k],2)+56*SuTable[i][2][k]*SuTable[i][4][k]+28*SuTable[i][1][k]*(10*pow(SuTable[i][2][k],2)+SuTable[i][5][k])+SuTable[i][7][k])+SuTable[i][8][k]);
        //   if(j==10) PTable[i][j][k]=1./Factorial(j)*PTable[i][0][k]*(pow(SuTable[i][0][k],j)\
        //                                                         +45 *pow(SuTable[i][0][k],8)*SuTable[i][1][k]+945*pow(SuTable[i][1][k],5)\
        //                                                         +120*pow(SuTable[i][0][k],7)*SuTable[i][2][k]+3150*pow(SuTable[i][1][k],3)*SuTable[i][3][k]+2100*pow(SuTable[i][2][k],2)*SuTable[i][3][k]\
        //                                                         +210*pow(SuTable[i][0][k],6)*(3.*pow(SuTable[i][1][k],2)+SuTable[i][3][k])+126*pow(SuTable[i][4][k],2)\
        //                                                         +252*pow(SuTable[i][0][k],5)*(10.*SuTable[i][1][k]*SuTable[i][2][k]+SuTable[i][4][k])+210*SuTable[i][3][k]*SuTable[i][5][k]+630*pow(SuTable[i][1][k],2)*(10*pow(SuTable[i][2][k],2)+SuTable[i][5][k])\
        //                                                         +210*pow(SuTable[i][0][k],4)*(15.*pow(SuTable[i][1][k],3)+10*pow(SuTable[i][2][k],2)+15*SuTable[i][1][k]*SuTable[i][3][k]+SuTable[i][5][k])+120*SuTable[i][2][k]*SuTable[i][6][k]\
        //                                                         +120*pow(SuTable[i][0][k],3)*(105*pow(SuTable[i][1][k],2)*SuTable[i][2][k]+35*SuTable[i][2][k]*SuTable[i][3][k]+21*SuTable[i][1][k]*SuTable[i][4][k]+SuTable[i][6][k])+45*SuTable[i][1][k]*(35*pow(SuTable[i][3][k],2)+56*SuTable[i][2][k]*SuTable[i][4][k]+SuTable[i][7][k])\
        //                                                         +45 *pow(SuTable[i][0][k],2)*(105*pow(SuTable[i][1][k],4)+210*pow(SuTable[i][1][k],2)*SuTable[i][3][k]+35*pow(SuTable[i][3][k],2)+56*SuTable[i][2][k]*SuTable[i][4][k]+28*SuTable[i][1][k]*(10*pow(SuTable[i][2][k],2)+SuTable[i][5][k])+SuTable[i][7][k])\
        //                                                         +10 *pow(SuTable[i][0][k],1)*(1260*pow(SuTable[i][1][k],3)*SuTable[i][2][k]+280*pow(SuTable[i][2][k],3)+378*pow(SuTable[i][1][k],2)*SuTable[i][4][k]+126*SuTable[i][3][k]*SuTable[i][4][k]+84*SuTable[i][2][k]*SuTable[i][5][k]+36*SuTable[i][1][k]*(35*SuTable[i][2][k]*SuTable[i][3][k]+SuTable[i][6][k])+SuTable[i][8][k])+SuTable[i][9][k]);
        // }

        ymin=0.0;
        ymax=ybin*(k+1);  //update current y value
        //QAG adaptive integration
        gsl_function Fun;
        Fun.function = &Function;
        Fun.params = nullptr;
        gsl_integration_workspace *workspace = gsl_integration_workspace_alloc(10000);
        res=0.0; err=0.0;
        gsl_integration_qag(&Fun, ymin, ymax, 0, 1e-3, 10000, 6, workspace, &res, &err);
        gsl_integration_workspace_free(workspace);
        SuTable[i][j][k]=res;
      }

      
      //Probility value
      if(j==0){
        ymin=0.0;
        ymax=std::log(pTbinmax/Qmed);
        //QAG adaptive integration
        gsl_function FBase;
        FBase.function = &FunctionBase;
        FBase.params = nullptr;
        gsl_integration_workspace *ws = gsl_integration_workspace_alloc(10000);
        res=0.0; err=0.0;
        gsl_integration_qag(&FBase, ymin, ymax, 0, 1e-3, 10000, 6, ws, &res, &err);
        gsl_integration_workspace_free(ws);
        P[i][j]=exp(-res);
      }
      if(j>0){
        for(int jj=0; jj<j; jj++){
          P[i][j]+=(double)(jj+1)/j*P[i][j-1-jj]*Su[i][jj]/Factorial(jj+1);
        }
      }
      // if(j>0){
      //   if(j==1)  P[i][j]=1./Factorial(j)*P[i][0]*(pow(Su[i][0],j));
      //   if(j==2)  P[i][j]=1./Factorial(j)*P[i][0]*(pow(Su[i][0],j)+Su[i][1]);
      //   if(j==3)  P[i][j]=1./Factorial(j)*P[i][0]*(pow(Su[i][0],j)\
      //                                           +3*pow(Su[i][0],1)*Su[i][1]+Su[i][2]);
      //   if(j==4)  P[i][j]=1./Factorial(j)*P[i][0]*(pow(Su[i][0],j)\
      //                                           +6*pow(Su[i][0],2)*Su[i][1]+3*pow(Su[i][1],2)\
      //                                           +4*pow(Su[i][0],1)*Su[i][2]+Su[i][3]);
      //   if(j==5)  P[i][j]=1./Factorial(j)*P[i][0]*(pow(Su[i][0],j)\
      //                                          +10*pow(Su[i][0],3)*Su[i][1]\
      //                                          +10*pow(Su[i][0],2)*Su[i][2]+10*Su[i][1]*Su[i][2]\
      //                                          +5 *pow(Su[i][0],1)*(3*pow(Su[i][1],2)+Su[i][3])+Su[i][4]);
      //   if(j==6)  P[i][j]=1./Factorial(j)*P[i][0]*(pow(Su[i][0],j)\
      //                                          +15*pow(Su[i][0],4)*Su[i][1]+15*pow(Su[i][1],3)\
      //                                          +20*pow(Su[i][0],3)*Su[i][2]+10*pow(Su[i][2],2)+15*Su[i][1]*Su[i][3]\
      //                                          +15*pow(Su[i][0],2)*(3*pow(Su[i][1],2)+Su[i][3])\
      //                                          +6 *pow(Su[i][0],1)*(10*Su[i][1]*Su[i][2]+Su[i][4])+Su[i][5]);
      //   if(j==7)  P[i][j]=1./Factorial(j)*P[i][0]*(pow(Su[i][0],j)\
      //                                          +21*pow(Su[i][0],5)*Su[i][1]\
      //                                          +35*pow(Su[i][0],4)*Su[i][2]+105*pow(Su[i][1],2)*Su[i][2]+35*Su[i][2]*Su[i][3]\
      //                                          +35*pow(Su[i][0],3)*(3*pow(Su[i][1],2)+Su[i][3])+21*Su[i][1]*Su[i][4]\
      //                                          +21*pow(Su[i][0],2)*(10*Su[i][1]*Su[i][2]+Su[i][4])\
      //                                          +7 *pow(Su[i][0],1)*(15*pow(Su[i][1],3)+10*pow(Su[i][2],2)+15*Su[i][1]*Su[i][3]+Su[i][5])+Su[i][6]);
      //   if(j==8)  P[i][j]=1./Factorial(j)*P[i][0]*(pow(Su[i][0],j)\
      //                                          +28*pow(Su[i][0],6)*Su[i][1]+105*pow(Su[i][1],4)\
      //                                          +56*pow(Su[i][0],5)*Su[i][2]+210*pow(Su[i][1],2)*Su[i][3]+35*pow(Su[i][3],2)\
      //                                          +70*pow(Su[i][0],4)*(3*pow(Su[i][1],2)+Su[i][3])+56*Su[i][2]*Su[i][4]\
      //                                          +56*pow(Su[i][0],3)*(10*Su[i][1]*Su[i][2]+Su[i][4])+28*Su[i][1]*(10*pow(Su[i][2],2)+Su[i][5])\
      //                                          +28*pow(Su[i][0],2)*(15*pow(Su[i][1],3)+10*pow(Su[i][2],2)+15*Su[i][1]*Su[i][3]+Su[i][5])\
      //                                          +8 *pow(Su[i][0],1)*(105*pow(Su[i][1],2)*Su[i][2]+35*Su[i][2]*Su[i][3]+21*Su[i][1]*Su[i][4]+Su[i][6])+Su[i][7]);
      //   if(j==9)  P[i][j]=1./Factorial(j)*P[i][0]*(pow(Su[i][0],j)\
      //                                         +36 *pow(Su[i][0],7)*Su[i][1]\
      //                                         +84 *pow(Su[i][0],6)*Su[i][2]+1260*pow(Su[i][1],3)*Su[i][2]+280*pow(Su[i][2],3)\
      //                                         +126*pow(Su[i][0],5)*(3*pow(Su[i][1],2)+Su[i][3])+378*pow(Su[i][1],2)*Su[i][4]+126*Su[i][3]*Su[i][4]\
      //                                         +126*pow(Su[i][0],4)*(10*Su[i][1]*Su[i][2]+Su[i][4])+84*Su[i][2]*Su[i][5]\
      //                                         +84 *pow(Su[i][0],3)*(15*pow(Su[i][1],3)+10*pow(Su[i][2],2)+15*Su[i][1]*Su[i][3]+Su[i][5])+36*Su[i][1]*(35*Su[i][2]*Su[i][3]+Su[i][6])\
      //                                         +36 *pow(Su[i][0],2)*(105*pow(Su[i][1],2)*Su[i][2]+35*Su[i][2]*Su[i][3]+21*Su[i][1]*Su[i][4]+Su[i][6])\
      //                                         +9  *pow(Su[i][0],1)*(105*pow(Su[i][1],4)+210*pow(Su[i][1],2)*Su[i][3]+35*pow(Su[i][3],2)+56*Su[i][2]*Su[i][4]+28*Su[i][1]*(10*pow(Su[i][2],2)+Su[i][5])+Su[i][7])+Su[i][8]);
      //   if(j==10) P[i][j]=1./Factorial(j)*P[i][0]*(pow(Su[i][0],j)\
      //                                         +45 *pow(Su[i][0],8)*Su[i][1]+945*pow(Su[i][1],5)\
      //                                         +120*pow(Su[i][0],7)*Su[i][2]+3150*pow(Su[i][1],3)*Su[i][3]+2100*pow(Su[i][2],2)*Su[i][3]\
      //                                         +210*pow(Su[i][0],6)*(3.*pow(Su[i][1],2)+Su[i][3])+126*pow(Su[i][4],2)\
      //                                         +252*pow(Su[i][0],5)*(10.*Su[i][1]*Su[i][2]+Su[i][4])+210*Su[i][3]*Su[i][5]+630*pow(Su[i][1],2)*(10*pow(Su[i][2],2)+Su[i][5])\
      //                                         +210*pow(Su[i][0],4)*(15.*pow(Su[i][1],3)+10*pow(Su[i][2],2)+15*Su[i][1]*Su[i][3]+Su[i][5])+120*Su[i][2]*Su[i][6]\
      //                                         +120*pow(Su[i][0],3)*(105*pow(Su[i][1],2)*Su[i][2]+35*Su[i][2]*Su[i][3]+21*Su[i][1]*Su[i][4]+Su[i][6])+45*Su[i][1]*(35*pow(Su[i][3],2)+56*Su[i][2]*Su[i][4]+Su[i][7])\
      //                                         +45 *pow(Su[i][0],2)*(105*pow(Su[i][1],4)+210*pow(Su[i][1],2)*Su[i][3]+35*pow(Su[i][3],2)+56*Su[i][2]*Su[i][4]+28*Su[i][1]*(10*pow(Su[i][2],2)+Su[i][5])+Su[i][7])\
      //                                         +10 *pow(Su[i][0],1)*(1260*pow(Su[i][1],3)*Su[i][2]+280*pow(Su[i][2],3)+378*pow(Su[i][1],2)*Su[i][4]+126*Su[i][3]*Su[i][4]+84*Su[i][2]*Su[i][5]+36*Su[i][1]*(35*Su[i][2]*Su[i][3]+Su[i][6])+Su[i][8])+Su[i][9]);
      // }

      ymin=0.0;
      ymax=std::log(pTbinmax/Qmed);
      //QAG adaptive integration
      gsl_function F;
      F.function = &Function;
      F.params = nullptr;
      gsl_integration_workspace *wsp = gsl_integration_workspace_alloc(10000);
      res=0.0; err=0.0;
      gsl_integration_qag(&F, ymin, ymax, 0, 1e-3, 10000, 6, wsp, &res, &err);
      gsl_integration_workspace_free(wsp);
      Su[i][j]=res;

      PSum+=P[i][j];
      nMul+=(j+1)*P[i][j];
      OutputFile << P[i][j] << " ";
      // std::cout << "pT=" << pT+pTbin/2. << ",  n=" << j+1 << ",  Pn=" << P[i][j] << ",  PSum" << PSum << std::endl;
    }
    OutputFile << nMul << std::endl;
    // std::cout << "pT=" << pT+pTbin/2. << ",  PSum=" << PSum << ",  nMul=" << nMul << std::endl;
  }

  OutputFile.close();
}



//####################################################################################################
int main(){

  CalProbability();

  return 0;
}
