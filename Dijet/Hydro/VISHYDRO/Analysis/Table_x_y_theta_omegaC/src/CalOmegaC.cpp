#include <iostream>
#include <fstream>
#include <iomanip>
#include <cmath>
#include "OmegaC.h"


namespace CalOmegaC{

  double Function(double range, void *para){
    OmegaC *mypara = (OmegaC*) para;
    mypara->setOmegaC(range);
    return mypara->getOmegaC();
  }


  //------------------------------------------------------------
  void CalculateOmegaC(){

    OmegaC *MyOmegaC = new OmegaC();
    double qhat0=1.0;
    MyOmegaC->setqhat0(qhat0);

    //Range of tau
    double taumin=0.6;
    double taumax=20.;
    //Range of the integration: from rangemin to rangemax
    double rangemin=taumin;
    double rangemax=taumax;

    int count=0;
    std::ofstream OutputFile;
    OutputFile.open("../Output/x_y_theta_omegaC.dat", std::ofstream::out);
    int xnum=301,ynum=301,thetanum=36;
    double xmin=-15.,xmax=15.;
    double ymin=-15.,ymax=15.;
    double xbin=(xmax-xmin)/(xnum-1);
    double ybin=(ymax-ymin)/(ynum-1);
    double thetabin=2.*M_PI/thetanum;  //divide thetabin
    //Known conditions: ix,iy,itheta
    for(int ix=0; ix<xnum; ix++){
      for(int iy=0; iy<ynum; iy++){
        for(int itheta=0; itheta<thetanum; itheta++){
          
          double x=ix*xbin+xmin;
          double y=iy*ybin+ymin;
          double theta=itheta*thetabin;  //attention thetabin
          // double theta=(itheta+0.5)*thetabin;  //attention thetabin
          MyOmegaC->setxytheta(x,y,theta);
        
          //Calculate the gluon frequency (OmagaC)
          double res=0.0, err=0.0;
          MyOmegaC->setupQAG(&Function, rangemin, rangemax, MyOmegaC);
          MyOmegaC->calculateQAG(res,err);
       
          OutputFile << std::fixed << std::setprecision(4) << x << "   " << y << "   " << std::setprecision(16) << theta << "   " << res << std::endl;
          OutputFile.unsetf(std::ios::fixed);

          count+=1;
          std::cout << count << "  " << x << "  " << y << "  " << theta << "  " << " OmegaC = " << res << "   " << err << std::endl;
          std::cout << std::endl;

          MyOmegaC->cleanupQAG();
        }
      }
    }

    OutputFile.close();
  }

}



//####################################################################################################
//Main function
int main(){

  CalOmegaC::CalculateOmegaC();

  return 0;
}
