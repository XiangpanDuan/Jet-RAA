#include <iostream>
#include <fstream>
#include <sstream>
#include <iomanip>
#include <cmath>
#include "TrilinearInterpolation.h"


namespace CalOmegaC{

  double OmegaC(double range, void *params){
    TrilinearInterpolation *TriInter = (TrilinearInterpolation*) params;
    TriInter->setOmegaC(range);
    return TriInter->getOmegaC();
  }


  //------------------------------------------------------------
  void calcTrilinearInterpolation(){

    double qhat0=1.0;  //adjusted in the later calculations
    TrilinearInterpolation *MyTrilinear = new TrilinearInterpolation(qhat0);
    // // test
    // MyTrilinear->TriInterpolation(0.25,0.2,0.6);
    // MyTrilinear->TriInterpolation(0.25,0.2,6.6);

    //Range of tau
    double taumin=0.6;
    double taumax=20.;
    //Range of the integration: from rangemin to rangemax
    double rangemin=taumin;
    double rangemax=taumax;

    int count=0;
    std::ofstream OutputFile;
    OutputFile.open("../OutputFile/x_y_theta_omegaC.dat", std::ofstream::out);
    int xnum=301,ynum=301,thetanum=36;
    double xmin=-15.,xmax=15.;
    double ymin=-15.,ymax=15.;
    double xbin=(xmax-xmin)/(xnum-1);
    double ybin=(ymax-ymin)/(ynum-1);
    double thetabin=2.*M_PI/thetanum;  //divide thetabin
    //Known conditions: ix,iy,itheta
    for(int ix=0; ix<xnum; ix++){
      // if(ix%10!=0) continue;  //control OutputFile data format
	    for(int iy=0; iy<ynum; iy++){
        // if(iy%10!=0) continue;  //control OutputFile data format
	  	  for(int itheta=0; itheta<thetanum; itheta++){
          
          double x=ix*xbin+xmin;
          double y=iy*ybin+ymin;
          double theta=itheta*thetabin;  //attention thetabin
          // double theta=(itheta+0.5)*thetabin;  //attention thetabin
          MyTrilinear->setInputData(x,y,theta);  //set known table (x,y,theta)

        
          //Calculate the gluon frequency (OmagaC)
          double res=0.0, err=0.0;
          MyTrilinear->setupQAG(&OmegaC, rangemin, rangemax, MyTrilinear);
          MyTrilinear->calculateQAG(res,err);
       
          OutputFile << std::fixed << std::setprecision(4) << x << "   " << y << "   " << std::setprecision(16) << theta << "   " << res << std::endl; //control radix point number
          OutputFile.unsetf(std::ios::fixed);

          count=count+1;
          // if(itheta==0) std::cout << count << "  " << ix << "  " << " OmegaC = " << res << "   " << err << std::endl;
          std::cout << count << "  " << x << "  " << y << "  " << theta << "  " << " OmegaC = " << res << "   " << err << std::endl;
          std::cout << std::endl;

          MyTrilinear->cleanupQAG();
        }
      }
    }

    OutputFile.close();
  }

}



//####################################################################################################
//Main function
int main(){

  CalOmegaC::calcTrilinearInterpolation();

  return 0;
}
