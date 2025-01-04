#include <iostream>
#include <fstream>
#include <sstream>
#include <iomanip>
#include <cmath>
#include "TrilinearInterpolation.h"

using namespace std;

double bx,by;

double OmegaC(double range[], size_t dim, void *p){
  TrilinearInterpolation *TriInter=(TrilinearInterpolation*)p;
  TriInter->setOmegaC(range[0]);
  return TriInter->getOmegaC();
}


//*************************************************************************************************************************************************************************
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
  int dim=1;
  double rangemin[]={taumin};
  double rangemax[]={taumax};

  int count=0;
  ofstream output;
  output.open("../Output/x_y_theta_omegaC_5_test301.dat", ofstream::out);
  int xnum=301,ynum=301,thetanum=36;
  double xmin=-15.,xmax=15.;
  double ymin=-15.,ymax=15.;
  double xbin=(xmax-xmin)/(xnum-1);
  double ybin=(ymax-ymin)/(ynum-1);
  double thetabin=2.*M_PI/thetanum; //divide thetabin
  //Known conditions: ix,iy,itheta
  for(int ix=0; ix<xnum; ix++){
    // if(ix%10!=0) continue; //control output data format
    for(int iy=0; iy<ynum; iy++){
      // if(iy%10!=0) continue; //control output data format
      for(int itheta=0; itheta<thetanum; itheta++){
        
        double x=ix*xbin+xmin;
        double y=iy*ybin+ymin;
        double theta=itheta*thetabin; //attention thetabin
        // double theta=(itheta+0.5)*thetabin; //same set with Lin Chen
        MyTrilinear->setInputData(x,y,theta); //set known table (x,y,theta)

        size_t calls = 100000;
        MyTrilinear->setCalls(calls);
        double res,err;
        //Calculate the gluon frequency (OmagaC)
        MyTrilinear->setupMC(&OmegaC, dim, rangemin, rangemax, MyTrilinear);
        MyTrilinear->calculate(res,err);
       
        output << fixed << setprecision(4) << x << "   " << y << "   " << setprecision(16) << theta << "   " << res << endl; //control radix point number
        output.unsetf(ios::fixed);

        count=count+1;
        // if(itheta==0) cout << count << "  " << ix << "  " << " OmegaC = " << res << "   " << err << endl;
        cout << count << "  " << x << "  " << y << "  " << theta << "  " << " OmegaC = " << res << "   " << err << endl;
        cout << endl;

        MyTrilinear->cleanupMC();
      }
    }
  }
  output.close();

}



//*************************************************************************************************************************************************************************
int main(){

  calcTrilinearInterpolation();

  return 0;
}
