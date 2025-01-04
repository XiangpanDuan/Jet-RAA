#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <cmath>
#include <iomanip>
#include <chrono>
#include "DiJetLOELoss.h"


namespace CalDiJetLOFixedELoss{

  //Verify TAB in Glauber model
  double MyFun_VerifyTAB(double range[], size_t dim, void *params){
    DiJetLOELoss *jetp = (DiJetLOELoss*) params;
    jetp->setGlauber(range[0],range[1]);
    return jetp->getGlauber();
  }


  //------------------------------------------------------------
  //Cross section calculation: d^2sigma/dpTdy
  void CalSigma(){

    //Record start time
    auto start = std::chrono::high_resolution_clock::now();

    const int nEloss=Input::nEloss;

    //Open Multi-Processing
    #pragma omp parallel for num_threads(nEloss)
    for(int j=0; j<nEloss; j++){

      Particle *Parton  = new Particle(Input::name);
      DiJetLOELoss *Jet = new DiJetLOELoss(Input::A,Input::B,Input::Cen,Input::Ecm,Parton);  //PbPb 5020GeV
    
      //Set initial conditions
      Jet->setpTScale(std::pow(2.,Input::scale));  //pT scale to control the pT errer bar
      Jet->setNf(Input::nf);
      Jet->setLambdaQCD(Input::nloop);
      //Set energy loss
      Jet->setRsize(Input::Rsize);
      Jet->setalphaSmed(Input::alphaSmed);
      Jet->setQmed(Input::Qmed);

      //Range of x,y,theta in Glauber model
      double xmin=Input::xmin;
      double xmax=Input::xmax;
      double ymin=Input::ymin;
      double ymax=Input::ymax;
      double thetamin=Input::thetamin;
      double thetamax=Input::thetamax;
        
      size_t dim=3;
      double rangemin[]={xmin, ymin, thetamin};
      double rangemax[]={xmax, ymax, thetamax};
        
      size_t calls=Input::calls;
      Jet->setCalls(calls);

      double res,err;
      Jet->setupMC(&MyFun_VerifyTAB, dim, rangemin, rangemax, Jet);
      Jet->calculateMC(res,err);
      Jet->cleanupMC();

      std::cout << "Determine whether the output is normalized: " << res << std::endl;

      delete Jet;

      //Record end time
      auto end = std::chrono::high_resolution_clock::now();
      std::chrono::duration<double> duration = end - start;
      std::cout << "Function execution time: " << duration.count() << " seconds" << std::endl;
    }
    
  }

}



//####################################################################################################
//Main function
int main(){

  std::cout << "######################################################################" << std::endl;
  std::cout << "# Function Test                                                      #" << std::endl;
  std::cout << "######################################################################" << std::endl;

  if(MODE==0 || GLB==0){
    std::cerr << "Error in Test.cpp!!!" << std::endl;
    exit(EXIT_FAILURE);
  }


  CalDiJetLOFixedELoss::CalSigma();
    
  return 0;
}


