#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <cmath>
#include <iomanip>
#include <chrono>
#include "DiJetLO.h"


namespace CalDiJetLO{

  double SigQuark(double range[], size_t dim, void *params){
    DiJetLO *jetp = (DiJetLO*) params;
    // std::string type="quark";
    jetp->setpT(range[0]);
    jetp->setLOParameters(range[1],range[2]);

    return jetp->SigmaLOQuark();
  }

  double SigGluon(double range[], size_t dim, void *params){
    DiJetLO *jetp = (DiJetLO*) params;
    // std::string type="gluon";
    jetp->setpT(range[0]);
    jetp->setLOParameters(range[1],range[2]);
    
    return jetp->SigmaLOGluon();
  }


  //------------------------------------------------------------
  //Cross section calculation: d^2sigma/dpTdy
  void CalSigma(){

    //Record start time
    auto start = std::chrono::high_resolution_clock::now();

    //Output files
    std::stringstream ss[3];
    std::string OutputFile[3];
    std::ofstream OutStream[3];
    for(int i=0; i<3; i++){
      ss[i] << "../Output/pp/Test_dsigma_dpTdy_" << Input::Ecm << "GeV_type" << i << ".dat";
      OutputFile[i]=ss[i].str();
      OutStream[i].open(OutputFile[i], std::ofstream::out);
    }


    Particle *Parton = new Particle(Input::name);
    DiJetLO  *Jet    = new DiJetLO(Input::Ecm,Parton);  //pp 5020GeV
    
    //Set initial conditions
    Jet->setpTScale(std::pow(2.,Input::scale));
    Jet->setNf(Input::nf);
    Jet->setLambdaQCD(Input::nloop);

    //Range of jet rapidity
    double rapmin3=-Input::rap3;
    double rapmax3= Input::rap3;
    double rapmin4=-Input::rap4;
    double rapmax4= Input::rap4;
    // double drap3=rapmax3-rapmin3;  //triggered jet rapidity
    double drap4=rapmax4-rapmin4;     //observed  jet rapidity

    //Final observed jet momentum
    // const int pTnum=21;
    // double pTRange[pTnum+1]={50.,60.,70.,80.,90.,100.,112.,125.,141.,158.,177.,199.,223.,251.,281.,316.,354.,398.,501.,630.,999.,1200.};
    // double pT,pTbin,pTmin,pTmax;
    int    pTnum=Input::pTnum;
    double pTmin=Input::pTmin;
    double pTmax=Input::pTmax;
    double pTbin=(pTmax-pTmin)/pTnum;
    for(int k=0; k<pTnum; k++){
      double pT=pTmin+pTbin/2.+pTbin*k;
      double pTbinmin=pT-pTbin/2.;
      double pTbinmax=pT+pTbin/2.;
      // pTbinmin=pTRange[i];
      // pTbinmax=pTRange[i+1];
      // pTbin=pTbinmax-pTbinmin;
      // pT=(pTbinmin+pTbinmax)/2.;

      //Range of the integration: from rangemin to rangemax
      //pp collisions
      size_t dim=3;
      double rangemin[]={pTbinmin, rapmin3, rapmin4};
      double rangemax[]={pTbinmax, rapmax3, rapmax4};
      std::cout << pT << "  " << dim << std::endl;
      
      size_t calls = Input::calls;
      Jet->setCalls(calls);

      double res_q,err_q;
      double res_g,err_g;

      //Quark cross section
      Jet->setupMC(&SigQuark, dim, rangemin, rangemax, Jet);
      Jet->calculateMC(res_q,err_q);
      Jet->cleanupMC();
      // std::cout << "quark d^2sigma/dpTdy = " <<  res_q/(pTbin*drap4)*0.001 << " nb with err = " << err_q/(pTbin*drap4)*0.001 << " nb." << std::endl;
      OutStream[1] << pT << "   " << res_q/(pTbin*drap4)*0.001 << "   " << err_q/(pTbin*drap4)*0.001 << std::endl;

      //Gluon cross section
      Jet->setupMC(&SigGluon, dim, rangemin, rangemax, Jet);
      Jet->calculateMC(res_g,err_g);
      Jet->cleanupMC();
      // std::cout << "gluon d^2sigma/dpTdy = " <<  res_g/(pTbin*drap4)*0.001 << " nb with err = " << err_g/(pTbin*drap4)*0.001 << " nb." << std::endl;
      OutStream[2] << pT << "   " << res_g/(pTbin*drap4)*0.001 << "   " << err_g/(pTbin*drap4)*0.001 << std::endl;

      //Total cross section
      std::cout << "total d^2sigma/dpTdy = " <<  (res_q+res_g)/(pTbin*drap4)*0.001 << " nb with err = " << (err_q+err_g)/(pTbin*drap4)*0.001 << " nb." << std::endl;
      OutStream[0] << pT << "   " << (res_q+res_g)/(pTbin*drap4)*0.001 << "   " << (err_q+err_g)/(pTbin*drap4)*0.001 << std::endl;


      //Record end time
      auto end=std::chrono::high_resolution_clock::now();
      std::chrono::duration<double> duration=end-start;
      std::cout << "Function execution time: " << duration.count() << " seconds" << std::endl;
    }
    delete Jet;

    for(int i=0; i<3; i++){
      OutStream[i].close();
    }
    
  }

}



//####################################################################################################
//Main function
int main(){

  std::cout << "######################################################################" << std::endl;
  std::cout << "# Leading-Order Cross Section Calculation in pp collisions           #" << std::endl;
  std::cout << "######################################################################" << std::endl;

  if(MODE!=0){
    std::cerr << "Error in pp collisions!!!" << std::endl;
    exit(EXIT_FAILURE);
  }
    
    
  CalDiJetLO::CalSigma();

  return 0;
}