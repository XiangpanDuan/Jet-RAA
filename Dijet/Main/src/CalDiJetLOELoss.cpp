#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <cmath>
#include <iomanip>
#include <chrono>
#include "DiJetLOELoss.h"


namespace CalDiJetLOELoss{

  double SigQuark(double *range, size_t dim, void *params){
    DiJetLOELoss *jetp = (DiJetLOELoss*) params;
    std::string type="quark";
    //Glauber model
    jetp->setGlauber(range[3],range[4]);
    //Energy loss in BDMPS
    jetp->setEnergyLoss(type,range,dim);
    //Kinematic parameters in leading order
    jetp->setParameters(range[1],range[2]);

    return jetp->SigmaLOELossQuark();
  }

  double SigGluon(double *range, size_t dim, void *params){
    DiJetLOELoss *jetp = (DiJetLOELoss*) params;
    std::string type="gluon";
    //Glauber model
    jetp->setGlauber(range[3],range[4]);
    //Energy loss in BDMPS
    jetp->setEnergyLoss(type,range,dim);
    //Kinematic parameters in leading order
    jetp->setParameters(range[1],range[2]);

    return jetp->SigmaLOELossGluon();
  }


  //------------------------------------------------------------
  //Cross section calculation: d^2sigma/dpTdy
  void CalSigma(){

    //Record start time
    auto start = std::chrono::high_resolution_clock::now();

    const int nEloss=Input::nEloss;

    //Output files
    std::stringstream ss[3][nEloss];
    std::string OutputFile[3][nEloss];
    std::ofstream OutStream[3][nEloss];
    for(int i=0; i<3; i++){
      for(int j=0; j<nEloss; j++){
        if(nEloss==1){
          if(MULTI==0 && TAUS==0) ss[i][j] << "../Output/PbPb/dsigma_dpTdy_" << Input::Ecm << "GeV_qhat" << Input::qhat0 << "_single_type" << i << ".dat";
          if(MULTI==1 && TAUS==0) ss[i][j] << "../Output/PbPb/dsigma_dpTdy_" << Input::Ecm << "GeV_qhat" << Input::qhat0 << "_multi_Qmed"  << Input::Qmed << "_type" << i << ".dat";
          if(MULTI==0 && TAUS==1) ss[i][j] << "../Output/PbPb/dsigma_dpTdy_" << Input::Ecm << "GeV_qhat" << Input::qhat0 << "_tau_single_type" << i << ".dat";
          if(MULTI==1 && TAUS==1) ss[i][j] << "../Output/PbPb/dsigma_dpTdy_" << Input::Ecm << "GeV_qhat" << Input::qhat0 << "_tau_multi_Qmed"  << Input::Qmed << "_type" << i << ".dat";
        }
        if(nEloss!=1){
          if(MULTI==2 && TAUS==0) ss[i][j] << "../Output/PbPb/Probability/dsigma_dpTdy_" << Input::Ecm << "GeV_qhat" << Input::qhat0 << "_probability_type" << i << j << ".dat";
          if(MULTI==2 && TAUS==1) ss[i][j] << "../Output/PbPb/Probability/dsigma_dpTdy_" << Input::Ecm << "GeV_qhat" << Input::qhat0 << "_tau_probability_type" << i << j << ".dat";
        }
        OutputFile[i][j]=ss[i][j].str();
        OutStream[i][j].open(OutputFile[i][j], std::ofstream::out);
      }
    }

    //Open Multi-Processing
    #pragma omp parallel for num_threads(nEloss)
    for(int j=0; j<nEloss; j++){

      Particle *Parton  = new Particle(Input::name);
      DiJetLOELoss *Jet = new DiJetLOELoss(Input::A,Input::B,Input::Cen,Input::Ecm,Parton);
    
      //Set initial conditions
      Jet->setpTScale(std::pow(2.,Input::scale));  //pT scale to control the pT errer bar
      Jet->setNf(Input::nf);
      Jet->setLambdaQCD(Input::nloop);
      //Set energy loss
      Jet->setRsize(Input::Rsize);
      Jet->setalphaSmed(Input::alphaSmed);
      Jet->setQmed(Input::Qmed);
      if(GLB==1)  Jet->setqhat0(Input::qhat0);
      if(TAUS==0) Jet->setOmegaCTable();  //energy loss table
      if(TAUS==1) Jet->setTempTable();    //temperature table

      // //Test tau0
      // const int pTnum=100;
      // double pT,pTbin;
      // for(int i=0; i<pTnum; i++){
      //   pTbin=10.;
      //   pT=55.+pTbin*i;
      //   Jet->setpT(pT);
      //   Jet->setTauStart();
      //   Jet->setOmegaCTau(0.,0.,0.);
      // }

      //Range of jet rapidity
      double rapmin3=-Input::rap3;
      double rapmax3= Input::rap3;
      double rapmin4=-Input::rap4;
      double rapmax4= Input::rap4;
      // double drap3=rapmax3-rapmin3;  //triggered jet rapidity
      double drap4=rapmax4-rapmin4;     //observed  jet rapidity
      //Range of jet energy loss
      double elossmin=Input::elossmin;
      double elossmax=Input::elossmax;
      //Range of x,y,theta in Glauber model
      double xmin=Input::xmin;
      double xmax=Input::xmax;
      double ymin=Input::ymin;
      double ymax=Input::ymax;
      double thetamin=Input::thetamin;
      double thetamax=Input::thetamax;
    
      //Final observed jet momentum
      // const int pTnum=20;
      // double pTRange[pTnum+1]={50.,60.,70.,80.,90.,100.,112.,125.,141.,158.,177.,199.,223.,251.,281.,316.,354.,398.,501.,630.,999.};
      // double pT,pTbin,pTmin,pTmax;
      int    pTnum=Input::pTnum;
      double pTmin=Input::pTmin;
      double pTmax=Input::pTmax;
      double pTbin=(pTmax-pTmin)/pTnum;
      for(int k=0; k<pTnum; k++){
        double pT=pTmin+pTbin/2.+pTbin*k;
        double pTbinmin=pT-pTbin/2.;
        double pTbinmax=pT+pTbin/2.;
        // double pTbinmin=pTRange[k];
        // double pTbinmax=pTRange[k+1];
        // double pTbin=pTbinmax-pTbinmin;
        // double pT=(pTbinmin+pTbinmax)/2.;
      
        //Range of the integration: from rangemin to rangemax
        size_t dim=j+7;
        double *rangemin = new double[dim];
        double *rangemax = new double[dim];
        for(int l=0; l<dim; l++){
          if(l==0) {rangemin[l]=pTbinmin; rangemax[l]=pTbinmax;}
          if(l==1) {rangemin[l]=rapmin3;  rangemax[l]=rapmax3;}
          if(l==2) {rangemin[l]=rapmin4;  rangemax[l]=rapmax4;}
          if(l==3) {rangemin[l]=xmin;     rangemax[l]=xmax;}
          if(l==4) {rangemin[l]=ymin;     rangemax[l]=ymax;}
          if(l==5) {rangemin[l]=thetamin; rangemax[l]=thetamax;}
          if(l>=6) {rangemin[l]=elossmin; rangemax[l]=elossmax;}
        }
        std::cout << pT << "  " << dim << std::endl;
        
        size_t calls = Input::calls+20000*(dim-7);
        Jet->setCalls(calls);

        double resq,errq;
        double resg,errg;

        //Quark cross section
        Jet->setupMC(&SigQuark, dim, rangemin, rangemax, Jet);
        Jet->calculateMC(resq,errq);
        Jet->cleanupMC();
        // std::cout << "quark d^2sigma/dpTdy = " <<  resq/(pTbin*drap4)*0.001 << " nb with err = " << errq/(pTbin*drap4)*0.001 << " nb." << std::endl;
        OutStream[1][j] << pT << "   " << resq/(pTbin*drap4)*0.001 << "   " << errq/(pTbin*drap4)*0.001 << std::endl;

        //Gluon cross section
        Jet->setupMC(&SigGluon, dim, rangemin, rangemax, Jet);
        Jet->calculateMC(resg,errg);
        Jet->cleanupMC();
        // std::cout << "gluon d^2sigma/dpTdy = " <<  resg/(pTbin*drap4)*0.001 << " nb with err = " << errg/(pTbin*drap4)*0.001 << " nb." << std::endl;
        OutStream[2][j] << pT << "   " << resg/(pTbin*drap4)*0.001 << "   " << errg/(pTbin*drap4)*0.001 << std::endl;

        //Total cross section
        std::cout << "total d^2sigma/dpTdy = " <<  (resq+resg)/(pTbin*drap4)*0.001 << " nb with err = " << (errq+errg)/(pTbin*drap4)*0.001 << " nb." << std::endl;
        OutStream[0][j] << pT << "   " << (resq+resg)/(pTbin*drap4)*0.001 << "   " << (errq+errg)/(pTbin*drap4)*0.001 << std::endl;
      
        delete [] rangemin;
        delete [] rangemax;
      }
      delete Jet;

      //Record end time
      auto end = std::chrono::high_resolution_clock::now();
      std::chrono::duration<double> duration = end - start;
      std::cout << "Function execution time: " << duration.count() << " seconds" << std::endl;
    }

    for(int i=0; i<3; i++){
      for(int j=0; j<nEloss; j++){
        OutStream[i][j].close();
      }
    }
    
  }

}



//####################################################################################################
//Main function
int main(){

  std::cout << "######################################################################" << std::endl;
  std::cout << "# Leading-Order Cross Section Calculation with BDMPS Energy Loss     #" << std::endl;
  std::cout << "######################################################################" << std::endl;

  if(MODE==0 || GLB==0){
    std::cerr << "Error in AA collisions with BDMPS Energy Loss!!!" << std::endl;
    exit(EXIT_FAILURE);
  }


  CalDiJetLOELoss::CalSigma();
    
  return 0;
}


