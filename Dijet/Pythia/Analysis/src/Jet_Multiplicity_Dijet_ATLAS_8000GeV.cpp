#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <cstdlib>
#include <vector>
#include <cmath>
#include "TFile.h"
#include "TH1D.h"
#include "fastjet/PseudoJet.hh"
#include "fastjet/ClusterSequence.hh"
#include "fastjet/Selector.hh"


//####################################################################################################
//Jet cutoff
int    mode=0;                            //0:parton; 1:hadron
double jet_ptmin=50.00;                   //pt lowlimit for full jets
double particle_ptmin=0.5;                //pt lowlimit for input particles
double constituent_ptmin=0.5;             //pt lowlimit for constituent particles with pT > 0.5, 2.0, 5.0 GeV
double Rsize=0.4;                         //jet cone size
double particle_etamax=2.5;               //eta upperlimit for input particles
double jet_etamax=particle_etamax-Rsize;  //eta upperlimit for full jets
int    jet_nhardest=2;                    //n hardest jets

//Histograms
TH1D *NumEvent = new TH1D("NumEvent","NumEvent",1,0.,1.);
TH1D *NumJet   = new TH1D("NumJet",  "NumJet",  1,0.,1.);
TH1D *pTDis    = new TH1D("pTDis",   "pTDis",   200,0.,2000.);
TH1D *EtaDis   = new TH1D("EtaDis",  "EtaDis",  100,-5.,5.);
TH1D *PhiDis   = new TH1D("PhiDis",  "PhiDis",  36,-M_PI,M_PI);
TH1D *SigmaDis = new TH1D("SigmaDis","SigmaDis",200,0.,2000.);
//BinSets
const int ntype=3;
const int npTbin=12;
TH1D *Multi[ntype];
TH1D *MultipT[ntype][npTbin];
double pTrange[13]={0.,50.,100.,200.,300.,400.,500.,600.,800.,1000.,1200.,1500.,1800.};
// double pTrange[25]={0.,50.,60.,70.,80.,90.,100.,112.,125.,141.,158.,177.,199.,223.,251.,281.,316.,354.,398.,501.,630.,999.,1200.,1500.,1800.};


//####################################################################################################
//Function declaration
double GetCharge(const int particle_id);
void WriteHistograms(const std::string OutputString);
void DeleteHistograms();


//####################################################################################################
int main(int argc, char **argv)
{
	if(argc!=1){
        std::cerr << "Usage error!" << std::endl;
        exit(EXIT_FAILURE);
    }

	//------------------------------------------------------------
	//Histograms set
	int multibin=100;
	double multimin=0., multimax=100.;
	for(int i=0; i<ntype; i++){
		Multi[i] = new TH1D(Form("Multi%d",i),Form("Multi%d",i),2000,0.,2000.);
		for(int j=0; j<npTbin; j++){
			MultipT[i][j] = new TH1D(Form("Multi%d_pT%d",i,j),Form("Multi%d_pT%d",i,j),multibin,multimin,multimax);
		}
	}
	
	//------------------------------------------------------------
	//Data input set
	double ptdiff=0.0;
	if     (jet_ptmin==50.)   ptdiff=5.;
	else if(jet_ptmin==200.)  ptdiff=30.;
	else if(jet_ptmin==500.)  ptdiff=50.;
	else if(jet_ptmin==1000.) ptdiff=50.;
	std::stringstream ssin,ssout;
	if(mode==0) ssout << "../Output/ATLAS_DiJet_Parton_8000GeV_" << jet_ptmin << "GeV_R0p4_pTmin" << particle_ptmin << "GeV_ConpTmin" << constituent_ptmin << "GeV_Tune21_Scale_pTmin2GeV";
	// if(mode==0) ssout << "../Output/ATLAS_DiJet_Parton_8000GeV_" << jet_ptmin << "GeV_R0p4_pTmin" << particle_ptmin << "GeV_ConpTmin" << constituent_ptmin << "GeV_Tune21";
	if(mode==1) ssout << "../Output/ATLAS_DiJet_Hadron_8000GeV_" << jet_ptmin << "GeV_R0p4_pTmin" << particle_ptmin << "GeV_ConpTmin" << constituent_ptmin << "GeV_Tune21";
	std::string OutputString=ssout.str();
	if(mode==0) ssin  << "../../Simulation/Output/ATLAS_DiJet_Parton_8000GeV_" << jet_ptmin-ptdiff << "GeV_Tune21_Scale_pTmin2GeV.dat";
	// if(mode==0) ssin  << "../../Simulation/Output/ATLAS_DiJet_Parton_8000GeV_" << jet_ptmin-ptdiff << "GeV_Tune21.dat";
	if(mode==1) ssin  << "../../Simulation/Output/ATLAS_DiJet_Hadron_8000GeV_" << jet_ptmin-ptdiff << "GeV_Tune21.dat";
	std::string InputString=ssin.str();
	std::ifstream InputFile;
	InputFile.open(InputString);
	if(!InputFile){
        std::cerr << "Error: One could not open the input file!" << std::endl;
        exit(EXIT_FAILURE);
    }

	//------------------------------------------------------------
    //Event loop
	int nEvents=200000;
	std::cout << "Total events: " << nEvents << std::endl;
	for(int iEvent=1; iEvent<=nEvents; iEvent++){
		//Write current working event
		if(iEvent%5000==0) std::cout << "Working on event #" << iEvent << std::endl;
		
		//Pythia variables
		int    ievent,nTracks;
		double sigma,sigmaerr;
		int    Hard_status[4],Hard_id[4];
		double Hard_px[4],Hard_py[4],Hard_pz[4],Hard_e[4],Hard_m[4],Hard_eta[4],Hard_phi[4];
		int    index,id;
		double px,py,pz,energy,mass;

		//Particle loop
		std::string str;
		InputFile >> str;
		if(str=="#"){
			InputFile >> ievent >> nTracks >> sigma >> sigmaerr;
			sigma*=1.e+6;  //transform unit from mb to nb
			//Hard process
			for(int i=0; i<4; i++){
				std::string sstr;
				InputFile >> sstr;
				if(sstr=="##") {InputFile >> Hard_status[i];}
				if(Hard_status[i]==-21 || Hard_status[i]==-23){
					InputFile >> Hard_id[i] >> Hard_px[i] >> Hard_py[i] >> Hard_pz[i] >> Hard_e[i] >> Hard_m[i] >> Hard_eta[i] >> Hard_phi[i];
				}
				else{
					std::cerr << "Hard process is error!" << std::endl;
					exit(EXIT_FAILURE);
				}
			}

			//Read in input particles
			std::vector<fastjet::PseudoJet> input_particles;
			for(int iTrack=0; iTrack<nTracks; iTrack++){
				InputFile >> index >> id >> px >> py >> pz >> energy >> mass;
				// InputFile >> index >> id >> charge >> px >> py >> pz >> energy >> mass;
    			//Push event onto back of input_particles vector
                fastjet::PseudoJet particle(px,py,pz,energy);
                particle.set_user_index(id);
				if(std::abs(particle.eta())<particle_etamax && particle.pt()>particle_ptmin){  //all particles
   			    // if(std::abs(particle.eta())<particle_etamax && particle.pt()>particle_ptmin && GetCharge(id)!=0){  //charged particles
					input_particles.push_back(particle);
				}
			}
			

			//------------------------------------------------------------
			//Create a jet definition for the clustering
			fastjet::JetDefinition jet_def(fastjet::antikt_algorithm, Rsize);
			//Run the jet clustering with the above jet definition
    		fastjet::ClusterSequence clust_seq(input_particles, jet_def);
			//Get the resulting jets ordered in pt
			std::vector<fastjet::PseudoJet> inclusive_jets = sorted_by_pt(clust_seq.inclusive_jets(jet_ptmin));
			//Selector specifying the jet range
			fastjet::Selector jet_selector = fastjet::SelectorAbsEtaMax(jet_etamax) * fastjet::SelectorNHardest(jet_nhardest);  //select jet_etamax > |eta| and 2 jets with largest transverse momenta
			inclusive_jets  = jet_selector(inclusive_jets);


			//------------------------------------------------------------
			//Jet analysis start
			if(inclusive_jets.size()!=2) continue;                             //ATLAS condition
			if((inclusive_jets[0].pt()/inclusive_jets[1].pt())>1.5) continue;  //ATLAS condition
			NumEvent->Fill(0.5,1.0);
			for(unsigned int iJet=0; iJet<inclusive_jets.size(); iJet++){
				//Find back-to-back jets and label as quark or gluon jet
				int lable=0;
				for(int i=2; i<4; i++){
					if(std::abs(Hard_phi[i]-inclusive_jets[iJet].phi_std())<(1./3.*M_PI) || std::abs(Hard_phi[i]-inclusive_jets[iJet].phi_std())>(5./3.*M_PI)) lable=i;
				}
				if(lable==0) break;

				NumJet->Fill(0.5,1.0);
	  			pTDis->Fill(inclusive_jets[iJet].pt(),1./10.);
				EtaDis->Fill(inclusive_jets[iJet].eta(),1./0.1);
	  			PhiDis->Fill(inclusive_jets[iJet].phi_std(),1./(2.*M_PI/36.));
				SigmaDis->Fill(inclusive_jets[iJet].pt(),1./(10.*2.*jet_etamax));

				//pT range
				int ptbin=0;  //invalid the minimum pt bin
				for(int j=0; j<npTbin; j++){
					if(inclusive_jets[iJet].pt()>=pTrange[j] && inclusive_jets[iJet].pt()<pTrange[j+1]) {ptbin=j; break;}
				}

				int sumcharge=0;
				std::vector<fastjet::PseudoJet> constituents=inclusive_jets[iJet].constituents();
				for(unsigned int iCon=0; iCon<constituents.size(); iCon++){
					int ischarge=0;
					if(mode==0 && constituents[iCon].user_index()!=22)            ischarge=1;  //parton
					if(mode==1 && GetCharge(constituents[iCon].user_index())!=0.) ischarge=1;  //hadron
					if(ischarge==1 && constituents[iCon].pt()>constituent_ptmin){
						sumcharge+=1;
					}
				}
				// if(sumcharge==0) continue;  //No 0
				//Multiplicity
				MultipT[0][ptbin]->Fill(sumcharge,1.0);                                                                 //all jet
				if(std::abs(Hard_id[lable])>=1 && std::abs(Hard_id[lable])<=6) MultipT[1][ptbin]->Fill(sumcharge,1.0);  //quark jet
				if(std::abs(Hard_id[lable])==21) MultipT[2][ptbin]->Fill(sumcharge,1.0);                                //gluon jet
			}
        
		}//End of particle loop

		if(iEvent==nEvents) SigmaDis->Scale(sigma/nEvents);  //cross section

	}//End of event loop
    
	
	//------------------------------------------------------------
	//Histogram normalization
	for(int i=0; i<ntype; i++){
		Multi[i]->SetBins(npTbin,pTrange);
		for(int j=1; j<npTbin; j++){  //j!=0: delete the minimum pt bin
			Multi[i]->SetBinContent(j+1,MultipT[i][j]->GetMean());
			Multi[i]->SetBinError(j+1,MultipT[i][j]->GetMeanError());
		}
	}


	//------------------------------------------------------------
	//Save data and free memory
	WriteHistograms(OutputString);
    DeleteHistograms();
	InputFile.close();
	
    return 0;
}


//####################################################################################################
double GetCharge(const int particle_id)
{
	double charge=0.0;
	
	//QUARKS
	if      (particle_id== 1)    charge=-1/3.;   // d
	else if (particle_id== 2)    charge= 2/3.;   // u
	else if (particle_id== 3)    charge=-1/3.;   // s
	else if (particle_id== 4)    charge= 2/3.;   // c
	else if (particle_id== 5)    charge=-1/3.;   // b
	else if (particle_id== 6)    charge= 2/3.;   // t
	else if (particle_id==-1)    charge= 1/3.;
	else if (particle_id==-2)    charge=-2/3.;
	else if (particle_id==-3)    charge= 1/3.;
	else if (particle_id==-4)    charge=-2/3.;
	else if (particle_id==-5)    charge= 1/3.;
	else if (particle_id==-6)    charge=-2/3.;
	//LEPTONS
	else if (particle_id== 11)   charge=-1;      // e-
	else if (particle_id==-11)   charge= 1;
	//LIGHT I = 1 MESONS
	else if (particle_id== 211)  charge= 1;      // π+
	else if (particle_id== 213)  charge= 1;      // ρ(770)+
	else if (particle_id==-211)  charge=-1;
	else if (particle_id==-213)  charge=-1;
	//STRANGE MESONS
	else if (particle_id== 321)  charge= 1;      // Κ+
	else if (particle_id== 323)  charge= 1;      // Κ*(892)+
	else if (particle_id==-321)  charge=-1;
	else if (particle_id==-323)  charge=-1;
	//CHARMED MESONS
	else if (particle_id== 411)  charge= 1;      // D+
	else if (particle_id== 413)  charge= 1;      // D*(2010)+
	else if (particle_id== 431)  charge= 1;      // Ds+
	else if (particle_id== 433)  charge= 1;      // Ds*+
	else if (particle_id==-411)  charge=-1;
	else if (particle_id==-413)  charge=-1;
	else if (particle_id==-431)  charge=-1;
	else if (particle_id==-433)  charge=-1;
	//BOTTOM MESONS
	else if (particle_id== 521)  charge= 1;      // B+
	else if (particle_id== 523)  charge= 1;      // B*+
	else if (particle_id==-521)  charge=-1;
	else if (particle_id==-523)  charge=-1;
	//LIGHT BARYONS
	else if (particle_id== 2212) charge= 1;      // p+
	else if (particle_id== 2224) charge= 2;      // Δ++
	else if (particle_id== 2214) charge= 1;      // Δ+
	else if (particle_id== 1114) charge=-1;      // Δ-
	else if (particle_id==-2212) charge=-1;
	else if (particle_id==-2224) charge=-2;
	else if (particle_id==-2214) charge=-1;
	else if (particle_id==-1114) charge=1;
	//STRANGE BARYONS
	else if (particle_id== 3222) charge= 1;      // Σ+
	else if (particle_id== 3112) charge=-1;      // Σ-
	else if (particle_id== 3224) charge= 1;      // Σ*+
	else if (particle_id== 3114) charge=-1;      // Σ*-
	else if (particle_id== 3312) charge=-1;      // Ξ-
	else if (particle_id== 3314) charge=-1;      // Ξ*-
	else if (particle_id== 3334) charge=-1;      // Ω-
	else if (particle_id==-3222) charge=-1;
	else if (particle_id==-3112) charge= 1;
	else if (particle_id==-3224) charge=-1;
	else if (particle_id==-3114) charge= 1;
	else if (particle_id==-3312) charge= 1;
	else if (particle_id==-3314) charge= 1;
	else if (particle_id==-3334) charge= 1;
	//CHARMED BARYONS
	else if (particle_id== 4122) charge= 1;      // Λc+
	else if (particle_id== 4222) charge= 2;      // Σc++
	else if (particle_id== 4212) charge= 1;      // Σc+
	else if (particle_id== 4224) charge= 2;      // Σc*++
	else if (particle_id== 4232) charge= 1;      // Ξc+
	else if (particle_id== 4322) charge= 1;      // Ξ'c+
	else if (particle_id== 4324) charge= 1;      // Ξc*+
	else if (particle_id==-4122) charge=-1;
	else if (particle_id==-4222) charge=-2;
	else if (particle_id==-4212) charge=-1;
	else if (particle_id==-4224) charge=-2;
	else if (particle_id==-4232) charge=-1;
	else if (particle_id==-4322) charge=-1;
	else if (particle_id==-4324) charge=-1;
	
	return charge;
}


//####################################################################################################
void WriteHistograms(const std::string OutputString)
{
    std::stringstream ss;
    ss << OutputString << ".root";
    std::string Output=ss.str();
    TFile *File = new TFile(Output.c_str(),"RECREATE");
	
	NumJet->Write();
	NumEvent->Write();
	pTDis->Write();
	EtaDis->Write();
    PhiDis->Write();
	SigmaDis->Write();
	for(int i=0; i<ntype; i++){
		Multi[i]->Write();
		for(int j=0; j<npTbin; j++){
			MultipT[i][j]->Write();
		}
	}
	
	File->Write();
	File->Close();
	delete File;
}


//####################################################################################################
void DeleteHistograms()
{	
	delete NumJet;
	delete NumEvent;
	delete pTDis;
    delete EtaDis;
    delete PhiDis;
	delete SigmaDis;
	for(int i=0; i<ntype; i++){
		delete Multi[i];
		for(int j=0; j<npTbin; j++){
			delete MultipT[i][j];
		}
	}

}

