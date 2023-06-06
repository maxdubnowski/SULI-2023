#define FlatTreeAnalyzer_cxx
#include "FlatTreeAnalyzer.h"

#include <TH1D.h>
#include <TFile.h>
#include <TString.h>
#include <TMath.h>
#include <TVector3.h>
#include <TLorentzVector.h>

#include <iomanip>
#include <sstream>
#include <iostream>
#include <vector>
#include <iterator>

using namespace std;

//Function to divide by the bin width and to get xsecs
void Reweight(TH1D* h);

//----------------------------------------//

void FlatTreeAnalyzer::Loop() {

	//----------------------------------------//	

	if (fChain == 0) return;
	Long64_t nentries = fChain->GetEntriesFast();
	Long64_t nbytes = 0, nb = 0;

	double Units = 1E38; // so that the extracted cross-section is in 10^{-38} cm^{2}
	double A = 40.; // so that we can have xsecs per nucleus

	int NNeut = 6; //Number of Neutrons

	int NInte = 6; // Interaction processes: All, QE, MEC, RES, DIS, COH
	std::vector<TString> InteractionLabels = {"","QE","MEC","RES","DIS","COH"};

	//----------------------------------------//	

        // Output file

	TString FileNameAndPath = "OutputFiles/FlatTreeAnalyzerOutput_"+fOutputFile+".root";
	TFile* file = new TFile(FileNameAndPath,"recreate");

	std::cout << std::endl << "------------------------------------------------" << std::endl << std::endl;
	std::cout << "File " << FileNameAndPath << " to be created" << std::endl << std::endl;
	
	//----------------------------------------//

	// Plot declaration

	TH1D* TrueMuonCosThetaPlot[NInte][NNeut];
	TH1D* TrueDeltaPtPlot[NInte][NNeut];
	TH1D* TrueNeutronMultiplicityPlot[NInte][NNeut];

	// Loop over the interaction processes

	for (int inte = 0; inte < NInte; inte++) {
	  

	  //--------------------------------------------------//
	  for (int neut =0; neut < NNeut; neut++){
	    TrueMuonCosThetaPlot[inte][neut] = new TH1D(InteractionLabels[inte]+"TrueMuonCosThetaPlot_Neutrons"+to_string(neut),";cos(#theta_{#mu})",10,-1.,1.);
	    TrueDeltaPtPlot[inte][neut] = new TH1D(InteractionLabels[inte]+"TrueDeltaPtPlot_Neutrons"+to_string(neut),";#delta p_{t}",20,0.,1.);
	    TrueNeutronMultiplicityPlot[inte][neut] = new TH1D(InteractionLabels[inte]+"TrueNeutronMultiplicityPlot_Neutrons"+to_string(neut),";Number of Neutrons",6,-0.5,5.5);
	  }

	  //--------------------------------------------------//

	} // End of the loop over the interaction processes							

	//----------------------------------------//

	// Counters

	int CounterEventsPassedSelection = 0;
	int CounterQEEventsPassedSelection = 0;
	int CounterMECEventsPassedSelection = 0;
	int CounterRESEventsPassedSelection = 0;
	int CounterDISEventsPassedSelection = 0;
	int CounterCOHEventsPassedSelection = 0;	

	//----------------------------------------//
	
	// Loop over the events

	for (Long64_t jentry=0; jentry<nentries;jentry++) {

	  //----------------------------------------//	
	
	  Long64_t ientry = LoadTree(jentry);
	  if (ientry < 0) break; nb = fChain->GetEntry(jentry); nbytes += nb;
	  if (jentry%1000 == 0) std::cout << jentry/1000 << " k " << std::setprecision(3) << double(jentry)/nentries*100. << " %"<< std::endl;

	  //----------------------------------------//	
		
	  double weight = fScaleFactor*Units*A*Weight;	

	  //----------------------------------------//	

	  // Signal definition

	  if (PDGLep != 13) { continue; } // make sure that we have only a muon in the final state
	  if (cc != 1) { continue; } // make sure that we have only CC interactions		

	  int ProtonTagging = 0, ChargedPionTagging = 0, NeutralPionTagging = 0, MuonTagging = 0, NeutronTagging = 0;
	  vector <int> ProtonID; ProtonID.clear();
	  vector <int> MuonID; MuonID.clear();
	  vector <int> NeutronID; NeutronID.clear();

	   

	  // Example selection with CC1p0pi (units in GeV/c)
	  // Loop over final state particles
	  
	  

	  for (int i = 0; i < nfsp; i++) {
	   
	    double pf = TMath::Sqrt( px[i]*px[i] + py[i]*py[i] + pz[i]*pz[i]);

	    if (pdg[i] == 2112 ){
	      NeutronTagging++;
	      NeutronID.push_back(i);
	      
	    }


	    if (pdg[i] == 13 && pf > 0.1) {
	      MuonTagging ++;
	      MuonID.push_back(i);
	    }


	    if (pdg[i] == 2212 && (pf > 0.3 || pf < 1.) ) {
	      ProtonTagging ++;
	      ProtonID.push_back(i);
	    }


	    if (fabs(pdg[i]) == 211 && pf > 0.07)  {
	      ChargedPionTagging ++;
	    }


	    if (pdg[i] == 111)  {
	      NeutralPionTagging ++;
	    }
	  } // End of the loop over the final state particles

	  // If the signal definition is not satisfied, continue

	  if ( ProtonTagging != 1 || ChargedPionTagging != 0 || NeutralPionTagging != 0 || MuonTagging !=1) { continue; }


	 
	  // -- Creating the Transverse Missing Momentum -- //
	  TVector3 proton_vector(px[ProtonID[0]], py[ProtonID[0]], pz[ProtonID[0]]);
	  TVector3 muon_vector(px[MuonID[0]], py[MuonID[0]], pz[MuonID[0]]);
	  TVector3 deltapt_vector(muon_vector.X()+proton_vector.X(),muon_vector.Y()+proton_vector.Y(),0 );
	
	  double transP = deltapt_vector.Mag();
	  //----------------------------------------//	

	  // https://arxiv.org/pdf/2106.15809.pdf




	  CounterEventsPassedSelection++;
	
	  // Classify the events based on the interaction type

	  int genie_mode = -1.;
	  if (TMath::Abs(Mode) == 1) { CounterQEEventsPassedSelection++; genie_mode = 1; } // QE
	  else if (TMath::Abs(Mode) == 2) { CounterMECEventsPassedSelection++; genie_mode = 2; } // MEC
	  else if (
		   TMath::Abs(Mode) == 11 || TMath::Abs(Mode) == 12 || TMath::Abs(Mode) == 13 ||
		   TMath::Abs(Mode) == 17 || TMath::Abs(Mode) == 22 || TMath::Abs(Mode) == 23
		   ) { CounterRESEventsPassedSelection++; genie_mode = 3; } // RES
	  else if (TMath::Abs(Mode) == 21 || TMath::Abs(Mode) == 26) { CounterDISEventsPassedSelection++; genie_mode = 4; } // DIS
	  else if (TMath::Abs(Mode) == 16) { CounterCOHEventsPassedSelection++; genie_mode = 5;} // COH
	  else { continue; }  

	  // Feb 8 2022: Only case that is not covered is 15 = diffractive

	  //----------------------------------------//

	  // filling in the histo regardless of interaction mode
	  /*
	    TrueMuonCosThetaPlot[0][neut]->Fill(CosLep,weight);
	    TrueDeltaPtPlot[0][neut]->Fill(transP,weight);
	    TrueNeutronMultiplicityPlot[0][neut]->Fill(NeutronTagging,weight);
	    //----------------------------------------//
	    
	    // filling in the histo based on the interaction mode
	    
	    TrueMuonCosThetaPlot[genie_mode][neut]->Fill(CosLep,weight);
	    TrueDeltaPtPlot[genie_mode][neut]->Fill(transP,weight);
	    TrueNeutronMultiplicityPlot[genie_mode][neut]->Fill(NeutronTagging,weight);
	    //----------------------------------------//
	    */

	  if (NeutronTagging ==0){
	   TrueMuonCosThetaPlot[0][NeutronTagging]->Fill(CosLep,weight);
	    TrueDeltaPtPlot[0][NeutronTagging]->Fill(transP,weight);
	    TrueNeutronMultiplicityPlot[0][NeutronTagging]->Fill(NeutronTagging,weight);
	    //----------------------------------------//
	    
	    // filling in the histo based on the interaction mode
	    
	    TrueMuonCosThetaPlot[genie_mode][NeutronTagging]->Fill(CosLep,weight);
	    TrueDeltaPtPlot[genie_mode][NeutronTagging]->Fill(transP,weight);
	    TrueNeutronMultiplicityPlot[genie_mode][NeutronTagging]->Fill(NeutronTagging,weight);
	  }

	  else if (NeutronTagging == 1){

	    TrueMuonCosThetaPlot[0][NeutronTagging]->Fill(CosLep,weight);
	    TrueDeltaPtPlot[0][NeutronTagging]->Fill(transP,weight);
	    TrueNeutronMultiplicityPlot[0][NeutronTagging]->Fill(NeutronTagging,weight);
	    //----------------------------------------//
	    
	    // filling in the histo based on the interaction mode
	    
	    TrueMuonCosThetaPlot[genie_mode][NeutronTagging]->Fill(CosLep,weight);
	    TrueDeltaPtPlot[genie_mode][NeutronTagging]->Fill(transP,weight);
	    TrueNeutronMultiplicityPlot[genie_mode][NeutronTagging]->Fill(NeutronTagging,weight);
	  }

	  else if (NeutronTagging == 2){

	    TrueMuonCosThetaPlot[0][NeutronTagging]->Fill(CosLep,weight);
	    TrueDeltaPtPlot[0][NeutronTagging]->Fill(transP,weight);
	    TrueNeutronMultiplicityPlot[0][NeutronTagging]->Fill(NeutronTagging,weight);
	    //----------------------------------------//
	    
	    // filling in the histo based on the interaction mode
	    
	    TrueMuonCosThetaPlot[genie_mode][NeutronTagging]->Fill(CosLep,weight);
	    TrueDeltaPtPlot[genie_mode][NeutronTagging]->Fill(transP,weight);
	    TrueNeutronMultiplicityPlot[genie_mode][NeutronTagging]->Fill(NeutronTagging,weight);
	  }

	  else if (NeutronTagging == 3){
	    TrueMuonCosThetaPlot[0][NeutronTagging]->Fill(CosLep,weight);
	    TrueDeltaPtPlot[0][NeutronTagging]->Fill(transP,weight);
	    TrueNeutronMultiplicityPlot[0][NeutronTagging]->Fill(NeutronTagging,weight);
	    //----------------------------------------//
	    
	    // filling in the histo based on the interaction mode
	    
	    TrueMuonCosThetaPlot[genie_mode][NeutronTagging]->Fill(CosLep,weight);
	    TrueDeltaPtPlot[genie_mode][NeutronTagging]->Fill(transP,weight);
	    TrueNeutronMultiplicityPlot[genie_mode][NeutronTagging]->Fill(NeutronTagging,weight);
	  }

	  else if (NeutronTagging == 4){
	    TrueMuonCosThetaPlot[0][NeutronTagging]->Fill(CosLep,weight);
	    TrueDeltaPtPlot[0][NeutronTagging]->Fill(transP,weight);
	    TrueNeutronMultiplicityPlot[0][NeutronTagging]->Fill(NeutronTagging,weight);
	    //----------------------------------------//
	    
	    // filling in the histo based on the interaction mode
	    
	    TrueMuonCosThetaPlot[genie_mode][NeutronTagging]->Fill(CosLep,weight);
	    TrueDeltaPtPlot[genie_mode][NeutronTagging]->Fill(transP,weight);
	    TrueNeutronMultiplicityPlot[genie_mode][NeutronTagging]->Fill(NeutronTagging,weight);
	  }
	  
	  else if (NeutronTagging >=5){
	    TrueMuonCosThetaPlot[0][NeutronTagging]->Fill(CosLep,weight);
	    TrueDeltaPtPlot[0][NeutronTagging]->Fill(transP,weight);
	    TrueNeutronMultiplicityPlot[0][NeutronTagging]->Fill(NeutronTagging,weight);
	    //----------------------------------------//
	    
	    // filling in the histo based on the interaction mode
	    
	    TrueMuonCosThetaPlot[genie_mode][NeutronTagging]->Fill(CosLep,weight);
	    TrueDeltaPtPlot[genie_mode][NeutronTagging]->Fill(transP,weight);
	    TrueNeutronMultiplicityPlot[genie_mode][NeutronTagging]->Fill(NeutronTagging,weight);
	  }

	} // End of the loop over the events

	//----------------------------------------//	

	std::cout << "Percetage of events passing the selection cuts = " << 
	double(CounterEventsPassedSelection)/ double(nentries)*100. << " %" << std::endl; std::cout << std::endl;

	std::cout << "Success percetage in selecting QE events = " << 
	double(CounterQEEventsPassedSelection)/ double(CounterEventsPassedSelection)*100. << " %" << std::endl; std::cout << std::endl;

	std::cout << "Success percetage in selecting MEC events = " << 
	double(CounterMECEventsPassedSelection)/ double(CounterEventsPassedSelection)*100. << " %" << std::endl; std::cout << std::endl;

	std::cout << "Success percetage in selecting RES events = " << 
	double(CounterRESEventsPassedSelection)/ double(CounterEventsPassedSelection)*100. << " %" << std::endl; std::cout << std::endl;

	std::cout << "Success percetage in selecting DIS events = " << 
	double(CounterDISEventsPassedSelection)/ double(CounterEventsPassedSelection)*100. << " %" << std::endl; std::cout << std::endl;

	std::cout << "Success percetage in selecting COH events = " << 
	double(CounterCOHEventsPassedSelection)/ double(CounterEventsPassedSelection)*100. << " %" << std::endl; std::cout << std::endl;	

	//----------------------------------------//	
	//----------------------------------------//	

	// Division by bin width to get the cross sections	
	// Loop over the interaction processes

	for (int inte = 0; inte < NInte; inte++) {
	  for (int neut =0; neut < NNeut; neut++){
		//----------------------------------------//
	
		Reweight(TrueMuonCosThetaPlot[inte][neut]);
		Reweight(TrueDeltaPtPlot[inte][neut]);
		Reweight(TrueNeutronMultiplicityPlot[inte][neut]);
		//----------------------------------------//
	  }
	} // End of the loop over the interaction processes		

	//----------------------------------------//		
		
	file->cd();
	file->Write();
	fFile->Close();

	std::cout << std::endl;
	std::cout << "File " << FileNameAndPath +" has been created created " << std::endl; 
	std::cout << std::endl;

	std::cout << std::endl << "------------------------------------------------" << std::endl << std::endl;

	//----------------------------------------//		

} // End of the program

//----------------------------------------//		

void Reweight(TH1D* h) {

  int NBins = h->GetXaxis()->GetNbins();

  for (int i = 0; i < NBins; i++) {

    double CurrentEntry = h->GetBinContent(i+1);
    double NewEntry = CurrentEntry / h->GetBinWidth(i+1);

    double CurrentError = h->GetBinError(i+1);
    double NewError = CurrentError / h->GetBinWidth(i+1);

    h->SetBinContent(i+1,NewEntry); 
    h->SetBinError(i+1,NewError); 
    //h->SetBinError(i+1,0.000001); 

  }

}

//----------------------------------------//		
