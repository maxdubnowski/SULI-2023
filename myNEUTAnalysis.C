#define myNEUTAnalysis_cxx
#include "myNEUTAnalysis.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <TH1D.h>
#include <TFile.h>
#include <iomanip>
#include <TString.h>
#include <TMath.h>
#include <TVector3.h>
#include <TLorentzVector.h>
#include <TRandom.h>

#include <sstream>
#include <iostream>
#include <vector>
#include <iterator>

#include "../myClasses/Tools.h"
#include "../myClasses/STV_Tools.h"
#include "../../../../../Secondary_Code/mySimFunctions.cpp"

using namespace std;

void myNEUTAnalysis::Loop() {

	if (fChain == 0) return;
	Long64_t nentries = fChain->GetEntriesFast();
	Long64_t nbytes = 0, nb = 0;
	
	// ---------------------------------------------------------------------------------------------------------------------------------	

	TString FileNameAndPath = "OutputFiles/STVAnalysis_"+GeneratorName+".root";
	TFile* file = new TFile(FileNameAndPath,"recreate");
	std::cout << std::endl << std::endl;
	
	// ---------------------------------------------------------------------------------------------------------------------------------	

	TH1D* TrueMuonMomentumPlot = new TH1D("TrueMuonMomentumPlot",LabelXAxisMuonMomentum,NBinsMuonMomentum,ArrayNBinsMuonMomentum);
	TH1D* TrueMuonPhiPlot = new TH1D("TrueMuonPhiPlot",LabelXAxisMuonPhi,NBinsMuonPhi,ArrayNBinsMuonPhi);
	TH1D* TrueMuonCosThetaPlot = new TH1D("TrueMuonCosThetaPlot",LabelXAxisMuonCosTheta,NBinsMuonCosTheta,ArrayNBinsMuonCosTheta);
	TH1D* TrueMuonCosThetaSingleBinPlot = new TH1D("TrueMuonCosThetaSingleBinPlot",LabelXAxisMuonCosTheta,1,-1.,1.);

	TH1D* TrueProtonMomentumPlot = new TH1D("TrueProtonMomentumPlot",LabelXAxisProtonMomentum,NBinsProtonMomentum,ArrayNBinsProtonMomentum);
	TH1D* TrueProtonPhiPlot = new TH1D("TrueProtonPhiPlot",LabelXAxisProtonPhi,NBinsProtonPhi,ArrayNBinsProtonPhi);
	TH1D* TrueProtonCosThetaPlot = new TH1D("TrueProtonCosThetaPlot",LabelXAxisProtonCosTheta,NBinsProtonCosTheta,ArrayNBinsProtonCosTheta);

	TH1D* TrueECalPlot = new TH1D("TrueECalPlot",LabelXAxisECal,NBinsECal,ArrayNBinsECal);
	TH1D* TrueECalLowPTPlot = new TH1D("TrueECalLowPTPlot",LabelXAxisECal,NBinsECal,ArrayNBinsECal);
	TH1D* TrueECalMidPTPlot = new TH1D("TrueECalMidPTPlot",LabelXAxisECal,NBinsECal,ArrayNBinsECal);
	TH1D* TrueECalHighPTPlot = new TH1D("TrueECalHighPTPlot",LabelXAxisECal,NBinsECal,ArrayNBinsECal);
	TH1D* TrueEQEPlot = new TH1D("TrueEQEPlot",LabelXAxisEQE,NBinsEQE,ArrayNBinsEQE);	
	TH1D* TrueQ2Plot = new TH1D("TrueQ2Plot",LabelXAxisQ2,NBinsQ2,ArrayNBinsQ2);
	
	TH1D* TrueDeltaPTPlot = new TH1D("TrueDeltaPTPlot",LabelXAxisDeltaPT,NBinsDeltaPT,ArrayNBinsDeltaPT);
	TH1D* TrueDeltaAlphaTPlot = new TH1D("TrueDeltaAlphaTPlot",LabelXAxisDeltaAlphaT,NBinsDeltaAlphaT,ArrayNBinsDeltaAlphaT);
	TH1D* TrueDeltaPhiTPlot = new TH1D("TrueDeltaPhiTPlot",LabelXAxisDeltaPhiT,NBinsDeltaPhiT,ArrayNBinsDeltaPhiT);	

	TH1D* TrueCCQEMuonMomentumPlot = new TH1D("TrueCCQEMuonMomentumPlot",LabelXAxisMuonMomentum,CCQENBinsMuonMomentum,CCQEArrayNBinsMuonMomentum);
	TH1D* TrueCCQEMuonPhiPlot = new TH1D("TrueCCQEMuonPhiPlot",LabelXAxisMuonPhi,CCQENBinsMuonPhi,CCQEArrayNBinsMuonPhi);
	TH1D* TrueCCQEMuonCosThetaPlot = new TH1D("TrueCCQEMuonCosThetaPlot",LabelXAxisMuonCosTheta,CCQENBinsMuonCosTheta,CCQEArrayNBinsMuonCosTheta);	

	TH1D* TrueCCQEProtonMomentumPlot = new TH1D("TrueCCQEProtonMomentumPlot",LabelXAxisProtonMomentum,CCQENBinsProtonMomentum,CCQEArrayNBinsProtonMomentum);
	TH1D* TrueCCQEProtonPhiPlot = new TH1D("TrueCCQEProtonPhiPlot",LabelXAxisProtonPhi,CCQENBinsProtonPhi,CCQEArrayNBinsProtonPhi);
	TH1D* TrueCCQEProtonCosThetaPlot = new TH1D("TrueCCQEProtonCosThetaPlot",LabelXAxisProtonCosTheta,CCQENBinsProtonCosTheta,CCQEArrayNBinsProtonCosTheta);

	TH1D* TrueCCQEECalPlot = new TH1D("TrueCCQEECalPlot",LabelXAxisECal,CCQENBinsECal,CCQEArrayNBinsECal);
	TH1D* TrueCCQEQ2Plot = new TH1D("TrueCCQEQ2Plot",LabelXAxisQ2,CCQENBinsQ2,CCQEArrayNBinsQ2);	

	TH1D* TruekMissPlot = new TH1D("TruekMissPlot",LabelXAxiskMiss,NBinskMiss,ArrayNBinskMiss);
	TH1D* TruePMissMinusPlot = new TH1D("TruePMissMinusPlot",LabelXAxisPMissMinus,NBinsPMissMinus,ArrayNBinsPMissMinus);
	TH1D* TruePMissPlot = new TH1D("TruePMissPlot",LabelXAxisPMiss,NBinsPMiss,ArrayNBinsPMiss);

	TH1D* TrueDeltaPLPlot = new TH1D("TrueDeltaPLPlot",LabelXAxisDeltaPL,NBinsDeltaPL,ArrayNBinsDeltaPL);
	TH1D* TrueDeltaPnPlot = new TH1D("TrueDeltaPnPlot",LabelXAxisDeltaPn,NBinsDeltaPn,ArrayNBinsDeltaPn);
	TH1D* TrueDeltaPtxPlot = new TH1D("TrueDeltaPtxPlot",LabelXAxisDeltaPtx,NBinsDeltaPtx,ArrayNBinsDeltaPtx);
	TH1D* TrueDeltaPtyPlot = new TH1D("TrueDeltaPtyPlot",LabelXAxisDeltaPty,NBinsDeltaPty,ArrayNBinsDeltaPty);
	TH1D* TrueAPlot = new TH1D("TrueAPlot",LabelXAxisA,NBinsA,ArrayNBinsA);

	// For now and until box opening
	int NBins2DAnalysis = 4;
		
	TH2D* TrueCosThetaMuPmuPlot = new TH2D("TrueCosThetaMuPmuPlot",LabelXAxisMuonCosTheta+LabelXAxisMuonMomentum
			,NBins2DAnalysis,ArrayNBinsMuonCosTheta[0],ArrayNBinsMuonCosTheta[NBinsMuonCosTheta]
			,NBins2DAnalysis,ArrayNBinsMuonMomentum[0],ArrayNBinsMuonMomentum[NBinsMuonMomentum]);
			
	TH2D* TrueCosThetaPPpPlot = new TH2D("TrueCosThetaPPpPlot",LabelXAxisProtonCosTheta+LabelXAxisProtonMomentum
			,NBins2DAnalysis,ArrayNBinsProtonCosTheta[0],ArrayNBinsProtonCosTheta[NBinsProtonCosTheta]
			,NBins2DAnalysis,ArrayNBinsProtonMomentum[0],ArrayNBinsProtonMomentum[NBinsProtonMomentum]);

	// ---------------------------------------------------------------------------------------------------------------------------------

	// Counters STVlike analysis

	int CounterSTVlikeEventsPassedSelection = 0;
	int CounterSTVlikeQEEventsPassedSelection = 0;
	int CounterSTVlikeMECEventsPassedSelection = 0;
	int CounterSTVlikeRESEventsPassedSelection = 0;
	int CounterSTVlikeDISEventsPassedSelection = 0;

	// ----------------------------------------------------------------------------------------------------------------------------------------

	Tools tools;

	// ---------------------------------------------------------------------------------------------------------------------------------	
	
	for (Long64_t jentry=0; jentry<nentries;jentry++) {
	
		Long64_t ientry = LoadTree(jentry);
		if (ientry < 0) break; nb = fChain->GetEntry(jentry); nbytes += nb;
		if (jentry%1000 == 0) std::cout << jentry/1000 << " k " << std::setprecision(3) << double(jentry)/nentries*100. << " %"<< std::endl;

		// ---------------------------------------------------------------------------------------------------------------------------------
		
//		double weight = 1.;	
		double A = 40.;
		double weight = fScaleFactor*Units*A*Weight;	

		// ---------------------------------------------------------------------------------------------------------------------------------

		if (PDGLep != 13) { continue; } // make sure that we have only CC interactions

		// ---------------------------------------------------------------------------------------------------------------------------------

		int ProtonTagging = 0, ChargedPionTagging = 0, NeutralPionTagging = 0, ElectronTagging = 0, GammaTagging = 0, MuonTagging = 0, TrueHeavierMesonCounter = 0;
		vector <int> ProtonID; ProtonID.clear();
		vector <int> MuonID; MuonID.clear();		

		for (int i = 0; i < nfsp; i++) {
		
			double pf = TMath::Sqrt( px[i]*px[i] + py[i]*py[i] + pz[i]*pz[i]);

			if (pdg[i] == MuonPdg && pf > ArrayNBinsMuonMomentum[0]) {

				MuonTagging ++;
				MuonID.push_back(i);

			}

			if (pdg[i] == ProtonPdg && pf > ArrayNBinsProtonMomentum[0]) {

				ProtonTagging ++;
				ProtonID.push_back(i);

			}

			if (fabs(pdg[i]) == AbsChargedPionPdg && pf > ChargedPionMomentumThres)  {

				ChargedPionTagging ++;

			}

			if (pdg[i] == NeutralPionPdg)  {

				NeutralPionTagging ++;

			}

			if ( pdg[i] != NeutralPionPdg && fabs(pdg[i]) != AbsChargedPionPdg && tools.is_meson_or_antimeson(pdg[i]) ) { TrueHeavierMesonCounter++; }

		} // End of the loop over the final state particles

		if ( ProtonTagging != 1 || ChargedPionTagging != 0 || NeutralPionTagging != 0 || MuonTagging !=1) { continue; }
		if ( TrueHeavierMesonCounter != 0 ) { continue; }

		// ----------------------------------------------------------------------------------------------------------------------------------

		// Muon

		TVector3 Muon3Vector(px[MuonID.at(0)],py[MuonID.at(0)],pz[MuonID.at(0)]); // GeV
		TLorentzVector Muon4Vector(px[MuonID.at(0)],py[MuonID.at(0)],pz[MuonID.at(0)],ELep); // GeV
		double MuonMomentum = Muon3Vector.Mag(); // GeV / c
		double MuonCosTheta = CosLep;
		double MuonPhi = Muon4Vector.Phi() * 180. / TMath::Pi(); // deg

		// -------------------------------------------------------------------------------------------------------------------------------

		// Proton

		int ProtonIndex = ProtonID.at(0);
		TVector3 Proton3Vector(px[ProtonIndex],py[ProtonIndex],pz[ProtonIndex]); // GeV
		TLorentzVector Proton4Vector(px[ProtonIndex],py[ProtonIndex],pz[ProtonIndex],E[ProtonIndex]); // GeV
		double ProtonMomentum = Proton3Vector.Mag(); // GeV / c

		double ProtonCosTheta = Proton3Vector.CosTheta();
		double ProtonPhi = Proton4Vector.Phi() * 180. / TMath::Pi(); // deg
		double ProtonEnergy = Proton4Vector.E(); // GeV
		double ProtonKE = ProtonEnergy - ProtonMass_GeV; // GeV

		// ---------------------------------------------------------------------------------------------------------------------------------

		// Relative Proton - Muon Angles

		double DeltaThetaProtonMuon = Proton3Vector.Angle(Muon3Vector) * 180. / TMath::Pi();
		if (DeltaThetaProtonMuon < 0.) { DeltaThetaProtonMuon += 180.; }
		if (DeltaThetaProtonMuon > 180.) { DeltaThetaProtonMuon -= 180.; }

		double DeltaPhiProtonMuon = Muon3Vector.DeltaPhi(Proton3Vector)* 180. / TMath::Pi();
		if (DeltaPhiProtonMuon < 0.) { DeltaPhiProtonMuon += 360.; }
		if (DeltaPhiProtonMuon > 360.) { DeltaPhiProtonMuon -= 360.; }

		// ---------------------------------------------------------------------------------------------------------------------------------

		STV_Tools stv_tool(Muon3Vector,Proton3Vector,ELep,ProtonEnergy);

		double PTmissMomentum = stv_tool.ReturnPt();
		double TrueDeltaAlphaT = stv_tool.ReturnDeltaAlphaT();
		double TrueDeltaPhiT = stv_tool.ReturnDeltaPhiT();
		double ECal = stv_tool.ReturnECal();
		double EQE = stv_tool.ReturnEQE();
		double TrueQ2 = stv_tool.ReturnQ2();	

		double TruekMiss = stv_tool.ReturnkMiss();
		double TruePMissMinus = stv_tool.ReturnPMissMinus();
		double TrueMissMomentum = stv_tool.ReturnPMiss();

		double TruePL = stv_tool.ReturnPL();
		double TruePn = stv_tool.ReturnPn();
		double TruePtx = stv_tool.ReturnPtx();
		double TruePty = stv_tool.ReturnPty();
		double TrueA = stv_tool.ReturnA();

		// ---------------------------------------------------------------------------------------------------------------------------------

		if (MuonMomentum > ArrayNBinsMuonMomentum[0]
		    && ProtonMomentum > ArrayNBinsProtonMomentum[0]
		) {

			if (		
			    // Same evenst fill all the STVLike plots 
			    
			    PTmissMomentum > ArrayNBinsDeltaPT[0] //&& PTmissMomentum < ArrayNBinsDeltaPT[NBinsDeltaPT]
			    && TrueDeltaAlphaT > ArrayNBinsDeltaAlphaT[0] && TrueDeltaAlphaT < ArrayNBinsDeltaAlphaT[NBinsDeltaAlphaT]
		 	    && TrueDeltaPhiT > ArrayNBinsDeltaPhiT[0] && TrueDeltaPhiT < ArrayNBinsDeltaPhiT[NBinsDeltaPhiT]			    
			     
			    && MuonMomentum < ArrayNBinsMuonMomentum[NBinsMuonMomentum]  
			    && ProtonMomentum < ArrayNBinsProtonMomentum[NBinsProtonMomentum]
			    && MuonCosTheta > ArrayNBinsMuonCosTheta[0]
			    && MuonCosTheta < ArrayNBinsMuonCosTheta[NBinsMuonCosTheta]
			    && ProtonCosTheta > ArrayNBinsProtonCosTheta[0]
			    && ProtonCosTheta < ArrayNBinsProtonCosTheta[NBinsProtonCosTheta]
			    
//			    && ECal > ArrayNBinsECal[0] && ECal < ArrayNBinsECal[NBinsECal]
//			    && EQE > ArrayNBinsEQE[0] && EQE < ArrayNBinsEQE[NBinsEQE]
//			    && TrueQ2 > ArrayNBinsQ2[0] && TrueQ2 < ArrayNBinsQ2[NBinsQ2]			    

			) {
			
				TrueDeltaPTPlot->Fill(PTmissMomentum,weight);
				TrueDeltaAlphaTPlot->Fill(TrueDeltaAlphaT,weight);
				TrueDeltaPhiTPlot->Fill(TrueDeltaPhiT,weight);			

				TrueMuonMomentumPlot->Fill(MuonMomentum,weight);
				TrueMuonPhiPlot->Fill(MuonPhi,weight);
				TrueMuonCosThetaPlot->Fill(MuonCosTheta,weight);
				TrueMuonCosThetaSingleBinPlot->Fill(MuonCosTheta,weight);

				TrueProtonMomentumPlot->Fill(ProtonMomentum,weight);
				TrueProtonPhiPlot->Fill(ProtonPhi,weight);
				TrueProtonCosThetaPlot->Fill(ProtonCosTheta,weight);

				TrueECalPlot->Fill(ECal,weight);
				if (PTmissMomentum > LowPT[0] && PTmissMomentum < HighPT[0]) { TrueECalLowPTPlot->Fill(ECal,weight); }
				if (PTmissMomentum > LowPT[1] && PTmissMomentum < HighPT[1]) { TrueECalMidPTPlot->Fill(ECal,weight); }
				if (PTmissMomentum > LowPT[2] && PTmissMomentum < HighPT[2]) { TrueECalHighPTPlot->Fill(ECal,weight); }
				TrueEQEPlot->Fill(EQE,weight);				
				TrueQ2Plot->Fill(TrueQ2,weight);

				TruekMissPlot->Fill(TruekMiss,weight);
				TruePMissMinusPlot->Fill(TruePMissMinus,weight);
				TruePMissPlot->Fill(TrueMissMomentum,weight);

				TrueDeltaPLPlot->Fill(TruePL,weight);
				TrueDeltaPnPlot->Fill(TruePn,weight);
				TrueDeltaPtxPlot->Fill(TruePtx,weight);
				TrueDeltaPtyPlot->Fill(TruePty,weight);
				TrueAPlot->Fill(TrueA,weight);

				// 2D Analysis
		
				TrueCosThetaMuPmuPlot->Fill(MuonCosTheta,MuonMomentum,weight);
				TrueCosThetaPPpPlot->Fill(ProtonCosTheta,ProtonMomentum,weight);

				CounterSTVlikeEventsPassedSelection++;
				if (Mode == 1) { CounterSTVlikeQEEventsPassedSelection++; } // QE
				if (Mode == 2) { CounterSTVlikeMECEventsPassedSelection++; } // MEC
				if (Mode == 11 || Mode == 12 || Mode == 13) { CounterSTVlikeRESEventsPassedSelection++; } // RES
				if (Mode == 26) { CounterSTVlikeDISEventsPassedSelection++; } // DIS

				if (		
					// CCQElike measurement
						
					PTmissMomentum < 0.35
					&& TMath::Abs(DeltaPhiProtonMuon - 180) < 35
					&& TMath::Abs(DeltaThetaProtonMuon - 90) < 55

					&& MuonMomentum > CCQEArrayNBinsMuonMomentum[0]  
					&& ProtonMomentum > CCQEArrayNBinsProtonMomentum[0]					
					&& MuonMomentum < CCQEArrayNBinsMuonMomentum[CCQENBinsMuonMomentum]  
					&& ProtonMomentum < CCQEArrayNBinsProtonMomentum[CCQENBinsProtonMomentum]
					
					&& MuonCosTheta > CCQEArrayNBinsMuonCosTheta[0]
					&& MuonCosTheta < CCQEArrayNBinsMuonCosTheta[CCQENBinsMuonCosTheta]
					&& ProtonCosTheta > CCQEArrayNBinsProtonCosTheta[0]
					&& ProtonCosTheta < CCQEArrayNBinsProtonCosTheta[CCQENBinsProtonCosTheta]

				) {

					TrueCCQEMuonMomentumPlot->Fill(MuonMomentum,weight);
					TrueCCQEMuonPhiPlot->Fill(MuonPhi,weight);
					TrueCCQEMuonCosThetaPlot->Fill(MuonCosTheta,weight);

					TrueCCQEProtonMomentumPlot->Fill(ProtonMomentum,weight);
					TrueCCQEProtonPhiPlot->Fill(ProtonPhi,weight);
					TrueCCQEProtonCosThetaPlot->Fill(ProtonCosTheta,weight);	

					TrueCCQEECalPlot->Fill(ECal,weight);
					TrueCCQEQ2Plot->Fill(TrueQ2,weight);								

				}				

			}
			
		} // End of the angle selection cuts && the demand that we fill the plots with the same events

		// ----------------------------------------------------------------------------------------------------------------------------------	

	} // End of the loop over the events
	
	// ----------------------------------------------------------------------------------------------------------------------------------------

	// STVlike analysis
	
	std::cout << std::endl << "----------------------- STVlike analysis -------------------------" << std::endl << std::endl;

	std::cout << "Percetage of events passing the selection cuts = " << 
	double(CounterSTVlikeEventsPassedSelection)/ double(nentries)*100. << " %" << std::endl; std::cout << std::endl;

	std::cout << "Success percetage in selecting QE events = " << 
	double(CounterSTVlikeQEEventsPassedSelection)/ double(CounterSTVlikeEventsPassedSelection)*100. << " %" << std::endl; std::cout << std::endl;

	std::cout << "Success percetage in selecting MEC events = " << 
	double(CounterSTVlikeMECEventsPassedSelection)/ double(CounterSTVlikeEventsPassedSelection)*100. << " %" << std::endl; std::cout << std::endl;

	std::cout << "Success percetage in selecting RES events = " << 
	double(CounterSTVlikeRESEventsPassedSelection)/ double(CounterSTVlikeEventsPassedSelection)*100. << " %" << std::endl; std::cout << std::endl;

	std::cout << "Success percetage in selecting DIS events = " << 
	double(CounterSTVlikeDISEventsPassedSelection)/ double(CounterSTVlikeEventsPassedSelection)*100. << " %" << std::endl; std::cout << std::endl;

	// ---------------------------------------------------------------------------------------------------------------------------------------
	// -------------------------------------------------------------------------------------------------------------------------------------------

	double ScalingFactor = 1.;

	// Division by bin width to get the cross sections
	
	Reweight(TrueMuonMomentumPlot,ScalingFactor);
	Reweight(TrueMuonPhiPlot,ScalingFactor);
	Reweight(TrueMuonCosThetaPlot,ScalingFactor);
	Reweight(TrueMuonCosThetaSingleBinPlot,ScalingFactor);

	Reweight(TrueProtonMomentumPlot,ScalingFactor);
	Reweight(TrueProtonPhiPlot,ScalingFactor);
	Reweight(TrueProtonCosThetaPlot,ScalingFactor);

	Reweight(TrueECalPlot,ScalingFactor);
	Reweight(TrueEQEPlot,ScalingFactor);	
	Reweight(TrueQ2Plot,ScalingFactor);

	Reweight(TrueCCQEMuonMomentumPlot,ScalingFactor);
	Reweight(TrueCCQEMuonPhiPlot,ScalingFactor);
	Reweight(TrueCCQEMuonCosThetaPlot,ScalingFactor);

	Reweight(TrueCCQEProtonMomentumPlot,ScalingFactor);
	Reweight(TrueCCQEProtonPhiPlot,ScalingFactor);
	Reweight(TrueCCQEProtonCosThetaPlot,ScalingFactor);

	Reweight(TrueCCQEECalPlot,ScalingFactor);
	Reweight(TrueCCQEQ2Plot,ScalingFactor);	
	
	Reweight(TrueDeltaPTPlot,ScalingFactor);
	Reweight(TrueDeltaAlphaTPlot,ScalingFactor);
	Reweight(TrueDeltaPhiTPlot,ScalingFactor);

	Reweight(TruekMissPlot,ScalingFactor);
	Reweight(TruePMissPlot,ScalingFactor);
	Reweight(TruePMissMinusPlot,ScalingFactor);

	Reweight(TrueDeltaPLPlot,ScalingFactor);
	Reweight(TrueDeltaPnPlot,ScalingFactor);
	Reweight(TrueDeltaPtxPlot,ScalingFactor);
	Reweight(TrueDeltaPtyPlot,ScalingFactor);
	Reweight(TrueAPlot,ScalingFactor);

	Reweight2D(TrueCosThetaMuPmuPlot,ScalingFactor);
	Reweight2D(TrueCosThetaPPpPlot,ScalingFactor);
	
	// --------------------------------------------------------------------------------------------------------------------------------------------	
	
	file->cd();
	file->Write();

	std::cout << std::endl;
	std::cout << "File " << FileNameAndPath +" has been created created " << std::endl; 
	std::cout << std::endl;

} // End of the program
