{

	vector<TString> WhichSample; vector<TString> WhichName;

	//----------------------------------------//

	WhichSample.push_back("NEUT"); WhichName.push_back("NEUT");
	WhichSample.push_back("bnb.ub.num.neut_5_4_0_1.flat"); WhichName.push_back("NEUTv5401_LFG");
	WhichSample.push_back("bnb.ub.num.neut_5_4_0_1_RFG.flat"); WhichName.push_back("NEUTv5401_RFG");	
//	WhichSample.push_back("bnb.ub.num.neut_5_4_0_1_EffSF.flat"); WhichName.push_back("NEUTv5401_EffSF");
	WhichSample.push_back("NuWroCard_CC_Ar_uBFlux_flat"); WhichName.push_back("NuWrov190201_LFG");
	WhichSample.push_back("NuWroCard_CC_Ar_uBFlux_RFG_flat"); WhichName.push_back("NuWrov190201_RFG");		

	//----------------------------------------//

	gROOT->ProcessLine(".L ../myClasses/Util.C+");
	gROOT->ProcessLine(".L ../myClasses/STV_Tools.cxx+");
	gROOT->ProcessLine(".L ../myClasses/Tools.cxx+");	

	gROOT->ProcessLine(".L myNEUTAnalysis.cxx+");

	for (int i =0;i < (int)(WhichSample.size()); i++) {

		gROOT->ProcessLine("myNEUTAnalysis(\""+WhichSample[i]+"\",\""+WhichName[i]+"\").Loop()");

	}
	//gROOT->ProcessLine(".q");
};
