{

	vector<TString> WhichSample; vector<TString> WhichName;

	//----------------------------------------//

	WhichSample.push_back("GENIE_v3_0_6_G18_10a_02a.flat"); WhichName.push_back("GENIE");
	WhichSample.push_back("NEUT.flat"); WhichName.push_back("NEUT.flat");
	WhichSample.push_back("NuWro.flat"); WhichName.push_back("NuWro");
	WhichSample.push_back("GiBUU.flat"); WhichName.push_back("GiBUU");				

	//----------------------------------------//

	gROOT->ProcessLine(".L FlatTreeAnalyzer.cxx+");

	for (int i =0;i < (int)(WhichSample.size()); i++) {

		gROOT->ProcessLine("FlatTreeAnalyzer(\""+WhichSample[i]+"\",\""+WhichName[i]+"\").Loop()");

	}
	//gROOT->ProcessLine(".q");
};
