{
	gROOT->ProcessLine(".L ../myClasses/STV_Tools.cxx+");	

	gROOT->ProcessLine(".L myNEUTAnalysis.C+");
	gROOT->ProcessLine("myNEUTAnalysis().Loop()");
	//gROOT->ProcessLine(".q");
};
