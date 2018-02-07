void mark_run(){

	std::cout<<"Loading genie_CC and plotfunctions"<<std::endl;
	gROOT->ProcessLine(".L genie_CC.c");
	gROOT->ProcessLine(".L plotfunctions.c");

	gROOT->ProcessLine("run_all_genie_study()");
	gROOT->ProcessLine("stack_study()");


}
