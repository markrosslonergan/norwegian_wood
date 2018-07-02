void mark_run2(){

	std::cout<<"Loading genie_CC and plotfunctions"<<std::endl;
	gROOT->ProcessLine(".L plotfunctions.c");
	gROOT->ProcessLine(".L genie_mu.c");

	gROOT->ProcessLine("run_all_genie_mu()");
	gROOT->ProcessLine("stack_mu()");


}
