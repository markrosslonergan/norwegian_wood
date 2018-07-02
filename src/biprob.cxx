#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <cmath>
#include <vector>
#include <unistd.h>
#include <getopt.h>
#include <cstring>

#include "TFile.h"
#include "TTree.h"
#include "TH1F.h"
#include "TString.h"
#include "TNtuple.h"
#include "TChain.h"
#include "TMath.h"
#include "TSystem.h"
#include "TMatrixT.h"
#include "TRandom.h"
#include "TError.h"
#include "TCanvas.h"
#include "TH2F.h"
#include "TGraph.h"

#include "params.h"
#include "SBNconfig.h"
#include "SBNchi.h"
#include "SBNdet.h"
#include "SBNspec.h"
#include "SBNosc.h"
#include "SBgeN.h"
#include "SBNfit.h"
#include "SBNfit3pN.h"
#include "prob.h"
#include "SBNprob.h"

#define no_argument 0
#define required_argument 1
#define optional_argument 2

using namespace sbn;









/*************************************************************
 *************************************************************
 *		BEGIN example.cxx
 ************************************************************
 ************************************************************/
int main(int argc, char* argv[])
{

	double T34=0;
	double T24=0;

	std::string xml = "default.xml";
	int iarg = 0;
	opterr=1;
	int index; 
	int test_mode=0;
	std::string filename = "default.root";
	/*************************************************************
	 *************************************************************
	 *		Command Line Argument Reading
	 ************************************************************
	 ************************************************************/

	const struct option longopts[] = 
	{
		{"xml", 		required_argument, 	0, 'x'},
		{"test",		required_argument,	0, 't'},
		{"file",		required_argument,	0, 'f'},
		{0,			no_argument, 		0,  0},
	};


	while(iarg != -1)
	{
		iarg = getopt_long(argc,argv, "x:t:f:", longopts, &index);

		switch(iarg)
		{
			case 'x':
				xml = optarg;
				break;
			case 'f':
				filename = optarg;//`strtof(optarg,NULL);
				break;

			case 't':
				T34 = strtod(optarg,NULL);
				break;
			case '?':
			case 'h':
				std::cout<<"Allowed arguments:"<<std::endl;
				std::cout<<"\t-x\t--xml\t\tInput .xml file for SBNconfig"<<std::endl;
				return 0;
		}

	}



	std::ofstream dunestream,icarusstream,sbndstream;
	dunestream.open ("DUNE_bi_prob.dat");

//	std::vector<double> angles = {33.6, 42.1, 8.5, 0,0,0};
//	std::vector<double> angles_oct = {33.6, 49.6, 8.5, 0,0,0};

	std::vector<double> angles = {33.6, 42.1, 8.5, 14,3,0};
	std::vector<double> angles_oct = {33.6, 49.6, 8.5, 14,3,0};
	std::vector<double> phases = {0,0,0};


	std::vector<double> phases_180 = {180,0,0};
	//http://lbne2-docdb.fnal.gov/cgi-bin/RetrieveFile?docid=10688&filename=DUNE-CDR-physics-volume.pdf&version=10
	std::vector<double> mass_splittings = {7.5*pow(10,-5), 2.457*pow(10,-3),2};
	std::vector<double> mass_splittings_inv = {7.5*pow(10,-5), -2.449*pow(10,-3),2};


//plot with
//plot 'DUNE_bi_prob.dat' u 2:3 w l , 'DUNE_bi_prob.dat' u 4:5 w l, 'DUNE_bi_prob.dat' u 6:7 w l, 'DUNE_bi_prob.dat' u 8:9 w l
	
	for(double d14 = 0; d14 < 360; d14 +=7.5){
	for(double dcp = 0; dcp < 360; dcp +=7.5){
		dunestream<<dcp<<" ";
		for(int h =0; h<2; h++){
			for(int o =0; o<2; o++){
				phases.at(0) = dcp;
				phases.at(1) = d14;
	
				std::vector<double>mass;
				std::vector<double>ang;
				if(h==0)mass=mass_splittings;	
				if(h==1)mass=mass_splittings_inv;	

				if(o==0)ang=angles;
				if(o==1)ang=angles_oct;

				SBNprob myprob(4, ang, phases, mass);
			
				dunestream<<myprob.probabilityMatterExact(1,0,1,3.0,1300)<<" "<<myprob.probabilityMatterExact(1,0,-1,3.0,1300)<<" ";
				//dunestream<<myprob.probabilityGlobes(1,0,1,3.0,1300);//" "<<myantiprob.probabilityGlobes(1,0,-1,3.0,1300)<<" ";



				//dunestream<<myprob.probabilityGlobes(1,0,1,3.0,1300)<<" "<<myantiprob.probabilityGlobes(1,0,-1,3.0,1300)<<" ";
			}
		}
	dunestream<<std::endl;
	}
	}
	//plotProbabilityMatter(int a, int b, double Emin, double Emax, double L, double percen, double n);
	//myprob.plotProbabilityMatter(2, 1, -3, 2, L_DUNE, 0.15, 1500,&dunestream); //  mu->e
	/*myprob.plotProbabilityMatter(2, 2, -3, 2, L_DUNE, 0.15, 1500,&dunestream); //  mu->mu
	  myprob.plotProbabilityMatter(2, 3, -3, 2, L_DUNE, 0.15, 1500,&dunestream); //  mu -> tau
	  myprob.plotProbabilityMatter(2, 4, -3, 2, L_DUNE, 0.15, 1500,&dunestream); //  mu -> s
	  myprob.plotProbabilityMatter(1, 4, -3, 2, L_DUNE, 0.15, 1500,&dunestream); //  e ->s

	  myprob.plotProbabilityMatter(2, 1, -4, 1, L_SBND, 0.15, 1500,&sbndstream); //  mu->e
	  myprob.plotProbabilityMatter(2, 2, -4, 1, L_SBND, 0.15, 1500,&sbndstream); //  mu->mu
	  myprob.plotProbabilityMatter(2, 3, -4, 1, L_SBND, 0.15, 1500,&sbndstream); //  mu -> tau
	  myprob.plotProbabilityMatter(2, 4, -4, 1, L_SBND, 0.15, 1500,&sbndstream); //  mu -> s
	  myprob.plotProbabilityMatter(1, 4, -4, 1, L_SBND, 0.15, 1500,&sbndstream); //  e ->s

	  myprob.plotProbabilityMatter(2, 1, -4, 1, L_ICARUS, 0.15, 1500,&icarusstream); //  mu->e
	  myprob.plotProbabilityMatter(2, 2, -4, 1, L_ICARUS, 0.15, 1500,&icarusstream); //  mu->mu
	  myprob.plotProbabilityMatter(2, 3, -4, 1, L_ICARUS, 0.15, 1500,&icarusstream); //  mu -> tau
	  myprob.plotProbabilityMatter(2, 4, -4, 1, L_ICARUS, 0.15, 1500,&icarusstream); //  mu -> s
	  myprob.plotProbabilityMatter(1, 4, -4, 1, L_ICARUS, 0.15, 1500,&icarusstream); //  e ->s

	  myprob.setAntiNeutrinoMode(true);
	//AntiNeutrino
	myprob.plotProbabilityMatter(-2, 1, -3, 2, L_DUNE, 0.15, 1500,&dunestream); //  mu->e
	myprob.plotProbabilityMatter(-2, 2, -3, 2, L_DUNE, 0.15, 1500,&dunestream); //  mu->mu
	myprob.plotProbabilityMatter(-2, 3, -3, 2, L_DUNE, 0.15, 1500,&dunestream); //  mu -> tau
	myprob.plotProbabilityMatter(-2, 4, -3, 2, L_DUNE, 0.15, 1500,&dunestream); //  mu -> s
	myprob.plotProbabilityMatter(-1, 4, -3, 2, L_DUNE, 0.15, 1500,&dunestream); //  e ->s
	 */

	dunestream.close();
	//	sbndstream.close();
	//	icarusstream.close();


	/*
	   if(false){//this is giff mode!
	   myprob.t34 = T34*3.14159/180.0;
	   myprob.init();



	   std::cout<<"# t34: "<<T34<<std::endl;
	//int plotProbabilityMatter(int a, int b, double Emin, double Emax, double L, double percen, double n);
	//		myprob.plotProbabilityMatter(1, 0, -1, 2, 1300, 0.1, 1500);
	//		myprob.plotProbabilityMatter(1, 1, -1, 2, 1300, 0.1, 1500);
	//		myprob.plotProbabilityMatter(1, 2, -1, 2, 1300, 0.1, 1500);
	//		myprob.plotProbabilityMatter(1, 3, -1, 2, 1300, 0.1, 1500);
	}
	 */
}


