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



#include "genDune.h"

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


	std::string xml = "../../xml/dune.xml";
	int iarg = 0;
	opterr=1;
	int index; 
	int test_mode=0;
	std::string filename = "default.root";
	std::string which_mode = "default";
	/*************************************************************
	 *************************************************************
	 *		Command Line Argument Reading
	 ************************************************************
	 ************************************************************/

	const struct option longopts[] = 
	{
		{"xml", 		required_argument, 	0, 'x'},
		{"test",		required_argument,	0, 't'},
		{"mode",		required_argument,	0, 'm'},
		{"file",		required_argument,	0, 'f'},
		{0,			no_argument, 		0,  0},
	};


	while(iarg != -1)
	{
		iarg = getopt_long(argc,argv, "x:m:t:f:", longopts, &index);

		switch(iarg)
		{
			case 'x':
				xml = optarg;
				break;
			case 'f':
				filename = optarg;//`strtof(optarg,NULL);
				break;
			case 'm':
				which_mode = optarg;//`strtof(optarg,NULL);
				break;

			case 't':
				test_mode = strtof(optarg,NULL);
				break;
			case '?':
			case 'h':
				std::cout<<"Allowed arguments:"<<std::endl;
				std::cout<<"\t-x\t--xml\t\tInput .xml file for SBNconfig"<<std::endl;
				return 0;
		}

	}

	//SBNspec * dune_spec = new SBNspec("~/work/pheno/DUNE+SBN/yeonjae/sb_macros/DUNE_bf"  , xml);
	//dune_spec->writeOut("dune_test.root");

	if(which_mode == "default"){
	//*************************************************************************//
	//*************************************************************************//
	std::cout<<" Starting Default Mode: "<<std::endl;
	//*************************************************************************//
	//*************************************************************************//
		std::vector<double> angles = {30, 44, 8, 0,0,0};
		std::vector<double> phases = {0,0,0};
		std::vector<double> mass_splittings = {7.5*pow(10,-5), 2.552*pow(10,-3),0};
		genDUNE bkg_only(xml);

		bkg_only.prob = new SBNprob(4, angles, phases, mass_splittings);
		bkg_only.preCalculateProbs();

		bkg_only.doMC("three_neutrino");
		
		return 0;



	} else if(which_mode == "dcp"){
		//*************************************************************************//
		//*************************************************************************//
		std::cout<<" Starting Default Mode: "<<std::endl;
		//*************************************************************************//
		//*************************************************************************//



		std::vector<double> angles = {30, 44, 8, 0,0,0};
		std::vector<double> phases = {0,0,0};
		std::vector<double> mass_splittings = {7.5*pow(10,-5), 2.552*pow(10,-3),0};
		genDUNE dcp0(xml);
		genDUNE dcp180 = dcp0;

		dcp0.prob = new SBNprob(4, angles, phases, mass_splittings);
		dcp0.preCalculateProbs();

		phases.at(0) = 180;
		dcp180.prob = new SBNprob(4, angles, phases, mass_splittings);
		dcp180.preCalculateProbs();


		dcp0.doMC("gentest0");
		dcp180.doMC("gentest180");


		//SBNspec dune_0("~/work/pheno/DUNE+SBN/build/src/dune_dcp0",xml);
		//dune_0.compressVector(); 
		//SBNspec dune_180("~/work/pheno/DUNE+SBN/build/src/dune_dcp180",xml); 
		//dune_180.compressVector(); 

		SBNchi mychi0(dcp0);
		SBNchi mychi180(dcp180);

		std::ostringstream out;
		out << std::fixed;
		out << std::setprecision(2);

		std::ofstream dunestream;
		dunestream.open ("DUNE_data.dat");

		for(double dcp = 0.0; dcp <= 360; dcp += 5.0){

			genDUNE testGen(xml);
			phases.at(0) = dcp; 

			testGen.prob = new SBNprob(4, angles, phases, mass_splittings);
			testGen.preCalculateProbs();

			testGen.doMC("gentest");
			//testGen.writeOut("dat/out_"+std::to_string(dcp)+".root");
			testGen.compressVector();
			double chi0 = mychi0.CalcChi(&testGen);
			double chi180 = mychi180.CalcChi(&testGen);

			dunestream<<"ANK: "<<dcp<<" "<<chi0<<" "<<chi180<<" "<<std::min(chi0,chi180)<<std::endl;
		}

		dunestream.close();

		return 0;

	}
}




