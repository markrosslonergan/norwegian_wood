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


	genDUNE testGen(xml);

	std::ostringstream out;
	out << std::fixed;
	out << std::setprecision(2);

	testGen.doMC("gentest");
	testGen.writeOut("out.root");

return 0;
}







