#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <cmath>
#include <vector>
#include <unistd.h>
#include <getopt.h>
#include <cstring>
#include <algorithm>
#include <iterator>

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


	template <typename T>
std::string to_string_prec(const T a_value, const int n = 6)
{
	std::ostringstream out;
	out <<std::fixed<< std::setprecision(n) << a_value;
	return out.str();
}






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


	TRandom3 *rangen= new TRandom3();
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
	//these might have to be updated
	std::string dune_xml = "../../xml/dune.xml";
	std::string sbn_xml = "../../xml/sbn_osc.xml";

	std::vector<double> angles = {33.6, 42.1, 8.5, 10,10,0};
	std::vector<double> phases = {0,0,0};
	std::vector<double> mass_splittings = {7.5*pow(10,-5), 2.457*pow(10,-3), 1.0};
	std::vector<double> mass_splittings_inv = {7.5*pow(10,-5), -2.449*pow(10,-3), 1.00};

	std::vector<double> theta23 = {40,41,42,43,44,45,46,47,48,49,50,51,52};
	std::vector<std::string> order_names = {"NO","IO"};
	std::vector<double> order_vals = {2.457*pow(10,-3), -2.449*pow(10,-3)};

	//construct a signalModel (this is because the way SBN took in its parameters is different. must change that in future)
	neutrinoModel signalModel(mass_splittings.back(), sin(angles.at(3)*3.1415/180.0), sin(angles.at(3)*3.14159/180));
	
	//so currently the chi^2 = chi_dune^2 + chi_sbn^2
	//each chi_dune^2 and chi_sbn^2 needs a covariance matrix, a signal and a background 

	//Lets get covariance matricies for both DUNE and SBN programme
	//this might have to be updated 
	TFile *fcov_dune_in = new TFile("../../covar/covariance_matrices_xcheck_1408x1408.root","read");
	TFile *fcov_sbn_in = new TFile("../../covar/SBN_covariance_matrices_xcheck_690x690.root","read");

	TMatrixT<double> * frac_covar_DUNE = (TMatrixT<double>*)fcov_dune_in->Get("TMatrixT<double>;1");
	TMatrixT<double> * frac_covar_SBN = (TMatrixT<double>*)fcov_sbn_in->Get("TMatrixT<double>;7");

	//First lets generate the DUNE spectra		
	//The dune spectra is calculated on the FLY with full matter effects, so takes a while.     
	genDUNE * dune_sig_spec = new genDUNE(dune_xml);	
	dune_sig_spec->prob = new SBNprob(4,angles,phases, mass_splittings);
	dune_sig_spec->preCalculateProbs();		
	dune_sig_spec->doMC();

	//And for background, load up an example precomputed background
	std::string dune_bkg_name = order_names.at(0)+"_DCP_"+to_string_prec(0.0,3)+"_T23_"+to_string_prec(theta23.at(6),3);
	SBNspec * dune_bkg_spec = new SBNspec(("precomp/"+dune_bkg_name+".SBNspec").c_str(),dune_xml);	


	//Now create a background only spectra from preloaded. This should point to whipping star!
	std::string whipping_star_location = "/home/mark/work/SBNfit/whipping_star/";

	SBNspec* sbn_bkg_spec = new SBNspec((whipping_star_location+"data/precomp/SBN_bkg_all").c_str(), sbn_xml);
	sbn_bkg_spec->compressVector();

	//Now create a oscillation spectra, constructed the same.
	SBNosc * sbn_sig_spec = new SBNosc((whipping_star_location+"data/precomp/SBN_bkg_all").c_str(),sbn_xml);
	sbn_sig_spec->compressVector();
	sbn_sig_spec->load_model(signalModel);
	sbn_sig_spec->OscillateThis();
	sbn_sig_spec->compressVector();
	

	SBNchi dune_chi(*dune_bkg_spec, *frac_covar_DUNE);
	SBNchi sbn_chi(*sbn_bkg_spec, *frac_covar_SBN);

	double ans_dune = dune_chi.CalcChi(dune_sig_spec);
	double ans_sbn = sbn_chi.CalcChi(sbn_sig_spec);
	
	std::cout<<ans_dune<<" "<<ans_sbn<<" "<<ans_dune+ans_sbn<<std::endl;
	
	
	//truth->compareSBNspecs(sterile,"DUNE_compare.root");




}

