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

std::vector<double> element_convert_to_params (double Ue4, double Um4, double phase14, double phase24, double phase34, bool simple);


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

	std::vector<double> theta23 = {38,40,42,44,45,46,47,49,51,53};
	std::vector<std::string> order_names = {"NO","IO"};
	std::vector<double> order_vals = {2.457*pow(10,-3), -2.449*pow(10,-3)};
    
    bool want_simple = true;
    
    
    
    //will read txt file and um...yeah...
    std::ifstream in("/a/data/westside/yjwa/NW/norwegian_wood/build/src/process_t14_t24_t34_d14.txt");
    
    std::string str;
    std::string delimiter = "_";
    
    std::ofstream paramchistream;
    paramchistream.open ("test.txt");//yj

    
    TFile *fcov_dune_in = new TFile("../../covar/covariance_matrices_xcheck_1408x1408.root","read");
    TFile *fcov_sbn_in = new TFile("../../covar/SBN_covariance_matrices_xcheck_690x690.root","read");
    
    TMatrixT<double> * frac_covar_DUNE = (TMatrixT<double>*)fcov_dune_in->Get("TMatrixT<double>;1");
    TMatrixT<double> * frac_covar_SBN = (TMatrixT<double>*)fcov_sbn_in->Get("TMatrixT<double>;7");
    
    std::string whipping_star_location = "/a/data/westside/yjwa/NW/whipping_star/";//yj
    
    //SBNspec* sbn_bkg_spec = new SBNspec((whipping_star_location+"data/precomp/SBN_bkg_all").c_str(), sbn_xml);
    //std::cout << "Combined fit :: lodaing SBN bkg ...  " << dune_bkg_name << std::endl;
    //sbn_bkg_spec->compressVector();
    
    
    //neutrinoModel bkgModel(0.1,0.0,0.0);//not really sure if this returns 3+0 osc.

    
    SBNosc * sbn_bkg_spec = new SBNosc((whipping_star_location+"data/precomp/SBN_bkg_all").c_str(),sbn_xml);
    sbn_bkg_spec->setBothMode();
    sbn_bkg_spec->compressVector();
    //sbn_bkg_spec->load_model(bkgModel);
    //sbn_bkg_spec->OscillateThis();
    //sbn_bkg_spec->compressVector();
    
    

    int line_count = 0;
    while (std::getline(in, str)) {
      std::cout << " this shoud be working " << std::endl;
        // output the line
        //std::cout << str << std::endl;
        
        size_t pos = 0;
        std::string sub;
        
        int which_param = 0;
        
        std::vector<std::string> lineinfile;
	int after_process_num = 0;        
	std::string precalc_name;
        while( (pos = str.find(delimiter)) != std::string::npos) {
            
            sub = str.substr(0,pos);
            
            
            //std::cout << sub << std::endl;
            lineinfile.push_back(sub);
            str.erase(0, pos+delimiter.length());
            if (after_process_num == 0){
                
	      precalc_name = str;
	    }
	    after_process_num++;

        }
        lineinfile.push_back(str);
        //std::cout << str << std::endl;
        line_count++;
        
        
        //std::cout << lineinfile.at(0) << " " << lineinfile.at(1) << " " << lineinfile.at(2) << std::endl;
        //if (lineinfile.at(2) > 2.00) {
	//  std::cout << "This mass " <<  lineinfile.at(2) << " is out of bound, want to skip the evts" << std::endl;
	//  continue;
        //}
        /*
         std::vector<double> test_sterile_angles = element_convert_to_params (lineinfile.at(0), lineinfile.at(1), 0, 0, 0, want_simple);
        
        //std::cout << "t14 : " << test_sterile_angles.at(0) << " , t24: " << test_sterile_angles.at(1) << " , t34: " << test_sterile_angles.at(2) <<std::endl;
        angles.at(3) = test_sterile_angles.at(0);
        */
	double t14;
	double t24;

	t14 = std::stod(lineinfile.at(2));
	t24 = std::stod(lineinfile.at(4));

	//std::cout << "trying to call signal model:" << std::endl;

        neutrinoModel signalModel(1.7, std::sin(t14*3.14159/180.), std::cos(t14*3.14159/180.)*std::sin(t24*3.14159/180.));//m4, ue4. um4
        
	//std::cout << "called signal model : " << std::endl;
        
        
        
        //Now create a oscillation spectra, constructed the same.
        SBNosc * sbn_sig_spec = new SBNosc((whipping_star_location+"data/precomp/SBN_bkg_all").c_str(),sbn_xml);

	//std::cout << " this shoud be working ... " << std::endl;
        sbn_sig_spec->compressVector();

	//std::cout << " this shoud be working .... " << std::endl;
        sbn_sig_spec->load_model(signalModel);
	//std::cout << " this shoud be working ..... " << std::endl;
        sbn_sig_spec->OscillateThis();
	//std::cout << " this shoud be working ...... " << std::endl;
        sbn_sig_spec->compressVector();
	//std::cout << " this shoud be working ....... " << std::endl;
        
        
        
        
        //	SBNchi dune_chi(*dune_bkg_spec, *frac_covar_DUNE);
        SBNchi sbn_chi(*sbn_bkg_spec, *frac_covar_SBN);
	//std::cout << " this shoud be working ........ " << std::endl;
        //	double ans_dune = dune_chi.CalcChi(dune_sig_spec);
        double ans_sbn = sbn_chi.CalcChi(sbn_sig_spec);
	//std::cout << " this shoud be working ......... " << std::endl;
        //std::cout<<ans_dune<<" "<<ans_sbn<<" "<<ans_dune+ans_sbn<<std::endl;
        
        //std::cout<<ans_sbn<<std::endl;
        
	paramchistream<<precalc_name<<"_SBN_"<<ans_sbn<<std::endl;
	std::cout << precalc_name << "_SBN_" << ans_sbn << std::endl;


        
	//	std::cout << " this shoud be working .......... " << std::endl;
        
        
    }
    std::cout << "How many lines in 90CL : " << line_count << std::endl;

    std::cout << "convert element to params test" << std::endl;
    
    
    
    
    

	//construct a signalModel (this is because the way SBN took in its parameters is different. must change that in future)
	//neutrinoModel signalModel(mass_splittings.back(), sin(angles.at(3)*3.1415/180.0), sin(angles.at(3)*3.14159/180));
	//mass_splittings.back() == m41, angles.at(3) = theta41
    
    
    
	//so currently the chi^2 = chi_dune^2 + chi_sbn^2
	//each chi_dune^2 and chi_sbn^2 needs a covariance matrix, a signal and a background 

	//Lets get covariance matricies for both DUNE and SBN programme
	//this might have to be updated 
	//TFile *fcov_dune_in = new TFile("../../covar/covariance_matrices_xcheck_1408x1408.root","read");
	//TFile *fcov_sbn_in = new TFile("../../covar/SBN_covariance_matrices_xcheck_690x690.root","read");

	//TMatrixT<double> * frac_covar_DUNE = (TMatrixT<double>*)fcov_dune_in->Get("TMatrixT<double>;1");
	//TMatrixT<double> * frac_covar_SBN = (TMatrixT<double>*)fcov_sbn_in->Get("TMatrixT<double>;7");

	//First lets generate the DUNE spectra		
	//The dune spectra is calculated on the FLY with full matter effects, so takes a while.     
	
    //for now just SBN for CL coverage test
    /*
    genDUNE * dune_sig_spec = new genDUNE(dune_xml);
	dune_sig_spec->prob = new SBNprob(4,angles,phases, mass_splittings);
	dune_sig_spec->preCalculateProbs();		
	dune_sig_spec->doMC();
    

	//And for background, load up an example precomputed background
	std::string dune_bkg_name = order_names.at(0)+"_DCP_"+to_string_prec(0.0,3)+"_T23_"+to_string_prec(theta23.at(6),3);
	SBNspec * dune_bkg_spec = new SBNspec(("precomp/"+dune_bkg_name+".SBNspec").c_str(),dune_xml);	
    std::cout << "Combined fit :: lodaing DUNE bkg ...  " << dune_bkg_name << std::endl;

    
    
	//Now create a background only spectra from preloaded. This should point to whipping star!
	//std::string whipping_star_location = "/home/mark/work/SBNfit/whipping_star/";
    std::string whipping_star_location = "/Users/yeon-jaejwa/sandbox/ws/whipping_star/";//yj

	SBNspec* sbn_bkg_spec = new SBNspec((whipping_star_location+"data/precomp/SBN_bkg_all").c_str(), sbn_xml);
    
    
    //std::cout << "Combined fit :: lodaing SBN bkg ...  " << dune_bkg_name << std::endl;
    
	sbn_bkg_spec->compressVector();

	//Now create a oscillation spectra, constructed the same.
	SBNosc * sbn_sig_spec = new SBNosc((whipping_star_location+"data/precomp/SBN_bkg_all").c_str(),sbn_xml);
	sbn_sig_spec->compressVector();
	sbn_sig_spec->load_model(signalModel);
	sbn_sig_spec->OscillateThis();
	sbn_sig_spec->compressVector();
    
    
	

//	SBNchi dune_chi(*dune_bkg_spec, *frac_covar_DUNE);
	SBNchi sbn_chi(*sbn_bkg_spec, *frac_covar_SBN);

//	double ans_dune = dune_chi.CalcChi(dune_sig_spec);
	double ans_sbn = sbn_chi.CalcChi(sbn_sig_spec);
	
	//std::cout<<ans_dune<<" "<<ans_sbn<<" "<<ans_dune+ans_sbn<<std::endl;
    
    std::cout<<ans_sbn<<std::endl;
	
	
	//truth->compareSBNspecs(sterile,"DUNE_compare.root");*/




}


std::vector<double> element_convert_to_params (double Ue4, double Um4, double phase14, double phase24, double phase34, bool simple){
    
    /*
     Following PMNS matrix parametrization #17 and #20 (they have same 4th column) from
     https://drive.google.com/file/d/0B6i555wZRmqDWWUtaUhxMUxELXc/view
     
     Ue4 = e^(-i*delta14) * sin(theta14)
     Um4 = cos(theta14) * cos(theta24)
     Ut4 = cos(theta14) * cos(theta24) * e^(-i*delta34) * sin(theta34)
     Us4 = cos(theta14) * cos(theta24) * cos(theta34)
     
     |Us4|^2 = 1 - |Ue4|^2 - |Um4|^2 -|Ut4|^2
     
     unknowns are theta14, 24, 34, Ut4, Us4
     
     */
    
    
    
    double theta14 = -1000;
    double theta24 = -1000;
    double theta34 = -1000;
    
    // for the simpler case first,
    // Ut4 = 0 , phase14, 24, 34 = 0
    
    if (simple){
        theta14 = std::asin(Ue4);
        theta24 = std::asin(Um4/std::cos(theta14));
        double Us4sq = 1 - Ue4*Ue4 - Um4*Um4;
        double Us4 = std::sqrt(Us4sq);
        theta34 = std::acos(Us4/std::cos(theta14)/std::cos(theta24));
    }
    
    //radians to degrees
    theta14 *= 180./3.14159;
    theta24 *= 180./3.14159;
    theta34 *= 180./3.14159;
    
    std::vector<double> thetas;
    thetas.push_back(theta14);
    thetas.push_back(theta24);
    thetas.push_back(theta34);

    
    return thetas;
    
}
