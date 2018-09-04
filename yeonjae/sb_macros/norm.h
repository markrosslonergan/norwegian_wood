#ifndef NORM_H_
#define NORM_H_

//#include "TRandom.h"
#include "TRandom3.h"
#include "TString.h"
#include <vector>
#include <string>
#include "TLorentzVector.h"//yj

#include "params.h"

struct gst_file{

	int nu_type;
	int beam_type;
	TString filename;
	double norm;
	std::vector<std::string> hists;
	std::vector<bool> in_use;	

	std::string tag;

	gst_file(TString innam, int inbeam , int intype, double innorm) : filename(innam), beam_type(inbeam), nu_type(intype), norm(innorm){ tag = "dune";}; 
	gst_file(TString innam, int inbeam , int intype, double innorm, std::string intag) : filename(innam), beam_type(inbeam), nu_type(intype), norm(innorm), tag(intag){}; 
	
	int setHistLocations(std::vector<string> strin){
		std::string horncurrent;
		if(beam_type==0) horncurrent = "nu";
		if(beam_type==1) horncurrent = "nubar";

		for(auto &s: strin){

			if(s!=""){
				hists.push_back( horncurrent+"_"+tag+"_"+s);
				in_use.push_back(true);
			}else{
				hists.push_back("");
				in_use.push_back(false);
			}
		}

		return 0;
	};

};




double get_normalization(TString name){

	//NEUTRINO MODE
	if(name == "gntp.0.numu50k_gst" || name == "gntp.0.numuflux_nuebeam50k_gst"){
		return 1.10634;
	}else if(name == "gntp.0.nutau20k_gst"){
		return 2.76584; 
	}else if(name == "gntp.0.numubar10k_gst"  || name == "gntp.0.numubarflux_nuebarbeam10k_gst"){
		return 0.181760;
	}else if(name == "gntp.0.nutaubar20k_gst"){
		return 0.0908800;
	}else if(name == "gntp.0.nue20k_gst"){
		return 0.0323551;
	}else if(name == "gntp.0.nuebar20k_gst"){
		return 0.00404113;

		//ANTI NEUTRINO MODE 
	}else if(name == "gntp.0.RHC_FD_numuflux_numubeam20k_gst"){
		return 0.256093;
	} else if( name == "gntp.0.RHC_FD_numuflux_nutaubeam10k_gst"){
		return 0.512185;
	} else if(name =="gntp.0.RHC_FD_numuflux_nuebeam10k_gst"){
		return 0.512185;

	}else if(name == "gntp.0.RHC_FD_numubarflux_numubarbeam50k_gst"){
		return 0.395868;
	}else if(name == "gntp.0.RHC_FD_numubarflux_nutaubarbeam20k_gst"){
		return 0.989670;
	} else if(name == "gntp.0.RHC_FD_numubarflux_nuebarbeam20k_gst"){
		return 0.989670;
	}else if(name == "gntp.0.RHC_FD_nuebarflux_nuebarbeam10k_gst"){
		return 0.0213827;
	}else if(name == "gntp.0.RHC_FD_nueflux_nuebeam10k_gst"){
		return 0.0215175;

	}else{
		std::cout<<"get_normalization || ERROR, passed in an non-existant filename/norm not calculated for it"<<std::endl;
		exit(EXIT_FAILURE);

	}	



}

double get_mulike_nutau_kNN_eff(double energy){

    double bkg_nutau_effi_ratio[40] = {1.,
        1.,
        0.325581395348837,
        0.397050482132729,
        0.425581395348838,
        0.455324357405141,
        0.48062015503876,
        0.475631951466128,
        0.53359173126615,
        0.516015796401931,
        0.552501761804088,
        0.544920440636476,
        0.52093023255814,
        0.551573187414501,
        0.55098389982111,
        0.625116279069768,
        0.605306256141501,
        0.622851365015167,
        0.674091057975762,
        0.721288014311271,
        0.824127906976744,
        0.787226657410621,
        0.837209302325581,
        0.901610017889087,
        0.900436046511628,
        0.887949260042283,
        0.945736434108528,
        0.913728432108027,
        0.850712678169541,
        0.877414268821443,
        0.976744186046511,
        0.913728432108027,
        0.960990247561891,
        0.850712678169541,
        0.801804928844152,
        0.871556350626118,
        0.96124031007752,
        0.993023255813955,
        1.,
        0.837209302325581};

	int pos = floor(energy/0.25);
	if (pos > 39) return 1.0;
	
	return bkg_nutau_effi_ratio[pos];

    
}

double get_mulike_NC_kNN_eff(double energy){

	 double bkg_NC_effi_ratio[40] = {
    1.,
    1.,
    0.48837209302324,
    0.488372093023265,
    0.651162790697676,
    0.651162790697676,
    0.488372093023265,
    0.488372093023265,
    0.488372093023265,
    0.390697674418608,
    0.325581395348838,
    0.781395348837217,
    0.488372093023257,
    0.390697674418608,
    0.390697674418608,
    0.488372093023257,
    0.488372093023257,
    0.418604651162796,
    0.418604651162796,
    0.418604651162796,
    0.732558139534889,
    0.488372093023257,
    0.418604651162796,
    0.418604651162796,
    0.488372093023259,
    0.813953488372078,
    0.651162790697683,
    0.434108527131789,
    0.293023255813956,
    0.390697674418608,
    0.48837209302325,
    0.558139534883728,
    0.558139534883728,
    0.355179704016919,
    0.325581395348842,
    0.325581395348842,
    0.177589852008459,
    0.355179704016919,
    0.390697674418608,
    0.217054263565894
        };

	int pos = floor(energy/0.25);
	if (pos > 39) return 1.0;
	
	return bkg_NC_effi_ratio[pos];

    
}

 
double get_elike_nutau_kNN_eff(double energy){

	double bkg_nutau_effi_ratio[40] = {
        0.185177101265301,
        0.0569097745791012,
        0.142302952010498,
        0.271548091251609,
        0.365215909373156,
        0.486457828608005,
        0.463130167051557,
        0.49821981228389,
        0.508886165054105,
        0.523416189301326,
        0.560997511737675,
        0.518695827424634,
        0.545142155498725,
        0.530004233403559,
        0.562003158877504,
        0.509813049542501,
        0.509139483299998,
        0.486598707120586,
        0.48577769000742,
        0.466012371584107,
        0.442657261327112,
        0.444240085377149,
        0.483963356356026,
        0.42227539095115,
        0.43934426885497,
        0.375047996507843,
        0.407865995461402,
        0.363710583961687,
        0.390076613505933,
        0.357954552368202,
        0.396268301879026,
        0.387585758062179,
        0.385305674134454,
        0.306145466796258,
        0.380393505599226,
        0.374360284108998,
        0.413174986327187,
        0.372564212374476,
        0.349898958801367,
        0.334581841047651};
   	int pos = floor(energy/0.25);
	if (pos > 39) return 0.3345;
	
	return bkg_nutau_effi_ratio[pos];

    
}

double get_elike_NC_kNN_eff(double energy){
 
    double bkg_NC_effi_ratio[40] = {
        0.18047465517046,
        0.101151907725579,
        0.232938787430535,
        0.337010457724242,
        0.648029150535554,
        0.501602340213943,
        0.645180823288289,
        0.931730203198994,
        0.761877907924052,
        0.776838758789639,
        0.372765819839769,
        0.84718572080561,
        0.65233019873389,
        0.593016077460186,
        0.507954119678036,
        0.6212656967605,
        0.57345478472973,
        0.620662846465288,
        0.507703967080246,
        0.559117027037734,
        0.348976890933601,
        0.349160950475233,
        0.372580987352924,
        0.232986701371153,
        0.279438769960226,
        0.207053277186504,
        0.155284049293359,
        0.310483034741583,
        0.20695622854475,
        0.232985115308552,
        0.310737887138773,
        0.132925008629442,
        0.0932994577334251,
        0.186445714386841,
        0.349048037607936,
        0.399259898770444,
        0.18651271829387,
        0.310834997254636,
        0.103361131972826,
        0.116469686089845};
        
	int pos = floor(energy/0.25);
	if (pos > 39) return bkg_NC_effi_ratio[39];
	
	return bkg_NC_effi_ratio[pos];

    
}
       
    //double bkg_nue_effi_ratio[40] = {
    //        };
    
   


#endif
