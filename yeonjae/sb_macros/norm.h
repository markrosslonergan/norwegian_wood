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

	gst_file(TString innam, int inbeam , int intype, double innorm) : filename(innam), beam_type(inbeam), nu_type(intype), norm(innorm){}; 

	int setHistLocations(std::vector<string> strin){
		std::string horncurrent;
		if(beam_type==0) horncurrent = "nu";
		if(beam_type==1) horncurrent = "nubar";

		for(auto &s: strin){

			if(s!=""){
				hists.push_back( horncurrent+"_dune_"+s);
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


#endif
