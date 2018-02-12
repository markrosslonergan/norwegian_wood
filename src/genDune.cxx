#include "genDune.h"

bool genDUNE::eventSelection(int file){
	return true;
}


int genDUNE::fillHistograms(int file, int uni, double wei){
	//	nu_dune_elike_fullosc;24	nu_dune_elike_fullosc
	// 	nu_dune_elike_intrinsic;24	nu_dune_elike_intrinsic
	// 	nu_dune_elike_mumisid;24	nu_dune_elike_mumisid
	//	nu_dune_elike_taumisid;24	nu_dune_elike_taumisid
	//	nu_dune_elike_ncmisid;18	nu_dune_elike_ncmisid
	//	nu_dune_mulike_intrinsic;16	nu_dune_mulike_intrinsic
	//	nu_dune_mulike_taumisid;16	nu_dune_mulike_taumisid
	//	nu_dune_mulike_ncmisid;16	nu_dune_mulike_ncmisid

	double Enu_true = *vmapD[file]["Etrue"];
	double Enu_reco =  *vmapD[file]["Ereco"];
	double weight = *vmapD[file]["Weight"];
	int nutype = *vmapI[file]["NuType"];

	//write a map for osc patterns 

	
	double oscprob_far =1.0; 
	double oscprob_near =1.0; 

	if(oscillation_patterns.at(file).first >0){	
		oscprob_far = prob->probabilityMatterExact(oscillation_patterns.at(file).first, oscillation_patterns.at(file).second ,Enu_true, 1300);
		oscprob_near =prob->probabilityMatterExact(oscillation_patterns.at(file).first, oscillation_patterns.at(file).second ,Enu_true, 1);
	}

	//std::cout<<Enu_true<<" "<<Enu_reco<<" on file: "<<multisim_name.at(file)<<" "<<nutype<<" PROB: "<<oscprob<<std::endl;
	//std::cout<<"Map Hist: "<<map_hist[multisim_name.at(file)]<<std::endl;

	hist.at(map_hist[multisim_name.at(file)]).Fill(Enu_reco, oscprob_far*weight*far_detector_weight);
	hist.at(map_hist[near_detector_names.at(file)]).Fill(Enu_reco, oscprob_near*weight*near_detector_weight);

	return 0;
}

int genDUNE::tidyHistograms(){
	return 0;
}







std::pair<int,int> genDUNE::getOscPattern(std::string name){
	//This is a silly wy fo doing this..
	std::pair<int, int> ans;

	for(std::string &m:mode_names){
		for(std::string &d:detector_names){
			for(int i=0; i<num_channels; i++){
				for(int j=0; j<num_subchannels.at(i); j++){
					if(name == m+"_"+d+"_"+channel_names.at(i)+"_"+subchannel_names.at(i).at(j)){

						int pattern = subchannel_osc_patterns.at(i).at(j); 

						if(pattern == 0){
							//dont oscillate!
							ans.first = -99; ans.second =-99;
						}else if(pattern >0){
							//Neutrino Mode
							std::string s_pattern = std::to_string(pattern);
							ans.first = (int)s_pattern.at(0) - '0';
							ans.second = (int)s_pattern.at(1) - '0';
						}else{
							//Antineutrino.
							std::string s_pattern = std::to_string(pattern);
							ans.first = (int)s_pattern.at(1) -'0';
							ans.second = (int)s_pattern.at(2) -'0';
							prob->useAntiNeutrino = true;
							prob->init();
						}

					}

				}
			}
		}
	}	
	ans.first = ans.first-1;
	ans.second= ans.second-1;
	return ans;
}
