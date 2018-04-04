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
		//std::cout<<oscillation_patterns.at(file).first<<" "<<oscillation_patterns.at(file).second<<std::endl;	

		oscprob_far = interpolate_prob_far(oscillation_patterns.at(file).first,oscillation_patterns.at(file).second, Enu_true);
		if(oscillation_patterns.at(file).first != oscillation_patterns.at(file).second){
			oscprob_near =0;
		}

//		oscprob_near = interpolate_prob_near(oscillation_patterns.at(file).first,oscillation_patterns.at(file).second, Enu_true);
		//oscprob_far = prob->probabilityMatterExact(oscillation_patterns.at(file).first, oscillation_patterns.at(file).second ,Enu_true, 1300);
		//oscprob_near = prob->probabilityMatterExact(oscillation_patterns.at(file).first, oscillation_patterns.at(file).second ,Enu_true, 1);
	}

	//std::cout<<Enu_true<<" "<<Enu_reco<<" on file: "<<multisim_name.at(file)<<" "<<nutype<<" PROB: "<<oscprob_far<<std::endl;
	//std::cout<<"Map Hist: "<<map_hist[multisim_name.at(file)]<<std::endl;


	hist.at(map_hist[multisim_name.at(file)]).Fill(Enu_reco, oscprob_far*weight*far_detector_weight);
	hist.at(map_hist[near_detector_names.at(file)]).Fill(Enu_reco, oscprob_near*weight*near_detector_weight);

	return 0;
}

int genDUNE::tidyHistograms(){
	return 0;
}

double genDUNE::interpolate_prob_far(int a, int b, double enu){
	int rnd = std::floor(enu/0.01);

	double p1 = precalc_prob_far.at(a).at(b).at(rnd);
	double p2 = precalc_prob_far.at(a).at(b).at(rnd+1);
	
	double e1 = (double)rnd*0.01;
	double e2 = (double)(rnd+1.0)*0.01;


	return lin_interp(enu, p1,e1,p2,e2    );
}



double genDUNE::interpolate_prob_near(int a, int b, double enu){
	int rnd = std::floor(enu/0.01);

	double p1 = precalc_prob_near.at(a).at(b).at(rnd);
	double p2 = precalc_prob_near.at(a).at(b).at(rnd+1);
	
	double e1 = (double)rnd*0.01;
	double e2 = (double)(rnd+1.0)*0.01;


	return this->lin_interp(enu, p1,e1,p2,e2    );
}

double genDUNE::lin_interp(double ein, double p1, double e1, double p2, double e2){

return p1+(ein-e1)*(p2-p1)/(e2-e1);

}

int genDUNE::preCalculateProbs(){
	
	std::cout<<"Starting precalculating probs"<<std::endl;
	for(int i = 0; i < 4; i++){
		std::vector<std::vector<double>> tmp_near;
		std::vector<std::vector<double>> tmp_far;

		for(int j = 0; j < 4; j++){
			std::vector<double> tmpen_near;
			std::vector<double> tmpen_far;

			for(double en = 0.001; en < 50; en+= 0.01){
				tmpen_far.push_back(prob->probabilityMatterExact(i,j,en,1300));
				tmpen_near.push_back(prob->probabilityVacuumExact(i,j,en,1.0));
				//std::cout<<"P"<<i<<j<<" "<<en<<" "<<tmpen_far.back()<<" "<<tmpen_near.back()<<std::endl;
			}

		tmp_near.push_back(tmpen_near);
		tmp_far.push_back(tmpen_far);

		} 
	precalc_prob_far.push_back(tmp_far);	
	precalc_prob_near.push_back(tmp_near);	

	}
	std::cout<<"Done precalculating probs"<<std::endl;

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
