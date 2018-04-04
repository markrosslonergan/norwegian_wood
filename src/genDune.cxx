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

	if(oscillation_patterns.at(file).first != 0 && oscillation_patterns.at(file).second != 0){
		//std::cout<<oscillation_patterns.at(file).first<<" "<<oscillation_patterns.at(file).second<<std::endl;	

		oscprob_far = interpolate_prob_far(oscillation_patterns.at(file).first, oscillation_patterns.at(file).second, Enu_true);

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
	double p1,p2;

	if(a<0 && b<0 ){
		p1 = precalc_prob_farbar.at(-a-1).at(-b-1).at(rnd);
		p2 = precalc_prob_farbar.at(-a-1).at(-b-1).at(rnd+1);
	}else{
		p1 = precalc_prob_far.at(a-1).at(b-1).at(rnd);
		p2 = precalc_prob_far.at(a-1).at(b-1).at(rnd+1);
	}


	double e1 = (double)rnd*0.01;
	double e2 = (double)(rnd+1.0)*0.01;


	return lin_interp(enu, p1,e1,p2,e2 );
}



double genDUNE::interpolate_prob_near(int a, int b, double enu){
	int rnd = std::floor(enu/0.01);
	double p1,p2;

	if(a<0){
		p1 = precalc_prob_nearbar.at(-a-1).at(-b-1).at(rnd);
		p2 = precalc_prob_nearbar.at(-a-1).at(-b-1).at(rnd+1);
	}else{
		p1 = precalc_prob_near.at(a-1).at(b-1).at(rnd);
		p2 = precalc_prob_near.at(a-1).at(b-1).at(rnd+1);
	}

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
	std::cout<<"Starting precalculating probs BAR"<<std::endl;
	prob->useAntiNeutrino = true;
	prob->init();
	for(int i = 0; i < 4; i++){
		std::vector<std::vector<double>> tmp_nearbar;
		std::vector<std::vector<double>> tmp_farbar;

		for(int j = 0; j < 4; j++){
			std::vector<double> tmpen_nearbar;
			std::vector<double> tmpen_farbar;

			for(double en = 0.001; en < 50; en+= 0.01){
				tmpen_farbar.push_back(prob->probabilityMatterExact(i,j,en,1300));
				tmpen_nearbar.push_back(prob->probabilityVacuumExact(i,j,en,1.0));
				//std::cout<<"P"<<i<<j<<" "<<en<<" "<<tmpen_farbar.back()<<" "<<tmpen_nearbar.back()<<std::endl;
			}

			tmp_nearbar.push_back(tmpen_nearbar);
			tmp_farbar.push_back(tmpen_farbar);


		} 
		precalc_prob_farbar.push_back(tmp_farbar);	
		precalc_prob_nearbar.push_back(tmp_nearbar);	

	}
	std::cout<<"Done precalculating probs"<<std::endl;

	return 0;
}




int genDUNE::getOscPattern(){

	for(std::string &m:mode_names){
		for(std::string &d:detector_names){
			for(int i=0; i<num_channels; i++){
				for(int j=0; j<num_subchannels.at(i); j++){
					std::string name = m+"_"+d+"_"+channel_names.at(i)+"_"+subchannel_names.at(i).at(j);


					int pattern = subchannel_osc_patterns.at(i).at(j); 


					if(pattern == 0){


						std::pair<int, int> oscpair; oscpair.first = 0; oscpair.second = 0;
						osc_pattern_map[name] = oscpair;

					}else if(pattern >0){
						std::string s_pattern = std::to_string(pattern);

						int init_flavor = (int)s_pattern.at(0) - '0';
						int final_flavor = (int)s_pattern.at(1) - '0';

						std::pair<int, int> oscpair; oscpair.first = init_flavor; oscpair.second = final_flavor;
						osc_pattern_map[name] = oscpair;
					}else{
						std::string s_pattern = std::to_string(-pattern);

						int init_flavor = (int)s_pattern.at(0) - '0';
						int final_flavor = (int)s_pattern.at(1) - '0';

						std::pair<int, int> oscpair; oscpair.first = -init_flavor; oscpair.second = -final_flavor;
						osc_pattern_map[name] = oscpair;



					}




				}
			}
		}
	}	
	return 0;
}
