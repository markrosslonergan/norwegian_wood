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

	double Enu_true = *vmapD[file]["Etrue"];//yj vampD: vector of maps string::Double*
	double Enu_reco =  *vmapD[file]["Ereco"];
	double weight = *vmapD[file]["Weight"];
	int nutype = *vmapI[file]["NuType"];//yj vmapI: vector of maps string::Int*

	//write a map for osc patterns 

	double oscprob_far =1.0; 
	double oscprob_near =1.0; 


	//FILE is for input file, nothing to do with the order of fullnames. remember that dammin. the oscillation patterns must be matched with multisim_name(file) maps!

	//Map lookup of strings is slow
	//std::pair<int,int> oscillation_pattern = osc_pattern_map.at(multisim_name.at(file));
	std::pair<int,int> * oscillation_pattern = &osc_pattern_vec.at(file);  
	
	
	//yj what is std::pair? just two inputs.

	if(oscillation_pattern->first != 0 && oscillation_pattern->second != 0){
		//std::cout<<oscillation_patterns.at(file).first<<" "<<oscillation_patterns.at(file).second<<std::endl;	

		oscprob_far = interpolate_prob_far(oscillation_pattern->first, oscillation_pattern->second, Enu_true);
		oscprob_near = interpolate_prob_near(oscillation_pattern->first, oscillation_pattern->second, Enu_true);
		//std::cout<<"INV: "<<oscillation_pattern.first<<" "<<oscillation_pattern.second<<" "<<multisim_name.at(file)<<std::endl;

		//if(oscillation_patterns.at(file).first != oscillation_patterns.at(file).second){
		//	oscprob_near =0;
		//}
		//oscprob_near = interpolate_prob_near(oscillation_patterns.at(file).first,oscillation_patterns.at(file).second, Enu_true);
		//oscprob_far = prob->probabilityMatterExact(oscillation_patterns.at(file).first, oscillation_patterns.at(file).second ,Enu_true, 1300);
		//oscprob_near = prob->probabilityMatterExact(oscillation_patterns.at(file).first, oscillation_patterns.at(file).second ,Enu_true, 1);
	}


	//std::cout<<Enu_true<<" "<<Enu_reco<<" on file: "<<multisim_name.at(file)<<" "<<nutype<<" PROB: "<<oscprob_far<<std::endl;
	//std::cout<<"Map Hist: "<<map_hist[multisim_name.at(file)]<<std::endl;

	int m1 = map_hist_vec.at(file);
	int m2 = map_hist_vec_near.at(file);

	hist.at(m1).Fill(Enu_reco, oscprob_far*weight*far_detector_weight);
	hist.at(m2).Fill(Enu_reco, oscprob_near*weight*near_detector_weight);

//	hist.at(map_hist[multisim_name.at(file)]).Fill(Enu_reco, oscprob_far*weight*far_detector_weight);
//	hist.at(map_hist[near_detector_name_map.at(multisim_name.at(file))]).Fill(Enu_reco, oscprob_near*weight*near_detector_weight);

	return 0;
}

int genDUNE::tidyHistograms(){
	return 0;
}

double genDUNE::interpolate_prob_far(int a, int b, double enu){
	int rnd = std::floor(enu/0.01);
	double p1,p2;

	if(a<0 && b<0 ){
		//		std::cout<<"INV2: "<<-a-1<<" "<<-b-1<<" anti "<<std::endl;
		p1 = precalc_prob_farbar.at(-a-1).at(-b-1).at(rnd);
		p2 = precalc_prob_farbar.at(-a-1).at(-b-1).at(rnd+1);
	}else{
		//		std::cout<<"INV2: "<<a-1<<" "<<b-1<<" nuetrino "<<std::endl;
		p1 = precalc_prob_far.at(a-1).at(b-1).at(rnd);
		p2 = precalc_prob_far.at(a-1).at(b-1).at(rnd+1);
	}


	double e1 = (double)rnd*0.01;
	double e2 = (double)(rnd+1.0)*0.01;


	return lin_interp(enu, p1,e1,p2,e2 );
}



double genDUNE::interpolate_prob_near(int a, int b, double enu){//yj
	int rnd = std::floor(enu/0.01); //energy step is 0.01
	double p1,p2;

	if(a<0 && b < 0){
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



int genDUNE::loadPreCalculatedProbs(std::string dir){

	std::ifstream s1;   s1.open (dir+"prob_precalc_nu_near.dat");
	std::ifstream s2;   s2.open (dir+"prob_precalc_nu_far.dat");
	std::ifstream s3;   s3.open (dir+"prob_precalc_nubar_near.dat");
	std::ifstream s4;   s4.open (dir+"prob_precalc_nubar_far.dat");

	precalc_prob_far.resize(4);
	precalc_prob_near.resize(4);
	precalc_prob_farbar.resize(4);
	precalc_prob_nearbar.resize(4);

	for(int i=0; i<4; i++){
		precalc_prob_far.at(i).resize(4);
		precalc_prob_near.at(i).resize(4);
		precalc_prob_farbar.at(i).resize(4);
		precalc_prob_nearbar.at(i).resize(4);
	}

	if(s1.good())  // same as: if (myfile.good())
	{
		std::string line;
		while (getline(s1, line )){
			{
				//std::cout<<line<<std::endl;
				std::stringstream ss(line);
				std::istream_iterator<std::string> begin(ss);
				std::istream_iterator<std::string> end;
				std::vector<std::string> vstrings(begin, end);
				//std::copy(vstrings.begin(), vstrings.end(), std::ostream_iterator<std::string>(std::cout, "\n"));
				int a = std::stoi(vstrings.at(0));
				int b = std::stoi(vstrings.at(1));
				for(int t=2; t<vstrings.size();t++){
					precalc_prob_near.at(a).at(b).push_back(std::stod(vstrings.at(t)));
				}
			}
		}
	}  


	if(s2.good())  // same as: if (myfile.good())
	{
		std::string line;
		while (getline(s2, line )){
			{
				//std::cout<<line<<std::endl;
				std::stringstream ss(line);
				std::istream_iterator<std::string> begin(ss);
				std::istream_iterator<std::string> end;
				std::vector<std::string> vstrings(begin, end);
				//std::copy(vstrings.begin(), vstrings.end(), std::ostream_iterator<std::string>(std::cout, "\n"));
				int a = std::stoi(vstrings.at(0));
				int b = std::stoi(vstrings.at(1));
				for(int t=2; t<vstrings.size();t++){
					precalc_prob_far.at(a).at(b).push_back(std::stod(vstrings.at(t)));
				}
			}
		}
	}  


	if(s3.good())  // same as: if (myfile.good())
	{
		std::string line;
		while (getline(s3, line )){
			{
				//std::cout<<line<<std::endl;
				std::stringstream ss(line);
				std::istream_iterator<std::string> begin(ss);
				std::istream_iterator<std::string> end;
				std::vector<std::string> vstrings(begin, end);
				//std::copy(vstrings.begin(), vstrings.end(), std::ostream_iterator<std::string>(std::cout, "\n"));
				int a = std::stoi(vstrings.at(0));
				int b = std::stoi(vstrings.at(1));
				for(int t=2; t<vstrings.size();t++){
					precalc_prob_nearbar.at(a).at(b).push_back(std::stod(vstrings.at(t)));
				}
			}
		}
	}  


	if(s4.good())  // same as: if (myfile.good())
	{
		std::string line;
		while (getline(s4, line )){
			{
				//std::cout<<line<<std::endl;
				std::stringstream ss(line);
				std::istream_iterator<std::string> begin(ss);
				std::istream_iterator<std::string> end;
				std::vector<std::string> vstrings(begin, end);
				//std::copy(vstrings.begin(), vstrings.end(), std::ostream_iterator<std::string>(std::cout, "\n"));
				int a = std::stoi(vstrings.at(0));
				int b = std::stoi(vstrings.at(1));
				for(int t=2; t<vstrings.size();t++){
					precalc_prob_farbar.at(a).at(b).push_back(std::stod(vstrings.at(t)));
				}
			}
		}
	}  


	s1.close(); s2.close(); s3.close(); s4.close();
	return 0;   
}

int genDUNE::preCalculateProbs(){

	std::ofstream s1;   s1.open ("prob_precalc_nu_near.dat");
	std::ofstream s2;   s2.open ("prob_precalc_nu_far.dat");
	std::ofstream s3;   s3.open ("prob_precalc_nubar_near.dat");
	std::ofstream s4;   s4.open ("prob_precalc_nubar_far.dat");

	//yj i,j, are the number of flavors.
	std::cout<<"Starting precalculating probs"<<std::endl;
	for(int i = 0; i < 4; i++){
		std::vector<std::vector<double>> tmp_near;
		std::vector<std::vector<double>> tmp_far;

		for(int j = 0; j < 4; j++){
			std::vector<double> tmpen_near;
			std::vector<double> tmpen_far;

			s1<<i<<" "<<j<<" "<<" ";
			s2<<i<<" "<<j<<" "<<" ";
			for(double en = 0.001; en < 50; en+= 0.01){
				tmpen_far.push_back(prob->probabilityMatterExact(i,j,1,en,1300));//prob?
				tmpen_near.push_back(prob->probabilityVacuumExact(i,j,1,en,1.0));//yj near detector is too far. .525 km?
				s1<<tmpen_far.back()<<" ";
				s2<<tmpen_near.back()<<" ";
			}
			s1<<"\n";
			s2<<"\n";

			tmp_near.push_back(tmpen_near);
			tmp_far.push_back(tmpen_far);


		} 
		precalc_prob_far.push_back(tmp_far);	
		precalc_prob_near.push_back(tmp_near);	

	}
	std::cout<<"Starting precalculating probs BAR"<<std::endl;
	//prob->useAntiNeutrino = true;
	//prob->init();

	for(int i = 0; i < 4; i++){
		std::vector<std::vector<double>> tmp_nearbar;
		std::vector<std::vector<double>> tmp_farbar;

		for(int j = 0; j < 4; j++){
			std::vector<double> tmpen_nearbar;
			std::vector<double> tmpen_farbar;

			s3<<i<<" "<<j<<" "<<" ";
			s4<<i<<" "<<j<<" "<<" ";

			for(double en = 0.001; en < 50; en+= 0.01){
				tmpen_farbar.push_back(prob->probabilityMatterExact(i,j,-1,en,1300));
				tmpen_nearbar.push_back(prob->probabilityVacuumExact(i,j,-1,en,1.0));//yj near detector is too far. .525 km?

				s3<<tmpen_farbar.back()<<" ";
				s4<<tmpen_nearbar.back()<<" ";
			}
			s3<<"\n";
			s4<<"\n";



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
