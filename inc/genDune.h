#ifndef GENDUNE_H
#define GENDUNE_H

#include <string>
#include <vector>
#include "SBgeN.h"
#include "SBNprob.h"



class genDUNE : public sbn::SBgeN {
    
    bool debug = false;
    public:
		sbn::SBNprob * prob;//yj SBNprob prob
		genDUNE(std::string xmlname):SBgeN(xmlname){
            
            if(debug) std::cout<<"Searching for memory leak 1"<<std::endl;

            
            
			prob = new sbn::SBNprob(4);
		
			near_detector_weight = 5.0;
			far_detector_weight = 1.0;
            
            if(debug) std::cout<<"Searching for memory leak 2"<<std::endl;


			for(int i = 0; i < multisim_name.size(); i++){
				std::string s = multisim_name.at(i);
				
                if(debug) std::cout<<s<<std::endl;
				
                std::string delimiter = "_";
				size_t pos = 0;
				std::string token;
				std::vector<std::string> breakdown_filename;

	
				if(debug) std::cout<<"Multisim: "<<s<<std::endl;
				while (true) {
					pos = s.find(delimiter);
					token = s.substr(0, pos);
					s.erase(0, pos + delimiter.length());
					breakdown_filename.push_back(token);
					if(pos == std::string::npos) break;
				}
				if(breakdown_filename.at(0) == "nu"){
					std::string near_name = "nu_near_"+breakdown_filename.at(2)+"_"+breakdown_filename.at(3);
					if(debug) std::cout<<"Mapped to: "<<near_name<<" @ "<<map_hist[near_name]<<std::endl;
					near_detector_name_map[ multisim_name.at(i)] = near_name;
					near_map_hist[ multisim_name.at(i)] = map_hist[near_name] ;
				} else if (breakdown_filename.at(0)=="nubar"){
					std::string near_name = "nubar_near_"+breakdown_filename.at(2)+"_"+breakdown_filename.at(3);
					if(debug) std::cout<<"Mapped to: "<<near_name<<" @ "<<map_hist[near_name]<<std::endl;
					near_detector_name_map[multisim_name.at(i)] = near_name;
					near_map_hist[multisim_name.at(i)] = map_hist[near_name]  ;			
				}


				if (map_hist.count(multisim_name.at(i))!=1){
					if(debug) std::cout<<"ERROR: one of the input multisim_tree's is not in the xml config: "<<multisim_name.at(i)<<std::endl;
					exit(EXIT_FAILURE);
				}


			}
            
            if(debug) std::cout<<"Searching for memory leak 3"<<std::endl;



			this->getOscPattern();
            
            if(debug) std::cout<<"Searching for memory leak 4"<<std::endl;


			for(int i=0; i< fullnames.size(); i++){
				oscillation_patterns.push_back(osc_pattern_map[fullnames.at(i)]);
				if(debug) std::cout<<"FRT: "<<i<<" "<<fullnames.at(i)<<" "<<oscillation_patterns.back().first<<" "<<oscillation_patterns.back().second<<std::endl;
			}
            
            if(debug) std::cout<<"Searching for memory leak 5"<<std::endl;

	

		}

		std::map<std::string, std::pair<int,int>> osc_pattern_map;

		std::vector<std::vector<std::vector<double>>> precalc_prob_farbar;
		std::vector<std::vector<std::vector<double>>> precalc_prob_nearbar;
	

		std::vector<std::vector<std::vector<double>>> precalc_prob_far;
		std::vector<std::vector<std::vector<double>>> precalc_prob_near;
	
		double interpolate_prob_far(int a, int b, double enu);
		double interpolate_prob_near(int a, int b, double enu);
		int preCalculateProbs();

		std::map<std::string, std::string> near_detector_name_map;
		std::map<std::string, int> near_map_hist;

		std::vector<std::string> near_detector_names;
		std::vector<std::pair<int,int>> oscillation_patterns;

		double near_detector_weight;
		double far_detector_weight;

		bool eventSelection(int file);
		int fillHistograms(int file, int uni, double wei);
		int tidyHistograms();
		~genDUNE(){};
	
		int getOscPattern();
		double lin_interp(double ein, double p1, double e1, double p2, double e2);

};
#endif
