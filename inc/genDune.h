#ifndef GENDUNE_h
#define GENDUNE_H

#include <string>
#include <vector>
#include "SBgeN.h"
#include "SBNprob.h"

class genDUNE : public sbn::SBgeN {
	public:
		sbn::SBNprob * prob;
		genDUNE(std::string xmlname):SBgeN(xmlname){
			prob = new sbn::SBNprob(4);
		
			near_detector_weight = 5.0;
			far_detector_weight = 1.0;

			for(int i = 0; i < multisim_name.size(); i++){
				std::string s = multisim_name.at(i);
				std::cout<<s<<std::endl;
				std::string delimiter = "_";
				size_t pos = 0;
				std::string token;
				std::vector<std::string> breakdown_filename;

				while (true) {
					pos = s.find(delimiter);
					token = s.substr(0, pos);
					s.erase(0, pos + delimiter.length());
					breakdown_filename.push_back(token);
					if(pos == std::string::npos) break;
				}
				std::string near_name = "nu_near_"+breakdown_filename.at(2)+"_"+breakdown_filename.at(3);
				near_detector_names.push_back(near_name);
			}
		}

		std::vector<std::string> near_detector_names;

		double near_detector_weight;
		double far_detector_weight;

		bool eventSelection(int file);
		int fillHistograms(int file, int uni, double wei);
		int tidyHistograms();
		~genDUNE(){};


};
#endif
