#ifndef GENDUNE_h
#define GENDUNE_H

#include <string>
#include <vector>
#include "SBgeN.h"

class genDUNE : public sbn::SBgeN {
	public:
		genDUNE(std::string xmlname):SBgeN(xmlname){};

		bool eventSelection(int file);
		int fillHistograms(int file, int uni, double wei);
		int tidyHistograms();
		~genDUNE(){};


};
#endif
