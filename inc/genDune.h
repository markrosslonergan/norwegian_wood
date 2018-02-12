#ifndef GENDUNE_h
#define GENDUNE_H

#include <string>
#include <vector>
#include "SBgeN.h"
#include "SBNprob.h"

class genDUNE : public sbn::SBgeN {
	public:
		sbn::SBNprob * prob;
		genDUNE(std::string xmlname):SBgeN(xmlname){prob = new sbn::SBNprob(4);};
		

		bool eventSelection(int file);
		int fillHistograms(int file, int uni, double wei);
		int tidyHistograms();
		~genDUNE(){};


};
#endif
