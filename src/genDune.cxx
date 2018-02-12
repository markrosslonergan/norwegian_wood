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

	double oscprob = prob->probabilityMatterExact(2, 2 ,Enu_true, 1300);
	//std::cout<<Enu_true<<" "<<Enu_reco<<" on file: "<<multisim_name.at(file)<<" "<<nutype<<" PROB: "<<oscprob<<std::endl;
	//std::cout<<"Map Hist: "<<map_hist[multisim_name.at(file)]<<std::endl;
	
	//std::cout<<multisim_name.at(file)<<" "<<near_detector_names.at(file)<<std::endl;
	hist.at(map_hist[multisim_name.at(file)]).Fill(Enu_reco, oscprob*weight*far_detector_weight);
	hist.at(map_hist[near_detector_names.at(file)]).Fill(Enu_reco, oscprob*weight*near_detector_weight);

	return 0;
}

int genDUNE::tidyHistograms(){
	return 0;
}







