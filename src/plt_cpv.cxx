#include "plotting_tools.h"

using namespace sbn;


int calc_cpv_3p1( std::ofstream* dunestream,  std::string outfile, std::string xml, double t14, double t24, double t34, TMatrixT<double> *m){
	bool is_verbose = false;

	if(is_verbose)std::cout<<"Opening stream"<<std::endl;
	dunestream->open(outfile.c_str());
	if(is_verbose)std::cout<<"Opened stream"<<std::endl;

	std::vector<double> angles = {33.6, 42.1, 8.5, 0,0,0};
	std::vector<double> angles_oct = {33.6, 49.6, 8.5, 0,0,0};
	std::vector<double> phases = {0,0,0};
	std::vector<double> phases_180 = {180,0,0};
	//http://lbne2-docdb.fnal.gov/cgi-bin/RetrieveFile?docid=10688&filename=DUNE-CDR-physics-volume.pdf&version=10
	std::vector<double> mass_splittings = {7.5*pow(10,-5), 2.457*pow(10,-3),0};
	std::vector<double> mass_splittings_inv = {7.5*pow(10,-5), -2.449*pow(10,-3),0};

	std::vector<std::string> order_names = {"NO","IO"};
	std::vector<double> order_vals = {2.457*pow(10,-3), -2.449*pow(10,-3)};

	std::vector<double> theta23 = {38,40,42,44,45,46,47,49,51,53};
	std::vector<double> theta34 = {0,10,15,20,25};

	std::vector<double> min_chi;

	std::vector<double> v_dcp = {0.0,180.0};
	std::vector<double> v_d14 = {0.0,180.0};

	for(double tru_dcp = 0; tru_dcp<=360; tru_dcp+=15){
		for(double tru_d14 = 0; tru_d14<360; tru_d14+=45){
			std::vector<double> chi_all;

			if(is_verbose) std::cout<<"True_dcp "<<tru_dcp<<" True_d14 "<<tru_d14<<std::endl;

			double temp_chi1 = 9999.;
			double temp_chi2 = 9999.;
			double temp_chi3 = 9999.;
			double temp_chi4 = 9999.;

			int tru_i23 = 2;

			std::string truth_name = order_names.at(0)+"_DCP_"+to_string_prec(tru_dcp,3)+"_T23_"+to_string_prec(theta23.at(3),3)+"_T14_"+to_string_prec(t14,3)+"_T24_"+to_string_prec(t24,3)+"_T34_"+to_string_prec(t34,3)+"_D14_"+to_string_prec(tru_d14,6);

			SBNspec * truth = new SBNspec(("/a/data/westside/yjwa/NW/DUNE_SBN_condor/condor_tests/"+truth_name+".SBNspec").c_str(),xml, is_verbose);
			truth->compressVector();
			std::cout << "assume truth : " <<truth_name << std::endl;
			SBNchi mychi(*truth,*m);


			for(int ord = 0; ord < order_names.size(); ord++){
				for(int i23 =0; i23 < theta23.size(); i23++){
					for(int i34 =0; i34 < theta34.size(); i34++){
						std::string name1 = order_names.at(ord)+"_DCP_"+to_string_prec(0.0,3)+"_T23_"+to_string_prec(theta23.at(i23),3)+"_T14_"+to_string_prec(t14,3)+"_T24_"+to_string_prec(t24,3)+"_T34_"+to_string_prec(theta34.at(i34),3)+"_D14_"+to_string_prec(0.0,6);
						std::string name2 = order_names.at(ord)+"_DCP_"+to_string_prec(0.0,3)+"_T23_"+to_string_prec(theta23.at(i23),3)+"_T14_"+to_string_prec(t14,3)+"_T24_"+to_string_prec(t24,3)+"_T34_"+to_string_prec(theta34.at(i34),3)+"_D14_"+to_string_prec( 180.0,6);
						std::string name3 = order_names.at(ord)+"_DCP_"+to_string_prec(180.0,3)+"_T23_"+to_string_prec(theta23.at(i23),3)+"_T14_"+to_string_prec(t14,3)+"_T24_"+to_string_prec(t24,3)+"_T34_"+to_string_prec(theta34.at(i34),3)+"_D14_"+to_string_prec(0.0,6);
						std::string name4 = order_names.at(ord)+"_DCP_"+to_string_prec(180.0,3)+"_T23_"+to_string_prec(theta23.at(i23),3)+"_T14_"+to_string_prec(t14,3)+"_T24_"+to_string_prec(t24,3)+"_T34_"+to_string_prec(theta34.at(i34),3)+"_D14_"+to_string_prec(180.0,6);

						SBNspec * test1 = new SBNspec(("/a/data/westside/yjwa/NW/DUNE_SBN_condor/condor_tests/"+name1+".SBNspec").c_str(),xml, is_verbose);
						SBNspec * test2 = new SBNspec(("/a/data/westside/yjwa/NW/DUNE_SBN_condor/condor_tests/"+name2+".SBNspec").c_str(),xml, is_verbose);
						SBNspec * test3 = new SBNspec(("/a/data/westside/yjwa/NW/DUNE_SBN_condor/condor_tests/"+name3+".SBNspec").c_str(),xml, is_verbose);
						SBNspec * test4 = new SBNspec(("/a/data/westside/yjwa/NW/DUNE_SBN_condor/condor_tests/"+name4+".SBNspec").c_str(),xml, is_verbose);

						test1->compressVector();
						test2->compressVector();
						test3->compressVector();
						test4->compressVector();

						double chi1 = mychi.CalcChi(test1);
						double chi2 = mychi.CalcChi(test2);
						double chi3 = mychi.CalcChi(test3);
						double chi4 = mychi.CalcChi(test4);


						if (chi1<temp_chi1){
							temp_chi1 = chi1;
						}

						if (chi2<temp_chi2){
							temp_chi2 = chi2;
						}

						if (chi3<temp_chi3){
							temp_chi3 = chi3;
						}

						if (chi4<temp_chi4){
							temp_chi4 = chi4;
						}
						delete test1;
						delete test2;
						delete test3;
						delete test4;
					}
				}
			}


			(*dunestream)<<"order_true_dcp "<<tru_dcp<<" "<<tru_d14<<" "<<temp_chi1<<" "<<temp_chi2<<" "<<temp_chi3<<" "<<temp_chi4<<" "<<t14<<" "<<t24<<" "<<t34<<" "<<std::endl;
		}
	}


	return 0;
}
