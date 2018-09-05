#include "plotting_tools.h"

using namespace sbn;


int calc_neutrino_ordering( std::ofstream* dunestream,  std::string outfile, std::string xml){
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

	TFile *fin = new TFile("/a/data/westside/yjwa/NW/norwegian_wood/covar/covariance_matrices_xcheck_1408x1408.root","read");
	TMatrixT<double> * m = (TMatrixT<double>*)fin->Get("TMatrixT<double>;1");

	std::vector<std::string> order_names = {"NO","IO"};
	std::vector<double> order_vals = {2.457*pow(10,-3), -2.449*pow(10,-3)};

	std::vector<double> theta23 = {38,40,42,44,45,46,47,49,51,53};

	std::vector<double> min_chi;

	for(double tru_dcp = 0; tru_dcp<=360; tru_dcp+=15){
		std::vector<double> chi_all;

		if(is_verbose) std::cout<<"True_dcp "<<tru_dcp<<std::endl;

		double temp_chi = 9999.;
		double temp_chi_sub = 9999.;

		int tru_i23 = 2;
		std::string truth_name = order_names.at(1)+"_DCP_"+to_string_prec(tru_dcp, 3)+"_T23_"+to_string_prec(theta23.at(tru_i23),3);
		if(is_verbose) std::cout<<"Loading precomputed "<<truth_name<<std::endl; 
		SBNspec * truth = new SBNspec(("precomp/"+truth_name+".SBNspec").c_str(),xml, is_verbose);
		truth->compressVector();
		std::cout << "assume truth : " <<truth_name << std::endl;
		SBNchi mychi(*truth,*m);

		for(double dcp=0; dcp <=360; dcp+=15){
			for(int i23 =0; i23 < theta23.size(); i23++){
				std::string name = order_names.at(1)+"_DCP_"+to_string_prec(dcp,3)+"_T23_"+to_string_prec(theta23.at(i23),3);  

				SBNspec * test = new SBNspec(("precomp/"+name+".SBNspec").c_str(),xml, is_verbose);

				double chi = mychi.CalcChi(test);
				if (chi<temp_chi){
					temp_chi = chi;
				}
				delete test;
			}
		}

		for(double dcp=0; dcp <=360; dcp+=15){
			for(int i23 =0; i23 < theta23.size(); i23++){
				std::string name_sub = order_names.at(0)+"_DCP_"+to_string_prec(dcp,3)+"_T23_"+to_string_prec(theta23.at(i23),3);

				SBNspec * test_sub = new SBNspec(("precomp/"+name_sub+".SBNspec").c_str(),xml, is_verbose);

				double chi_sub = mychi.CalcChi(test_sub);


				if (chi_sub<temp_chi_sub){
					temp_chi_sub = chi_sub;
				}
				delete test_sub;
			}
		}


		double delta_chi = temp_chi_sub- temp_chi;

		(*dunestream)<<"order true_dcp "<<tru_dcp<<" "<<delta_chi<<std::endl;
		//std::cout<<"order true_dcp "<<tru_dcp<<" "<<delta_chi <<std::endl;                
	}


	return 0;
}	

int calc_neutrino_ordering_3p1( std::ofstream* dunestream,  std::string outfile, std::string xml, double t14, double t24, double t34, double d14){
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

	TFile *fin = new TFile("/a/data/westside/yjwa/NW/norwegian_wood/covar/covariance_matrices_xcheck_1408x1408.root","read");
	TMatrixT<double> * m = (TMatrixT<double>*)fin->Get("TMatrixT<double>;1");

	std::vector<std::string> order_names = {"NO","IO"};
	std::vector<double> order_vals = {2.457*pow(10,-3), -2.449*pow(10,-3)};

	std::vector<double> theta23 = {38,40,42,44,45,46,47,49,51,53};

	std::vector<double> min_chi;

	for(double tru_dcp = 0; tru_dcp<=360; tru_dcp+=15){
		std::vector<double> chi_all;

		if(is_verbose) std::cout<<"True_dcp "<<tru_dcp<<std::endl;

		double temp_chi = 9999.;
		double temp_chi_sub = 9999.;

		int tru_i23 = 2;

		std::string truth_name = order_names.at(1)+"_DCP_"+to_string_prec(tru_dcp,3)+"_T23_"+to_string_prec(theta23.at(3),3)+"_T14_"+to_string_prec(t14,3)+"_T24_"+to_string_prec(t24,3)+"_T34_"+to_string_prec(t34,3)+"_D14_"+to_string_prec(d14,6);
                     SBNspec * truth = new SBNspec(("/a/data/westside/yjwa/NW/DUNE_SBN_condor/condor_tests/"+truth_name+".SBNspec").c_str(),xml, is_verbose);
		truth->compressVector();
		std::cout << "assume truth : " <<truth_name << std::endl;
		SBNchi mychi(*truth,*m);

		for(double dcp=0; dcp <=360; dcp+=15){
			for(int i23 =0; i23 < theta23.size(); i23++){
				std::string name = order_names.at(0)+"_DCP_"+to_string_prec(dcp,3)+"_T23_"+to_string_prec(theta23.at(i23),3);  

				SBNspec * test = new SBNspec(("/a/data/westside/markross/norwegian_wood/build/src/precomp/"+name+".SBNspec").c_str(),xml, is_verbose);

				double chi = mychi.CalcChi(test);
				if (chi<temp_chi){
					temp_chi = chi;
				}
				delete test;
			}
		}

		for(double dcp=0; dcp <=360; dcp+=15){
			for(int i23 =0; i23 < theta23.size(); i23++){
				std::string name_sub = order_names.at(1)+"_DCP_"+to_string_prec(dcp,3)+"_T23_"+to_string_prec(theta23.at(i23),3);

				SBNspec * test_sub = new SBNspec(("/a/data/westside/markross/norwegian_wood/build/src/precomp/"+name_sub+".SBNspec").c_str(),xml, is_verbose);

				double chi_sub = mychi.CalcChi(test_sub);


				if (chi_sub<temp_chi_sub){
					temp_chi_sub = chi_sub;
				}
				delete test_sub;
			}
		}


		double delta_chi = temp_chi_sub- temp_chi;

		(*dunestream)<<"order_true_dcp "<<tru_dcp<<" "<<delta_chi<<" "<<t14<<" "<<t24<<" "<<t34<<" "<<d14<<std::endl;
		//std::cout<<"order true_dcp "<<tru_dcp<<" "<<delta_chi <<std::endl;                
	}


	return 0;
}	
