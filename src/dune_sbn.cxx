#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <cmath>
#include <vector>
#include <unistd.h>
#include <getopt.h>
#include <cstring>
#include <algorithm>
#include <iterator>

#include "TFile.h"
#include "TTree.h"
#include "TH1F.h"
#include "TString.h"
#include "TNtuple.h"
#include "TChain.h"
#include "TMath.h"
#include "TSystem.h"
#include "TMatrixT.h"
#include "TRandom.h"
#include "TError.h"
#include "TCanvas.h"
#include "TH2F.h"
#include "TGraph.h"

#include "params.h"
#include "SBNconfig.h"
#include "SBNchi.h"
#include "SBNdet.h"
#include "SBNspec.h"
#include "SBNosc.h"
#include "SBgeN.h"
#include "SBNfit.h"
#include "SBNfit3pN.h"



#include "genDune.h"

#define no_argument 0
#define required_argument 1
#define optional_argument 2

using namespace sbn;


	template <typename T>
std::string to_string_prec(const T a_value, const int n = 6)
{
	std::ostringstream out;
	out <<std::fixed<< std::setprecision(n) << a_value;
	return out.str();
}






/*************************************************************
 *************************************************************
 *		BEGIN example.cxx
 ************************************************************
 ************************************************************/
int main(int argc, char* argv[])
{


	std::string xml = "../../xml/dune.xml";
	int iarg = 0;
	opterr=1;
	int index; 
	int test_mode=0;
	std::string filename = "default.root";
	std::string which_mode = "default";


	TRandom3 *rangen= new TRandom3();
	/*************************************************************
	 *************************************************************
	 *		Command Line Argument Reading
	 ************************************************************
	 ************************************************************/

	const struct option longopts[] = 
	{
		{"xml", 		required_argument, 	0, 'x'},
		{"test",		required_argument,	0, 't'},
		{"mode",		required_argument,	0, 'm'},
		{"file",		required_argument,	0, 'f'},
		{0,			no_argument, 		0,  0},
	};


	while(iarg != -1)
	{
		iarg = getopt_long(argc,argv, "x:m:t:f:", longopts, &index);

		switch(iarg)
		{
			case 'x':
				xml = optarg;
				break;
			case 'f':
				filename = optarg;//`strtof(optarg,NULL);
				break;
			case 'm':
				which_mode = optarg;//`strtof(optarg,NULL);
				break;

			case 't':
				test_mode = strtof(optarg,NULL);
				break;
			case '?':
			case 'h':
				std::cout<<"Allowed arguments:"<<std::endl;
				std::cout<<"\t-x\t--xml\t\tInput .xml file for SBNconfig"<<std::endl;
				return 0;
		}

	}

	//SBNspec * dune_spec = new SBNspec("~/work/pheno/DUNE+SBN/yeonjae/sb_macros/DUNE_bf"  , xml);
	//dune_spec->writeOut("dune_test.root");

	if(which_mode == "default"){
		//*************************************************************************//
		//*************************************************************************//
		std::cout<<" Starting Default Mode: "<<std::endl;
		//*************************************************************************//
		//*************************************************************************//
		std::vector<double> angles = {30, 44, 8, 0,0,0};
		std::vector<double> phases = {0,0,0};
		std::vector<double> mass_splittings = {7.5*pow(10,-5), 2.552*pow(10,-3),0};
		genDUNE bkg_only(xml);

		bkg_only.prob = new SBNprob(4, angles, phases, mass_splittings);
		bkg_only.preCalculateProbs();

		bkg_only.doMC("three_neutrino");
		bkg_only.writeOut("three_neutrino.root");

		return 0;



	} else if(which_mode == "dcp"){
		//*************************************************************************//
		//*************************************************************************//
		std::cout<<" Starting Default Mode: "<<std::endl;
		//*************************************************************************//
		//*************************************************************************//
		genDUNE base(xml);


		//in degrees! t12, t23, t13 
		std::vector<double> angles = {33.62, 42.1304, 8.4750, 0,0,0};
		std::vector<double> angles_oct = {33.62, 49.6, 8.4750, 0,0,0};
		std::vector<double> phases = {0,0,0};
		std::vector<double> phases_180 = {180.0,0,0};
		//http://lbne2-docdb.fnal.gov/cgi-bin/RetrieveFile?docid=10688&filename=DUNE-CDR-physics-volume.pdf&version=10
		std::vector<double> mass_splittings = {7.5*pow(10,-5), 2.457*pow(10,-3),0};
		std::vector<double> mass_splittings_inv = {7.5*pow(10,-5), -2.449*pow(10,-3),0};


		TFile *fin = new TFile("/home/mark/work/pheno/DUNE+SBN/covar/covariance_matrices_xcheck_1408x1408.root","read");
		TMatrixT<double> * m = (TMatrixT<double>*)fin->Get("TMatrixT<double>;1");

		std::vector<genDUNE*> test_points;
		std::vector<SBNchi*> test_chis;
		std::vector<double> theta23_testpts = {40,42.1304,45,49.6,50};

		for(int i=0; i<theta23_testpts.size(); i++){
			std::vector<double> angles_test = {33.62, theta23_testpts.at(i), 8.4750,0,0,0};

			test_points.push_back(new genDUNE(xml));
			test_points.back()->prob = new SBNprob(4,angles_test,phases,mass_splittings);	
			test_points.back()->preCalculateProbs();	
			std::cout<<"DUNE_SBN:: "<<"Running MC for  NO dcp 0 theta23 "<<theta23_testpts.at(i)<<std::endl;
			test_points.back()->doMC("NO_dcp0_"+std::to_string(theta23_testpts.at(i)));
			test_chis.push_back(new SBNchi(*test_points.back(),*m));				


			test_points.push_back(new genDUNE(xml));
			test_points.back()->prob = new SBNprob(4,angles_test,phases_180,mass_splittings);	
			test_points.back()->preCalculateProbs();	
			std::cout<<"DUNE_SBN:: "<<"Running MC for  NO dcp 180 theta23 "<<theta23_testpts.at(i)<<std::endl;
			test_points.back()->doMC("NO_dcp180_"+std::to_string(theta23_testpts.at(i)));
			test_chis.push_back(new SBNchi(*test_points.back(),*m));				

			test_points.push_back(new genDUNE(xml));
			test_points.back()->prob = new SBNprob(4,angles_test,phases,mass_splittings_inv);	
			test_points.back()->preCalculateProbs();	
			std::cout<<"DUNE_SBN:: "<<"Running MC for  I0 dcp 0 theta23 "<<theta23_testpts.at(i)<<std::endl;
			test_points.back()->doMC("IO_dcp0_"+std::to_string(theta23_testpts.at(i)));
			test_chis.push_back(new SBNchi(*test_points.back(),*m));				

			test_points.push_back(new genDUNE(xml));
			test_points.back()->prob = new SBNprob(4,angles_test,phases_180,mass_splittings_inv);	
			test_points.back()->preCalculateProbs();	
			std::cout<<"DUNE_SBN:: "<<"Running MC for  IP dcp 180 theta23 "<<theta23_testpts.at(i)<<std::endl;
			test_points.back()->doMC("IO_dcp180_"+std::to_string(theta23_testpts.at(i)));
			test_chis.push_back(new SBNchi(*test_points.back(),*m));				
		}

		test_points.at(4)->compareSBNspecs(test_points.at(5),"comptets.root");

		//NO low octant
		/*	genDUNE dcp0   = base;
			genDUNE dcp180 = base;
			dcp0.prob = new SBNprob(4, angles, phases, mass_splittings);
			dcp0.preCalculateProbs();
			phases.at(0) = 180;
			dcp180.prob = new SBNprob(4, angles, phases, mass_splittings);
			dcp180.preCalculateProbs();

			std::cout<<"DUNE_SBN:: "<<"Running MC for NO LO 0"<<std::endl;
			dcp0.doMC("gentest0");
			std::cout<<"DUNE_SBN:: "<<"Running MC for NO LO 180"<<std::endl;
			dcp180.doMC("gentest180");
			std::cout<<"DUNE_SBN:: "<<"Constructing a chi^2 for NO LO 0"<<std::endl;
			SBNchi mychi0(dcp0,*m);
			std::cout<<"DUNE_SBN:: "<<"Constructing a chi^2 for NO LO 180"<<std::endl;
			SBNchi mychi180(dcp180,*m);


		//NO high octant
		genDUNE dcp0_o= base;
		genDUNE dcp180_o= base;
		phases.at(0) = 0;
		dcp0_o.prob = new SBNprob(4, angles_oct, phases, mass_splittings);
		dcp0_o.preCalculateProbs();
		phases.at(0) = 180;
		dcp180_o.prob = new SBNprob(4, angles_oct, phases, mass_splittings);
		dcp180_o.preCalculateProbs();

		std::cout<<"DUNE_SBN:: "<<"Running MC for NO HO 0"<<std::endl;
		dcp0_o.doMC("gentest0_o");
		std::cout<<"DUNE_SBN:: "<<"Running MC for NO HO 180"<<std::endl;
		dcp180_o.doMC("gentest180_o");
		std::cout<<"DUNE_SBN:: "<<"Constructing a chi^2 for NO HO 0"<<std::endl;
		SBNchi mychi0_o(dcp0_o,*m);
		std::cout<<"DUNE_SBN:: "<<"Constructing a chi^2 for NO HO 180"<<std::endl;
		SBNchi mychi180_o(dcp180_o,*m);


		//IO low octant
		genDUNE dcp0_i= base;
		genDUNE dcp180_i= base;
		phases.at(0) = 0;
		dcp0_i.prob = new SBNprob(4, angles, phases, mass_splittings_inv);
		dcp0_i.preCalculateProbs();
		phases.at(0) = 180;
		dcp180_i.prob = new SBNprob(4, angles, phases, mass_splittings_inv);
		dcp180_i.preCalculateProbs();

		std::cout<<"DUNE_SBN:: "<<"Running MC for IO LO 0"<<std::endl;
		dcp0_i.doMC("gentest0_i");
		std::cout<<"DUNE_SBN:: "<<"Running MC for IO LO 180"<<std::endl;
		dcp180_i.doMC("gentest180_i");
		std::cout<<"DUNE_SBN:: "<<"Constructing a chi^2 for IO LO 0"<<std::endl;
		SBNchi mychi0_i(dcp0_i,*m);
		std::cout<<"DUNE_SBN:: "<<"Constructing a chi^2 for IO LO 180"<<std::endl;
		SBNchi mychi180_i(dcp180_i,*m);



		//IO high octant
		genDUNE dcp0_i_o= base;
		genDUNE dcp180_i_o= base;
		phases.at(0) = 0;
		dcp0_i_o.prob = new SBNprob(4, angles_oct, phases, mass_splittings_inv);
		dcp0_i_o.preCalculateProbs();
		phases.at(0) = 180;
		dcp180_i_o.prob = new SBNprob(4, angles_oct, phases, mass_splittings_inv);
		dcp180_i_o.preCalculateProbs();

		std::cout<<"DUNE_SBN:: "<<"Running MC for IO HO 0"<<std::endl;
		dcp0_i_o.doMC("gentest0_i_o");
		std::cout<<"DUNE_SBN:: "<<"Running MC for IO HO 180"<<std::endl;
		dcp180_i_o.doMC("gentest180_i_o");
		std::cout<<"DUNE_SBN:: "<<"Constructing a chi^2 for IO HO 0"<<std::endl;
		SBNchi mychi0_i_o(dcp0_i_o,*m);
		std::cout<<"DUNE_SBN:: "<<"Constructing a chi^2 for IO HO 180"<<std::endl;
		SBNchi mychi180_i_o(dcp180_i_o,*m);
		*/



			std::ofstream dunestream;
		dunestream.open ("DUNE_dcp_IO.dat");
		bool print = true;

		for(double dcp = 0; dcp <= 360; dcp += 10.0){

			std::cout<<"DUNE_SBN:: "<<"On dcp "<<dcp<<std::endl;
			genDUNE *testGen = new genDUNE(xml);
			phases.at(0) = dcp; 

			testGen->prob = new SBNprob(4, angles, phases, mass_splittings);	
			testGen->preCalculateProbs();

			std::cout<<"DUNE_SBN:: "<<"Doing MC on dcp "<<dcp<<std::endl;
			testGen->doMC();
			testGen->compressVector();

			//SBNchi mychi(*testGen, *m);


			std::vector<double> chis;
			std::vector<double> chis_old;
			for(int i=0; i<test_points.size(); i++){
				test_points.at(i)->compressVector();
				chis.push_back( test_chis.at(i)->CalcChi(testGen));
				//TH1D hist_test = (TH1D) mychi.toyMC_varyInput(test_points.at(i), 1000);
				//TH1D hist_truth = (TH1D) mychi.toyMC_varyInput(testGen, 1000);

				//double ans = hist_test.GetMean()-hist_truth.GetMean();
				//if(ans<0)ans=0;
				//chis.push_back(ans);


			}
			/*std::cout<<"DUNE_SBN:: "<<"Calc chi^2  NOLO 180 on dcp "<<dcp<<std::endl;
			  double chi180 = mychi.CalcChi(&dcp180);
			  double chi0 = mychi.CalcChi(&dcp0);

			  std::cout<<"DUNE_SBN:: "<<"Calc chi^2  IOLO 180 on dcp "<<dcp<<std::endl;
			  double chi180_i = mychi.CalcChi(&dcp180_i);
			  double chi0_i = mychi.CalcChi(&dcp0_i);

			  std::cout<<"DUNE_SBN:: "<<"Calc chi^2  NOHO 180 on dcp "<<dcp<<std::endl;
			  double chi180_o = mychi.CalcChi(&dcp180_o);
			  double chi0_o = mychi.CalcChi(&dcp0_o);

			  std::cout<<"DUNE_SBN:: "<<"Calc chi^2  IOHO 180 on dcp "<<dcp<<std::endl;
			  double chi180_i_o = mychi.CalcChi(&dcp180_i_o);
			  double chi0_i_o = mychi.CalcChi(&dcp0_i_o);


			  std::cout<<"DUNE_SBN:: "<<"Calc chi^2  NOLO 180 on dcp "<<dcp<<std::endl;
			  double chi180 = mychi180.CalcChi(testGen);
			  double chi0 = mychi0.CalcChi(testGen);

			  std::cout<<"DUNE_SBN:: "<<"Calc chi^2  IOLO 180 on dcp "<<dcp<<std::endl;
			  double chi180_i = mychi180_i.CalcChi(testGen);
			  double chi0_i = mychi0_i.CalcChi(testGen);

			  std::cout<<"DUNE_SBN:: "<<"Calc chi^2  NOHO 180 on dcp "<<dcp<<std::endl;
			  double chi180_o = mychi180_o.CalcChi(testGen);
			  double chi0_o = mychi0_o.CalcChi(testGen);

			  std::cout<<"DUNE_SBN:: "<<"Calc chi^2  IOHO 180 on dcp "<<dcp<<std::endl;
			  double chi180_i_o = mychi180_i_o.CalcChi(testGen);
			  double chi0_i_o = mychi0_i_o.CalcChi(testGen);
			 */




			//std::vector<double> chis = {chi0, chi0_i, chi0_o, chi0_i_o, chi180, chi180_i, chi180_o, chi180_i_o};
			double chimin = 99999999;

			dunestream<<"ANK: "<<dcp;//" "<<chi0<<" "<<chi180<<" "<<std::min(chi0,chi180)<<std::endl;
			int l=0;
			int which = 0;
			for(auto &X: chis){
				dunestream<<" "<<X;
				if(X<chimin) {
					chimin = X;
					which = l;
				}
			}
			dunestream<<" MIN: "<<chimin<<std::endl;

			//Very important, root only likes so many files open.
			/*
			   if(chimin > 180 && print){
			   testGen->compareSBNspecs(test_points.at(which),std::to_string(dcp)+"_comp.root");
			   TFile *fchi = new TFile("chi.root","recreate");
			   fchi->cd();
			   TH2D cc =(TH2D)mychi.getChiogram();
			   cc.Write();
			   fchi->Close();
			   print = false;

			   }
			 */
			double vi=0;
			for(int i=0; i<testGen->num_bins_total_compressed; i++){
				//for(int j=0; j<testGen->num_bins_total_compressed; j++){
				//	std::cout<<"CHIOGRAM: "<<i<<" "<<j<<" "<<mychi.lastChi_vec.at(i).at(j)<<std::endl;
				//double vag = fabs(dcp0.compVec.at(i)-testGen->compVec.at(i))/sqrt(mychi.vMc.at(i).at(i));
				//std::cout<<"CHIO "<<i<<" "<<dcp0.compVec.at(i)<<" : "<<testGen->compVec.at(i)<<" +/- "<<sqrt(mychi.vMc.at(i).at(i))<<" "<<sqrt(testGen->compVec.at(i))<<"\t\t"<<vag<<std::endl;
				//vi+= vag;

				//}
			}	
			//std::cout<<"Total: "<<vi<<std::endl;

			delete testGen;
		}

		dunestream.close();

		return 0;

	}else if(which_mode == "mctest"){

		std::vector<double> angles = {33.62, 42.1304, 8.4750, 0,0,0};
		std::vector<double> angles_oct = {33.62, 49.6, 8.4750, 0,0,0};
		std::vector<double> phases = {0,0,0};
		std::vector<double> phases_180 = {180.0,0,0};
		//http://lbne2-docdb.fnal.gov/cgi-bin/RetrieveFile?docid=10688&filename=DUNE-CDR-physics-volume.pdf&version=10
		std::vector<double> mass_splittings = {7.5*pow(10,-5), 2.457*pow(10,-3),0};
		std::vector<double> mass_splittings_inv = {7.5*pow(10,-5), -2.449*pow(10,-3),0};

		TFile *fin = new TFile("/home/mark/work/pheno/DUNE+SBN/covar/covariance_matrices_xcheck_1408x1408.root","read");
		TMatrixT<double> * m = (TMatrixT<double>*)fin->Get("TMatrixT<double>;1");


		genDUNE* truth = new genDUNE(xml);
		phases.at(0)=10;
		truth->prob = new SBNprob(4,angles,phases,mass_splittings);
		truth->preCalculateProbs();		
		truth->doMC("MC_truth_cv");

		genDUNE* test0 = new genDUNE(xml);
		phases.at(0)=0;
		test0->prob = new SBNprob(4,angles,phases,mass_splittings);
		test0->preCalculateProbs();		
		test0->doMC("MC_test0_cv");

		genDUNE* test180 = new genDUNE(xml);	
		phases.at(0)= 180;
		test180->prob = new SBNprob(4,angles,phases,mass_splittings);
		test180->preCalculateProbs();		
		test180->doMC("MC_test180_cv");


		SBNchi mychi0(*truth,*m);	
		SBNchi mychi180(*truth,*m);	

		TFile *fp = new TFile("test.root","recreate");
		fp->cd();
		TCanvas *c = new TCanvas("please_work","please_work");
		c->Divide(2,1);

		TH1D hist_truth_0 = (TH1D) mychi0.toyMC_varyInput(truth, 2000);
		TH1D hist_truth_180 = (TH1D) mychi180.toyMC_varyInput(truth, 2000);

		TH1D hist_test_0 = (TH1D) mychi0.toyMC_varyInput(test0, 2000);
		TH1D hist_test_180 = (TH1D) mychi180.toyMC_varyInput(test180, 2000);

		hist_truth_0.Write();	
		hist_truth_180.Write();	
		hist_test_0.Write();	
		hist_test_180.Write();	


		c->cd(1);
		hist_truth_0.SetLineColor(kRed);
		hist_test_0.SetLineColor(kBlue);

		hist_truth_0.Draw("hist");
		hist_test_0.Draw("hist same");
		hist_truth_0.SetMaximum(std::max( hist_truth_0.GetMaximum(), hist_test_0.GetMaximum())*1.1);


		c->cd(2);
		hist_truth_180.SetLineColor(kRed);
		hist_test_180.SetLineColor(kBlue);

		hist_truth_180.Draw("hist");
		hist_test_180.Draw("hist same");
		hist_truth_180.SetMaximum(std::max( hist_truth_180.GetMaximum(), hist_test_180.GetMaximum())*1.1);


		c->Write();	

		fp->Close();


	}else if(which_mode=="gen"){

		std::vector<double> angles = {33.6, 42.1, 8.5, 0,0,0};
		std::vector<double> angles_oct = {33.6, 49.6, 8.5, 0,0,0};
		std::vector<double> phases = {0,0,0};
		std::vector<double> phases_180 = {180,0,0};
		//http://lbne2-docdb.fnal.gov/cgi-bin/RetrieveFile?docid=10688&filename=DUNE-CDR-physics-volume.pdf&version=10
		std::vector<double> mass_splittings = {7.5*pow(10,-5), 2.457*pow(10,-3),0};
		std::vector<double> mass_splittings_inv = {7.5*pow(10,-5), -2.449*pow(10,-3),0};


		TFile *fin = new TFile("/home/mark/work/pheno/DUNE+SBN/covar/covariance_matrices_xcheck_1408x1408.root","read");
		TMatrixT<double> * m = (TMatrixT<double>*)fin->Get("TMatrixT<double>;1");

		std::vector<std::string> order_names = {"NO","IO"};
		std::vector<double> order_vals = {2.457*pow(10,-3), -2.449*pow(10,-3)};

		std::vector<double> theta23 = {40,41,42,43,44,45,46,47,48,49,50,51,52};
		std::vector<double> theta23_NO = {42.1-2,42.1-1, 42.1, 42.1+1,42.1+2};
		std::vector<double> theta23_IO = {49.6-2,49.6-1, 49.6, 49.6+1,49.6+2};

		for(double dcp = 0; dcp<=360; dcp+=5){
			for(int ord = 0; ord<2; ord++){
				for(int i23 =0; i23 < theta23.size(); i23++){

					std::string name = order_names.at(ord)+"_DCP_"+to_string_prec(dcp,3)+"_T23_"+to_string_prec(theta23.at(i23),3);
					std::cout<<"GEN: on value "<<name<<std::endl;

					genDUNE* testpt = new genDUNE(xml);
					phases.at(0) = dcp;
					angles.at(1) = theta23.at(i23);
					mass_splittings.at(1) = order_vals.at(ord);

					testpt->prob = new SBNprob(4,angles,phases, mass_splittings);
					testpt->preCalculateProbs();		

					testpt->doMC("precomp/"+name);

					delete testpt;
				}
			}	
		}
	}else if(which_mode =="cpv"){
		std::ofstream dunestream;
		dunestream.open ("DUNE_cpv2.dat");

		std::vector<double> angles = {33.6, 42.1, 8.5, 0,0,0};
		std::vector<double> angles_oct = {33.6, 49.6, 8.5, 0,0,0};
		std::vector<double> phases = {0,0,0};
		std::vector<double> phases_180 = {180,0,0};
		//http://lbne2-docdb.fnal.gov/cgi-bin/RetrieveFile?docid=10688&filename=DUNE-CDR-physics-volume.pdf&version=10
		std::vector<double> mass_splittings = {7.5*pow(10,-5), 2.457*pow(10,-3),0};
		std::vector<double> mass_splittings_inv = {7.5*pow(10,-5), -2.449*pow(10,-3),0};

		TFile *fin = new TFile("/home/mark/work/pheno/DUNE+SBN/covar/covariance_matrices_xcheck_1408x1408.root","read");
		TMatrixT<double> * m = (TMatrixT<double>*)fin->Get("TMatrixT<double>;1");

		std::vector<std::string> order_names = {"NO","IO"};
		std::vector<double> order_vals = {2.457*pow(10,-3), -2.449*pow(10,-3)};

		std::vector<double> theta23 = {40,41,42,43,44,45,46,47,48,49,50,51,52};
		std::vector<double> theta23_NO = {42.1-2,42.1-1, 42.1, 42.1+1,42.1+2};
		std::vector<double> theta23_IO = {49.6-2,49.6-1, 49.6, 49.6+1,49.6+2};


		for(double tru_dcp = 0; tru_dcp<=360; tru_dcp+=5){

			std::string truth_name = order_names.at(0)+"_DCP_"+to_string_prec(tru_dcp,3)+"_T23_"+to_string_prec(theta23.at(10),3);
			//std::string truth_name = order_names.at(1)+"_DCP_"+to_string_prec(tru_dcp,3)+"_T23_"+to_string_prec(theta23.at(8),3);
			SBNspec * truth = new SBNspec(("precomp/"+truth_name+".SBNspec").c_str(),xml);
			truth->compressVector();	

			SBNchi mychi(*truth,*m);	
			std::vector<double> chi_all;
			std::vector<double> chi_0pi;

			//Here, vary the truth by poissonian_noise and build a distributuon for this. This is p-d/d then

			for(double dcp = 0; dcp<=360; dcp+=5){
				for(int ord = 0; ord<2; ord++){
					for(int i23 =0; i23 < theta23.size(); i23++){

						std::string name = order_names.at(ord)+"_DCP_"+to_string_prec(dcp,3)+"_T23_"+to_string_prec(theta23.at(i23),3);

						SBNspec * test = new SBNspec(("precomp/"+name+".SBNspec").c_str(),xml);	

						double chi = mychi.CalcChi(test);

						if(dcp<1 || (dcp < 181 && dcp > 179) ){
							chi_0pi.push_back(chi);
						}	

						chi_all.push_back(chi);
						delete test;
					}
				}
			}	

			double min_chi_all = 1e20; //std::min_element(chi_all.begin(), chi_all.end());	
			double min_chi_0pi = 1e20; //std::min_element(chi_0pi.begin(), chi_0pi.end());	
			for(auto &X: chi_all){
				if( X< min_chi_all) min_chi_all =X;
			}		
			for(auto &X: chi_0pi){
				if( X< min_chi_0pi) min_chi_0pi =X;
			}		


			if(min_chi_all>min_chi_0pi) std::cout<<"ERROR WARBING WARRRRBBBBING! global min is bigger than 0-pi min"<<std::endl;

			dunestream<<"CPV true_dcp "<<tru_dcp<<" "<<min_chi_0pi<<" "<<min_chi_all<<" "<<min_chi_0pi-min_chi_all<<std::endl;

		}
	}if(which_mode =="detail"){
		//std::ofstream dunestream;
		//dunestream.open ("DUNE_cpv_IO.dat");

		std::vector<double> angles = {33.6, 42.1, 8.5, 0,0,0};
		std::vector<double> angles_oct = {33.6, 49.6, 8.5, 0,0,0};
		std::vector<double> phases = {0,0,0};
		std::vector<double> phases_180 = {180,0,0};
		//http://lbne2-docdb.fnal.gov/cgi-bin/RetrieveFile?docid=10688&filename=DUNE-CDR-physics-volume.pdf&version=10
		std::vector<double> mass_splittings = {7.5*pow(10,-5), 2.457*pow(10,-3),0};
		std::vector<double> mass_splittings_inv = {7.5*pow(10,-5), -2.449*pow(10,-3),0};

		TFile *fin = new TFile("/home/mark/work/pheno/DUNE+SBN/covar/covariance_matrices_xcheck_1408x1408.root","read");
		TMatrixT<double> * m = (TMatrixT<double>*)fin->Get("TMatrixT<double>;1");

		std::vector<std::string> order_names = {"NO","IO"};
		std::vector<double> order_vals = {2.457*pow(10,-3), -2.449*pow(10,-3)};

		std::vector<double> theta23_NO = {42.1-2,42.1-1, 42.1, 42.1+1,42.1+2};
		std::vector<double> theta23_IO = {49.6-2,49.6-1, 49.6, 49.6+1,49.6+2};

		std::vector<double> theta23 = {42.1-2,42.1-1, 42.1, 42.1+1,42.1+2,49.6-2,49.6-1, 49.6, 49.6+1,49.6+2};

		//pick one true dcp
		double tru_dcp = 90;

		std::string truth_name = order_names.at(0)+"_DCP_"+to_string_prec(tru_dcp,3)+"_T23_"+to_string_prec(theta23_NO.at(1),3);
		SBNspec * truth = new SBNspec(("precomp/"+truth_name+".SBNspec").c_str(),xml);
		truth->compressVector();	

		SBNchi mychi(*truth,*m);	
		std::vector<double> chi_all;
		std::vector<double> chi_0pi;





		//Here, vary the truth by poissonian_noise and build a distributuon for this. This is p-d/d then
		std::vector<SBNspec> MCspecs; 
		TFile *fout= new TFile("DUNE_detail.root","recreate");	
		fout->cd();

		TH1D hist_truth("MCtruth","MCtruth",100,0,1000);
		TH1D hist_test0("MCtest0","MCtest0",100,0,1000);
		TH1D hist_test180("MCtest180","MCtest180",100,0,1000);
		TH1D hist_test("MCtest","MCtest",100,0,1000);

		int NMC=100;

		std::vector<double> idcp = {0,180};

		for(int i=0; i<=NMC; i++){
			SBNspec tmp = *truth;
			tmp.poissonScale(rangen);
			SBNchi tmpChi(tmp,*m);
			double thischi = tmpChi.CalcChi(truth);
			hist_truth.Fill(thischi);


			for(int ord = 0; ord<2; ord++){
				for(int i23 =0; i23 < theta23_NO.size(); i23++){
					std::vector<double> the;
					if(ord==0) the=theta23_NO;
					if(ord==1) the=theta23_IO;

					std::string name0 = order_names.at(ord)+"_DCP_"+to_string_prec(0.0,3)+"_T23_"+to_string_prec(the.at(i23),3);
					std::string name180 = order_names.at(ord)+"_DCP_"+to_string_prec(180.0,3)+"_T23_"+to_string_prec(the.at(i23),3);

					SBNspec * test0 = new SBNspec(("precomp/"+name0+".SBNspec").c_str(),xml);	
					SBNspec * test180 = new SBNspec(("precomp/"+name180+".SBNspec").c_str(),xml);	

					double chi0 = tmpChi.CalcChi(test0);
					double chi180 = tmpChi.CalcChi(test180);

					hist_test0.Fill(chi0);
					hist_test180.Fill(chi180);
					hist_test.Fill(std::min(chi0,chi180));

					delete test0;
					delete test180;
				}
			}	
		}

		fout->cd();
		hist_truth.Write();
		hist_test.Write();
		hist_test0.Write();
		hist_test180.Write();
		fout->Close();
	}else if(which_mode =="compare"){




		std::vector<double> angles = {33.6, 42.1, 8.5, 0,0,0};
		std::vector<double> angles_oct = {33.6, 49.6, 8.5, 0,0,0};
		std::vector<double> phases = {0,0,0};
		std::vector<double> phases_180 = {180,0,0};
		std::vector<double> mass_splittings = {7.5*pow(10,-5), 2.457*pow(10,-3),0};
		std::vector<double> mass_splittings_inv = {7.5*pow(10,-5), -2.449*pow(10,-3),0};

		TFile *fin = new TFile("/home/mark/work/pheno/DUNE+SBN/covar/covariance_matrices_xcheck_1408x1408.root","read");
		TMatrixT<double> * m = (TMatrixT<double>*)fin->Get("TMatrixT<double>;1");

		std::vector<std::string> order_names = {"NO","IO"};
		std::vector<double> order_vals = {2.457*pow(10,-3), -2.449*pow(10,-3)};

		std::vector<double> theta23 = {40,41,42,43,44,45,46,47,48,49,50,51,52};
		std::vector<double> theta23_NO = {42.1-2,42.1-1, 42.1, 42.1+1,42.1+2};
		std::vector<double> theta23_IO = {49.6-2,49.6-1, 49.6, 49.6+1,49.6+2};


		std::string truth_name = order_names.at(0)+"_DCP_"+to_string_prec(0.0,3)+"_T23_"+to_string_prec(theta23.at(1),3);
		SBNspec * truth = new SBNspec(("precomp/"+truth_name+".SBNspec").c_str(),xml);
		truth->compressVector();	


		truth->writeOut("3p1_cp0.root");

		std::string reco_name = order_names.at(0)+"_DCP_"+to_string_prec(90.0,3)+"_T23_"+to_string_prec(theta23.at(1),3);
		SBNspec * reco = new SBNspec(("precomp/"+reco_name+".SBNspec").c_str(),xml);
		reco->compressVector();	



		//Attempt a sterile oscillation
		std::vector<double> angles3p1 = {33.6, 41, 8.5, 8.6,9.9,0};
		std::vector<double> phases3p1 = {0,0,0};
		std::vector<double> mass_splittings3p1 = {7.5*pow(10,-5), 2.457*pow(10,-3),0.92};

		genDUNE * sterile = new genDUNE(xml);	
		sterile->prob = new SBNprob(4,angles3p1,phases3p1, mass_splittings3p1);
		sterile->preCalculateProbs();		
		sterile->doMC();

		SBNchi mychi(*truth,*m);	


		std::cout<<"Compare Chi: "<<mychi.CalcChi(sterile)<<std::endl;
		truth->compareSBNspecs(sterile,"DUNE_compare.root");

		//		std::cout<<"Compare Chi: "<<mychi.CalcChi(reco)<<std::endl;
		//		truth->compareSBNspecs(reco,"DUNE_compare.root");


		truth->writeSpec("comp");

		TFile *fchi = new TFile("DUNE_compare.root","update");
		fchi->cd();
		TH2D cc =(TH2D)mychi.getChiogram();
		cc.Write();
		fchi->Close();




	}else if(which_mode =="scan"){
		std::ofstream dunestream;
		dunestream.open ("DUNE_sterilescan.dat");

		TFile *fin = new TFile("/home/mark/work/pheno/DUNE+SBN/covar/covariance_matrices_xcheck_1408x1408.root","read");
		TMatrixT<double> * m = (TMatrixT<double>*)fin->Get("TMatrixT<double>;1");
		std::vector<double> angles = {33.6, 42.1, 8.5, 0,0,0};
		std::vector<double> phases = {0,0,0};
		std::vector<double> mass_splittings = {7.5*pow(10,-5), 2.457*pow(10,-3),0};

	
		std::vector<std::string> order_names = {"NO","IO"};
		std::vector<double> order_vals = {2.457*pow(10,-3), -2.449*pow(10,-3)};
		std::vector<double> theta23 = {40,41,42,43,44,45,46,47,48,49,50,51,52};

		std::vector<double> angles3p1 = {33.6, 41, 8.5, 8.6,9.9,0};
		std::vector<double> phases3p1 = {0,0,0};
		std::vector<double> mass_splittings3p1 = {7.5*pow(10,-5), 2.457*pow(10,-3),0.92};

		genDUNE * sterile = new genDUNE(xml);	
		sterile->prob = new SBNprob(4,angles3p1,phases3p1, mass_splittings3p1);
		sterile->preCalculateProbs();		
		sterile->doMC();

		SBNchi mychi(*sterile,*m);	
		std::vector<double> chi_all;

		//Here, vary the truth by poissonian_noise and build a distributuon for this. This is p-d/d then
		double chimin =999;
		for(double dcp = 0; dcp<=360; dcp+=5){
			for(int ord = 0; ord<2; ord++){
				for(int i23 =0; i23 < theta23.size(); i23++){

					std::string name = order_names.at(ord)+"_DCP_"+to_string_prec(dcp,3)+"_T23_"+to_string_prec(theta23.at(i23),3);

					SBNspec * test = new SBNspec(("precomp/"+name+".SBNspec").c_str(),xml);	

					double chi = mychi.CalcChi(test);

					chi_all.push_back(chi);
					dunestream<<"STERILESCAN "<<dcp<<" "<<ord<<" "<<theta23.at(i23)<<" "<<chi<<std::endl;
					if(chi<chimin)chimin=chi;
					delete test;
				}
			}
		}	
		dunestream<<"STERILESCAN MIN: "<<chimin<<std::endl;

	}



}

