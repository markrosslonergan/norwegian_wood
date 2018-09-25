
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
#include "plotting_tools.h"

#define no_argument 0
#define required_argument 1
#define optional_argument 2

using namespace sbn;

#include <limits>

std::fstream& GotoLine(std::fstream& file, unsigned int num){
  file.seekg(std::ios::beg);
  for(int i=0; i < num - 1; ++i){
    file.ignore(std::numeric_limits<std::streamsize>::max(),'\n');
  }
  return file;
}




/*************************************************************
 *************************************************************
 *		BEGIN example.cxx
 ************************************************************
 ************************************************************/
int main(int argc, char* argv[])
{


  std::string xml = "/a/data/westside/yjwa/NW/norwegian_wood/xml/dune.xml";
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
      {0,			no_argument, 		0, 0},
      {"process",			required_argument, 	0,  'p'},
      {"cpv4", required_argument, 0, 'c'},
    
      {"massorder", required_argument, 0, 's'},
      {"order", required_argument, 0, 'o'},
      {"quickcpv",required_argument, 0, 'q'}
    };

    
  int which_process = -1 ;

  while(iarg != -1)
    {
      iarg = getopt_long(argc,argv, "x:m:t:f:p:c:o:q:s:", longopts, &index);

      //if(0 <iarg && iarg< 1000){
        
      //     process = iarg;
            
      //std::cout << "Process is " << process << std::endl;
      //}
        
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
	case 'p':
	  //whiche_mode = "gen4";
	  which_process = atoi(optarg)+1;
	  which_mode = "process_gen4";

	  std::cout<<"optarg : " <<which_process << std::endl;
	  std::cout<<"optarg : " <<which_process << " mode : " << which_mode<< std::endl;
	  break;
	case 'o':
	  //whiche_mode = "gen4";
	  which_process = atoi(optarg)+1;
	  which_mode = "process_order";

	  std::cout<<"optarg : " <<which_process << std::endl;
	  std::cout<<"optarg : " <<which_process << " mode : " << which_mode<< std::endl;
	  break;
	case 'q':
	  //whiche_mode = "gen4";
	  which_process = atoi(optarg)+1;
	  which_mode = "process_cpv3p1";

	  std::cout<<"optarg : " <<which_process << std::endl;
	  std::cout<<"optarg : " <<which_process << " mode : " << which_mode<< std::endl;
	  break;



	case 'c':
	  which_process = atoi(optarg)+1;
	  which_mode = "process_cpv4";
	  std::cout<<"optarg : " <<which_process << std::endl;
	  std::cout<<"optarg : " <<which_process << " mode : " << which_mode<< std::endl;
	case 's':
	  which_process = atoi(optarg)+1;
	  which_mode = "process_order4";
	  std::cout<<"optarg : " <<which_process << std::endl;
	  std::cout<<"optarg : " <<which_process << " mode : " << which_mode<< std::endl;

	  break;
	case 't':
	  test_mode = strtof(optarg,NULL);
	  break;
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
    genDUNE bkg_only(xml); //yj genDUNE!

    if(false){
	  bkg_only.loadPreCalculatedProbs("./");
    }else{
  	  bkg_only.prob = new SBNprob(4, angles, phases, mass_splittings);
  	  bkg_only.preCalculateProbs();
    }



    bkg_only.doMC("three_neutrino");
    bkg_only.writeOut("three_neutrino.root");

    return 0;



  }


  else if(which_mode =="order"){
   
    std::ofstream dunestream;
    dunestream.open ("DUNE_order_NO.dat");
        
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
     
     
      double temp_chi = 9999.;
      double temp_chi_sub = 9999.;
      for(int tru_i23 =0; tru_i23 < theta23.size(); tru_i23++){
	std::string truth_name = order_names.at(0)+"_DCP_"+to_string_prec(tru_dcp , 3)+"_T23_"+to_string_prec(theta23.at(tru_i23),3);
	SBNspec * truth = new SBNspec(("precomp/"+truth_name+".SBNspec").c_str(),xml);
	truth->compressVector();
	std::cout << "assume truth : " <<truth_name << std::endl;
	SBNchi mychi(*truth,*m);
	
	for(double dcp = 0; dcp<360; dcp+=15){
	for(int i23 =0; i23 < theta23.size(); i23++){
	  std::string name = order_names.at(1)+"_DCP_"+to_string_prec(dcp,3)+"_T23_"+to_string_prec(theta23.at(i23),3);               
	  std::string name_sub = order_names.at(0)+"_DCP_"+to_string_prec(dcp,3)+"_T23_"+to_string_prec(theta23.at(i23),3);

	  SBNspec * test = new SBNspec(("precomp/"+name+".SBNspec").c_str(),xml);
	  SBNspec * test_sub = new SBNspec(("precomp/"+name_sub+".SBNspec").c_str(),xml);
	  //std::cout<< "test : " << name << std:: endl;

	  double chi = mychi.CalcChi(test);
	  double chi_sub = mychi.CalcChi(test_sub);
	  double delta_chi = chi;
	  // - chi_sub;

	  //std::cout << "delta_chi"<<delta_chi <<" , temp min : "<< temp_chi << std::endl;
  
	  if (chi<temp_chi){
	    temp_chi = chi;
	  }
	  if (chi_sub < temp_chi_sub){

	    temp_chi_sub = chi_sub;
	  }

	  chi_all.push_back(delta_chi);
	  delete test;
	  delete test_sub;
	}
	}
      }
      //std::cout<<"order true_dcp "<<tru_dcp<<" min : "<< temp_chi <<std::endl;
    
      std::cout<<"order true_dcp "<<tru_dcp<<" chi: "<<temp_chi <<" , chi_sub : "<< temp_chi_sub << " , delta_chi : " << temp_chi-temp_chi_sub <<std::endl;
      dunestream<<"order true_dcp "<<tru_dcp<<" "<< temp_chi-temp_chi_sub <<std::endl;
    }    
  }


  else if(which_mode =="order4"){
    std::ofstream dunestream;
    dunestream.open ("DUNE_order4_NO.dat");
 
    TFile *fin = new TFile("/a/data/westside/yjwa/NW/norwegian_wood/covar/covariance_matrices_xcheck_1408x1408.root","read");
    TMatrixT<double> * m = (TMatrixT<double>*)fin->Get("TMatrixT<double>;1");
        
    std::vector<std::string> order_names = {"NO","IO"};
    std::vector<double> order_vals = {2.457*pow(10,-3), -2.449*pow(10,-3)};
    std::vector<double> theta23 = {38,40,42,44,45,46,47,49,51,53};
    std::vector<double> theta34 = {0, 10, 15, 20,25};


    double t14 = 8.316;
    double t24 = 6.885;


    std::vector<double> min_chi;
    for(double tru_dcp2 = 90; tru_dcp2<360; tru_dcp2+=45){            
      for(double tru_dcp = 0; tru_dcp<=360; tru_dcp+=15){

	std::vector<double> chi_all;

	double temp_chi = 9999.;
	double temp_chi_sub = 9999.;

	for(int tru_i34 =0; tru_i34 < theta34.size(); tru_i34++){
	  std::string truth_name = order_names.at(1)+"_DCP_"+to_string_prec(tru_dcp,3)+"_T23_"+to_string_prec(theta23.at(3),3)+"_T14_"+to_string_prec(t14,3)+"_T24_"+to_string_prec(t24,3)+"_T34_"+to_string_prec(theta34.at(tru_i34),3)+"_D14_"+to_string_prec(tru_dcp2,6);

	  SBNspec * truth = new SBNspec(("/a/data/westside/yjwa/NW/DUNE_SBN_condor/condor_tests/"+truth_name+".SBNspec").c_str(),xml);
	  truth->compressVector();
	  std::cout << "assume truth : " <<truth_name << std::endl;
	  SBNchi mychi(*truth,*m);
	
	  for(double dcp = 0; dcp<360; dcp+=15){
	    for(double dcp2 = 0; dcp2<360; dcp2+=45){
	      for(int i23=0; i23<theta23.size(); i23++){

		for(int i34=0; i34<theta34.size(); i34++){
		  std::string name = order_names.at(0)+"_DCP_"+to_string_prec(dcp,3)+"_T23_"+to_string_prec(theta23.at(i23),3)+"_T14_"+to_string_prec(t14,3)+"_T24_"+to_string_prec(t24,3)+"_T34_"+to_string_prec(theta34.at(i34),3)+"_D14_"+to_string_prec(dcp2,6);

		  std::string name_sub = order_names.at(1)+"_DCP_"+to_string_prec(dcp,3)+"_T23_"+to_string_prec(theta23.at(i23),3)+"_T14_"+to_string_prec(t14,3)+"_T24_"+to_string_prec(t24,3)+"_T34_"+to_string_prec(theta34.at(i34),3)+"_D14_"+to_string_prec(dcp2,6);


		  SBNspec * test = new SBNspec(("/a/data/westside/yjwa/NW/DUNE_SBN_condor/condor_tests/"+name+".SBNspec").c_str(),xml);
		  SBNspec * test_sub = new SBNspec(("/a/data/westside/yjwa/NW/DUNE_SBN_condor/condor_tests/"+name_sub+".SBNspec").c_str(),xml);

		  double chi = mychi.CalcChi(test);
		  double chi_sub = mychi.CalcChi(test_sub);

		  if (chi<temp_chi){
		    temp_chi = chi;
		  }
		  if (chi_sub < temp_chi_sub){
		    temp_chi_sub = chi_sub;
		  }
		  //chi_all.push_back(delta_chi);
		  delete test;
		  delete test_sub;
		}
	      }
	    }
	  }
	  //std::cout<<"order true_dcp "<<tru_dcp<<" min : "<< temp_chi <<std::endl;
      
	  std::cout<<"order true_dcp "<<tru_dcp<<" chi: "<<temp_chi <<" , chi_sub : "<< temp_chi_sub << " , delta_chi : " << temp_chi-temp_chi_sub <<std::endl;
	  dunestream<<"order true_dcp "<<tru_dcp<< "tru_d14 "<<tru_dcp2<<" "<< temp_chi-temp_chi_sub <<std::endl;
	    
	}
      }
    }
  }	    


  else if(which_mode =="process_order4"){
   
      std::string process_line;
      std::string delimiter = "_";
      std::fstream datlist("/a/data/westside/yjwa/NW/DUNE_SBN_condor/condor_order/IOlist.txt");
      std::fstream file("/a/data/westside/yjwa/NW/norwegian_wood/process_t14_t24_t34_d14.txt");
      
      GotoLine(file,which_process+620);
      file >> process_line;      
      std::cout << process_line;
      
      std::vector<std::string> lineinfile;
      
      size_t pos = 0;
      std::string sub;
      std::string precalc_name;
      std::string line;
      
      int after_process_num = 0;
      std::cout << " ? " << std::endl;
      while( (pos = process_line.find(delimiter)) != std::string::npos) {
          sub = process_line.substr(0,pos);
          lineinfile.push_back(sub);
          process_line.erase(0, pos+delimiter.length());
	  std::cout << sub << std::endl;    
          if (after_process_num == 0){
              precalc_name = process_line;
          }
          after_process_num++;
      }
      lineinfile.push_back(process_line);
      std::cout << " ?? " << std::endl;
      bool donejob = false;
            
      double t14 = std::stod(lineinfile.at(2));
 std::cout << " ??? " << std::endl;
          
 double t24 = std::stod(lineinfile.at(4));
    std::cout << " ???? " << std::endl;
        double t34 = std::stod(lineinfile.at(6));
 std::cout << " ????? " << std::endl;
       
    double d14 = std::stod(lineinfile.at(8));
      
      while( std::getline (datlist, line)){
          if(line.find(precalc_name) != std::string::npos){
	    //	    donejob = true;
          }
      }
      
      if (donejob){
          std::cout << "This is donejob" <<std::endl;
          exit(0);
      }
      //for (int i = 0 ; i< lineinfile.size() ; i++){
      //  std::cout << "=== " << i << "-th element in process line is "  <<lineinfile.at(i) << std::endl;
          
      //}
      TFile *fin = new TFile("/a/data/westside/yjwa/NW/norwegian_wood/covar/covariance_matrices_xcheck_1408x1408.root","read");
      TMatrixT<double> * m = (TMatrixT<double>*)fin->Get("TMatrixT<double>;1");
        
      std::vector<std::string> order_names = {"NO","IO"};
      std::vector<double> order_vals = {2.457*pow(10,-3), -2.449*pow(10,-3)};
      std::vector<double> theta23 = {38,40,42,44,45,46,47,49,51,53};
      std::vector<double> theta34 = {0, 10, 15, 20,25};

      std::ofstream dunestream;
      dunestream.open(order_names.at(1)+"_T23_"+to_string_prec(theta23.at(3),3)+"_T14_"+to_string_prec(t14,3)+"_T24_"+to_string_prec(t24,3)+"_T34_"+to_string_prec(t34,3)+"_D14_"+to_string_prec(d14,6));
	
      for(double tru_dcp = 0; tru_dcp<=360; tru_dcp+=15){
	std::vector<double> chi_all;
	double temp_chi = 9999.;
	double temp_chi_sub = 9999.;
	std::string truth_name = order_names.at(1)+"_DCP_"+to_string_prec(tru_dcp,3)+"_T23_"+to_string_prec(theta23.at(3),3)+"_T14_"+to_string_prec(t14,3)+"_T24_"+to_string_prec(t24,3)+"_T34_"+to_string_prec(t34,3)+"_D14_"+to_string_prec(d14,6);
	
	SBNspec * truth = new SBNspec(("/a/data/westside/yjwa/NW/DUNE_SBN_condor/condor_tests/"+truth_name+".SBNspec").c_str(),xml);
	truth->compressVector();
	std::cout << "assume truth : " <<truth_name << std::endl;
	SBNchi mychi(*truth,*m);
	
	for(double dcp = 0; dcp<360; dcp+=15){
	  //  for(double dcp2 = 0; dcp2<360; dcp2+=45){
	    for(int i23=0; i23<theta23.size(); i23++){
	      //  for(int i34=0; i34<theta34.size(); i34++){
		std::string name = order_names.at(0)+"_DCP_"+to_string_prec(dcp,3)+"_T23_"+to_string_prec(theta23.at(i23),3)+"_T14_"+to_string_prec(t14,3)+"_T24_"+to_string_prec(t24,3)+"_T34_"+to_string_prec(t34,3)+"_D14_"+to_string_prec(d14,6);

		std::string name_sub = order_names.at(1)+"_DCP_"+to_string_prec(dcp,3)+"_T23_"+to_string_prec(theta23.at(i23),3)+"_T14_"+to_string_prec(t14,3)+"_T24_"+to_string_prec(t24,3)+"_T34_"+to_string_prec(t34,3)+"_D14_"+to_string_prec(d14,6);


	      SBNspec * test = new SBNspec(("/a/data/westside/yjwa/NW/DUNE_SBN_condor/condor_tests/"+name+".SBNspec").c_str(),xml);
	      SBNspec * test_sub = new SBNspec(("/a/data/westside/yjwa/NW/DUNE_SBN_condor/condor_tests/"+name_sub+".SBNspec").c_str(),xml);

	      double chi = mychi.CalcChi(test);
	      double chi_sub = mychi.CalcChi(test_sub);

	      if (chi<temp_chi){
		temp_chi = chi;
	      }
	      if (chi_sub < temp_chi_sub){
		temp_chi_sub = chi_sub;
	      }
	      //chi_all.push_back(delta_chi);
	      delete test;
	      delete test_sub;
	      }
	    
	  
	//std::cout<<"order true_dcp "<<tru_dcp<<" min : "<< temp_chi <<std::endl;
	}
	std::cout<<"order true_dcp "<<tru_dcp<<" chi: "<<temp_chi <<" , chi_sub : "<< temp_chi_sub << " , delta_chi : " << temp_chi-temp_chi_sub <<std::endl;
	dunestream<<"order true_dcp "<<tru_dcp<< " tru_d14 "<<d14<<" "<< temp_chi-temp_chi_sub <<std::endl;
	    
      }
  }
  	    


	
  else if(which_mode=="process_gen4"){//yj, attemp precompute 3+1
        
        
        
    std::string process_line;
    std::string delimiter = "_";

    std::fstream rootlist("/a/data/westside/yjwa/NW/DUNE_SBN_condor/condor_dm1ev/rootlist.txt");
        
    std::fstream file("/a/data/westside/yjwa/NW/norwegian_wood/process_t14_t24_t34_d14.txt");

    if (file.fail()){
      std::cout << "file is null " << std::endl; 
    }

    //std::fstream file("/Users/yeon-jaejwa/sandbox/NW/norwegian_wood/build/src/process_t14_t24_t34_d14.txt");
      
     
    //file.seekg(which_process-1);
        
        
    //getline.
        
    GotoLine(file,which_process);//310 condor jobs
        
        
        
        
        
    file >> process_line;
        
    std::cout << process_line << std::endl;
        
    //std::cin.get();
        
 
    std::vector<std::string> lineinfile;
        
    size_t pos = 0;
    std::string sub;
        
    std::string precalc_name;

    std::string line;
        
    int after_process_num = 0;

    std::cout << process_line << std::endl;

    while( (pos = process_line.find(delimiter)) != std::string::npos) {
            
            
            
            
      sub = process_line.substr(0,pos);
            
            

      //std::cout << stod(sub) << std::endl;
      lineinfile.push_back(sub);
      process_line.erase(0, pos+delimiter.length());
            
            
      if (after_process_num == 0){
                
	precalc_name = process_line;
      }
      after_process_num++;
            
    }
    lineinfile.push_back(process_line);
    //std::cout << str << std::endl;
    //line_count++;
        
    bool donejob = false;

    while( std::getline (rootlist, line)){
      if(line.find(precalc_name) != std::string::npos){

	//donejob = true;
      }
	    
    }

    if (donejob){
      std::cout << "This is donejob" <<std::endl;
      exit(0);
    }
        
    for (int i = 0 ; i< lineinfile.size() ; i++){
        
    std::cout << "=== " << i << "-th element in process line is "  <<lineinfile.at(i) << std::endl;
            
            
    }
        
    std::cout << precalc_name << std::endl;
        
    std::vector<double> angles = {34.5, 0, 8.45, 0,0,0};//6 mixing angles
    //std::vector<double> angles_oct = {33.6, 49.6, 8.5, 0,0,0};//dif octant
    std::vector<double> phases = {0,0,0};
    //std::vector<double> phases_180 = {180,0,0};
    //http://lbne2-docdb.fnal.gov/cgi-bin/RetrieveFile?docid=10688&filename=DUNE-CDR-physics-volume.pdf&version=10
    std::vector<double> mass_splittings = {7.55*pow(10,-5), 2.50*pow(10,-3),0};
    std::vector<double> mass_splittings_inv = {7.55*pow(10,-5), -2.42*pow(10,-3),0};
      
        
    TFile *fin = new TFile("/a/data/westside/yjwa/NW/norwegian_wood/covar/covariance_matrices_xcheck_1408x1408.root","read");
    TMatrixT<double> * m = (TMatrixT<double>*)fin->Get("TMatrixT<double>;1");
        
    std::vector<std::string> order_names = {"NO","IO"};
    std::vector<double> order_vals = {2.457*pow(10,-3), -2.449*pow(10,-3)};
        
    //std::vector<double> theta23 = {40,41,42,43,44,45,46,47,48,49,50,51,52};
    //std::vector<double> theta
    std::vector<double> theta23 = {40.8, 41.8, 42.8, 43.8, 45., 45.8, 46.8, 47.8, 48.8, 49.8, 50.8};
        
    double t14=10000;
    double t24=10000;
    double t34=10000;
    double d14=10000;

    //char * str_t14 = lineinfile.at(2);

    t14 = std::stod(lineinfile.at(2));
    t24 = std::stod(lineinfile.at(4));
    t34 = std::stod(lineinfile.at(6));
    d14 = std::stod(lineinfile.at(8));


    //for(double dcp = 0; dcp<=360; dcp+=15){//
    for(double dcp =0; dcp <360; dcp += 30){
      //for (double dcp = 0; dcp<=360; dcp+=180){//yj
      for(int ord = 0; ord<2; ord++){
	
	if (ord == 1){
	  //angles.at(0) = ;
	  angles.at(2) = 8.53;
	}

	for(int i23 =0; i23 < theta23.size(); i23++){
          
          
	  //    for(int i14 =0; i14 < theta14.size(); i14++){
	  //for(int i24 =0; i24 < theta24.size(); i24++){
	  //  for(int i34 =0; i34 < theta34.size(); i34++){
                                
                                
	  std::string name = order_names.at(ord)+"_DCP_"+to_string_prec(dcp,3)+"_T23_"+to_string_prec(theta23.at(i23),3)+"_"+precalc_name;
	  std::cout<<"GEN: on value "<<name<<std::endl;
	  std::cout<<dcp<<" "<<ord<<std::endl;  
	  genDUNE* testpt = new genDUNE(xml);//yj i was here

	  //testpt->loadPreCalculatedProbs("/a/data/westside/markross/norwegian_wood/build/src/");


 
	  std::cout << "gen dune worked " <<std::endl; 

	  phases.at(0) = dcp;
	  angles.at(1) = theta23.at(i23);
                                
	  angles.at(3) = t14;
	  angles.at(4) = t24;
	  angles.at(5) = t34;

	  phases.at(1) = d14;
                                
	  mass_splittings.at(1) = order_vals.at(ord);
	  mass_splittings.at(2) = 1.0;
                    
	  testpt->prob = new SBNprob(4,angles,phases, mass_splittings);
	  
	  std::cout << "sbnprob worked" <<std::endl;

	  testpt->preCalculateProbs();

	  std::cout << "pre calc" << std::endl; 
                    
	  testpt->doMC(name);

	  std::cout << "do MC worked  " << std::endl;
          
	  delete testpt;
	}
      }
    }
      
  } 
    
    
  else if(which_mode=="gen"){//yj, precompute 3+0

    std::vector<double> angles = {33.6, 42.1, 8.5, 0,0,0};//6 mixing angles
    std::vector<double> angles_oct = {33.6, 49.6, 8.5, 0,0,0};//dif octant
    std::vector<double> phases = {0,0,0};
    std::vector<double> phases_180 = {180,0,0};
    //http://lbne2-docdb.fnal.gov/cgi-bin/RetrieveFile?docid=10688&filename=DUNE-CDR-physics-volume.pdf&version=10
    std::vector<double> mass_splittings = {7.5*pow(10,-5), 2.457*pow(10,-3),0};
    std::vector<double> mass_splittings_inv = {7.5*pow(10,-5), -2.449*pow(10,-3),0};


    TFile *fin = new TFile("/a/data/westside/yjwa/NW/norwegian_wood/covar/covariance_matrices_xcheck_1408x1408.root","read");
    TMatrixT<double> * m = (TMatrixT<double>*)fin->Get("TMatrixT<double>;1");

    std::vector<std::string> order_names = {"NO","IO"};
    std::vector<double> order_vals = {2.457*pow(10,-3), -2.449*pow(10,-3)};

    //std::vector<double> theta23 = {40,41,42,43,44,45,46,47,48,49,50,51,52};
    std::vector<double> theta23 = {38,40,42,44,45,46,47,49,51,53};
    //std::vector<double> theta23_NO = {42.1-2,42.1-1, 42.1, 42.1+1,42.1+2};
    //std::vector<double> theta23_IO = {49.6-2,49.6-1, 49.6, 49.6+1,49.6+2};

    for(double dcp = 0; dcp<=360; dcp+=15){//
      //for (double dcp = 0; dcp<=360; dcp+=180){//yj
      for(int ord = 0; ord<2; ord++){
	for(int i23 =0; i23 < theta23.size(); i23++){
	  //for (int i23=0; i23<1;i23++){//yj
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

  }
    
  else if(which_mode=="gen4"){//yj, attemp precompute 3+1
        
        
        
    std::vector<double> angles = {33.6, 42.1, 8.5, 0,0,0};//6 mixing angles
    std::vector<double> angles_oct = {33.6, 49.6, 8.5, 0,0,0};//dif octant
    std::vector<double> phases = {0,0,0};
    std::vector<double> phases_180 = {180,0,0};
    //http://lbne2-docdb.fnal.gov/cgi-bin/RetrieveFile?docid=10688&filename=DUNE-CDR-physics-volume.pdf&version=10
    std::vector<double> mass_splittings = {7.5*pow(10,-5), 2.457*pow(10,-3),0};
    std::vector<double> mass_splittings_inv = {7.5*pow(10,-5), -2.449*pow(10,-3),0};
        
        
    TFile *fin = new TFile("/a/data/westside/yjwa/NW/norwegian_wood/covar/covariance_matrices_xcheck_1408x1408.root","read");
    TMatrixT<double> * m = (TMatrixT<double>*)fin->Get("TMatrixT<double>;1");
        
    std::vector<std::string> order_names = {"NO","IO"};
    std::vector<double> order_vals = {2.457*pow(10,-3), -2.449*pow(10,-3)};
        
    std::vector<double> theta23 = {40,41,42,43,44,45,46,47,48,49,50,51,52};
    //std::vector<double> theta
        
    std::vector<double> theta14 = {0,10};
    std::vector<double> theta24 = {0,10};
    std::vector<double> theta34 = {0};//rather coarse for now for the purpose of test.
        
    std::vector<double> theta23_NO = {42.1-2,42.1-1, 42.1, 42.1+1,42.1+2};
    std::vector<double> theta23_IO = {49.6-2,49.6-1, 49.6, 49.6+1,49.6+2};
        
    for(double dcp = 0; dcp<=360; dcp+=5){//
      //for (double dcp = 0; dcp<=360; dcp+=180){//yj
      for(int ord = 0; ord<2; ord++){
	for(int i23 =0; i23 < theta23.size(); i23++){
                    
	  for(int i14 =0; i14 < theta14.size(); i14++){
	    for(int i24 =0; i24 < theta24.size(); i24++){
	      for(int i34 =0; i34 < theta34.size(); i34++){
                                
                                
		std::string name = order_names.at(ord)+"_DCP_"+to_string_prec(dcp,3)+"_T23_"+to_string_prec(theta23.at(i23),3)+"_T14_"+to_string_prec(theta14.at(i14),3)+"_T24_"+to_string_prec(theta24.at(i24),3)+"_T34_"+to_string_prec(theta34.at(i34),3);
		std::cout<<"GEN: on value "<<name<<std::endl;
                    
		genDUNE* testpt = new genDUNE(xml);
		phases.at(0) = dcp;
		angles.at(1) = theta23.at(i23);
                                
		angles.at(3) = theta14.at(i14);
		angles.at(4) = theta24.at(i24);
		angles.at(5) = theta34.at(i34);

                                
		mass_splittings.at(1) = order_vals.at(ord);
		mass_splittings.at(2) = 1.;
                    
		testpt->prob = new SBNprob(4,angles,phases, mass_splittings);
		testpt->preCalculateProbs();
                    
		testpt->doMC("precomp4/"+name);
                    
		delete testpt;
	      }
	    }
	  }
	}
      }
    }
        
  }
    
    
  else if(which_mode=="globalfit"){//yj, produce process number and grid point
        
    std::ofstream globalfitstream;
    globalfitstream.open ("process_t14_t24_t34_d14_4ue2um2.txt");//y
    std::ofstream globaldatstream;
    globaldatstream.open ("IO_dat_4ue2um2.txt");
        
        
    std::vector<double> angles = {33.6, 42.1, 8.5, 0,0,0};//6 mixing angles
    std::vector<double> angles_oct = {33.6, 49.6, 8.5, 0,0,0};//dif octant
    std::vector<double> phases = {0,0,0};
    std::vector<double> phases_180 = {180,0,0};
    //http://lbne2-docdb.fnal.gov/cgi-bin/RetrieveFile?docid=10688&filename=DUNE-CDR-physics-volume.pdf&version=10
    std::vector<double> mass_splittings = {7.5*pow(10,-5), 2.457*pow(10,-3),0};
    std::vector<double> mass_splittings_inv = {7.5*pow(10,-5), -2.449*pow(10,-3),0};
        
        
    TFile *fin = new TFile("/a/data/westside/yjwa/NW/norwegian_wood/covar/covariance_matrices_xcheck_1408x1408.root","read");
    TMatrixT<double> * m = (TMatrixT<double>*)fin->Get("TMatrixT<double>;1");
        
    std::vector<std::string> order_names = {"NO","IO"};
    std::vector<double> order_vals = {2.457*pow(10,-3), -2.449*pow(10,-3)};
        
    //std::vector<double> theta23 = {40,41,42,43,44,45,46,47,48,49,50,51,52};
    std::vector<double> theta23 = {38,40,42,44,45,46,47,49,51,53};

    //std::vector<double> theta
        
    std::vector<double> theta14 = {5.736, 7.026, 8.316, 10.308, 12.3};
    std::vector<double> theta24 = {3.654, 5.2695, 6.885, 9.9075, 12.93};
    std::vector<double> theta34 = {0, 10, 15, 20,25};
    std::vector<double> delta14 = {0, 45, 90, 135, 180, 225, 270, 315};
    //std::vector<double> theta23_NO = {42.1-2,42.1-1, 42.1, 42.1+1,42.1+2};
    //std::vector<double> theta23_IO = {49.6-2,49.6-1, 49.6, 49.6+1,49.6+2};
        
    int process_num = 0;
        
    for(int i14 =0; i14 < theta14.size(); i14++){
      for(int i24 =0; i24 < theta24.size(); i24++){
	for(int i34 =0; i34 < theta34.size(); i34++){
                    
	  for (int d14 = 0; d14 < delta14.size(); d14++){
                        
	    double sin2theta = 4*std::sin( theta14.at(i14)*3.14159/180. )*std::sin( theta14.at(i14)*3.14159/180. )*std::cos( theta14.at(i14)*3.14159/180. )*std::cos(theta14.at(i14)*3.14159/180.)*std::sin(theta24.at(i24)*3.14159/180.)*std::sin(theta24.at(i24)*3.14159/180.);
	    std::string name = "_T14_"+to_string_prec(theta14.at(i14),3)+"_T24_"+to_string_prec(theta24.at(i24),3)+"_T34_"+to_string_prec(theta34.at(i34),3)+"_D14_"+to_string_prec(delta14.at(d14))+"_4UE2UM2_"+to_string_prec(sin2theta,6);
                        
	    process_num++;
                        
	    globalfitstream<<process_num<<name<<std::endl;

	    if(sin2theta<0.002217){
	      globaldatstream <<"cpv4_IO_T14_"+to_string_prec(theta14.at(i14),3)+"_T24_"+to_string_prec(theta24.at(i24),3)+"_T34_"+to_string_prec(theta34.at(i34),3)+"_D14_"+to_string_prec(delta14.at(d14),6)+".dat"<<std::endl;
	    }
	  }
	}
      }
    }
  }

  else if(which_mode =="process_cpv3p1"){
      
      std::string process_line;
      std::string delimiter = "_";
      
      std::fstream datlist("/a/data/westside/markross/DUNE_SBN_condor/order/datlistIO.txt");
     
 
      std::fstream file("/a/data/westside/markross/norwegian_wood/process_t14_t24_t34_d14.txt");
      //std::fstream file("/a/data/westside/markross/DUNE_SBN_condor/order/makeup.list");
     
      //file.seekg(which_process-1); 
      //
      int working_line = which_process;      
      
      //getline.
      
      //std::cout<<which_process<<std::endl; 
      GotoLine(file, working_line);//320+320+90 condor job s 160
      
      file >> process_line;
      
      //std::cout<<"Anything?"<<std::endl;     
      std::cout << process_line;
      //std::cout<<"Anything?"<<std::endl;     
 
      //std::cin.get();
      
      
      
      std::vector<std::string> lineinfile;
      
      size_t pos = 0;
      std::string sub;
      
      std::string precalc_name;
      
      std::string line;
      
      int after_process_num = 0;
      
      while( (pos = process_line.find(delimiter)) != std::string::npos) {
          
          
          sub = process_line.substr(0,pos);
          
	 // std::cout << stod(sub) << std::endl;

          lineinfile.push_back(sub);
          process_line.erase(0, pos+delimiter.length());
          
          
          if (after_process_num == 0){
              
              precalc_name = process_line;
          }
          after_process_num++;
          
      }
      lineinfile.push_back(process_line);
      //std::cout << str << std::endl;
      //line_count++;
      
      bool donejob = false;
      
      while( std::getline (datlist, line)){
          if(line.find(precalc_name) != std::string::npos){
              
	    //donejob = true;
          }
          
      }
      
      if (donejob){
          std::cout << "This is donejob" <<std::endl;
          exit(0);
      }
      
      //for (int i = 0 ; i< lineinfile.size() ; i++){
      //  std::cout << "=== " << i << "-th element in process line is "  <<lineinfile.at(i) << std::endl;
      //}

      //std::ofstream dunestream;
      //dunestream.open ("DUNE_cpv3+1to3+0_NO_cpv4_testpoint.dat");
     


 
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
      std::vector<double> theta23_NO = {42.1-2,42.1-1, 42.1, 42.1+1,42.1+2};
      std::vector<double> theta23_IO = {49.6-2,49.6-1, 49.6, 49.6+1,49.6+2};
      
      double t14=10000;
      double t24=10000;
      double t34=10000;
      double d14=10000;
      
      //char * str_t14 = lineinfile.at(2);
	std::cout<<"Getting t14..etc.."<<std::endl;     
 
      t14 = std::stod(lineinfile.at(2));
      t24 = std::stod(lineinfile.at(4));
      t34 = std::stod(lineinfile.at(6));
      d14 = std::stod(lineinfile.at(8));

 
	std::ofstream dunestream;
	std::cout<<"Starting on CPV true plot"<<std::endl;
	std::string out_name  = "cpv3p1_"+order_names.at(0)+"_T14_"+to_string_prec(t14,3)+"_T24_"+to_string_prec(t24,3)+"_T34_"+to_string_prec(t34,3)+"."+to_string_prec(working_line,2)+".dat";
	calc_cpv_3p1( &dunestream, out_name, xml, t14, t24, t34, m);

  }



     
  else if(which_mode =="process_order"){
      
      std::string process_line;
      std::string delimiter = "_";
      
      std::fstream datlist("/a/data/westside/markross/DUNE_SBN_condor/order/datlistIO.txt");
     
 
      std::fstream file("/a/data/westside/markross/norwegian_wood/process_t14_t24_t34_d14.txt");
      //std::fstream file("/a/data/westside/markross/DUNE_SBN_condor/order/makeup.list");
     
      //file.seekg(which_process-1);
      
      
      //getline.
      
      //std::cout<<which_process<<std::endl; 
      GotoLine(file, which_process+960);//320+320+90 condor job s 160
      
      file >> process_line;
      
      //std::cout<<"Anything?"<<std::endl;     
      std::cout << process_line;
      //std::cout<<"Anything?"<<std::endl;     
 
      //std::cin.get();
      
      
      
      std::vector<std::string> lineinfile;
      
      size_t pos = 0;
      std::string sub;
      
      std::string precalc_name;
      
      std::string line;
      
      int after_process_num = 0;
      
      while( (pos = process_line.find(delimiter)) != std::string::npos) {
          
          
          sub = process_line.substr(0,pos);
          
	 // std::cout << stod(sub) << std::endl;

          lineinfile.push_back(sub);
          process_line.erase(0, pos+delimiter.length());
          
          
          if (after_process_num == 0){
              
              precalc_name = process_line;
          }
          after_process_num++;
          
      }
      lineinfile.push_back(process_line);
      //std::cout << str << std::endl;
      //line_count++;
      
      bool donejob = false;
      
      while( std::getline (datlist, line)){
          if(line.find(precalc_name) != std::string::npos){
              
	    //donejob = true;
          }
          
      }
      
      if (donejob){
          std::cout << "This is donejob" <<std::endl;
          exit(0);
      }
      
      //for (int i = 0 ; i< lineinfile.size() ; i++){
      //  std::cout << "=== " << i << "-th element in process line is "  <<lineinfile.at(i) << std::endl;
      //}

      //std::ofstream dunestream;
      //dunestream.open ("DUNE_cpv3+1to3+0_NO_cpv4_testpoint.dat");
     


 
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
      std::vector<double> theta23_NO = {42.1-2,42.1-1, 42.1, 42.1+1,42.1+2};
      std::vector<double> theta23_IO = {49.6-2,49.6-1, 49.6, 49.6+1,49.6+2};
      
      double t14=10000;
      double t24=10000;
      double t34=10000;
      double d14=10000;
      
      //char * str_t14 = lineinfile.at(2);
	std::cout<<"Getting t14..etc.."<<std::endl;     
 
      t14 = std::stod(lineinfile.at(2));
      t24 = std::stod(lineinfile.at(4));
      t34 = std::stod(lineinfile.at(6));
      d14 = std::stod(lineinfile.at(8));

 
	std::ofstream dunestream;
	std::cout<<"Starting on mass ordering plot"<<std::endl;
	std::string out_name  = "order_"+order_names.at(1)+"_T14_"+to_string_prec(t14,3)+"_T24_"+to_string_prec(t24,3)+"_T34_"+to_string_prec(t34,3)+"_D14_"+to_string_prec(d14,6)+".dat";
	calc_neutrino_ordering_3p1( &dunestream, out_name, xml, t14, t24, t34, d14 );

  }



  else if(which_mode =="process_cpv4"){
      
      std::string process_line;
      std::string delimiter = "_";
      
      std::fstream datlist("/a/data/westside/yjwa/NW/DUNE_SBN_condor/condor_2/datlistIO.txt");
      
      std::fstream file("/a/data/westside/yjwa/NW/norwegian_wood/process_t14_t24_t34_d14.txt");
      //std::fstream file("/Users/yeon-jaejwa/sandbox/NW/norwegian_wood/build/src/process_t14_t24_t34_d14.txt");
     
      //file.seekg(which_process-1);
      
      
      //getline.
      
      GotoLine(file,which_process+600);//320+320+90 condor job s 160
      
      
      
      
      
      file >> process_line;
      
      std::cout << process_line;
      
      //std::cin.get();
      
      
      
      std::vector<std::string> lineinfile;
      
      size_t pos = 0;
      std::string sub;
      
      std::string precalc_name;
      
      std::string line;
      
      int after_process_num = 0;
      
      while( (pos = process_line.find(delimiter)) != std::string::npos) {
          
          
          
          
          sub = process_line.substr(0,pos);
          
          
          
          //std::cout << stod(sub) << std::endl;
          lineinfile.push_back(sub);
          process_line.erase(0, pos+delimiter.length());
          
          
          if (after_process_num == 0){
              
              precalc_name = process_line;
          }
          after_process_num++;
          
      }
      lineinfile.push_back(process_line);
      //std::cout << str << std::endl;
      //line_count++;
      
      bool donejob = false;
      
      while( std::getline (datlist, line)){
          if(line.find(precalc_name) != std::string::npos){
              
	    //donejob = true;
          }
          
      }
      
      if (donejob){
          std::cout << "This is donejob" <<std::endl;
          exit(0);
      }
      
      //for (int i = 0 ; i< lineinfile.size() ; i++){
          
      //  std::cout << "=== " << i << "-th element in process line is "  <<lineinfile.at(i) << std::endl;
          
          
      //}

      //std::ofstream dunestream;
      //dunestream.open ("DUNE_cpv3+1to3+0_NO_cpv4_testpoint.dat");
      
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
      std::vector<double> theta23_NO = {42.1-2,42.1-1, 42.1, 42.1+1,42.1+2};
      std::vector<double> theta23_IO = {49.6-2,49.6-1, 49.6, 49.6+1,49.6+2};
      
      double t14=10000;
      double t24=10000;
      double t34=10000;
      double d14=10000;
      
      //char * str_t14 = lineinfile.at(2);
      
      t14 = std::stod(lineinfile.at(2));
      t24 = std::stod(lineinfile.at(4));
      t34 = std::stod(lineinfile.at(6));
      d14 = std::stod(lineinfile.at(8));


                      
    std::ofstream dunestream;
    dunestream.open("cpv4_"+order_names.at(1)+"_T14_"+to_string_prec(t14,3)+"_T24_"+to_string_prec(t24,3)+"_T34_"+to_string_prec(t34,3)+"_D14_"+to_string_prec(d14,6)+".dat");
    std::cout << "cpv4_"+order_names.at(1)+"_T14_"+to_string_prec(t14,3)+"_T24_"+to_string_prec(t24,3)+"_T34_"+to_string_prec(t34,3)+"_D14_"+to_string_prec(d14,6)+".dat" << std::endl;
                      
                      
    for (double tru_dcp = 0; tru_dcp<360; tru_dcp+=15){
                          std::string truth_name = order_names.at(1)+"_DCP_"+to_string_prec(tru_dcp,3)+"_T23_"+to_string_prec(theta23.at(3),3)+"_T14_"+to_string_prec(t14,3)+"_T24_"+to_string_prec(t24,3)+"_T34_"+to_string_prec(t34,3)+"_D14_"+to_string_prec(d14,6);
                          SBNspec * truth = new SBNspec(("/a/data/westside/yjwa/NW/DUNE_SBN_condor/condor_tests/"+truth_name+".SBNspec").c_str(),xml);
                          //		std::ofstream dunestream;
                          //dunestream.open("cpv4_"+truth_name+".dat");
                          truth->compressVector();
                          SBNchi mychi(*truth,*m);
                          //Here, vary the truth by poissonian_noise and build a distributuon for this. This is p-d/d then
                          
                          std::vector<double> chi_all;
                          std::vector<double> chi_0pi;
                          
                          std::vector<double> chi_all_sub;
                          std::vector<double> chi_0pi_sub;
                          
                          //min of (chi2 at 0 - chi2 at true, chi2 at 180 - chi2 at true)
                          
                          
                          for(double dcp = 0; dcp<=360; dcp+=15){
                              
                              //std::cout << "tru dcp : "<<tru_dcp << " , dcp : " << dcp << std::endl;
                              for(int ord = 0; ord<2; ord++){
                                  for(int i23 =0; i23 < theta23.size(); i23++){
                                      
                                      std::string name = order_names.at(ord)+"_DCP_"+to_string_prec(dcp,3)+"_T23_"+to_string_prec(theta23.at(i23),3);
                                      
                                      SBNspec * test = new SBNspec(("/a/data/westside/yjwa/NW/norwegian_wood/build/src/precomp/"+name+".SBNspec").c_str(),xml);
                                      
                                      std::string name_sub = order_names.at(ord)+"_DCP_"+to_string_prec(tru_dcp,3)+"_T23_"+to_string_prec(theta23.at(i23),3);
                                      
                                      SBNspec * test_sub = new SBNspec(("/a/data/westside/yjwa/NW/norwegian_wood/build/src/precomp/"+name_sub+".SBNspec").c_str(),xml);
                                      
                                      
                                      
                                      double chi = mychi.CalcChi(test);
                                      
                                      
                                      double chi_sub = mychi.CalcChi(test_sub);
                                      
                                      if(dcp<1 || (dcp < 181 && dcp > 179) ){
                                          chi_0pi.push_back(chi);
                                          
                                          chi_0pi_sub.push_back(chi-chi_sub);
                                          std::cout << "chi_0pi"<< chi << " , chi_0pi-sub : " << chi_sub << std::endl;
                                      }
                                      
                                      chi_all.push_back(chi);
                                      
                                      chi_all_sub.push_back(chi-chi_sub);
                                      
                                      
                                      delete test;
                                      delete test_sub;
                                  }
                              }
                          }
                          double min_chi_all = 1e20; //std::min_element(chi_all.begin(), chi_all.end());                                                                           
                          double min_chi_0pi = 1e20; //std::min_element(chi_0pi.begin(), chi_0pi.end());                                                                    
                          
                          double min_chi_all_sub = 1e20; //std::min_element(chi_all.begin(), chi_all.end());                                                                
                          double min_chi_0pi_sub = 1e20; //std::min_element(chi_0pi.begin(), chi_0pi.end());                                                                
                          
                          for(auto &X: chi_all){
                              if( X< min_chi_all) min_chi_all =X;
                          }
                          for(auto &X: chi_0pi){
                              if( X< min_chi_0pi) min_chi_0pi =X;
                          }
                          
                          for(auto &X: chi_all_sub){
                              if( X< min_chi_all_sub) min_chi_all_sub =X;
                          }
                          for(auto &X: chi_0pi_sub){
                              if( X< min_chi_0pi_sub) min_chi_0pi_sub =X;
                          }
                          
                          
                          if(min_chi_all>min_chi_0pi) std::cout<<"ERROR WARBING WARRRRBBBBING! global min is bigger than 0-pi min"<<std::endl;
                          
                          dunestream<<"CPV true_dcp "<<tru_dcp<<" "<<min_chi_0pi<<" "<<min_chi_all<<" "<<min_chi_0pi-min_chi_all<<" "<< min_chi_0pi_sub<<" "<<min_chi_all_sub <<  std::endl;
                      }
                      
                      //one loop for dcp ends here

  }


  else if(which_mode =="bestfit_cpv4"){
      
    std::ofstream dunestream;
    dunestream.open ("DUNE_cpv3+1to3+1_NO_cpv4_testpoint.dat");
      
    //http://lbne2-docdb.fnal.gov/cgi-bin/RetrieveFile?docid=10688&filename=DUNE-CDR-physics-volume.pdf&version=1
      
    TFile *fin = new TFile("/a/data/westside/yjwa/NW/norwegian_wood/covar/covariance_matrices_xcheck_1408x1408.root","read");
    TMatrixT<double> * m = (TMatrixT<double>*)fin->Get("TMatrixT<double>;1");
      
    std::vector<std::string> order_names = {"NO","IO"};   
    std::vector<double> theta23 = {38,40,42,44,45,46,47,49,51,53};
     
    double t14=8.316;
    double t24=6.885;
    double t34=0.;
    double d14=0.;
                     
    for (double tru_dcp = 180; tru_dcp<360; tru_dcp+=15){
      for (double tru_dcp2 = 0; tru_dcp2<360; tru_dcp2+=45){

	std::string truth_name = order_names.at(0)+"_DCP_"+to_string_prec(tru_dcp,3)+"_T23_"+to_string_prec(theta23.at(3),3)+"_T14_"+to_string_prec(t14,3)+"_T24_"+to_string_prec(t24,3)+"_T34_"+to_string_prec(t34,3)+"_D14_"+to_string_prec(tru_dcp2,6);
	SBNspec * truth = new SBNspec(("/a/data/westside/yjwa/NW/DUNE_SBN_condor/condor_tests/"+truth_name+".SBNspec").c_str(),xml);
     
	truth->compressVector();
	SBNchi mychi(*truth,*m);
                          
	std::vector<double> chi_all;
	std::vector<double> chi_0pi;
                                                    
	for(double dcp = 0; dcp<360; dcp+=15){
	  for (double dcp2 = 0; dcp2 < 360; dcp2+=45){                     
	    //std::cout << "tru dcp : "<<tru_dcp << " , dcp : " << dcp << std::endl;
	    for(int ord = 0; ord<2; ord++){
	      for(int i23 =0; i23 < theta23.size(); i23++){
                                      
		std::string name = order_names.at(ord)+"_DCP_"+to_string_prec(dcp,3)+"_T23_"+to_string_prec(theta23.at(i23),3)+"_T14_"+to_string_prec(t14,3)+"_T24_"+to_string_prec(t24,3)+"_T34_"+to_string_prec(t34,3)+"_D14_"+to_string_prec(dcp2,6);
                                      
		SBNspec * test = new SBNspec(("/a/data/westside/yjwa/NW/DUNE_SBN_condor/condor_tests/"+name+".SBNspec").c_str(),xml);
                                      
		double chi = mychi.CalcChi(test);
                                      
		if(dcp<1 || (dcp < 181 && dcp > 179) ){
		  chi_0pi.push_back(chi);
		 
		}
		if(dcp2 < 1 || (dcp2 < 181 && dcp2 > 179)){
		  chi_0pi.push_back(chi);
		}
                                      
		chi_all.push_back(chi);
                                      
		                                      
		delete test;
	    
	      }
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
                          
      dunestream<<"CPV true_dcp "<<tru_dcp<<" "<<min_chi_0pi<<" "<<min_chi_all<<" "<<min_chi_0pi-min_chi_all<<" "<<  std::endl;
    }
  }
  }
                      
  //one loop for dcp ends here

      
    


    
  else if(which_mode =="cpv4"){
    //std::ofstream dunestream;
    //dunestream.open ("DUNE_cpv3+1to3+0_NO_cpv4_testpoint.dat");
        
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
    std::vector<double> theta23_NO = {42.1-2,42.1-1, 42.1, 42.1+1,42.1+2};
    std::vector<double> theta23_IO = {49.6-2,49.6-1, 49.6, 49.6+1,49.6+2};


    std::vector<double> theta14 = {5.736, 7.026, 8.316, 10.308, 12.3};
    std::vector<double> theta24 = {3.654, 5.2695, 6.885, 9.9075, 12.93};
    std::vector<double> theta34 = {0, 10, 15, 20, 25};
    std::vector<double> delta14 = {0, 45, 90, 135, 180, 225, 270, 315};

    for(double tru_d14 = 0; tru_d14<1; tru_d14++){
      for (double tru_t14 = 0; tru_t14<1; tru_t14++){
	for (double tru_t24 = 0; tru_t24<1; tru_t24++){
	  for (double tru_t34 = 0; tru_t34<1; tru_t34++){

	    std::ofstream dunestream;
	    dunestream.open("cpv4/cpv4_"+order_names.at(0)+"_T14_"+to_string_prec(theta14.at(tru_t14),3)+"_T24_"+to_string_prec(theta24.at(tru_t24),3)+"_T34_"+to_string_prec(theta34.at(tru_t34),3)+"_D14_"+to_string_prec(delta14.at(tru_d14),6)+".dat");


	    for (double tru_dcp = 0; tru_dcp<360; tru_dcp+=15){
	      std::string truth_name = order_names.at(0)+"_DCP_"+to_string_prec(tru_dcp,3)+"_T23_"+to_string_prec(theta23.at(3),3)+"_T14_"+to_string_prec(theta14.at(tru_t14),3)+"_T24_"+to_string_prec(theta24.at(tru_t24),3)+"_T34_"+to_string_prec(theta34.at(tru_t34),3)+"_D14_"+to_string_prec(delta14.at(tru_d14),6);
	      SBNspec * truth = new SBNspec(("/a/data/westside/yjwa/NW/DUNE_SBN_condor/condor_tests/"+truth_name+".SBNspec").c_str(),xml);
	      //		std::ofstream dunestream;
	      //dunestream.open("cpv4_"+truth_name+".dat");
	      truth->compressVector();
	      SBNchi mychi(*truth,*m);
	      //Here, vary the truth by poissonian_noise and build a distributuon for this. This is p-d/d then
            
	      std::vector<double> chi_all;
	      std::vector<double> chi_0pi;

	      std::vector<double> chi_all_sub;
	      std::vector<double> chi_0pi_sub;

	      //min of (chi2 at 0 - chi2 at true, chi2 at 180 - chi2 at true)


	      for(double dcp = 0; dcp<=360; dcp+=15){

		//std::cout << "tru dcp : "<<tru_dcp << " , dcp : " << dcp << std::endl;  
		for(int ord = 0; ord<2; ord++){
		  for(int i23 =0; i23 < theta23.size(); i23++){
                        
		    std::string name = order_names.at(ord)+"_DCP_"+to_string_prec(dcp,3)+"_T23_"+to_string_prec(theta23.at(i23),3);
                        
		    SBNspec * test = new SBNspec(("precomp/"+name+".SBNspec").c_str(),xml);
                        
		    std::string name_sub = order_names.at(ord)+"_DCP_"+to_string_prec(tru_dcp,3)+"_T23_"+to_string_prec(theta23.at(i23),3);
                        
		    SBNspec * test_sub = new SBNspec(("precomp/"+name_sub+".SBNspec").c_str(),xml);
                        
			

		    double chi = mychi.CalcChi(test);
                        

		    double chi_sub = mychi.CalcChi(test_sub);

		    if(dcp<1 || (dcp < 181 && dcp > 179) ){
		      chi_0pi.push_back(chi);
			
		      chi_0pi_sub.push_back(chi-chi_sub);
		      std::cout << "chi_0pi"<< chi << " , chi_0pi-sub : " << chi_sub << std::endl; 
		    }
                        
		    chi_all.push_back(chi);

		    chi_all_sub.push_back(chi-chi_sub);


		    delete test;
		    delete test_sub;
		  }
		}
	      }
	      double min_chi_all = 1e20; //std::min_element(chi_all.begin(), chi_all.end());                                                                           
	      double min_chi_0pi = 1e20; //std::min_element(chi_0pi.begin(), chi_0pi.end());                                                                    

	      double min_chi_all_sub = 1e20; //std::min_element(chi_all.begin(), chi_all.end());                                                                
	      double min_chi_0pi_sub = 1e20; //std::min_element(chi_0pi.begin(), chi_0pi.end());                                                                

	      for(auto &X: chi_all){
		if( X< min_chi_all) min_chi_all =X;
	      }
	      for(auto &X: chi_0pi){
		if( X< min_chi_0pi) min_chi_0pi =X;
	      }

	      for(auto &X: chi_all_sub){
		if( X< min_chi_all_sub) min_chi_all_sub =X;
	      }
	      for(auto &X: chi_0pi_sub){
		if( X< min_chi_0pi) min_chi_0pi_sub =X;
	      }


	      if(min_chi_all>min_chi_0pi) std::cout<<"ERROR WARBING WARRRRBBBBING! global min is bigger than 0-pi min"<<std::endl;

	      dunestream<<"CPV true_dcp "<<tru_dcp<<" "<<min_chi_0pi<<" "<<min_chi_all<<" "<<min_chi_0pi-min_chi_all<<" "<< min_chi_0pi_sub<< std::endl;
	    }
	
	    //one loop for dcp ends here
	
	  }
	}
      }	
    }
	
       
            
  }


    
  else if(which_mode =="cpv4_test"){
    std::ofstream dunestream;
    dunestream.open ("DUNE_cpv3+1to3+0_NO_cpv4_testpoint.dat");
        
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
    std::vector<double> theta23_NO = {42.1-2,42.1-1, 42.1, 42.1+1,42.1+2};
    std::vector<double> theta23_IO = {49.6-2,49.6-1, 49.6, 49.6+1,49.6+2};
        
        
    std::vector<double> theta14 = {5.736, 7.026, 8.316, 10.308, 12.3};
    std::vector<double> theta24 = {3.654, 5.2695, 6.885, 9.9075, 12.93};
    std::vector<double> theta34 = {0, 10, 15, 20, 25};
    std::vector<double> delta14 = {0, 45, 90, 135, 180, 225, 270, 315};
 
    for(double tru_dcp = 0; tru_dcp<=360; tru_dcp+=15){
      std::string truth_name = order_names.at(0)+"_DCP_"+to_string_prec(tru_dcp,3)+"_T23_"+to_string_prec(theta23.at(3),3)+"_T14_"+to_string_prec(theta14.at(4),3)+"_T24_"+to_string_prec(theta24.at(0),3)+"_T34_"+to_string_prec(theta34.at(2),3)+"_D14_"+to_string_prec(delta14.at(2),6);
      SBNspec * truth = new SBNspec(("/a/data/westside/yjwa/NW/DUNE_SBN_condor/condor_tests/"+truth_name+".SBNspec").c_str(),xml);
            
      truth->compressVector();
            
      SBNchi mychi(*truth,*m);
      std::vector<double> chi_all;
      std::vector<double> chi_0pi;
            
      std::vector<double> chi_all_sub;
      std::vector<double> delta_chi_0pi;
            
      std::string truth_name_sub = order_names.at(0)+"_DCP_"+to_string_prec(tru_dcp,3)+"_T23_"+to_string_prec(theta23.at(3),3);
      SBNspec * truth_sub = new SBNspec(("precomp/"+truth_name_sub+".SBNspec").c_str(),xml);
      truth_sub->compressVector();
            
      SBNchi mychi_sub(*truth_sub,*m);
            
      double chi_sub = mychi_sub.CalcChi(truth);
	
 
      for(double dcp = 0; dcp<360; dcp+=180){
	for(int ord = 0; ord<2; ord++){
	  //for(int i23 =0; i23 < theta23.size(); i23++){
                        
	  std::string name = order_names.at(ord)+"_DCP_"+to_string_prec(dcp,3)+"_T23_"+to_string_prec(theta23.at(3),3);
                        
	  SBNspec * test = new SBNspec(("precomp/"+name+".SBNspec").c_str(),xml);
	  test->compressVector();
	  SBNchi testchi(*test,*m);

	  double chi = testchi.CalcChi(truth);

	  std::cout<< "truth: "<<truth_name << ", truth3+0: " << truth_name_sub<< ", test: "<<name<<" , chi "<< chi << ", chi4to3 "<< chi_sub << " "<<chi-chi_sub<< std::endl;

        
	  //double chi_sub = chi - 0;
                        
	  if(dcp<1 || (dcp < 181 && dcp > 179) ){
	    chi_0pi.push_back(chi);
    
	  }
                        
	  chi_all.push_back(chi);

	  delete test;
                        
	  //}
	}
      }
            
      /*
	chisq true -> chisq(delta test _ CP ) - chisq()
      */
            
      double min_chi_all = 1e20; //std::min_element(chi_all.begin(), chi_all.end());
      double min_chi_0pi = 1e20; //std::min_element(chi_0pi.begin(), chi_0pi.end());

      for(auto &X: chi_all){
	if( X< min_chi_all) min_chi_all =X;
      }
      for(auto &X: chi_0pi){
	if( X< min_chi_0pi) min_chi_0pi =X;
      }

      dunestream<<"CPV true_dcp "<<tru_dcp<<" "<<min_chi_0pi<<" "<<min_chi_all<<" "<<min_chi_0pi-chi_sub<<" "<< std::endl;
    }
  }
  else if(which_mode =="cpv"){
    std::ofstream dunestream;
    dunestream.open ("DUNE_cpv_NO_4.dat");

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
    for(double tru_dcp = 0; tru_dcp<=360; tru_dcp+=15){

      //std::string truth_name = order_names.at(0)+"_DCP_"+to_string_prec(tru_dcp,3)+"_T23_"+to_string_prec(theta23.at(10),3);
      std::string truth_name = order_names.at(0)+"_DCP_"+to_string_prec(tru_dcp,3)+"_T23_"+to_string_prec(theta23.at(4),3);
      SBNspec * truth = new SBNspec(("precomp/"+truth_name+".SBNspec").c_str(),xml);
      truth->compressVector();	

      SBNchi mychi(*truth,*m);	
      std::vector<double> chi_all;
      std::vector<double> chi_0pi;

      //Here, vary the truth by poissonian_noise and build a distributuon for this. This is p-d/d then

      for(double dcp = 0; dcp<=360; dcp+=15){
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
  }




}
