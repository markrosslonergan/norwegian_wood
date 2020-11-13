#ifndef PLOTTINGTOOLS_H
#define PLOTTINGTOOLS_H

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

int calc_cpv_3p1( std::ofstream* dunestream,  std::string outfile, std::string xml, double t14, double t24, double t34, TMatrixT<double>*m);


int calc_neutrino_ordering( std::ofstream * dunestream,  std::string outfile, std::string xml);
int calc_neutrino_ordering_3p1( std::ofstream* dunestream,  std::string outfile, std::string xml, double t14, double t24, double t34, double d14);

template <typename T>
std::string to_string_prec(const T a_value, const int n = 6)
{
  std::ostringstream out;
  out <<std::fixed<< std::setprecision(n) << a_value;
  //what is std::fixed? // just returning the number with 6 digits.
  return out.str();
}




#endif

