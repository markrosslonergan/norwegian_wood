//genie_mu.c


//#include <iostream.h>

#include "TFile.h"
#include "TTree.h"
#include "TCanvas.h"
#include "TPad.h"
#include "TH1F.h"
#include "TH2F.h"

#include "THStack.h"
#include "TLegend.h"

#include "TMath.h"
#include "TLorentzVector.h"

#include "TRandom3.h"

#include "TRotation.h"

#include "TVector3.h"


#include "detector.h"



#include "TString.h"


#include "norm.h"

struct myphoton{
	TLorentzVector lorentz;
	int isPion;
	double convlength;
	double *convpoint;

};

struct mycharpion{
	TLorentzVector lorentz;
	double tracklength;

};

#include "functions.c"
#include "genie_tau_test.c" 

// genie_mu() : Similarly, for the ⌫μ CC sample, intrinsic beam ⌫μ CC events were assumed to be selected with an 80% recon- struction and identification e ciency. Potential back- ground contributions would result from NC ⇡± interac- tions where the ⇡± can be mis-identified as a muon. This background was mitigated by requiring that all contained muon-like tracks have a track length larger than 50 cm, and that all escaping tracks that have a track length of less than 1 m are rejected. This is the same methodology as what was followed in Ref. [23].


void bethe_test(){

	//cout << log(3) << endl;
	//cout << "bethe with betagamma 0.1: " << bethe(0.0995037) << endl;
	//cout << "bethe with betagamma 1: " << bethe(0.707107) << endl;
	//cout << "bethe with betagamma 10: " << bethe(0.995037) << endl;
	//cout << "bethe with betagamma 100: " << bethe(0.99995) << endl;
	TF1 *bethef = new TF1("bethe test", "bethe(x)", 0.01, 0.999);
	bethef->SetNpx(100000);
	TCanvas *c1 = new TCanvas("c1");
	c1->cd();
	//bethef->Draw();
	gPad->SetLogx();
	//c1->SaveAs("c1_1.pdf");

	double El = 10.;

	cout << "pion_track_length_test : "<< pion_track_length(El) << endl;

	TF1 f("CSDA", "CSDA_integrand(x)", 0.105, El);
	ROOT::Math::WrappedTF1 wf1(f);

	// Create the Integrator
	ROOT::Math::Integrator ig(ROOT::Math::IntegrationOneDim::kNONADAPTIVE);


	// Set parameters of the integration
	ig.SetFunction(wf1, false);
	ig.SetRelTolerance(0.001);
	double ans = ig.Integral(0.105,El);
	//std::cout << "integral result is " << ans <<std::endl;



}

void genie_mu(TString filename, int nu_mode, int nu_type){

	//want to try osc. + weak decay of tau.

	//given that gsttree gets the POT scaling
	double POT_norm = get_normalization(filename);
	std::cout<<"POT Normalization: "<<filename<<" "<<POT_norm<<std::endl;


	//cc
	//e->e
	//e->mu
	//e->tau
	//mu->e
	//mu->mu
	//mu->tau

	TRandom3 *rangen = new TRandom3(0);

	// TString filename = "gntp.0.numu50k_gst";

	// gntp.0.numuflux_nuebeam50k_gst
	// gntp.0.nutaubar20k_gst
	// gntp.0.nutau20k_gst
	// gntp.0.nue20k_gst
	// gntp.0.nuebar20k_gst
	// gntp.0.numu50k_gst
	// gntp.0.numubarflux_nuebarbeam10k_gst
	// gntp.0.numubar10k_gst


	TFile *f = new TFile(TString("../gst0to40/")+filename+TString(".root"));
	TTree *gstree = (TTree*)f->Get("gst");

	//TFile *fnew = new TFile(TString("../gst0to40/out_mu_100cm/")+filename+TString("_genie_mu_out.root"),"RECREATE");
	TFile *fnew = new TFile(TString("../gst0to40/out/")+filename+TString("_genie_mu_out.root"),"RECREATE");


	//std::vector<std::string> subchannels = {"nuefullosc","nuebarfullosc","intrinsic","ncmisid","numumisid","nutaumisid"};


	TFile * fntuple = new TFile("DUNE_ntuple.root","UPDATE");



	std::vector<TTree*> list_o_trees;

	std::vector<std::string> subchannels ;
	if(nu_mode == 0){ 
	subchannels = {"nu_dune_mulike_intrinsic",  "nu_dune_mulike_intrinsicbar", "nu_dune_mulike_taumisid", "nu_dune_mulike_taumisidbar","nu_dune_mulike_ncmisid"};

	} else if(nu_mode ==1){
	subchannels = {"nubar_dune_mulike_intrinsic",  "nubar_dune_mulike_intrinsicbar", "nubar_dune_mulike_taumisid", "nubar_dune_mulike_taumisidbar","nubar_dune_mulike_ncmisid"};

	}
	double m_Ereco=0;
	double m_Etrue=0;
	double m_l=0;
	double m_weight=0;
	int m_nutype=0;
	int m_ngen=0;



	for(auto s: subchannels){

		if(fntuple->GetListOfKeys()->Contains(s.c_str())){
			std::cout<<s<<" exists in file."<<std::endl;
			list_o_trees.push_back( (TTree*)fntuple->Get(s.c_str()) );	
			list_o_trees.back()->SetBranchAddress("Ereco",&m_Ereco );
			list_o_trees.back()->SetBranchAddress("Etrue",&m_Etrue );
			list_o_trees.back()->SetBranchAddress("L",&m_l);
			list_o_trees.back()->SetBranchAddress("Weight",&m_weight);
			list_o_trees.back()->SetBranchAddress("NuType",&m_nutype);
			list_o_trees.back()->SetBranchAddress("NGen",&m_ngen);


		}else{
			std::cout<<s<<" does not exist in file, creating it."<<std::endl;
			list_o_trees.push_back(  new TTree(s.c_str(), s.c_str())  );
			list_o_trees.back()->Branch("Ereco",&m_Ereco ,"Ereco/D");
			list_o_trees.back()->Branch("Etrue",&m_Etrue,"Etrue/D" );
			list_o_trees.back()->Branch("L",&m_l, "L/D" );
			list_o_trees.back()->Branch("Weight",&m_weight, "Weight/D");
			list_o_trees.back()->Branch("NuType",&m_nutype, "NuType/I");
			list_o_trees.back()->Branch("NGen",&m_ngen, "NGen/I");

		}
	}


	int neu;
	double energyf[100];
	double pf[100];
	int pdgf[100];
	double Q2;
	int nf;
	bool cc;
	bool nc;
	int fspl;
	double El;
	double pxl;
	double pyl;
	double pzl;
	int nfpi0;
	int Np;
	int nfpip;
	int nfpim;
	int Nph;
	int No;

	double Ev;

	double pxf[100];
	double pyf[100];
	double pzf[100];

	double Ep[100];
	double Epip[100];
	double Epim[100];
	int    pdgo[100];
	double Eo[100];
	double Eph[100];

	//double EM_thresh =0.03;
	double egamma_misidrate=0.06;
	double cc_efficiency_mu=0.9;
	double pion_mass = 0.13498;
	double convlength_thresh = 5.;//using 5cm for photon conversion length threshold.
	double NEUTRON_FACTOR = 0.6;


	gstree->SetBranchAddress("neu",&neu);
	gstree->SetBranchAddress("Ef",&energyf);
	gstree->SetBranchAddress("Ev",&Ev);
	gstree->SetBranchAddress("pf",&pf);
	gstree->SetBranchAddress("pxf",&pxf);
	gstree->SetBranchAddress("pyf",&pyf);
	gstree->SetBranchAddress("pzf",&pzf);
	gstree->SetBranchAddress("pdgf",&pdgf);
	gstree->SetBranchAddress("Q2",&Q2);
	gstree->SetBranchAddress("nf",&nf);
	gstree->SetBranchAddress("cc",&cc);
	gstree->SetBranchAddress("nc",&nc);
	gstree->SetBranchAddress("fspl",&fspl);//final state pdg lepton?
	gstree->SetBranchAddress("El",&El);
	gstree->SetBranchAddress("pxl",&pxl);
	gstree->SetBranchAddress("pyl",&pyl);
	gstree->SetBranchAddress("pzl",&pzl);

	gstree->SetBranchAddress("nfpip",&nfpip);
	gstree->SetBranchAddress("nfpim",&nfpim);



	SBN_detector DUNE(4,false);

	TH1F *hist_cc_reco = new TH1F("hist_cc_reco", "cc_reco", 80, 0., 20.);
	hist_cc_reco->Sumw2();

	TH1F *hist_nc_reco = new TH1F("hist_nc_reco", "nc_reco", 80, 0., 20.);
	hist_nc_reco->Sumw2();


	TH1F *hist_cc_reco_true = new TH1F("hist_cc_reco_true", "cc_reco_true", 80, 0., 20.);
	hist_cc_reco_true->Sumw2();

	TH1F *hist_cc_ehad = new TH1F("hist_cc_ehad", "cc_ehad", 80, 0., 20.);
	hist_cc_ehad->Sumw2();

	TH1F *hist_cc_ehad_true = new TH1F("hist_cc_ehad_true", "cc_ehad_true", 80, 0., 20.);
	hist_cc_ehad_true->Sumw2();

	TH1F *hist_cc_El_true = new TH1F("hist_cc_El_true", "cc_El_true", 80, 0., 20.);
	hist_cc_El_true->Sumw2();

	TH1F *hist_cc_El_smear = new TH1F("hist_cc_El_smear", "cc_El_smear", 80, 0., 20.);
	hist_cc_El_smear->Sumw2();

	TH1F *hist_nc_num_charpi = new TH1F("hist_nc_num_charpi", "nc_num_charpi",5,0.,5. );


	TH1F *hist_cc_reco_effi = new TH1F("hist_cc_reco_effi", "cc_reco_effi", 80, 0., 20.);
	hist_cc_reco_effi->Sumw2();

	TH1F *hist_nc_charpi = new TH1F("hist_nc_charpi", "nc_charpi", 80, 0., 5.);
	hist_nc_charpi->Sumw2();

	TH1F *hist_nc_charpi_trackl = new TH1F("hist_nc_charpi_trackl", "nc_charpi_trackl", 80, 0., 200.);

	TH2F *hist_nc_charpi_e_vs_trackl = new TH2F("hist_nc_charpi_e_vs_trackl", "nc_charpi_e_vs_trackl", 20, 0., 1.,20, 0., 1000.);
	//hist_nc_charpi->Sumw2();
	TH1F *hist_nc_mucharpi = new TH1F("hist_nc_mucharpi", "nc_mucharpi", 80, 0., 5.);
	hist_nc_mucharpi->Sumw2();

	for(int i=0; i < gstree->GetEntries(); i++){
		gstree->GetEntry(i);
		m_ngen = gstree->GetEntries();

		bool mode = true;
		if (neu<0){
			mode = false;
		}


		if (nc && (nfpim+nfpip==1)) {

			// if NC, pi+/- mimics muon and will look like one electron.
			// then sort out pi+/-


			double Enu_reco = 0;
			double Ehad=0;
			double Ehad_true=0.;

			double Pxhad =0;
			double Pyhad =0;
			double Pzhad =0;

			bool vis_vertex = false;

			std::vector<mycharpion> charpi;

			int charged_pion_count = 0.;

			double vertex_pos[3];

			DUNE.random_pos(rangen, vertex_pos);

			if(nf!=0){




				double kin_true =0;
				double kin_smeared = 0;

				//hadronic energy - pion energy
				for(int j=0; j<nf; j++){

					if(pdgf[j]==2212){
						kin_true=energyf[j]-MPROTON;
						kin_smeared=smear_energy_type(pdgf[j], kin_true, rangen);
						if (kin_smeared > p_thresh){
							Ehad += kin_smeared;
							Ehad_true += kin_true;
							Pxhad += pxf[j];
							Pyhad += pyf[j];
							Pzhad += pzf[j];
						}
					}
					//else if(abs(pdgf[j])==211){
					//    kin_true=energyf[j]-MPION;
					//    kin_smeared=smear_energy_type(pdgf[j], kin_true, 1, rangen);//what is '1'?-> fully contained.
					//    if (kin_smeared > pip_thresh) {
					//        Ehad += kin_smeared+MPION;

					//        Ehad_true += kin_true+MPION;

					//        Pxhad += pxf[j];
					//        Pyhad += pyf[j];
					//        Pzhad += pzf[j];


					//    }
					//}
					//else if(pdgf[j]==22){
					//   kin_true=energyf[j];
					//    kin_smeared=smear_energy(kin_true, EMsmear, rangen);
					//    Ehad += kin_smeared;
					//}
					else if(pdgf[j]==2112){// do i need to add neutron mass? or
						//cout << " neutron mass : " << sqrt(energyf[j]*energyf[j] - pxf[j]*pxf[j]-pyf[j]*pyf[j]-pzf[j]*pzf[j])<< endl;
						kin_true=energyf[j]-MNEUTRON;
						//  kin_smeared=smear_energy(kin_true, nsmear, rangen);
						kin_smeared=smear_energy_type(pdgf[j], kin_true, rangen);
						if (kin_smeared > n_thresh) {
							Ehad += NEUTRON_FACTOR*kin_smeared;

							Ehad_true += NEUTRON_FACTOR*kin_true;

							Pxhad += pxf[j];
							Pyhad += pyf[j];
							Pzhad += pzf[j];
						}//define NEUTRO FACTOR = 0.6 later.//mesons we add rest mass+kinetic energy, baryons we add kinetic energy
					}
					else if(pdgf[j]==321 || pdgf[j]==-321){
						//cout << " kaon321 mass : " << sqrt(energyf[j]*energyf[j] - pxf[j]*pxf[j]-pyf[j]*pyf[j]-pzf[j]*pzf[j])<< endl;
						kin_true=energyf[j]-MKAON321;
						// kin_smeared=smear_energy(kin_true, osmear, rangen);
						kin_smeared=smear_energy_type(pdgf[j], kin_true, rangen);
						Ehad += kin_smeared+MKAON321;

						Ehad_true += kin_true+MKAON321;

						Pxhad += pxf[j];
						Pyhad += pyf[j];
						Pzhad += pzf[j];
					}
					else if(pdgf[j]==311){
						//cout << " kaon311 mass : " << sqrt(energyf[j]*energyf[j] - pxf[j]*pxf[j]-pyf[j]*pyf[j]-pzf[j]*pzf[j])<< endl;
						kin_true=energyf[j]-MKAON311;
						// kin_smeared=smear_energy(kin_true, osmear, rangen);
						kin_smeared=smear_energy_type(pdgf[j], kin_true, rangen);
						Ehad += kin_smeared+MKAON311;

						Ehad_true += kin_true+MKAON311;

						Pxhad += pxf[j];
						Pyhad += pyf[j];
						Pzhad += pzf[j];
					}
					else if(pdgf[j]==3222){
						//cout << " sigma mass3222 : " << sqrt(energyf[j]*energyf[j] - pxf[j]*pxf[j]-pyf[j]*pyf[j]-pzf[j]*pzf[j])<< endl;
						kin_true=energyf[j]-MSIGMA3222;
						// kin_smeared=smear_energy(kin_true, osmear, rangen);
						kin_smeared=smear_energy_type(pdgf[j], kin_true, rangen);
						Ehad += kin_smeared+MSIGMA3222;
						Ehad_true += kin_true+MSIGMA3222;

						Pxhad += pxf[j];
						Pyhad += pyf[j];
						Pzhad += pzf[j];
					}
					else if(pdgf[j]==3112){
						//cout << " sigma mass3112 : " << sqrt(energyf[j]*energyf[j] - pxf[j]*pxf[j]-pyf[j]*pyf[j]-pzf[j]*pzf[j])<< endl;
						kin_true=energyf[j]-MSIGMA3112;
						// kin_smeared=smear_energy(kin_true, osmear, rangen);
						kin_smeared=smear_energy_type(pdgf[j], kin_true, rangen);
						Ehad += kin_smeared+MSIGMA3112;
						Ehad_true += kin_true+MSIGMA3112;

						Pxhad += pxf[j];
						Pyhad += pyf[j];
						Pzhad += pzf[j];
					}
					else if(pdgf[j]==3122){
						//cout << " sigma mass3122 : " << sqrt(energyf[j]*energyf[j] - pxf[j]*pxf[j]-pyf[j]*pyf[j]-pzf[j]*pzf[j])<< endl;
						kin_true=energyf[j]-MSIGMA3122;
						// kin_smeared=smear_energy(kin_true, osmear, rangen);
						kin_smeared=smear_energy_type(pdgf[j], kin_true, rangen);
						Ehad += kin_smeared+MSIGMA3122;
						Ehad_true += kin_true+MSIGMA3122;

						Pxhad += pxf[j];
						Pyhad += pyf[j];
						Pzhad += pzf[j];
					}
				}



				for(int k=0; k<nf; k++){

					if (abs(pdgf[k])==211){

						double track_length = pion_track_length(energyf[k]);
						double endpos[3] = {0,0,0};
						double temp_3vec[3] = {pxf[k],pyf[k],pzf[k]};
						get_endpoint(vertex_pos, track_length, temp_3vec, endpos);

						double charged_pion_smear = smear_energy_type(211, DUNE.is_fully_contained(vertex_pos, endpos), rangen);

						if(charged_pion_smear>MU_thresh){


							hist_nc_charpi_trackl->Fill(track_length);
							hist_nc_charpi->Fill(energyf[k]);

							hist_nc_charpi_e_vs_trackl->Fill(energyf[k], track_length);
							//std::cout<<energyf[j]<<" "<<track_length<<std::endl;

							double observable_L = 0.;

							if( DUNE.is_fully_contained(vertex_pos, endpos) ) {
								observable_L=track_length;
							}
							else {
								observable_L = DUNE.track_length_escape(vertex_pos,endpos);
							}


							if ( observable_L > 100.0 && DUNE.is_fiducial(vertex_pos) ){// !(DUNE.is_fully_contained(vertex_pos, endpos) && pion_track_length>50.0) &&
								charged_pion_count++;

								mycharpion temp_charpi;
								temp_charpi.lorentz.SetPxPyPzE(pxf[k],pyf[k],pzf[k],energyf[k]);
								temp_charpi.tracklength=observable_L;

								charpi.push_back(temp_charpi);

								hist_nc_mucharpi->Fill(energyf[k]);




								//Fill ttrees


							}

						}

					}
				}

			}

			hist_nc_num_charpi->Fill(charged_pion_count);

			//Check if we actually have a "visibe vertex"
			if(Ehad >= vertex_thresh)
			{
				vis_vertex = true;
			}

			double reco_E = 0.;

			if (charged_pion_count==1){
				double reco_E = charpi.at(0).lorentz.E();
				if (vis_vertex){
					reco_E+=Ehad;
				}
				hist_nc_reco->Fill(reco_E,cc_efficiency_mu);

				// fill ttrees
				m_Ereco = reco_E;
				m_Etrue = Ev;
				m_l = 1300 ;	
				m_weight = cc_efficiency_mu*POT_norm;
				m_nutype = neu;
				//std::cout<<"NCC: "<<m_Ereco<<" "<<m_Etrue<<" "<<m_l<<" "<<m_weight<<" "<<m_nutype<<std::endl;
				if( nu_mode == 0 && (filename == "gntp.0.numu50k_gst" || filename == "gntp.0.numubar10k_gst" )){
					list_o_trees.at(4)->Fill();
				}else if(nu_mode ==1 && (filename =="gntp.0.RHC_FD_numuflux_numubeam20k_gst" || filename == "gntp.0.RHC_FD_numubarflux_numubarbeam50k_gst")){
					list_o_trees.at(4)->Fill();
				}



			}

		}

		//CC
		if (cc) {

			//Is there a visible vertex and how much energy is there!
			double Enu_reco = 0;
			double Ehad=0;
			double Ehad_true=0.;

			double Pxhad =0;
			double Pyhad =0;
			double Pzhad =0;

			bool vis_vertex = false;


			if(nf!=0){
				double kin_true =0;
				double kin_smeared = 0;

				for(int j=0; j<nf; j++){

					if(pdgf[j]==2212){
						kin_true=energyf[j]-MPROTON;
						kin_smeared=smear_energy_type(pdgf[j], kin_true, rangen);
						if (kin_smeared > p_thresh){
							Ehad += kin_smeared;
							Ehad_true += kin_true;
							Pxhad += pxf[j];
							Pyhad += pyf[j];
							Pzhad += pzf[j];
						}
					}
					else if(abs(pdgf[j])==211){
						kin_true=energyf[j]-MPION;
						kin_smeared=smear_energy_type(pdgf[j], kin_true, 1, rangen);//what is '1'?-> fully contained.
						if (kin_smeared > pip_thresh) {
							Ehad += kin_smeared+MPION;

							Ehad_true += kin_true+MPION;

							Pxhad += pxf[j];
							Pyhad += pyf[j];
							Pzhad += pzf[j];
						}
					}
					//else if(pdgf[j]==22){
					//   kin_true=energyf[j];
					//    kin_smeared=smear_energy(kin_true, EMsmear, rangen);
					//    Ehad += kin_smeared;
					//}
					else if(pdgf[j]==2112){// do i need to add neutron mass? or
						//cout << " neutron mass : " << sqrt(energyf[j]*energyf[j] - pxf[j]*pxf[j]-pyf[j]*pyf[j]-pzf[j]*pzf[j])<< endl;
						kin_true=energyf[j]-MNEUTRON;
						//  kin_smeared=smear_energy(kin_true, nsmear, rangen);
						kin_smeared=smear_energy_type(pdgf[j], kin_true, rangen);
						if (kin_smeared > n_thresh) {
							Ehad += NEUTRON_FACTOR*kin_smeared;

							Ehad_true += NEUTRON_FACTOR*kin_true;

							Pxhad += pxf[j];
							Pyhad += pyf[j];
							Pzhad += pzf[j];
						}//define NEUTRO FACTOR = 0.6 later.//mesons we add rest mass+kinetic energy, baryons we add kinetic energy
					}
					else if(pdgf[j]==321 || pdgf[j]==-321){
						//cout << " kaon321 mass : " << sqrt(energyf[j]*energyf[j] - pxf[j]*pxf[j]-pyf[j]*pyf[j]-pzf[j]*pzf[j])<< endl;
						kin_true=energyf[j]-MKAON321;
						// kin_smeared=smear_energy(kin_true, osmear, rangen);
						kin_smeared=smear_energy_type(pdgf[j], kin_true, rangen);
						Ehad += kin_smeared+MKAON321;

						Ehad_true += kin_true+MKAON321;

						Pxhad += pxf[j];
						Pyhad += pyf[j];
						Pzhad += pzf[j];
					}
					else if(pdgf[j]==311){
						//cout << " kaon311 mass : " << sqrt(energyf[j]*energyf[j] - pxf[j]*pxf[j]-pyf[j]*pyf[j]-pzf[j]*pzf[j])<< endl;
						kin_true=energyf[j]-MKAON311;
						// kin_smeared=smear_energy(kin_true, osmear, rangen);
						kin_smeared=smear_energy_type(pdgf[j], kin_true, rangen);
						Ehad += kin_smeared+MKAON311;

						Ehad_true += kin_true+MKAON311;

						Pxhad += pxf[j];
						Pyhad += pyf[j];
						Pzhad += pzf[j];
					}
					else if(pdgf[j]==3222){
						//cout << " sigma mass3222 : " << sqrt(energyf[j]*energyf[j] - pxf[j]*pxf[j]-pyf[j]*pyf[j]-pzf[j]*pzf[j])<< endl;
						kin_true=energyf[j]-MSIGMA3222;
						// kin_smeared=smear_energy(kin_true, osmear, rangen);
						kin_smeared=smear_energy_type(pdgf[j], kin_true, rangen);
						Ehad += kin_smeared+MSIGMA3222;
						Ehad_true += kin_true+MSIGMA3222;

						Pxhad += pxf[j];
						Pyhad += pyf[j];
						Pzhad += pzf[j];
					}
					else if(pdgf[j]==3112){
						//cout << " sigma mass3112 : " << sqrt(energyf[j]*energyf[j] - pxf[j]*pxf[j]-pyf[j]*pyf[j]-pzf[j]*pzf[j])<< endl;
						kin_true=energyf[j]-MSIGMA3112;
						// kin_smeared=smear_energy(kin_true, osmear, rangen);
						kin_smeared=smear_energy_type(pdgf[j], kin_true, rangen);
						Ehad += kin_smeared+MSIGMA3112;
						Ehad_true += kin_true+MSIGMA3112;

						Pxhad += pxf[j];
						Pyhad += pyf[j];
						Pzhad += pzf[j];
					}
					else if(pdgf[j]==3122){
						//cout << " sigma mass3122 : " << sqrt(energyf[j]*energyf[j] - pxf[j]*pxf[j]-pyf[j]*pyf[j]-pzf[j]*pzf[j])<< endl;
						kin_true=energyf[j]-MSIGMA3122;
						// kin_smeared=smear_energy(kin_true, osmear, rangen);
						kin_smeared=smear_energy_type(pdgf[j], kin_true, rangen);
						Ehad += kin_smeared+MSIGMA3122;
						Ehad_true += kin_true+MSIGMA3122;

						Pxhad += pxf[j];
						Pyhad += pyf[j];
						Pzhad += pzf[j];
					}
				}
			}


			//Check if we actually have a "visibe vertex"
			if(Ehad >= vertex_thresh)
			{
				vis_vertex = true;
			}



			//Nrevertex. for stats.
			int Nrevertex = 1;
			for(int k=0;k<Nrevertex; k++){

				double vertex_weight = 1.0/((double) Nrevertex);


				double vertex_pos[3];

				DUNE.random_pos(rangen, vertex_pos);



				// intrinsic
				if (cc && abs(neu)==14) {

					//double pT = TMath::Sqrt(Pxtot*Pxtot + Pytot*Pytot);
					//double pL = Pztot;


					double Lmu = muon_track_length(El);
					double endpos[3] = {0,0,0};
					double temp_3vec[3] = {pxl,pyl,pzl};
					get_endpoint(vertex_pos, Lmu, temp_3vec, endpos);

					double El_smear = smear_energy_type(13, El, DUNE.is_fully_contained(vertex_pos, endpos), rangen);
					//double El_smear= El;

					double observable_L = 0.;

					if( DUNE.is_fully_contained(vertex_pos, endpos) ) {
						observable_L=Lmu;
					}
					else {
						observable_L = DUNE.track_length_escape(vertex_pos,endpos);
					}

					double reco_E = El_smear;

					if (vis_vertex){
						reco_E += Ehad;
					}

					if (El_smear > EM_thresh && DUNE.is_fiducial(vertex_pos) && observable_L>100.0 ) {// && vis_vertex


						//MissingEnergy.SetXYZM(-Pxtot,-Pytot,-Pztot,0.);

						//hist_PT->Fill(MissingEnergy.Et(),cc_efficiency*vertex_weight);

						//hist_lepton_PT->Fill(sqrt(pxl*pxl+pyl*pyl),cc_efficiency*vertex_weight);

						//hist_PT_vs_lepton_PT->Fill(MissingEnergy.Et(),sqrt(pxl*pxl+pyl*pyl),cc_efficiency*vertex_weight);


						//hist_cc_El->Fill(El_smear,cc_efficiency*vertex_weight);
						//hist_cc_ehad->Fill(Ehad,cc_efficiency*vertex_weight);
						//hist_cc_reco->Fill(reco_E,cc_efficiency*vertex_weight*prob_muflav(Ev,2));//prob_muflav(Ev,1)
						hist_cc_reco->Fill(reco_E,cc_efficiency_mu*vertex_weight*prob_muflav(Ev,2,mode));

						hist_cc_ehad->Fill(Ehad,cc_efficiency_mu*vertex_weight*prob_muflav(Ev,2,mode));
						hist_cc_ehad_true->Fill(Ehad_true,cc_efficiency_mu*vertex_weight*prob_muflav(Ev,2,mode));

						hist_cc_El_true->Fill(El,cc_efficiency_mu*vertex_weight*prob_muflav(Ev,2,mode));
						hist_cc_El_smear->Fill(El_smear,cc_efficiency_mu*vertex_weight*prob_muflav(Ev,2,mode));


						hist_cc_reco_true->Fill(El+Ehad_true,cc_efficiency_mu*vertex_weight*prob_muflav(Ev,2,mode));

						//Fill ttrees
						m_Ereco = reco_E;
						m_Etrue = Ev;
						m_l = 1300 ;	
						m_weight = cc_efficiency_mu*vertex_weight*POT_norm;
						m_nutype = neu;
						
						//std::cout<<"INT: "<<m_Ereco<<" "<<m_Etrue<<" "<<m_l<<" "<<m_weight<<" "<<m_nutype<<std::endl;
						if(nu_type >0){
							list_o_trees.at(0)->Fill();
						} else {

							list_o_trees.at(1)->Fill();
						}


						//hist_cc_El_true->Fill(El,cc_efficiency*vertex_weight);
						//hist_cc_ehad_true->Fill(Ehad_true,cc_efficiency*vertex_weight);
						//hist_cc_reco_true->Fill(El+Ehad_true,cc_efficiency*vertex_weight);

						TLorentzVector truelvec;

						truelvec.SetPxPyPzE(pxl,pyl,pzl,El);

						//hist_cc_ElCosTheta_true->Fill(truelvec.CosTheta(),cc_efficiency*vertex_weight);


						/*hist_PT_vs_PL->Fill(MissingEnergy.Pt(),Pztot);
						  hist_PT->Fill(MissingEnergy.Pt());
						  hist_PL->Fill(Pztot);
						  hist_PLsqOverPTsq->Fill(Pztot*Pztot/MissingEnergy.Pt()/MissingEnergy.Pt());
						  hist_PTsqOverPLsq->Fill(MissingEnergy.Pt()*MissingEnergy.Pt()/Pztot/Pztot);
						  hist_PTsqOverPL->Fill(MissingEnergy.Pt()*MissingEnergy.Pt()/Pztot);
						  hist_PTsq->Fill(MissingEnergy.Pt()*MissingEnergy.Pt());
						  hist_PTcube->Fill(MissingEnergy.Pt()*MissingEnergy.Pt()*MissingEnergy.Pt());
						  hist_PTquartic->Fill(MissingEnergy.Pt()*MissingEnergy.Pt()*MissingEnergy.Pt()*MissingEnergy.Pt());


						  hist_PT_vs_PL_weak->Fill(MissingEnergy.Pt(),Pztot);
						  hist_PT_weak->Fill(MissingEnergy.Pt());
						  hist_PL_weak->Fill(Pztot);
						  hist_PLsqOverPTsq_weak->Fill(Pztot*Pztot/MissingEnergy.Pt()/MissingEnergy.Pt());
						  hist_PTsqOverPLsq_weak->Fill(MissingEnergy.Pt()*MissingEnergy.Pt()/Pztot/Pztot);
						  hist_PTsqOverPL_weak->Fill(MissingEnergy.Pt()*MissingEnergy.Pt()/Pztot);
						  hist_PTsq_weak->Fill(MissingEnergy.Pt()*MissingEnergy.Pt());
						  hist_PTcube_weak->Fill(MissingEnergy.Pt()*MissingEnergy.Pt()*MissingEnergy.Pt());
						  hist_PTquartic_weak->Fill(MissingEnergy.Pt()*MissingEnergy.Pt()*MissingEnergy.Pt()*MissingEnergy.Pt());*/




					}




				}// end of numu cc.


				if (cc && abs(neu)==16) {

					double vertex_pos[3];
					DUNE.random_pos(rangen, vertex_pos);




					//double

					TLorentzVector daughter[3];
					daughter[0].SetPxPyPzE(0.0,0.0,0.0,0.0);
					daughter[1].SetPxPyPzE(0.0,0.0,0.0,0.0);
					daughter[2].SetPxPyPzE(0.0,0.0,0.0,0.0);



					double tau_mass = 1.76;

					TLorentzVector parentTau;
					parentTau.SetPxPyPzE(0.0,0.0,0.0,tau_mass);

					DecayIt(parentTau, rangen, daughter);

					TLorentzVector daughter12, daughter23;

					daughter12 = daughter[2]+daughter[1];
					daughter23 = daughter[1]+daughter[0];



					//tau vector in lab frame
					TLorentzVector labtau;
					labtau.SetPxPyPzE(pxl,pyl,pzl,El);

					daughter[0].Boost(pxl/El,pyl/El,pzl/El);//daughter 0 : muon
					daughter[1].Boost(pxl/El,pyl/El,pzl/El);
					daughter[2].Boost(pxl/El,pyl/El,pzl/El);

					TLorentzVector twoneutrinos = daughter[1]+daughter[2];

					double track_length = muon_track_length(daughter[0].E());
					double endpos[3] = {0,0,0};
					double temp_3vec[3] = {daughter[0].Px(),daughter[0].Py(),daughter[0].Pz()};
					get_endpoint(vertex_pos, track_length, temp_3vec, endpos);


					double El_smear;

					El_smear = smear_energy_type(13, daughter[0].E(),DUNE.is_fully_contained(vertex_pos, endpos), rangen);

					double reco_E = El_smear;

					if (vis_vertex){
						reco_E += Ehad;
					}

					hist_cc_reco_effi->Fill(reco_E,prob_muflav(Ev,3,mode));//prob_muflav(Ev,3)*

					double observable_L = 0.;

					if( DUNE.is_fully_contained(vertex_pos, endpos) ) {
						observable_L=track_length;
					}
					else {
						observable_L = DUNE.track_length_escape(vertex_pos,endpos);
					}


					if (El_smear > EM_thresh && DUNE.is_fiducial(vertex_pos) && observable_L>100.0) {//&& vis_vertex&&

						//Pxtot+=daughter[0].Px();
						//Pytot+=daughter[0].Py();
						//Pztot+=daughter[0].Pz();

						//double pT = TMath::Sqrt(Pxtot*Pxtot + Pytot*Pytot);
						//double pL = Pztot;

						hist_cc_reco->Fill(reco_E,prob_muflav(Ev,3,mode)*cc_efficiency_mu*vertex_weight*0.1739);


						//Fill ttrees
						m_Ereco = reco_E;
						m_Etrue = Ev;
						m_l = 1300 ;	
						m_weight = cc_efficiency_mu*vertex_weight*0.1739*POT_norm;
						m_nutype = neu;
						//std::cout<<"TAU: "<<m_Ereco<<" "<<m_Etrue<<" "<<m_l<<" "<<m_weight<<" "<<m_nutype<<std::endl;
						if(nu_type>0){

						list_o_trees.at(2)->Fill();
						}else{
						list_o_trees.at(3)->Fill();
						}



						/*hist_PT_vs_PL_weak->Fill(MissingEnergy.Pt(),Pztot,prob_muflav(Ev,3));
						  hist_PT_weak->Fill(MissingEnergy.Pt(),prob_muflav(Ev,3));
						  hist_PL_weak->Fill(Pztot,prob_muflav(Ev,3));
						  hist_PLsqOverPTsq_weak->Fill(Pztot*Pztot/MissingEnergy.Pt()/MissingEnergy.Pt(),prob_muflav(Ev,3));
						  hist_PTsqOverPLsq_weak->Fill(MissingEnergy.Pt()*MissingEnergy.Pt()/Pztot/Pztot,prob_muflav(Ev,3));
						  hist_PTsqOverPL_weak->Fill(MissingEnergy.Pt()*MissingEnergy.Pt()/Pztot,prob_muflav(Ev,3));
						  hist_PTsq_weak->Fill(MissingEnergy.Pt()*MissingEnergy.Pt(),prob_muflav(Ev,3));
						  hist_PTcube_weak->Fill(MissingEnergy.Pt()*MissingEnergy.Pt()*MissingEnergy.Pt(),prob_muflav(Ev,3));
						  hist_PTquartic_weak->Fill(MissingEnergy.Pt()*MissingEnergy.Pt()*MissingEnergy.Pt()*MissingEnergy.Pt(),prob_muflav(Ev,3));
						 */

						//Pxtot-=daughter[0].Px();
						//Pytot-=daughter[0].Py();
						//Pztot-=daughter[0].Pz();

					}


				}// end of nutau cc.
			}

		}
	}
	fnew->cd();
	fnew->Write();

	fntuple->cd();
	for(auto &t: list_o_trees){
		t->Write();
	}
	fntuple->Purge();
	fntuple->Close();
}

void run_all_genie_mu(){

	TString nutaubar = "gntp.0.nutaubar20k_gst";
	TString nutau = "gntp.0.nutau20k_gst";
	TString nue  = "gntp.0.nue20k_gst";
	TString nuebar = "gntp.0.nuebar20k_gst";
	TString numu  = "gntp.0.numu50k_gst";
	TString numubar = "gntp.0.numubar10k_gst";

	TString numu_nuebeam = "gntp.0.numuflux_nuebeam50k_gst";
	TString numubar_nuebarbeam = "gntp.0.numubarflux_nuebarbeam10k_gst";



	std::cout<<"Starting CC nue."<<std::endl;
	genie_mu(nue,0,1);
	std::cout<<"Starting CC numu."<<std::endl;
	genie_mu(numu,0,1);
	std::cout<<"Starting CC nutau."<<std::endl;
	genie_mu(nutau,0,1);

	std::cout<<"Starting CC nuebar."<<std::endl;
	genie_mu(nuebar,0,-1);
	std::cout<<"Starting CC numubar."<<std::endl;
	genie_mu(numubar,0,-1);
	std::cout<<"Starting CC nutaubar."<<std::endl;
	genie_mu(nutaubar,0,-1);

	std::cout<<"Starting wierd CC numu_nuebeam."<<std::endl;
	genie_mu(numu_nuebeam,0,1);
	std::cout<<"Starting wierd CC numubear_nuebarbeam."<<std::endl;
	genie_mu(numubar_nuebarbeam,0,-1);


	TString BHCnutaubar = "gntp.0.RHC_FD_numubarflux_nutaubarbeam20k_gst";
	TString BHCnutau = "gntp.0.RHC_FD_numuflux_nutaubeam10k_gst";

	TString BHCnue  = "gntp.0.RHC_FD_nueflux_nuebeam10k_gst";
	TString BHCnuebar = "gntp.0.RHC_FD_nuebarflux_nuebarbeam10k_gst";

	TString BHCnumu  = "gntp.0.RHC_FD_numuflux_numubeam20k_gst";
	TString BHCnumubar = "gntp.0.RHC_FD_numubarflux_numubarbeam50k_gst";

	TString BHCnumu_BHCnuebeam  = "gntp.0.RHC_FD_numuflux_nuebeam10k_gst";
	TString BHCnumubar_BHCnuebarbeam = "gntp.0.RHC_FD_numubarflux_nuebarbeam20k_gst";


	std::cout<<"Starting CC BHCnue."<<std::endl;
	genie_mu(BHCnue,1,1);
	std::cout<<"Starting CC BHCnumu."<<std::endl;
	genie_mu(BHCnumu,1,1);
	std::cout<<"Starting CC BHCnutau."<<std::endl;
	genie_mu(BHCnutau,1,1);

	std::cout<<"Starting CC BHCnuebar."<<std::endl;
	genie_mu(BHCnuebar,1,-1);
	std::cout<<"Starting CC BHCnumubar."<<std::endl;
	genie_mu(BHCnumubar,1,-1);
	std::cout<<"Starting CC BHCnutaubar."<<std::endl;
	genie_mu(BHCnutaubar,1,-1);

	std::cout<<"Starting wierd CC BHCnumu_BHCnuebeam."<<std::endl;
	genie_mu(BHCnumu_BHCnuebeam,1,1);
	std::cout<<"Starting wierd CC BHCnumubear_BHCnuebarbeam."<<std::endl;
	genie_mu(BHCnumubar_BHCnuebarbeam,1,-1);





}

