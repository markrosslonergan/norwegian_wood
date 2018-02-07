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

struct myphoton{
	TLorentzVector lorentz;
	int isPion;
	double convlength;
	double *convpoint;

};

#include "functions.c"
#include "genie_tau_test.c" // "genie_tau_test.c" contains a function which returns numu->nutau osc. probability.



//genie_study -> temporary checking out plots

//genie_CC_slimmed -> essencial plots

//genie_CC


void genie_study(TString filename){

	//want to try osc. + weak decay of tau.

	//cc
	//e->e
	//e->mu
	//e->tau
	//mu->e
	//mu->mu
	//mu->tau

	TRandom3 *rangen = new TRandom3(0);


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

	TFile *fnew = new TFile(TString("../gst0to40/out/")+filename+TString("_study_out.root"),"recreate");

	//Fullosc sample is taken from the numu flux ran through GENIE as a nue beam
	bool isfullosc=false;
	if(filename == "gntp.0.numuflux_nuebeam50k_gst" || filename == "gntp.0.numubarflux_nuebarbeam50k_gst"){
		isfullosc=true;
	}

	TFile * fntuple = new TFile("DUNE_ntuple.root","UPDATE");

	std::vector<std::string> subchannels = {"nu_dune_elike_fullosc", "nu_dune_elike_intrinsic" , "nu_dune_elike_mumisid", "nu_dune_elike_taumisid"};
	std::vector<TTree*> list_o_trees;

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
	int Npip;
	int Npim;
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
	double cc_efficiency=0.8;
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



	SBN_detector DUNE(4,false);

	TH1F *hist_cc_reco = new TH1F("hist_cc_reco", "cc_reco", 80, 0., 20.);
	hist_cc_reco->Sumw2();

	TH1F *hist_cc_ehad = new TH1F("hist_cc_ehad", "cc_ehad", 80, 0., 20.);
	hist_cc_ehad->Sumw2();

	TH1F *hist_cc_reco_effi = new TH1F("hist_cc_reco_effi", "cc_reco_effi", 80, 0., 20.);
	hist_cc_reco_effi->Sumw2();

	for(int i=0; i < gstree->GetEntries(); i++){
		gstree->GetEntry(i);

		if (nc) {
			continue;
		}

		bool mode = true;
		if (neu<0){
			mode = false;
		}

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
		int Nrevertex = 25;
		for(int k=0;k<Nrevertex; k++){

			double vertex_weight = 1.0/((double) Nrevertex);


			double vertex_pos[3];

			DUNE.random_pos(rangen, vertex_pos);

			//photon flow
			std::vector<myphoton> gammas;
			std::vector<myphoton> gammas_true;


			//radiation gammas
			if (nf!=0){
				for(int j=0; j<nf;j++ ){

					if(pdgf[j]==22){ //energyf[j]>EM_thresh

						myphoton temp_vec;
						temp_vec.lorentz.SetPxPyPzE(pxf[j],pyf[j],pzf[j],energyf[j]);


						double temp_energy_smeared = smear_energy(energyf[j], EMsmear,rangen);

						if (temp_energy_smeared > EM_thresh) {

							gammas_true.push_back(temp_vec);

							temp_vec.lorentz.SetE(temp_energy_smeared);
							gammas.push_back(temp_vec);

							gammas.back().isPion=0;

						}
					}
				}
			}

			//pi0 gammas
			if (nf!=0){

				for(int k=0;k<nf;k++){

					if(pdgf[k]==111){

						double boost[3] = {pxf[k]/energyf[k],pyf[k]/energyf[k],pzf[k]/energyf[k]};//boost vector for pi0

						double x=0;
						double y=0;
						double z=0;

						rangen->Sphere(x,y,z,1.);//generate random direction on a unit sphere.

						myphoton gamma1, gamma2;


						gamma1.lorentz.SetPxPyPzE(x*pion_mass/2,y*pion_mass/2,z*pion_mass/2,pion_mass/2);// assigning gamma1 along the random direction
						gamma2.lorentz.SetPxPyPzE(-x*pion_mass/2,-y*pion_mass/2,-z*pion_mass/2,pion_mass/2);// assigning gamma2 to the opposite directon of gamma1, each of these has a half of pion rest energy.


						gamma1.lorentz.Boost(pxf[k]/energyf[k],pyf[k]/energyf[k],pzf[k]/energyf[k]);//boost gamma1 and gamma2 to original pi0 frame.
						gamma2.lorentz.Boost(pxf[k]/energyf[k],pyf[k]/energyf[k],pzf[k]/energyf[k]);
						// check invariant amss of gamma1 gamm2 pair ==mpion

						double temp_energy_smeared1= smear_energy(gamma1.lorentz.E(), EMsmear,rangen);
						double temp_energy_smeared2= smear_energy(gamma2.lorentz.E(), EMsmear,rangen);

						if (temp_energy_smeared1 > EM_thresh) {
							gamma1.isPion=1;
							gammas_true.push_back(gamma1);

							gamma1.lorentz.SetE(temp_energy_smeared1);
							gammas.push_back(gamma1);
							//pion_gamma_count++;

						}
						if (temp_energy_smeared2 > EM_thresh) {
							gamma2.isPion=1;

							gammas_true.push_back(gamma2);

							gamma2.lorentz.SetE(temp_energy_smeared2);
							gammas.push_back(gamma2);
							//pion_gamma_count++;

						}

					}


				}


			}


			double convlength[gammas.size()];
			double convpoint[gammas.size()][3];


			for (int j=0; j< gammas.size(); j++){
				//convlength[j] = photon_conversion_length(gammas.at(j).lorentz.E(), rangen);
				convlength[j] = photon_conversion_length(gammas.at(j).lorentz.E(), rangen);
				gammas.at(j).convlength=convlength[j];

				double temp_3vec[3] = {gammas.at(j).lorentz.Px(),gammas.at(j).lorentz.Py(),gammas.at(j).lorentz.Pz()};


				int k = get_endpoint(vertex_pos, convlength[j], temp_3vec, convpoint[j]);
				gammas.at(j).convpoint=convpoint[j];

				//  cout<<"conversion length: " << convlength[j]<< ", point x: " <<convpoint[j][0] << ", y: "<<convpoint[j][1]<<", z:"<< convpoint[j][2]<<  endl;


			}



			std::vector<myphoton> background_photons;
			std::vector<myphoton> background_photons_true;

			int photon_count = 0;
			int background_count = 0;

			for (int m=0; m<gammas.size(); m++){
				if (vis_vertex){
					if (gammas.at(m).convlength > convlength_thresh){
						if(DUNE.is_fiducial(gammas.at(m).convpoint)) photon_count++;
					}
					else{
						if(DUNE.is_fiducial(gammas.at(m).convpoint)){

							background_count++;
							background_photons.push_back(gammas.at(m));
							background_photons_true.push_back(gammas_true.at(m));
						}
					}

				}
				else {
					if(DUNE.is_fiducial(gammas.at(m).convpoint)){
						background_count++;
						background_photons.push_back(gammas.at(m));
						background_photons_true.push_back(gammas_true.at(m));
					}
				}
			}
			// photon flow end.


			double Pxtot = Pxhad;
			double Pytot = Pyhad;
			double Pztot = Pzhad;

			//TLorentzVector MissingEnergy;

			if (cc && (abs(neu)==12 || abs(neu)==14) ){
				Pxtot+=pxl;
				Pytot+=pyl;
				Pztot+=pzl;
			}



			/************************************************************
			 *			Intrinsic Nue Section
			 *
			 ***********************************************************/

			// intrinsic
			if (cc && abs(neu)==12) {

				double pT = TMath::Sqrt(Pxtot*Pxtot + Pytot*Pytot);
				double pL = Pztot;

				double El_smear = smear_energy_type(11, El, rangen);

				double reco_E = El_smear;

				if (vis_vertex){
					reco_E += Ehad;
				}

				if (El_smear > EM_thresh && DUNE.is_fiducial(vertex_pos)) {// && vis_vertex


					//MissingEnergy.SetXYZM(-Pxtot,-Pytot,-Pztot,0.);

					//hist_PT->Fill(MissingEnergy.Et(),cc_efficiency*vertex_weight);

					//hist_lepton_PT->Fill(sqrt(pxl*pxl+pyl*pyl),cc_efficiency*vertex_weight);

					//hist_PT_vs_lepton_PT->Fill(MissingEnergy.Et(),sqrt(pxl*pxl+pyl*pyl),cc_efficiency*vertex_weight);


					//hist_cc_El->Fill(El_smear,cc_efficiency*vertex_weight);
					//hist_cc_ehad->Fill(Ehad,cc_efficiency*vertex_weight);
					hist_cc_reco->Fill(reco_E,cc_efficiency*vertex_weight*prob_muflav(Ev,1,mode));//prob_muflav(Ev,1)

					// fill ttrees
					m_Ereco = reco_E;
					m_Etrue = Ev;
					m_l = 1300 ;	
					m_weight = cc_efficiency*vertex_weight;
					m_nutype = neu;
					//std::cout<<"NCC: "<<m_Ereco<<" "<<m_Etrue<<" "<<m_l<<" "<<m_weight<<" "<<m_nutype<<std::endl;
			
					if(isfullosc){
						list_o_trees.at(0)->Fill();
					}else{
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




			}// end of nue cc.


			/************************************************************
			 *			Muon Section
			 *
			 ***********************************************************/

			if (cc && abs(neu)==14) {
				// hist_charlep->Fill(El);
				double prob = 0.;


				double pT = TMath::Sqrt(Pxtot*Pxtot + Pytot*Pytot);
				double pL = Pztot;

				double Lmu = muon_track_length(El);
				double endpos[3] = {0,0,0};
				double temp_3vec[3] = {pxl,pyl,pzl};
				get_endpoint(vertex_pos, Lmu, temp_3vec, endpos);

				double El_smear = smear_energy_type(13, DUNE.is_fully_contained(vertex_pos, endpos), rangen);


				//ehad+ e_mu -> for vis vertex check..

				double Ehad_mu = El_smear+Ehad;

				double Ehad_mu_true = El+Ehad_true;

				bool vis_vertex_mu = false;

				if(Ehad_mu >= vertex_thresh)
				{
					vis_vertex_mu = true;
				}

				if (El_smear > MU_thresh && DUNE.is_fiducial(vertex_pos)) {//vis_vertex_mu &&

					//numu cc events with a track of length >/ 1m were assumed to be identifiable as numu-induced CC events and were rejected. Those with a track length below 1m were accepted as potential mis-identified events, if nay photons in the event were accepted under the same conditions as in the NC single photon events, above.

					double observable_L = 0;

					if( DUNE.is_fully_contained(vertex_pos, endpos) ) {
						observable_L=Lmu;
					}
					else {
						observable_L = DUNE.track_length_escape(vertex_pos,endpos);
					}



					if(observable_L<100.0 && background_photons.size()==1 && photon_count ==0)
					{
						double reco_E = background_photons.at(0).lorentz.E();

						if (vis_vertex_mu){
							reco_E += Ehad_mu;
						}
						//   hist_cc_El->Fill(El_smear,cc_efficiency*vertex_weight*egamma_misidrate);
						//   hist_cc_ehad->Fill(Ehad_mu,cc_efficiency*vertex_weight*egamma_misidrate);
						hist_cc_reco->Fill(reco_E,cc_efficiency*vertex_weight*egamma_misidrate);

					// fill ttrees
					m_Ereco = reco_E;
					m_Etrue = Ev;
					m_l = 1300 ;	
					m_weight = cc_efficiency*vertex_weight*egamma_misidrate;
					m_nutype = neu;
					//std::cout<<"NCC: "<<m_Ereco<<" "<<m_Etrue<<" "<<m_l<<" "<<m_weight<<" "<<m_nutype<<std::endl;
	
					if(isfullosc){
					}else{
						list_o_trees.at(2)->Fill();
					}	


						//   hist_cc_El_true->Fill(El,cc_efficiency*vertex_weight*egamma_misidrate);
						//   hist_cc_ehad_true->Fill(Ehad_mu_true,cc_efficiency*vertex_weight*egamma_misidrate);
						//   hist_cc_reco_true->Fill(background_photons_true.at(0).lorentz.E()+Ehad_mu_true,cc_efficiency*vertex_weight*egamma_misidrate);





					}
				}
			}// end of numu cc.


			/************************************************************
			 *			Tau Section
			 *
			 ***********************************************************/


			if (cc && abs(neu)==16) {

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

				daughter[0].Boost(pxl/El,pyl/El,pzl/El);
				daughter[1].Boost(pxl/El,pyl/El,pzl/El);
				daughter[2].Boost(pxl/El,pyl/El,pzl/El);

				TLorentzVector twoneutrinos = daughter[1]+daughter[2];


				double El_smear;

				El_smear = smear_energy_type(11, daughter[0].E(), rangen);

				double reco_E = El_smear;

				if (vis_vertex){
					reco_E += Ehad;
				}

				hist_cc_reco_effi->Fill(reco_E,prob_muflav(Ev,3,mode));//prob_muflav(Ev,3)*

				if (El_smear > EM_thresh && DUNE.is_fiducial(vertex_pos)) {//&& vis_vertex&&

					Pxtot+=daughter[0].Px();
					Pytot+=daughter[0].Py();
					Pztot+=daughter[0].Pz();

					double pT = TMath::Sqrt(Pxtot*Pxtot + Pytot*Pytot);
					double pL = Pztot;

					hist_cc_reco->Fill(reco_E,prob_muflav(Ev,3,mode)*cc_efficiency*vertex_weight*0.1783);

					// fill ttrees
					m_Ereco = reco_E;
					m_Etrue = Ev;
					m_l = 1300 ;	
					m_weight = cc_efficiency*vertex_weight*0.1783;
					m_nutype = neu;
					//std::cout<<"NCC: "<<m_Ereco<<" "<<m_Etrue<<" "<<m_l<<" "<<m_weight<<" "<<m_nutype<<std::endl;
					if(!isfullosc)	list_o_trees.at(3)->Fill();

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

	fnew->cd();
	fnew->Write();

	fntuple->cd();
	for(auto &t: list_o_trees){
		t->Write();
	}
	fntuple->Purge();
	fntuple->Close();
}


void genie_CC_slimmed(){


	TRandom3 *rangen = new TRandom3(0);

	TString filename = "gntp.0.nutau20k_gst";


	TFile *f = new TFile(filename+TString(".root"));
	TTree *gstree = (TTree*)f->Get("gst");

	TFile *fnew = new TFile(filename+TString("_out_CC_slimmed.root"),"recreate");

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
	int Npip;
	int Npim;
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
	double cc_efficiency=0.8;
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


	//TH1F *hist_PT_Weak = new TH1F("hist_PT_Weak", "PT_Weak", 80, 0., 2.);
	//hist_PT_Weak->Sumw2();
	//TH1F *hist_lepton_PT_Weak = new TH1F("hist_lepton_PT_Weak", "lepton_PT_Weak", 80, 0., 2.);
	//hist_lepton_PT_Weak->Sumw2();
	//TH2F *hist_PT_vs_lepton_PT_Weak = new TH2F("hist_PT_vs_lepton_PT_Weak", "PT_vs_lepton_PT_Weak", 80, 0., 2.,80.,0.,2.);
	//hist_PT_vs_lepton_PT_Weak->Sumw2();
	//TH1F *hist_twoneutrinos_PT_Weak = new TH1F("hist_twoneutrinos_PT_Weak", "twoneutrinos_PT_Weak", 80, 0., 2.);
	//hist_twoneutrinos_PT_Weak->Sumw2();

	TH1F *hist_PT_Polarized = new TH1F("hist_PT_Polarized", "PT_Polarized", 80, 0., 2.);
	hist_PT_Polarized->Sumw2();
	TH1F *hist_lepton_PT_Polarized = new TH1F("hist_lepton_PT_Polarized", "lepton_PT_Polarized", 80, 0., 2.);
	hist_lepton_PT_Polarized->Sumw2();
	TH2F *hist_PT_vs_lepton_PT_Polarized = new TH2F("hist_PT_vs_lepton_PT_Polarized", "PT_vs_lepton_PT_Polarized", 80, 0., 2.,80.,0.,2.);
	hist_PT_vs_lepton_PT_Polarized->Sumw2();
	TH1F *hist_twoneutrinos_PT_Polarized = new TH1F("hist_twoneutrinos_PT_Polarized", "twoneutrinos_PT_Polarized", 80, 0., 2.);
	hist_twoneutrinos_PT_Polarized->Sumw2();


	TH2F *hist_dalitz_DecayIt = new TH2F("hist_dalitz_DecayIt", "dalitz_DecayIt",200,0.,5.,200,0.,5.);
	hist_dalitz_DecayIt->Sumw2();
	TH2F *hist_dalitz_DecayItPolarized = new TH2F("hist_dalitz_DecayItPolarized", "dalitz_DecayItPolarized",200,0.,5.,200,0.,5.);
	hist_dalitz_DecayItPolarized->Sumw2();

	TH1F *hist_cc_ehad = new TH1F("hist_cc_ehad", "cc_ehad", 80, 0., 20.);
	hist_cc_ehad->Sumw2();
	TH1F *hist_cc_reco = new TH1F("hist_cc_reco", "cc_reco", 80, 0., 20.);
	hist_cc_reco->Sumw2();
	TH1F *hist_cc_El = new TH1F("hist_cc_El", "cc_El", 80, 0., 20.);
	hist_cc_El->Sumw2();

	TH1F *hist_cc_ehad_true = new TH1F("hist_cc_ehad_true", "cc_ehad_true", 80, 0., 20.);
	hist_cc_ehad_true->Sumw2();
	TH1F *hist_cc_reco_true = new TH1F("hist_cc_reco_true", "cc_reco_true", 80, 0., 20.);
	hist_cc_reco_true->Sumw2();
	TH1F *hist_cc_El_true = new TH1F("hist_cc_El_true", "cc_El_true", 40, 0., 20.);
	hist_cc_El_true->Sumw2();

	TH1F *hist_cc_ehad_pt05cut = new TH1F("hist_cc_ehad_pt05cut", "cc_ehad_pt05cut", 80, 0., 20.);
	hist_cc_ehad_pt05cut->Sumw2();
	TH1F *hist_cc_reco_pt05cut = new TH1F("hist_cc_reco_pt05cut", "cc_reco_pt05cut", 80, 0., 20.);
	hist_cc_reco_pt05cut->Sumw2();
	TH1F *hist_cc_El_pt05cut = new TH1F("hist_cc_El_pt05cut", "cc_El_pt05cut", 80, 0., 20.);
	hist_cc_El_pt05cut->Sumw2();

	TH1F *hist_cc_ehad_true_pt05cut = new TH1F("hist_cc_ehad_true_pt05cut", "cc_ehad_true_pt05cut", 80, 0., 20.);
	hist_cc_ehad_true_pt05cut->Sumw2();
	TH1F *hist_cc_reco_true_pt05cut = new TH1F("hist_cc_reco_true_pt05cut", "cc_reco_true_pt05cut", 80, 0., 20.);
	hist_cc_reco_true_pt05cut->Sumw2();
	TH1F *hist_cc_El_true_pt05cut = new TH1F("hist_cc_El_true_pt05cut", "cc_El_true_pt05cut", 40, 0., 20.);
	hist_cc_El_true_pt05cut->Sumw2();

	TH1F *hist_cc_ElCosTheta_true = new TH1F("hist_cc_ElCosTheta_true", "cc_ElCosTheta_true", 40, 0.6, 1.);
	hist_cc_ElCosTheta_true->Sumw2();

	TH1F *hist_cc_ElPhi_true = new TH1F("hist_cc_ElPhi_true", "cc_ElPhi_true", 80, -3.14, 3.14);
	hist_cc_ElPhi_true->Sumw2();
	TH2F *hist_cc_El_vs_CosTheta_true = new TH2F("hist_cc_El_vs_CosTheta_true", "cc_El_vs_CosTheta_true", 10, 0., 20., 10, 0.8, 1.0);
	hist_cc_El_vs_CosTheta_true->Sumw2();



	TH1F *hist_cc_ehad_WithSpin = new TH1F("hist_cc_ehad_WithSpin", "cc_ehad_WithSpin", 80, 0., 20.);
	hist_cc_ehad_WithSpin->Sumw2();
	TH1F *hist_cc_reco_WithSpin = new TH1F("hist_cc_reco_WithSpin", "cc_reco_WithSpin", 80, 0., 20.);
	hist_cc_reco_WithSpin->Sumw2();
	TH1F *hist_cc_El_WithSpin = new TH1F("hist_cc_El_WithSpin", "cc_El_WithSpin", 80, 0., 20.);
	hist_cc_El_WithSpin->Sumw2();
	TH1F *hist_cc_ehad_true_WithSpin = new TH1F("hist_cc_ehad_true_WithSpin", "cc_ehad_true_WithSpin", 80, 0., 20.);
	hist_cc_ehad_true_WithSpin->Sumw2();
	TH1F *hist_cc_reco_true_WithSpin = new TH1F("hist_cc_reco_true_WithSpin", "cc_reco_true_WithSpin", 80, 0., 20.);
	hist_cc_reco_true_WithSpin->Sumw2();
	TH1F *hist_cc_El_true_WithSpin = new TH1F("hist_cc_El_true_WithSpin", "cc_El_true_WithSpin", 40, 0., 20.);
	hist_cc_El_true_WithSpin->Sumw2();
	TH1F *hist_cc_ElCosTheta_true_WithSpin = new TH1F("hist_cc_ElCosTheta_true_WithSpin", "cc_ElCosTheta_true_WithSpin", 40, 0.6, 1.);
	hist_cc_ElCosTheta_true_WithSpin->Sumw2();
	TH1F *hist_cc_ElCosTheta_above2GeV_true_WithSpin = new TH1F("hist_cc_ElCosTheta_above2GeV_true_WithSpin", "cc_ElCosTheta_above2GeV_true_WithSpin", 40, 0.6, 1.);
	hist_cc_ElCosTheta_above2GeV_true_WithSpin->Sumw2();
	TH1F *hist_cc_ElPhi_true_WithSpin = new TH1F("hist_cc_ElPhi_true_WithSpin", "cc_ElPhi_true_WithSpin", 80, -3.14, 3.14);
	hist_cc_ElPhi_true_WithSpin->Sumw2();
	TH2F *hist_cc_El_vs_CosTheta_true_WithSpin = new TH2F("hist_cc_El_vs_CosTheta_true_WithSpin", "cc_El_vs_CosTheta_true_WithSpin", 10, 0., 20., 10, 0.8, 1.0);
	hist_cc_El_vs_CosTheta_true_WithSpin->Sumw2();


	TH1F *hist_cc_Ev_true = new TH1F("hist_cc_Ev_true", "cc_Ev_true", 40, 0., 40.);
	hist_cc_Ev_true->Sumw2();
	TH1F *hist_cc_tau_El_true = new TH1F("hist_cc_tau_El_true", "cc_tau_El_true", 80, 0., 40.);
	hist_cc_tau_El_true->Sumw2();
	TH1F *hist_cc_tau_ElCosTheta_true = new TH1F("hist_cc_tau_ElCosTheta_true", "cc_tau_ElCosTheta_true", 80, 0.8, 1.);
	hist_cc_tau_ElCosTheta_true->Sumw2();
	TH1F *hist_cc_tau_ElTheta_true = new TH1F("hist_cc_tau_ElTheta_true", "cc_tau_ElTheta_true", 80, 0., 1.);
	hist_cc_tau_ElTheta_true->Sumw2();
	TH1F *hist_cc_tau_ElPhi_true = new TH1F("hist_cc_tau_ElPhi_true", "cc_tau_ElPhi_true", 80, -3.14, 3.14);
	hist_cc_tau_ElPhi_true->Sumw2();


	TH1F *hist_tau_Ev_true = new TH1F("hist_tau_Ev_true", "tau_Ev_true", 40, 0., 40.);
	hist_tau_Ev_true->Sumw2();

	TH1F *hist_tau_Ev_osc = new TH1F("hist_tau_Ev_osc", "tau_Ev_osc", 40, 0., 40.);
	hist_tau_Ev_osc->Sumw2();


	TH1F *hist_nc_ehad = new TH1F("hist_nc_ehad", "nc_ehad", 80, 0., 20.);
	hist_nc_ehad->Sumw2();
	TH1F *hist_nc_reco = new TH1F("hist_nc_reco", "nc_reco", 80, 0., 20.);
	hist_nc_reco->Sumw2();
	TH1F *hist_nc_gamma = new TH1F("hist_nc_gamma", "nc_gamma", 80, 0., 20.);
	hist_nc_gamma->Sumw2();


	TH1F *hist_nc_ehad_true = new TH1F("hist_nc_ehad_true", "nc_ehad_true", 80, 0., 20.);
	hist_nc_ehad_true->Sumw2();
	TH1F *hist_nc_reco_true = new TH1F("hist_nc_reco_true", "nc_reco_true", 80, 0., 20.);
	hist_nc_reco_true->Sumw2();
	TH1F *hist_nc_gamma_true = new TH1F("hist_nc_gamma_true", "nc_gamma_true", 80, 0., 20.);
	hist_nc_gamma_true->Sumw2();

	TH2F *hist_lPT_vs_lPL = new TH2F("hist_lPT_vs_lPL", "lPT_vs_lPL", 20, 0., 4., 20, 0., 40.);
	TH2F *hist_PT_vs_PL = new TH2F("hist_PT_vs_PL", "PT_vs_PL", 20, 0., 40., 20, 0., 40.);
	TH1F *hist_PT = new TH1F("hist_PT", "PT", 40, 0., 2.);
	TH1F *hist_PL = new TH1F("hist_PL", "PL", 20, 0., 40.);

	TH1F *hist_PLsqOverPTsq = new TH1F("hist_PLsqOverPTsq", "PLsqOverPTsq", 20, 0., 400.);
	TH1F *hist_PTsqOverPLsq = new TH1F("hist_PTsqOverPLsq", "PTsqOverPLsq", 160, 0., 0.2);

	TH1F *hist_PTsqOverPL = new TH1F("hist_PTsqOverPL", "PTsqOverPL", 160, 0., 0.2);

	TH1F *hist_PTsq = new TH1F("hist_PTsq", "PTsq", 40, 0., 1.);
	TH1F *hist_PTcube = new TH1F("hist_PTcube", "PTcube", 40, 0., 0.5);
	TH1F *hist_PTquartic = new TH1F("hist_PTquartic", "PTquartic", 40, 0., 0.125);




	//weak decay
	TH2F *hist_lPT_vs_lPL_weak = new TH2F("hist_lPT_vs_lPL_weak", "lPT_vs_lPL_weak", 20, 0., 4., 20, 0., 40.);
	TH2F *hist_PT_vs_PL_weak = new TH2F("hist_PT_vs_PL_weak", "PT_vs_PL_weak", 20, 0., 40., 20, 0., 40.);
	TH1F *hist_PT_weak = new TH1F("hist_PT_weak", "PT_weak", 40, 0., 2.);
	TH1F *hist_PL_weak = new TH1F("hist_PL_weak", "PL_weak", 20, 0., 40.);

	TH1F *hist_PLsqOverPTsq_weak = new TH1F("hist_PLsqOverPTsq_weak", "PLsqOverPTsq_weak", 20, 0., 400.);
	TH1F *hist_PTsqOverPLsq_weak = new TH1F("hist_PTsqOverPLsq_weak", "PTsqOverPLsq_weak", 80, 0., 0.2);

	TH1F *hist_PTsqOverPL_weak = new TH1F("hist_PTsqOverPL_weak", "PTsqOverPL_weak", 40, 0., 0.2);

	TH1F *hist_PTsq_weak = new TH1F("hist_PTsq_weak", "PTsq_weak", 40, 0., 1.);
	TH1F *hist_PTcube_weak = new TH1F("hist_PTcube_weak", "PTcube_weak", 40, 0., 0.5);
	TH1F *hist_PTquartic_weak = new TH1F("hist_PTquartic_weak", "PTquartic_weak", 40, 0., 0.125);


	SBN_detector DUNE(4,false);

	//SBN_detector DUNE_NC(5,false);

	for(int i=0; i < gstree->GetEntries(); i++){
		gstree->GetEntry(i);

		if (nc) {

			continue;

		}

		bool mode = true;
		if (neu<0){
			mode = false;
		}

		//tau test plot
		hist_tau_Ev_true->Fill(Ev);
		hist_tau_Ev_osc->Fill(Ev,prob_muflav(Ev,3,mode));


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
					//cout << " proton mass : " << sqrt(energyf[j]*energyf[j] - pxf[j]*pxf[j]-pyf[j]*pyf[j]-pzf[j]*pzf[j])<< endl;
					kin_true=energyf[j]-MPROTON;
					//kin_smeared=smear_energy(kin_true, psmear, rangen);
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
					//cout << " pion mass : " << sqrt(energyf[j]*energyf[j] - pxf[j]*pxf[j]-pyf[j]*pyf[j]-pzf[j]*pzf[j])<< endl;
					kin_true=energyf[j]-MPION;
					//  kin_smeared=smear_energy(kin_true, pismear, rangen);
					kin_smeared=smear_energy_type(pdgf[j], kin_true, 1, rangen);
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
		//hist_vis_vertex->Fill(vis_vertex);


		//Nrevertex. for stats.
		int Nrevertex = 25;
		for(int k=0;k<Nrevertex; k++){

			double vertex_weight = 1.0/((double) Nrevertex);


			double vertex_pos[3];

			//double vertex_pos_NC[3];

			DUNE.random_pos(rangen, vertex_pos);
			//DUNE_NC.random_pos(rangen, vertex_pos_NC);



			//photon flow
			std::vector<myphoton> gammas;
			std::vector<myphoton> gammas_true;

			//       std::vector<myphoton> gammas_smeared;

			// TLorentzVector final_gamma[10];
			//radiation gammas
			if (nf!=0){
				for(int j=0; j<nf;j++ ){

					if(pdgf[j]==22){ //energyf[j]>EM_thresh

						myphoton temp_vec;
						temp_vec.lorentz.SetPxPyPzE(pxf[j],pyf[j],pzf[j],energyf[j]);


						double temp_energy_smeared = smear_energy(energyf[j], EMsmear,rangen);

						if (temp_energy_smeared > EM_thresh) {

							gammas_true.push_back(temp_vec);

							temp_vec.lorentz.SetE(temp_energy_smeared);
							gammas.push_back(temp_vec);

							gammas.back().isPion=0;

						}
					}
				}
			}

			//pi0 gammas
			if (nf!=0){

				for(int k=0;k<nf;k++){

					if(pdgf[k]==111){

						double boost[3] = {pxf[k]/energyf[k],pyf[k]/energyf[k],pzf[k]/energyf[k]};//boost vector for pi0

						double x=0;
						double y=0;
						double z=0;

						rangen->Sphere(x,y,z,1.);//generate random direction on a unit sphere.

						myphoton gamma1, gamma2;


						gamma1.lorentz.SetPxPyPzE(x*pion_mass/2,y*pion_mass/2,z*pion_mass/2,pion_mass/2);// assigning gamma1 along the random direction
						gamma2.lorentz.SetPxPyPzE(-x*pion_mass/2,-y*pion_mass/2,-z*pion_mass/2,pion_mass/2);// assigning gamma2 to the opposite directon of gamma1, each of these has a half of pion rest energy.


						gamma1.lorentz.Boost(pxf[k]/energyf[k],pyf[k]/energyf[k],pzf[k]/energyf[k]);//boost gamma1 and gamma2 to original pi0 frame.
						gamma2.lorentz.Boost(pxf[k]/energyf[k],pyf[k]/energyf[k],pzf[k]/energyf[k]);
						// check invariant amss of gamma1 gamm2 pair ==mpion

						double temp_energy_smeared1= smear_energy(gamma1.lorentz.E(), EMsmear,rangen);
						double temp_energy_smeared2= smear_energy(gamma2.lorentz.E(), EMsmear,rangen);

						if (temp_energy_smeared1 > EM_thresh) {
							gamma1.isPion=1;
							gammas_true.push_back(gamma1);

							gamma1.lorentz.SetE(temp_energy_smeared1);
							gammas.push_back(gamma1);
							//pion_gamma_count++;

						}
						if (temp_energy_smeared2 > EM_thresh) {
							gamma2.isPion=1;

							gammas_true.push_back(gamma2);

							gamma2.lorentz.SetE(temp_energy_smeared2);
							gammas.push_back(gamma2);
							//pion_gamma_count++;

						}

					}


				}


			}

			if (gammas.size()>0){
				//cout << "# of gammas over EM_thresh : " << gammas.size() << ", first gamma energy : "<< gammas.at(0).lorentz.E()<< endl;
			}// we have everyting in gammas.

			double convlength[gammas.size()];
			double convpoint[gammas.size()][3];

			// cout << "vertex point x: " <<vertex_pos[0] << ", y: "<<vertex_pos[1]<<", z:"<< vertex_pos[2]<<  endl;


			for (int j=0; j< gammas.size(); j++){
				//convlength[j] = photon_conversion_length(gammas.at(j).lorentz.E(), rangen);
				convlength[j] = photon_conversion_length(gammas.at(j).lorentz.E(), rangen);
				gammas.at(j).convlength=convlength[j];

				double temp_3vec[3] = {gammas.at(j).lorentz.Px(),gammas.at(j).lorentz.Py(),gammas.at(j).lorentz.Pz()};


				int k = get_endpoint(vertex_pos, convlength[j], temp_3vec, convpoint[j]);
				gammas.at(j).convpoint=convpoint[j];

				//  cout<<"conversion length: " << convlength[j]<< ", point x: " <<convpoint[j][0] << ", y: "<<convpoint[j][1]<<", z:"<< convpoint[j][2]<<  endl;


			}



			std::vector<myphoton> background_photons;
			std::vector<myphoton> background_photons_true;

			int photon_count = 0;
			int background_count = 0;

			for (int m=0; m<gammas.size(); m++){
				if (vis_vertex){
					if (gammas.at(m).convlength > convlength_thresh){
						if(DUNE.is_fiducial(gammas.at(m).convpoint)) photon_count++;
					}
					else{
						if(DUNE.is_fiducial(gammas.at(m).convpoint)){

							background_count++;
							background_photons.push_back(gammas.at(m));
							background_photons_true.push_back(gammas_true.at(m));
						}
					}

				}
				else {
					if(DUNE.is_fiducial(gammas.at(m).convpoint)){
						background_count++;
						background_photons.push_back(gammas.at(m));
						background_photons_true.push_back(gammas_true.at(m));
					}
				}
			}
			// photon flow end.





			double Pxtot = Pxhad;
			double Pytot = Pyhad;
			double Pztot = Pzhad;

			TLorentzVector MissingEnergy;

			if (cc && (abs(neu)==12 || abs(neu)==14) ){
				Pxtot+=pxl;
				Pytot+=pyl;
				Pztot+=pzl;
			}




			// intrinsic
			if (cc && abs(neu)==12) {
				// hist_charlep->Fill(El);

				//cout << " electron mass : " << sqrt(El*El-pxl*pxl-pyl*pyl-pzl*pzl)<< endl;

				//double El_smear = smear_energy(El, EMsmear, rangen);
				double El_smear = smear_energy_type(11, El, rangen);
				if (El_smear > EM_thresh && vis_vertex && DUNE.is_fiducial(vertex_pos)) {


					MissingEnergy.SetXYZM(-Pxtot,-Pytot,-Pztot,0.);

					//hist_PT->Fill(MissingEnergy.Et(),cc_efficiency*vertex_weight);

					//hist_lepton_PT->Fill(sqrt(pxl*pxl+pyl*pyl),cc_efficiency*vertex_weight);

					//hist_PT_vs_lepton_PT->Fill(MissingEnergy.Et(),sqrt(pxl*pxl+pyl*pyl),cc_efficiency*vertex_weight);


					//hist_cc_El->Fill(El_smear,cc_efficiency*vertex_weight);
					//hist_cc_ehad->Fill(Ehad,cc_efficiency*vertex_weight);
					//hist_cc_reco->Fill(El_smear+Ehad,cc_efficiency*vertex_weight);

					//hist_cc_El_true->Fill(El,cc_efficiency*vertex_weight);
					//hist_cc_ehad_true->Fill(Ehad_true,cc_efficiency*vertex_weight);
					//hist_cc_reco_true->Fill(El+Ehad_true,cc_efficiency*vertex_weight);

					TLorentzVector truelvec;

					truelvec.SetPxPyPzE(pxl,pyl,pzl,El);

					//hist_cc_ElCosTheta_true->Fill(truelvec.CosTheta(),cc_efficiency*vertex_weight);


					if (MissingEnergy.Et()<0.5){

						// hist_cc_El_pt05cut->Fill(El_smear,cc_efficiency*vertex_weight);
						// hist_cc_ehad_pt05cut->Fill(Ehad,cc_efficiency*vertex_weight);
						// hist_cc_reco_pt05cut->Fill(El_smear+Ehad,cc_efficiency*vertex_weight);

						//  hist_cc_El_true_pt05cut->Fill(El,cc_efficiency*vertex_weight);
						//  hist_cc_ehad_true_pt05cut->Fill(Ehad_true,cc_efficiency*vertex_weight);
						//  hist_cc_reco_true_pt05cut->Fill(El+Ehad_true,cc_efficiency*vertex_weight);


					}

					hist_PT_vs_PL->Fill(MissingEnergy.Pt(),Pztot);
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
					hist_PTquartic_weak->Fill(MissingEnergy.Pt()*MissingEnergy.Pt()*MissingEnergy.Pt()*MissingEnergy.Pt());




				}




			}// end of nue cc.

			if (cc && abs(neu)==14) {
				// hist_charlep->Fill(El);


				double Lmu = muon_track_length(El);
				double endpos[3] = {0,0,0};
				double temp_3vec[3] = {pxl,pyl,pzl};
				get_endpoint(vertex_pos, Lmu, temp_3vec, endpos);

				double El_smear = smear_energy_type(13, DUNE.is_fully_contained(vertex_pos, endpos), rangen);


				//ehad+ e_mu -> for vis vertex check..

				double Ehad_mu = El_smear+Ehad;

				double Ehad_mu_true = El+Ehad_true;

				bool vis_vertex_mu = false;

				if(Ehad_mu >= vertex_thresh)
				{
					vis_vertex_mu = true;
				}

				if (El_smear > MU_thresh && vis_vertex_mu && DUNE.is_fiducial(vertex_pos)) {

					//numu cc events with a track of length >/ 1m were assumed to be identifiable as numu-induced CC events and were rejected. Those with a track length below 1m were accepted as potential mis-identified events, if nay photons in the event were accepted under the same conditions as in the NC single photon events, above.

					double observable_L = 0;

					if( DUNE.is_fully_contained(vertex_pos, endpos) ) {
						observable_L=Lmu;
					}
					else {
						observable_L = DUNE.track_length_escape(vertex_pos,endpos);
					}

					if(observable_L<100.0 && background_photons.size()==1 && photon_count ==0)
					{
						//   hist_cc_El->Fill(El_smear,cc_efficiency*vertex_weight*egamma_misidrate);
						//   hist_cc_ehad->Fill(Ehad_mu,cc_efficiency*vertex_weight*egamma_misidrate);
						//   hist_cc_reco->Fill(background_photons.at(0).lorentz.E()+Ehad_mu,cc_efficiency*vertex_weight*egamma_misidrate);

						//   hist_cc_El_true->Fill(El,cc_efficiency*vertex_weight*egamma_misidrate);
						//   hist_cc_ehad_true->Fill(Ehad_mu_true,cc_efficiency*vertex_weight*egamma_misidrate);
						//   hist_cc_reco_true->Fill(background_photons_true.at(0).lorentz.E()+Ehad_mu_true,cc_efficiency*vertex_weight*egamma_misidrate);

					}
				}
			}// end of numu cc.



			if (cc && abs(neu)==16) {
				// hist_charlep->Fill(El);





				//decay tau -> nu tau , nu el, el.
				TLorentzVector daughter[3];
				daughter[0].SetPxPyPzE(0.0,0.0,0.0,0.0);
				daughter[1].SetPxPyPzE(0.0,0.0,0.0,0.0);
				daughter[2].SetPxPyPzE(0.0,0.0,0.0,0.0);

				TLorentzVector daughterWithSpin[3];
				daughterWithSpin[0].SetPxPyPzE(0.0,0.0,0.0,0.0);
				daughterWithSpin[1].SetPxPyPzE(0.0,0.0,0.0,0.0);
				daughterWithSpin[2].SetPxPyPzE(0.0,0.0,0.0,0.0);

				double thetaP = TMath::Pi()*El/Ev;

				TVector3 parent_polarization (0.,0.,-1.);

				if(neu == 16){
					parent_polarization.SetXYZ(-TMath::Sin(thetaP),0.0,TMath::Cos(thetaP));
				}
				else if(neu == -16){
					parent_polarization.SetXYZ(0.,0.,-1.);
				}

				double tau_mass = 1.76;

				TLorentzVector parentTau;
				parentTau.SetPxPyPzE(0.0,0.0,0.0,tau_mass);

				double p0[4] = {0,0,0,0};
				double p1[4] = {0,0,0,0};
				double p2[4] = {0,0,0,0};

				//threebodydecay(rangen, tau_mass, 0.00051, 0., 0., p0, p1, p2);

				DecayIt(parentTau, rangen, daughter);

				DecayItWithSpin(parentTau, rangen, daughterWithSpin, parent_polarization);

				TLorentzVector daughter12, daughter23;
				TLorentzVector daughterWithSpin12, daughterWithSpin23;

				daughter12 = daughter[2]+daughter[1];
				daughter23 = daughter[1]+daughter[0];

				daughterWithSpin12 = daughterWithSpin[2]+daughterWithSpin[1];
				daughterWithSpin23 = daughterWithSpin[1]+daughterWithSpin[0];



				hist_dalitz_DecayIt->Fill(daughter12.M()*daughter12.M(), daughter23.M()*daughter23.M() );
				hist_dalitz_DecayItPolarized->Fill(daughterWithSpin12.M()*daughterWithSpin12.M(), daughterWithSpin23.M()*daughterWithSpin23.M() );

				//tau vector in lab frame
				TLorentzVector labtau;
				labtau.SetPxPyPzE(pxl,pyl,pzl,El);

				hist_lPT_vs_lPL->Fill(labtau.Pt(),labtau.Pz());


				//  hist_cc_Ev_true->Fill(Ev,prob_muflav(Ev,3));
				//  hist_cc_tau_El_true->Fill(El,prob_muflav(Ev,3));
				//  hist_cc_tau_ElCosTheta_true->Fill(labtau.CosTheta(),prob_muflav(Ev,3));
				//  hist_cc_tau_ElTheta_true->Fill(labtau.Theta(),prob_muflav(Ev,3));
				//  hist_cc_tau_ElPhi_true->Fill(labtau.Phi(),prob_muflav(Ev,3));


				daughter[0].Boost(pxl/El,pyl/El,pzl/El);
				daughter[1].Boost(pxl/El,pyl/El,pzl/El);
				daughter[2].Boost(pxl/El,pyl/El,pzl/El);

				daughterWithSpin[0].Boost(pxl/El,pyl/El,pzl/El);
				daughterWithSpin[1].Boost(pxl/El,pyl/El,pzl/El);
				daughterWithSpin[2].Boost(pxl/El,pyl/El,pzl/El);


				TLorentzVector twoneutrinos = daughter[1]+daughter[2];
				TLorentzVector twoneutrinosPolarized = daughterWithSpin[1]+daughterWithSpin[2];


				//Fill DecayIt
				double El_smear;

				El_smear = smear_energy_type(11, daughter[0].E(), rangen);
				if (El_smear > EM_thresh && vis_vertex && DUNE.is_fiducial(vertex_pos)) {

					Pxtot+=daughter[0].Px();
					Pytot+=daughter[0].Py();
					Pztot+=daughter[0].Pz();

					MissingEnergy.SetXYZM(-Pxtot,-Pytot,-Pztot,0.);

					hist_PT_vs_PL_weak->Fill(MissingEnergy.Pt(),Pztot,prob_muflav(Ev,3,mode));
					hist_PT_weak->Fill(MissingEnergy.Pt(),prob_muflav(Ev,3,mode));
					hist_PL_weak->Fill(Pztot,prob_muflav(Ev,3,mode));
					hist_PLsqOverPTsq_weak->Fill(Pztot*Pztot/MissingEnergy.Pt()/MissingEnergy.Pt(),prob_muflav(Ev,3,mode));
					hist_PTsqOverPLsq_weak->Fill(MissingEnergy.Pt()*MissingEnergy.Pt()/Pztot/Pztot,prob_muflav(Ev,3,mode));
					hist_PTsqOverPL_weak->Fill(MissingEnergy.Pt()*MissingEnergy.Pt()/Pztot,prob_muflav(Ev,3,mode));
					hist_PTsq_weak->Fill(MissingEnergy.Pt()*MissingEnergy.Pt(),prob_muflav(Ev,3,mode));
					hist_PTcube_weak->Fill(MissingEnergy.Pt()*MissingEnergy.Pt()*MissingEnergy.Pt(),prob_muflav(Ev,3,mode));
					hist_PTquartic_weak->Fill(MissingEnergy.Pt()*MissingEnergy.Pt()*MissingEnergy.Pt()*MissingEnergy.Pt(),prob_muflav(Ev,3,mode));


					if (MissingEnergy.Et()<0.5){


						//_pt05cut

						//          hist_cc_El_pt05cut->Fill(El_smear,cc_efficiency*vertex_weight*0.1783*prob_muflav(Ev,3));
						//          hist_cc_ehad_pt05cut->Fill(Ehad,cc_efficiency*vertex_weight*0.1783*prob_muflav(Ev,3));
						//          hist_cc_reco_pt05cut->Fill(El_smear+Ehad,cc_efficiency*vertex_weight*0.1783*prob_muflav(Ev,3));

						//          hist_cc_El_true_pt05cut->Fill(daughter[0].E(),cc_efficiency*vertex_weight*0.1783*prob_muflav(Ev,3));
						//          hist_cc_ehad_true_pt05cut->Fill(Ehad_true,cc_efficiency*vertex_weight*0.1783*prob_muflav(Ev,3));
						//          hist_cc_reco_true_pt05cut->Fill(daughter[0].E()+Ehad_true,cc_efficiency*vertex_weight*0.1783*prob_muflav(Ev,3));


					}

					Pxtot-=daughter[0].Px();
					Pytot-=daughter[0].Py();
					Pztot-=daughter[0].Pz();

				}

				//Fill DecayItWithSpin
				El_smear = smear_energy_type(11, daughterWithSpin[0].E(), rangen);

				if (El_smear > EM_thresh && vis_vertex && DUNE.is_fiducial(vertex_pos)) {

					Pxtot+=daughterWithSpin[0].Px();
					Pytot+=daughterWithSpin[0].Py();
					Pztot+=daughterWithSpin[0].Pz();

					MissingEnergy.SetXYZM(-Pxtot,-Pytot,-Pztot,0.);
					//     hist_lepton_PT_WithSpin->Fill(daughterWithSpin[0].Pt(),prob_muflav(Ev,3));
					//     hist_PT_vs_lepton_PT_WithSpin->Fill(MissingEnergy.Et(),daughterWithSpin[0].Pt(),prob_muflav(Ev,3));


					//      hist_cc_El_WithSpin->Fill(El_smear,cc_efficiency*vertex_weight*0.1783*prob_muflav(Ev,3));
					//      hist_cc_ehad_WithSpin->Fill(Ehad,cc_efficiency*vertex_weight*0.1783*prob_muflav(Ev,3));
					//      hist_cc_reco_WithSpin->Fill(El_smear+Ehad,cc_efficiency*vertex_weight*0.1783*prob_muflav(Ev,3));

					//      hist_cc_El_true_WithSpin->Fill(daughterWithSpin[0].E(),cc_efficiency*vertex_weight*0.1783*prob_muflav(Ev,3));
					//      hist_cc_ehad_true_WithSpin->Fill(Ehad_true,cc_efficiency*vertex_weight*0.1783*prob_muflav(Ev,3));
					//      hist_cc_reco_true_WithSpin->Fill(daughterWithSpin[0].E()+Ehad_true,cc_efficiency*vertex_weight*0.1783*prob_muflav(Ev,3));

					//      hist_cc_ElCosTheta_true_WithSpin->Fill(daughterWithSpin[0].CosTheta(), cc_efficiency*vertex_weight*0.1783*prob_muflav(Ev,3));

					//       hist_cc_ElPhi_true_WithSpin->Fill(daughterWithSpin[0].Phi(), cc_efficiency*vertex_weight*0.1783*prob_muflav(Ev,3));
					//       hist_cc_El_vs_CosTheta_true_WithSpin->Fill(daughterWithSpin[0].E(), daughterWithSpin[0].CosTheta(), cc_efficiency*vertex_weight*0.1783*prob_muflav(Ev,3));

					//       hist_PT_WithSpin->Fill(MissingEnergy.Et(),cc_efficiency*vertex_weight*0.1783*prob_muflav(Ev,3));
					//       hist_twoneutrinos_PT_WithSpin->Fill(twoneutrinosWithSpin.Et(),cc_efficiency*vertex_weight*0.1783*prob_muflav(Ev,3));

					hist_PT_vs_PL->Fill(MissingEnergy.Pt(),Pztot,prob_muflav(Ev,3,mode));
					hist_PT->Fill(MissingEnergy.Pt(),prob_muflav(Ev,3,mode));
					hist_PL->Fill(Pztot,prob_muflav(Ev,3,mode));
					hist_PLsqOverPTsq->Fill(Pztot*Pztot/MissingEnergy.Pt()/MissingEnergy.Pt(),prob_muflav(Ev,3,mode));
					hist_PTsqOverPLsq->Fill(MissingEnergy.Pt()*MissingEnergy.Pt()/Pztot/Pztot,prob_muflav(Ev,3,mode));
					hist_PTsqOverPL->Fill(MissingEnergy.Pt()*MissingEnergy.Pt()/Pztot,prob_muflav(Ev,3,mode));
					hist_PTsq->Fill(MissingEnergy.Pt()*MissingEnergy.Pt(),prob_muflav(Ev,3,mode));
					hist_PTcube->Fill(MissingEnergy.Pt()*MissingEnergy.Pt()*MissingEnergy.Pt(),prob_muflav(Ev,3,mode));
					hist_PTquartic->Fill(MissingEnergy.Pt()*MissingEnergy.Pt()*MissingEnergy.Pt()*MissingEnergy.Pt(),prob_muflav(Ev,3,mode));




					Pxtot-=daughterWithSpin[0].Px();
					Pytot-=daughterWithSpin[0].Py();
					Pztot-=daughterWithSpin[0].Pz();

				}


			}// end of nutau cc.
		}
	}

	fnew->Write();

}


void genie_CC(){

	/*"The signal for nue apperance is an excess of charged-current(CC) nue and nuebar interactions over the expected background in the far detector. The background to nue appearance is composed of : (1) CC interactions of nue and nuebar intrinsic to the beam; (2) misidentified numu and numubar CC events; (3) neutral current (NC) backgrounds and (4) nutau and nutaubar CC events in which the tau's decay leptonically into electrons/positrons. NC and nutau backgrounds are due to interactions of higher-energy neutrinos but they contribute to backgrounds mainly at low energy, which is important for the sensitivity to CP violation."*/

	TRandom3 *rangen = new TRandom3(0);


	TFile *f = new TFile("gntp.0.gst_nue.root");
	TTree *gstree = (TTree*)f->Get("gst");

	TFile *fnew = new TFile("test_nue_CC_prob.root","recreate");

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
	int Npip;
	int Npim;
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
	double cc_efficiency=0.8;
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

	TH1F *hist_sizeofgammas = new TH1F("hist_sizeofgammas", "sizeofgammas", 10, 0., 10.);
	hist_sizeofgammas->Sumw2();

	TH1F *hist_PT = new TH1F("hist_PT", "PT", 80, 0., 2.);
	hist_PT->Sumw2();
	TH1F *hist_lepton_PT = new TH1F("hist_lepton_PT", "lepton_PT", 80, 0., 2.);
	hist_lepton_PT->Sumw2();
	TH2F *hist_PT_vs_lepton_PT = new TH2F("hist_PT_vs_lepton_PT", "PT_vs_lepton_EP", 80, 0., 2.,80.,0.,2.);
	hist_PT_vs_lepton_PT->Sumw2();
	TH1F *hist_twoneutrinos_PT = new TH1F("hist_twoneutrinos_PT", "twoneutrinos_PT", 80, 0., 2.);
	hist_twoneutrinos_PT->Sumw2();

	TH1F *hist_PT_WithSpin = new TH1F("hist_PT_WithSpin", "PT_WithSpin", 80, 0., 2.);
	hist_PT_WithSpin->Sumw2();
	TH1F *hist_lepton_PT_WithSpin = new TH1F("hist_lepton_PT_WithSpin", "lepton_PT_WithSpin", 80, 0., 2.);
	hist_lepton_PT_WithSpin->Sumw2();
	TH2F *hist_PT_vs_lepton_PT_WithSpin = new TH2F("hist_PT_vs_lepton_PT_WithSpin", "PT_vs_lepton_EP_WithSpin", 80, 0., 2.,80.,0.,2.);
	hist_PT_vs_lepton_PT_WithSpin->Sumw2();
	TH1F *hist_twoneutrinos_PT_WithSpin = new TH1F("hist_twoneutrinos_PT_WithSpin", "twoneutrinos_PT_WithSpin", 80, 0., 2.);
	hist_twoneutrinos_PT_WithSpin->Sumw2();

	TH1F *hist_PT_Uniform = new TH1F("hist_PT_Uniform", "PT_Uniform", 80, 0., 2.);
	hist_PT_Uniform->Sumw2();
	TH1F *hist_lepton_PT_Uniform = new TH1F("hist_lepton_PT_Uniform", "lepton_PT_Uniform", 80, 0., 2.);
	hist_lepton_PT_Uniform->Sumw2();
	TH2F *hist_PT_vs_lepton_PT_Uniform = new TH2F("hist_PT_vs_lepton_PT_Uniform", "PT_vs_lepton_EP_Uniform", 80, 0., 2.,80.,0.,2.);
	hist_PT_vs_lepton_PT_Uniform->Sumw2();
	TH1F *hist_twoneutrinos_PT_Uniform = new TH1F("hist_twoneutrinos_PT_Uniform", "twoneutrinos_PT_Uniform", 80, 0., 2.);
	hist_twoneutrinos_PT_Uniform->Sumw2();

	TH1F *hist_sizeofbackgroundgammas = new TH1F("hist_sizeofbackgroundgammas", "sizeofbackgroundgammas", 5, 0., 5.);
	hist_sizeofbackgroundgammas->Sumw2();
	TH1F *hist_sizeofisphoton = new TH1F("hist_sizeofisphoton", "sizeofisphoton", 5, 0., 5.);
	hist_sizeofisphoton->Sumw2();


	TH1F *hist_gamma1phi = new TH1F("hist_gamma1phi", "gamma1phi", 20., -5., 5.);
	hist_gamma1phi->Sumw2();
	TH1F *hist_gamma1costheta = new TH1F("hist_gamma1costheta", "gamma1costheta", 20., -1., 1.);
	hist_gamma1costheta->Sumw2();
	TH1F *hist_gamma1sintheta = new TH1F("hist_gamma1sintheta", "gamma1sintheta", 20., -1., 1.);
	hist_gamma1sintheta->Sumw2();

	TH1F *hist_angle_separation = new TH1F("hist_angle_separation", "angle_separation", 100., 0., 3.14);
	hist_angle_separation->Sumw2();

	TH1F *hist_angle_separation2 = new TH1F("hist_angle_separation2", "angle_separation2", 100., 0., 3.14);
	hist_angle_separation2->Sumw2();

	TH1F *hist_2gammainvmass = new TH1F("hist_2gammainvmass", "2gammainvmass", 20., 0., 0.2);
	hist_2gammainvmass->Sumw2();

	TH1F *hist_vis_vertex = new TH1F("hist_vis_vertex", "hist_vis_vertex", 2., 0., 2.);
	hist_vis_vertex->Sumw2();

	TH2F *hist_convlength_vs_E = new TH2F("hist_convlength_vs_E", "convlength_vs_E",20,0.,10.,20,0.,100.);
	hist_convlength_vs_E->Sumw2();

	TH2F *hist_dalitz = new TH2F("hist_dalitz", "dalitz",200,0.,5.,200,0.,5.);
	hist_dalitz->Sumw2();

	TH2F *hist_dalitz_DecayIt = new TH2F("hist_dalitz_DecayIt", "dalitz_DecayIt",200,0.,5.,200,0.,5.);
	hist_dalitz_DecayIt->Sumw2();

	TH2F *hist_dalitz_DecayItWithSpin = new TH2F("hist_dalitz_DecayItWithSpin", "dalitz_DecayItWithSpin",200,0.,5.,200,0.,5.);
	hist_dalitz_DecayItWithSpin->Sumw2();


	TH1F *hist_convlength = new TH1F("hist_convlength", "convlength",20,0.,100.);
	hist_convlength->Sumw2();

	TH1F *hist_cc_ehad = new TH1F("hist_cc_ehad", "cc_ehad", 80, 0., 20.);
	hist_cc_ehad->Sumw2();
	TH1F *hist_cc_reco = new TH1F("hist_cc_reco", "cc_reco", 80, 0., 20.);
	hist_cc_reco->Sumw2();
	TH1F *hist_cc_El = new TH1F("hist_cc_El", "cc_El", 80, 0., 20.);
	hist_cc_El->Sumw2();

	TH1F *hist_cc_ehad_true = new TH1F("hist_cc_ehad_true", "cc_ehad_true", 80, 0., 20.);
	hist_cc_ehad_true->Sumw2();
	TH1F *hist_cc_reco_true = new TH1F("hist_cc_reco_true", "cc_reco_true", 80, 0., 20.);
	hist_cc_reco_true->Sumw2();
	TH1F *hist_cc_El_true = new TH1F("hist_cc_El_true", "cc_El_true", 40, 0., 20.);
	hist_cc_El_true->Sumw2();







	TH1F *hist_cc_ehad_pt05cut = new TH1F("hist_cc_ehad_pt05cut", "cc_ehad_pt05cut", 80, 0., 20.);
	hist_cc_ehad_pt05cut->Sumw2();
	TH1F *hist_cc_reco_pt05cut = new TH1F("hist_cc_reco_pt05cut", "cc_reco_pt05cut", 80, 0., 20.);
	hist_cc_reco_pt05cut->Sumw2();
	TH1F *hist_cc_El_pt05cut = new TH1F("hist_cc_El_pt05cut", "cc_El_pt05cut", 80, 0., 20.);
	hist_cc_El_pt05cut->Sumw2();

	TH1F *hist_cc_ehad_true_pt05cut = new TH1F("hist_cc_ehad_true_pt05cut", "cc_ehad_true_pt05cut", 80, 0., 20.);
	hist_cc_ehad_true_pt05cut->Sumw2();
	TH1F *hist_cc_reco_true_pt05cut = new TH1F("hist_cc_reco_true_pt05cut", "cc_reco_true_pt05cut", 80, 0., 20.);
	hist_cc_reco_true_pt05cut->Sumw2();
	TH1F *hist_cc_El_true_pt05cut = new TH1F("hist_cc_El_true_pt05cut", "cc_El_true_pt05cut", 40, 0., 20.);
	hist_cc_El_true_pt05cut->Sumw2();




	TH1F *hist_cc_ElCosTheta_true = new TH1F("hist_cc_ElCosTheta_true", "cc_ElCosTheta_true", 40, 0.6, 1.);
	hist_cc_ElCosTheta_true->Sumw2();



	TH1F *hist_cc_ElPhi_true = new TH1F("hist_cc_ElPhi_true", "cc_ElPhi_true", 80, -3.14, 3.14);
	hist_cc_ElPhi_true->Sumw2();
	TH2F *hist_cc_El_vs_CosTheta_true = new TH2F("hist_cc_El_vs_CosTheta_true", "cc_El_vs_CosTheta_true", 10, 0., 20., 10, 0.8, 1.0);
	hist_cc_El_vs_CosTheta_true->Sumw2();



	TH1F *hist_cc_ehad_WithSpin = new TH1F("hist_cc_ehad_WithSpin", "cc_ehad_WithSpin", 80, 0., 20.);
	hist_cc_ehad_WithSpin->Sumw2();
	TH1F *hist_cc_reco_WithSpin = new TH1F("hist_cc_reco_WithSpin", "cc_reco_WithSpin", 80, 0., 20.);
	hist_cc_reco_WithSpin->Sumw2();
	TH1F *hist_cc_El_WithSpin = new TH1F("hist_cc_El_WithSpin", "cc_El_WithSpin", 80, 0., 20.);
	hist_cc_El_WithSpin->Sumw2();
	TH1F *hist_cc_ehad_true_WithSpin = new TH1F("hist_cc_ehad_true_WithSpin", "cc_ehad_true_WithSpin", 80, 0., 20.);
	hist_cc_ehad_true_WithSpin->Sumw2();
	TH1F *hist_cc_reco_true_WithSpin = new TH1F("hist_cc_reco_true_WithSpin", "cc_reco_true_WithSpin", 80, 0., 20.);
	hist_cc_reco_true_WithSpin->Sumw2();
	TH1F *hist_cc_El_true_WithSpin = new TH1F("hist_cc_El_true_WithSpin", "cc_El_true_WithSpin", 40, 0., 20.);
	hist_cc_El_true_WithSpin->Sumw2();
	TH1F *hist_cc_ElCosTheta_true_WithSpin = new TH1F("hist_cc_ElCosTheta_true_WithSpin", "cc_ElCosTheta_true_WithSpin", 40, 0.6, 1.);
	hist_cc_ElCosTheta_true_WithSpin->Sumw2();
	TH1F *hist_cc_ElCosTheta_above2GeV_true_WithSpin = new TH1F("hist_cc_ElCosTheta_above2GeV_true_WithSpin", "cc_ElCosTheta_above2GeV_true_WithSpin", 40, 0.6, 1.);
	hist_cc_ElCosTheta_above2GeV_true_WithSpin->Sumw2();
	TH1F *hist_cc_ElPhi_true_WithSpin = new TH1F("hist_cc_ElPhi_true_WithSpin", "cc_ElPhi_true_WithSpin", 80, -3.14, 3.14);
	hist_cc_ElPhi_true_WithSpin->Sumw2();
	TH2F *hist_cc_El_vs_CosTheta_true_WithSpin = new TH2F("hist_cc_El_vs_CosTheta_true_WithSpin", "cc_El_vs_CosTheta_true_WithSpin", 10, 0., 20., 10, 0.8, 1.0);
	hist_cc_El_vs_CosTheta_true_WithSpin->Sumw2();


	TH1F *hist_cc_ehad_Uniform = new TH1F("hist_cc_ehad_Uniform", "cc_ehad_Uniform", 80, 0., 20.);
	hist_cc_ehad_Uniform->Sumw2();
	TH1F *hist_cc_reco_Uniform = new TH1F("hist_cc_reco_Uniform", "cc_reco_Uniform", 80, 0., 20.);
	hist_cc_reco_Uniform->Sumw2();
	TH1F *hist_cc_El_Uniform = new TH1F("hist_cc_El_Uniform", "cc_El_Uniform", 80, 0., 20.);
	hist_cc_El_Uniform->Sumw2();
	TH1F *hist_cc_ehad_true_Uniform = new TH1F("hist_cc_ehad_true_Uniform", "cc_ehad_true_Uniform", 80, 0., 20.);
	hist_cc_ehad_true_Uniform->Sumw2();
	TH1F *hist_cc_reco_true_Uniform = new TH1F("hist_cc_reco_true_Uniform", "cc_reco_true_Uniform", 80, 0., 20.);
	hist_cc_reco_true_Uniform->Sumw2();
	TH1F *hist_cc_El_true_Uniform = new TH1F("hist_cc_El_true_Uniform", "cc_El_true_Uniform", 40, 0., 20.);
	hist_cc_El_true_Uniform->Sumw2();
	TH1F *hist_cc_ElCosTheta_true_Uniform = new TH1F("hist_cc_ElCosTheta_true_Uniform", "cc_ElCosTheta_true_Uniform", 40, 0.6, 1.);
	hist_cc_ElCosTheta_true_Uniform->Sumw2();
	TH1F *hist_cc_ElPhi_true_Uniform = new TH1F("hist_cc_ElPhi_true_Uniform", "cc_ElPhi_true_Uniform", 80, -3.14, 3.14);
	hist_cc_ElPhi_true_Uniform->Sumw2();
	TH2F *hist_cc_El_vs_CosTheta_true_Uniform = new TH2F("hist_cc_El_vs_CosTheta_true_Uniform", "cc_El_vs_CosTheta_true_Uniform", 10, 0., 20., 10, 0.8, 1.0);
	hist_cc_El_vs_CosTheta_true_Uniform->Sumw2();



	TH1F *hist_cc_Ev_true = new TH1F("hist_cc_Ev_true", "cc_Ev_true", 40, 0., 40.);
	hist_cc_Ev_true->Sumw2();
	TH1F *hist_cc_tau_El_true = new TH1F("hist_cc_tau_El_true", "cc_tau_El_true", 80, 0., 40.);
	hist_cc_tau_El_true->Sumw2();
	TH1F *hist_cc_tau_ElCosTheta_true = new TH1F("hist_cc_tau_ElCosTheta_true", "cc_tau_ElCosTheta_true", 80, 0.8, 1.);
	hist_cc_tau_ElCosTheta_true->Sumw2();
	TH1F *hist_cc_tau_ElTheta_true = new TH1F("hist_cc_tau_ElTheta_true", "cc_tau_ElTheta_true", 80, 0., 1.);
	hist_cc_tau_ElTheta_true->Sumw2();
	TH1F *hist_cc_tau_ElPhi_true = new TH1F("hist_cc_tau_ElPhi_true", "cc_tau_ElPhi_true", 80, -3.14, 3.14);
	hist_cc_tau_ElPhi_true->Sumw2();


	TH1F *hist_tau_Ev_true = new TH1F("hist_tau_Ev_true", "tau_Ev_true", 40, 0., 40.);
	hist_tau_Ev_true->Sumw2();

	TH1F *hist_tau_Ev_osc = new TH1F("hist_tau_Ev_osc", "tau_Ev_osc", 40, 0., 40.);
	hist_tau_Ev_osc->Sumw2();


	TH1F *hist_nc_ehad = new TH1F("hist_nc_ehad", "nc_ehad", 80, 0., 20.);
	hist_nc_ehad->Sumw2();
	TH1F *hist_nc_reco = new TH1F("hist_nc_reco", "nc_reco", 80, 0., 20.);
	hist_nc_reco->Sumw2();
	TH1F *hist_nc_gamma = new TH1F("hist_nc_gamma", "nc_gamma", 80, 0., 20.);
	hist_nc_gamma->Sumw2();


	TH1F *hist_nc_ehad_true = new TH1F("hist_nc_ehad_true", "nc_ehad_true", 80, 0., 20.);
	hist_nc_ehad_true->Sumw2();
	TH1F *hist_nc_reco_true = new TH1F("hist_nc_reco_true", "nc_reco_true", 80, 0., 20.);
	hist_nc_reco_true->Sumw2();
	TH1F *hist_nc_gamma_true = new TH1F("hist_nc_gamma_true", "nc_gamma_true", 80, 0., 20.);
	hist_nc_gamma_true->Sumw2();


	SBN_detector DUNE(4,false);

	//SBN_detector DUNE_NC(5,false);

	for(int i=0; i < gstree->GetEntries(); i++){
		gstree->GetEntry(i);

		bool mode = true;
		if (neu<0){
			mode = false;
		}

		//Is there a visible vertex and how much energy is there!
		double Enu_reco = 0;
		double Ehad=0;
		double Ehad_true=0.;

		double Pxhad =0;
		double Pyhad =0;
		double Pzhad =0;

		bool vis_vertex = false;


		//tau test plot
		hist_tau_Ev_true->Fill(Ev);
		hist_tau_Ev_osc->Fill(Ev,prob_muflav(Ev,3,mode));


		if(nf!=0){
			double kin_true =0;
			double kin_smeared = 0;

			for(int j=0; j<nf; j++){

				if(pdgf[j]==2212){
					//cout << " proton mass : " << sqrt(energyf[j]*energyf[j] - pxf[j]*pxf[j]-pyf[j]*pyf[j]-pzf[j]*pzf[j])<< endl;
					kin_true=energyf[j]-MPROTON;
					//kin_smeared=smear_energy(kin_true, psmear, rangen);
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
					//cout << " pion mass : " << sqrt(energyf[j]*energyf[j] - pxf[j]*pxf[j]-pyf[j]*pyf[j]-pzf[j]*pzf[j])<< endl;
					kin_true=energyf[j]-MPION;
					//  kin_smeared=smear_energy(kin_true, pismear, rangen);
					kin_smeared=smear_energy_type(pdgf[j], kin_true, 1, rangen);
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
		hist_vis_vertex->Fill(vis_vertex);


		//Nrevertex. for stats.
		int Nrevertex = 25;
		for(int k=0;k<Nrevertex; k++){

			double vertex_weight = 1.0/((double) Nrevertex);


			double vertex_pos[3];

			//double vertex_pos_NC[3];

			DUNE.random_pos(rangen, vertex_pos);
			//DUNE_NC.random_pos(rangen, vertex_pos_NC);



			//photon flow
			std::vector<myphoton> gammas;
			std::vector<myphoton> gammas_true;

			//       std::vector<myphoton> gammas_smeared;

			// TLorentzVector final_gamma[10];
			//radiation gammas
			if (nf!=0){
				for(int j=0; j<nf;j++ ){

					if(pdgf[j]==22){ //energyf[j]>EM_thresh

						myphoton temp_vec;
						temp_vec.lorentz.SetPxPyPzE(pxf[j],pyf[j],pzf[j],energyf[j]);


						double temp_energy_smeared = smear_energy(energyf[j], EMsmear,rangen);

						if (temp_energy_smeared > EM_thresh) {

							gammas_true.push_back(temp_vec);

							temp_vec.lorentz.SetE(temp_energy_smeared);
							gammas.push_back(temp_vec);

							gammas.back().isPion=0;

						}
					}
				}
			}

			//pi0 gammas
			if (nf!=0){

				for(int k=0;k<nf;k++){

					if(pdgf[k]==111){

						double boost[3] = {pxf[k]/energyf[k],pyf[k]/energyf[k],pzf[k]/energyf[k]};//boost vector for pi0

						double x=0;
						double y=0;
						double z=0;

						rangen->Sphere(x,y,z,1.);//generate random direction on a unit sphere.

						myphoton gamma1, gamma2;


						gamma1.lorentz.SetPxPyPzE(x*pion_mass/2,y*pion_mass/2,z*pion_mass/2,pion_mass/2);// assigning gamma1 along the random direction
						gamma2.lorentz.SetPxPyPzE(-x*pion_mass/2,-y*pion_mass/2,-z*pion_mass/2,pion_mass/2);// assigning gamma2 to the opposite directon of gamma1, each of these has a half of pion rest energy.

						hist_gamma1phi->Fill(gamma1.lorentz.Phi());
						hist_gamma1costheta->Fill(cos(gamma1.lorentz.Theta()));
						hist_gamma1sintheta->Fill(sin(gamma1.lorentz.Theta()));
						//cos(gamma1.Theta());
						//cout << "gamma1 phi : " << gamma1.lorentz.Phi() << endl;

						gamma1.lorentz.Boost(pxf[k]/energyf[k],pyf[k]/energyf[k],pzf[k]/energyf[k]);//boost gamma1 and gamma2 to original pi0 frame.
						gamma2.lorentz.Boost(pxf[k]/energyf[k],pyf[k]/energyf[k],pzf[k]/energyf[k]);
						// check invariant amss of gamma1 gamm2 pair ==mpion
						hist_2gammainvmass->Fill((gamma1.lorentz+gamma2.lorentz).M());
						double angle_separation = gamma1.lorentz.Angle(gamma2.lorentz.BoostVector());
						double angle_separation2 = gamma1.lorentz.BoostVector().Angle(gamma2.lorentz.BoostVector());

						hist_angle_separation->Fill(angle_separation);
						hist_angle_separation2->Fill(angle_separation2);

						double temp_energy_smeared1= smear_energy(gamma1.lorentz.E(), EMsmear,rangen);
						double temp_energy_smeared2= smear_energy(gamma2.lorentz.E(), EMsmear,rangen);

						if (temp_energy_smeared1 > EM_thresh) {
							gamma1.isPion=1;
							gammas_true.push_back(gamma1);

							gamma1.lorentz.SetE(temp_energy_smeared1);
							gammas.push_back(gamma1);
							//pion_gamma_count++;

						}
						if (temp_energy_smeared2 > EM_thresh) {
							gamma2.isPion=1;

							gammas_true.push_back(gamma2);

							gamma2.lorentz.SetE(temp_energy_smeared2);
							gammas.push_back(gamma2);
							//pion_gamma_count++;

						}

					}


				}


			}

			if (gammas.size()>0){
				//cout << "# of gammas over EM_thresh : " << gammas.size() << ", first gamma energy : "<< gammas.at(0).lorentz.E()<< endl;
			}// we have everyting in gammas.

			double convlength[gammas.size()];
			double convpoint[gammas.size()][3];

			// cout << "vertex point x: " <<vertex_pos[0] << ", y: "<<vertex_pos[1]<<", z:"<< vertex_pos[2]<<  endl;


			for (int j=0; j< gammas.size(); j++){
				//convlength[j] = photon_conversion_length(gammas.at(j).lorentz.E(), rangen);
				convlength[j] = photon_conversion_length(gammas.at(j).lorentz.E(), rangen);
				gammas.at(j).convlength=convlength[j];

				double temp_3vec[3] = {gammas.at(j).lorentz.Px(),gammas.at(j).lorentz.Py(),gammas.at(j).lorentz.Pz()};


				int k = get_endpoint(vertex_pos, convlength[j], temp_3vec, convpoint[j]);
				gammas.at(j).convpoint=convpoint[j];

				//  cout<<"conversion length: " << convlength[j]<< ", point x: " <<convpoint[j][0] << ", y: "<<convpoint[j][1]<<", z:"<< convpoint[j][2]<<  endl;

				hist_convlength_vs_E->Fill(gammas.at(j).lorentz.E(),gammas.at(j).convlength,vertex_weight);
				hist_convlength->Fill(gammas.at(j).convlength,vertex_weight);

			}



			std::vector<myphoton> background_photons;
			std::vector<myphoton> background_photons_true;

			int photon_count = 0;
			int background_count = 0;

			for (int m=0; m<gammas.size(); m++){
				if (vis_vertex){
					if (gammas.at(m).convlength > convlength_thresh){
						if(DUNE.is_fiducial(gammas.at(m).convpoint)) photon_count++;
					}
					else{
						if(DUNE.is_fiducial(gammas.at(m).convpoint)){

							background_count++;
							background_photons.push_back(gammas.at(m));
							background_photons_true.push_back(gammas_true.at(m));
						}
					}

				}
				else {
					if(DUNE.is_fiducial(gammas.at(m).convpoint)){
						background_count++;
						background_photons.push_back(gammas.at(m));
						background_photons_true.push_back(gammas_true.at(m));
					}
				}
			}
			// photon flow end.




			hist_sizeofbackgroundgammas->Fill(background_photons.size(),vertex_weight);

			double Pxtot = Pxhad;
			double Pytot = Pyhad;
			double Pztot = Pzhad;

			TLorentzVector MissingEnergy;

			if (cc && (abs(neu)==12 || abs(neu)==14) ){
				Pxtot+=pxl;
				Pytot+=pyl;
				Pztot+=pzl;
			}




			// intrinsic
			if (cc && abs(neu)==12) {
				// hist_charlep->Fill(El);

				//cout << " electron mass : " << sqrt(El*El-pxl*pxl-pyl*pyl-pzl*pzl)<< endl;

				//double El_smear = smear_energy(El, EMsmear, rangen);
				double El_smear = smear_energy_type(11, El, rangen);
				if (El_smear > EM_thresh && vis_vertex && DUNE.is_fiducial(vertex_pos)) {


					MissingEnergy.SetXYZM(-Pxtot,-Pytot,-Pztot,0.);

					hist_PT->Fill(MissingEnergy.Et(),cc_efficiency*vertex_weight);

					hist_lepton_PT->Fill(sqrt(pxl*pxl+pyl*pyl),cc_efficiency*vertex_weight);

					hist_PT_vs_lepton_PT->Fill(MissingEnergy.Et(),sqrt(pxl*pxl+pyl*pyl),cc_efficiency*vertex_weight);


					hist_cc_El->Fill(El_smear,cc_efficiency*vertex_weight);
					hist_cc_ehad->Fill(Ehad,cc_efficiency*vertex_weight);
					hist_cc_reco->Fill(El_smear+Ehad,cc_efficiency*vertex_weight);

					hist_cc_El_true->Fill(El,cc_efficiency*vertex_weight);
					hist_cc_ehad_true->Fill(Ehad_true,cc_efficiency*vertex_weight);
					hist_cc_reco_true->Fill(El+Ehad_true,cc_efficiency*vertex_weight);

					TLorentzVector truelvec;

					truelvec.SetPxPyPzE(pxl,pyl,pzl,El);

					hist_cc_ElCosTheta_true->Fill(truelvec.CosTheta(),cc_efficiency*vertex_weight);


					if (MissingEnergy.Et()<0.5){

						hist_cc_El_pt05cut->Fill(El_smear,cc_efficiency*vertex_weight);
						hist_cc_ehad_pt05cut->Fill(Ehad,cc_efficiency*vertex_weight);
						hist_cc_reco_pt05cut->Fill(El_smear+Ehad,cc_efficiency*vertex_weight);

						hist_cc_El_true_pt05cut->Fill(El,cc_efficiency*vertex_weight);
						hist_cc_ehad_true_pt05cut->Fill(Ehad_true,cc_efficiency*vertex_weight);
						hist_cc_reco_true_pt05cut->Fill(El+Ehad_true,cc_efficiency*vertex_weight);


					}


				}
			}// end of nue cc.

			if (cc && abs(neu)==14) {
				// hist_charlep->Fill(El);


				double Lmu = muon_track_length(El);
				double endpos[3] = {0,0,0};
				double temp_3vec[3] = {pxl,pyl,pzl};
				get_endpoint(vertex_pos, Lmu, temp_3vec, endpos);

				double El_smear = smear_energy_type(13, DUNE.is_fully_contained(vertex_pos, endpos), rangen);


				//ehad+ e_mu -> for vis vertex check..

				double Ehad_mu = El_smear+Ehad;

				double Ehad_mu_true = El+Ehad_true;

				bool vis_vertex_mu = false;

				if(Ehad_mu >= vertex_thresh)
				{
					vis_vertex_mu = true;
				}

				if (El_smear > MU_thresh && vis_vertex_mu && DUNE.is_fiducial(vertex_pos)) {

					//numu cc events with a track of length >/ 1m were assumed to be identifiable as numu-induced CC events and were rejected. Those with a track length below 1m were accepted as potential mis-identified events, if nay photons in the event were accepted under the same conditions as in the NC single photon events, above.

					double observable_L = 0;

					if( DUNE.is_fully_contained(vertex_pos, endpos) ) {
						observable_L=Lmu;
					}
					else {
						observable_L = DUNE.track_length_escape(vertex_pos,endpos);
					}

					if(observable_L<100.0 && background_photons.size()==1 && photon_count ==0)
					{
						hist_cc_El->Fill(El_smear,cc_efficiency*vertex_weight*egamma_misidrate);
						hist_cc_ehad->Fill(Ehad_mu,cc_efficiency*vertex_weight*egamma_misidrate);
						hist_cc_reco->Fill(background_photons.at(0).lorentz.E()+Ehad_mu,cc_efficiency*vertex_weight*egamma_misidrate);

						hist_cc_El_true->Fill(El,cc_efficiency*vertex_weight*egamma_misidrate);
						hist_cc_ehad_true->Fill(Ehad_mu_true,cc_efficiency*vertex_weight*egamma_misidrate);
						hist_cc_reco_true->Fill(background_photons_true.at(0).lorentz.E()+Ehad_mu_true,cc_efficiency*vertex_weight*egamma_misidrate);

					}
				}
			}// end of numu cc.



			if (cc && abs(neu)==16) {
				// hist_charlep->Fill(El);



				//decay tau -> nu tau , nu el, el.
				TLorentzVector daughter[3];
				daughter[0].SetPxPyPzE(0.0,0.0,0.0,0.0);
				daughter[1].SetPxPyPzE(0.0,0.0,0.0,0.0);
				daughter[2].SetPxPyPzE(0.0,0.0,0.0,0.0);

				TLorentzVector daughterWithSpin[3];
				daughterWithSpin[0].SetPxPyPzE(0.0,0.0,0.0,0.0);
				daughterWithSpin[1].SetPxPyPzE(0.0,0.0,0.0,0.0);
				daughterWithSpin[2].SetPxPyPzE(0.0,0.0,0.0,0.0);

				TVector3 parent_polarization (0.0,0.0,1.0);

				double CosThetaP = 1.0;

				if(neu == 16){
					if (El<3.3){
						CosThetaP=0.9;
					}
					else if (El<4.2){
						CosThetaP=0.8;
					}
					else if (El<4.7){
						CosThetaP=0.7;
					}
					else if (El<5.1){
						CosThetaP=0.6;
					}
					else if (El<5.4){
						CosThetaP=0.5;
					}
					else if (El<5.7){
						CosThetaP=0.4;
					}
					else if (El<6.0){
						CosThetaP=0.3;
					}
					else if (El<6.3){
						CosThetaP=0.2;
					}
					else if (El<6.45){
						CosThetaP=0.1;
					}
					else if (El<6.8){
						CosThetaP=0.;
					}
					else if (El<6.9){
						CosThetaP=-0.1;
					}
					else if (El<7.1){
						CosThetaP=-0.2;
					}
					else if (El<7.3){
						CosThetaP=-0.3;
					}
					else if (El<7.5){
						CosThetaP=-0.4;
					}
					else if (El<7.7){
						CosThetaP=-0.5;
					}
					else if (El<7.8){
						CosThetaP=-0.6;
					}
					else if (El<8.3){
						CosThetaP=-0.7;
					}
					else if (El<8.6){
						CosThetaP=-0.8;
					}
					//else if (El<){
					//    CosThetaP=-0.9;
					//}
					else{
						CosThetaP=-1.0;
					}
				}
				double SinThetaP = sqrt(1.-CosThetaP*CosThetaP);

				parent_polarization.SetX(-SinThetaP);
				parent_polarization.SetZ(CosThetaP);

				double tau_mass = 1.76;

				TLorentzVector parentTau;
				parentTau.SetPxPyPzE(0.0,0.0,0.0,tau_mass);

				//threebodydecay(tau_mass,rangen, daughter);

				double p0[4] = {0,0,0,0};
				double p1[4] = {0,0,0,0};
				double p2[4] = {0,0,0,0};

				threebodydecay(rangen, tau_mass, 0.00051, 0., 0., p0, p1, p2);

				TLorentzVector daughterUniform[3];
				daughterUniform[0].SetPxPyPzE(p0[1],p0[2],p0[3],p0[0]);
				daughterUniform[1].SetPxPyPzE(p1[1],p1[2],p1[3],p1[0]);
				daughterUniform[2].SetPxPyPzE(p2[1],p2[2],p2[3],p2[0]);

				DecayIt(parentTau, rangen, daughter);

				DecayItWithSpin(parentTau, rangen, daughterWithSpin, parent_polarization);

				TLorentzVector daughter12, daughter23;
				TLorentzVector daughterWithSpin12, daughterWithSpin23;
				TLorentzVector daughterUniform12, daughterUniform23;


				daughter12 = daughter[2]+daughter[1];
				daughter23 = daughter[1]+daughter[0];

				daughterWithSpin12 = daughterWithSpin[2]+daughterWithSpin[1];
				daughterWithSpin23 = daughterWithSpin[1]+daughterWithSpin[0];

				daughterUniform12 = daughterUniform[2]+daughterUniform[1];
				daughterUniform23 = daughterUniform[1]+daughterUniform[0];



				hist_dalitz_DecayIt->Fill(daughter12.M()*daughter12.M(), daughter23.M()*daughter23.M() );
				hist_dalitz_DecayItWithSpin->Fill(daughterWithSpin12.M()*daughterWithSpin12.M(), daughterWithSpin23.M()*daughterWithSpin23.M() );
				hist_dalitz->Fill(daughterUniform12.M()*daughterUniform12.M(), daughterUniform23.M()*daughterUniform23.M());

				//daughter[0].SetXYZM(p0[1],p0[2],p0[3],0.);
				//daughter[1].SetXYZM(p1[1],p1[2],p1[3],0.);
				//daughter[2].SetXYZM(p2[1],p2[2],p2[3],0.00051);

				//tau vector in lab frame
				TLorentzVector labtau;
				labtau.SetPxPyPzE(pxl,pyl,pzl,El);


				hist_cc_Ev_true->Fill(Ev,prob_muflav(Ev,3,mode));
				hist_cc_tau_El_true->Fill(El,prob_muflav(Ev,3,mode));
				hist_cc_tau_ElCosTheta_true->Fill(labtau.CosTheta(),prob_muflav(Ev,3,mode));
				hist_cc_tau_ElTheta_true->Fill(labtau.Theta(),prob_muflav(Ev,3,mode));
				hist_cc_tau_ElPhi_true->Fill(labtau.Phi(),prob_muflav(Ev,3,mode));


				daughter[0].Boost(pxl/El,pyl/El,pzl/El);
				daughter[1].Boost(pxl/El,pyl/El,pzl/El);
				daughter[2].Boost(pxl/El,pyl/El,pzl/El);

				daughterWithSpin[0].Boost(pxl/El,pyl/El,pzl/El);
				daughterWithSpin[1].Boost(pxl/El,pyl/El,pzl/El);
				daughterWithSpin[2].Boost(pxl/El,pyl/El,pzl/El);

				daughterUniform[0].Boost(pxl/El,pyl/El,pzl/El);
				daughterUniform[1].Boost(pxl/El,pyl/El,pzl/El);
				daughterUniform[2].Boost(pxl/El,pyl/El,pzl/El);

				TLorentzVector twoneutrinos = daughter[1]+daughter[2];
				TLorentzVector twoneutrinosWithSpin = daughterWithSpin[1]+daughterWithSpin[2];
				TLorentzVector twoneutrinosUniform = daughterUniform[1]+daughterUniform[2];

				//Fill DecayIt
				//double El_smear;

				double El_smear = smear_energy_type(11, daughter[0].E(), rangen);
				//if (El_smear > EM_thresh && vis_vertex && DUNE.is_fiducial(vertex_pos)) {

				Pxtot+=daughter[0].Px();
				Pytot+=daughter[0].Py();
				Pztot+=daughter[0].Pz();

				MissingEnergy.SetXYZM(-Pxtot,-Pytot,-Pztot,0.);
				hist_lepton_PT->Fill(daughter[0].Pt(),prob_muflav(Ev,3,mode));
				hist_PT_vs_lepton_PT->Fill(MissingEnergy.Et(),daughter[0].Pt(),prob_muflav(Ev,3,mode));


				hist_cc_El->Fill(El_smear,cc_efficiency*vertex_weight*0.1783*prob_muflav(Ev,3,mode));
				hist_cc_ehad->Fill(Ehad,cc_efficiency*vertex_weight*0.1783*prob_muflav(Ev,3,mode));
				hist_cc_reco->Fill(El_smear+Ehad,cc_efficiency*vertex_weight*0.1783*prob_muflav(Ev,3,mode));

				hist_cc_El_true->Fill(daughter[0].E(),cc_efficiency*vertex_weight*0.1783*prob_muflav(Ev,3,mode));
				hist_cc_ehad_true->Fill(Ehad_true,cc_efficiency*vertex_weight*0.1783*prob_muflav(Ev,3,mode));
				hist_cc_reco_true->Fill(daughter[0].E()+Ehad_true,cc_efficiency*vertex_weight*0.1783*prob_muflav(Ev,3,mode));

				hist_cc_ElCosTheta_true->Fill(daughter[0].CosTheta(), cc_efficiency*vertex_weight*0.1783*prob_muflav(Ev,3,mode));

				hist_cc_ElPhi_true->Fill(daughter[0].Phi(), cc_efficiency*vertex_weight*0.1783*prob_muflav(Ev,3,mode));
				hist_cc_El_vs_CosTheta_true->Fill(daughter[0].E(), daughter[0].CosTheta(), cc_efficiency*vertex_weight*0.1783*prob_muflav(Ev,3,mode));


				hist_PT->Fill(MissingEnergy.Et(),cc_efficiency*vertex_weight*0.1783*prob_muflav(Ev,3,mode));
				hist_twoneutrinos_PT->Fill(twoneutrinos.Et(),cc_efficiency*vertex_weight*0.1783*prob_muflav(Ev,3,mode));

				if (MissingEnergy.Et()<0.5){


					//_pt05cut

					hist_cc_El_pt05cut->Fill(El_smear,cc_efficiency*vertex_weight*0.1783*prob_muflav(Ev,3,mode));
					hist_cc_ehad_pt05cut->Fill(Ehad,cc_efficiency*vertex_weight*0.1783*prob_muflav(Ev,3,mode));
					hist_cc_reco_pt05cut->Fill(El_smear+Ehad,cc_efficiency*vertex_weight*0.1783*prob_muflav(Ev,3,mode));

					hist_cc_El_true_pt05cut->Fill(daughter[0].E(),cc_efficiency*vertex_weight*0.1783*prob_muflav(Ev,3,mode));
					hist_cc_ehad_true_pt05cut->Fill(Ehad_true,cc_efficiency*vertex_weight*0.1783*prob_muflav(Ev,3,mode));
					hist_cc_reco_true_pt05cut->Fill(daughter[0].E()+Ehad_true,cc_efficiency*vertex_weight*0.1783*prob_muflav(Ev,3,mode));





				}

				Pxtot-=daughter[0].Px();
				Pytot-=daughter[0].Py();
				Pztot-=daughter[0].Pz();

				//}

				//Fill DecayItWithSpin
				El_smear = smear_energy_type(11, daughterWithSpin[0].E(), rangen);
				//if (El_smear > EM_thresh && vis_vertex && DUNE.is_fiducial(vertex_pos)) {

				Pxtot+=daughterWithSpin[0].Px();
				Pytot+=daughterWithSpin[0].Py();
				Pztot+=daughterWithSpin[0].Pz();

				MissingEnergy.SetXYZM(-Pxtot,-Pytot,-Pztot,0.);
				hist_lepton_PT_WithSpin->Fill(daughterWithSpin[0].Pt(),prob_muflav(Ev,3,mode));
				hist_PT_vs_lepton_PT_WithSpin->Fill(MissingEnergy.Et(),daughterWithSpin[0].Pt(),prob_muflav(Ev,3,mode));


				hist_cc_El_WithSpin->Fill(El_smear,cc_efficiency*vertex_weight*0.1783*prob_muflav(Ev,3,mode));
				hist_cc_ehad_WithSpin->Fill(Ehad,cc_efficiency*vertex_weight*0.1783*prob_muflav(Ev,3,mode));
				hist_cc_reco_WithSpin->Fill(El_smear+Ehad,cc_efficiency*vertex_weight*0.1783*prob_muflav(Ev,3,mode));

				hist_cc_El_true_WithSpin->Fill(daughterWithSpin[0].E(),cc_efficiency*vertex_weight*0.1783*prob_muflav(Ev,3,mode));
				hist_cc_ehad_true_WithSpin->Fill(Ehad_true,cc_efficiency*vertex_weight*0.1783*prob_muflav(Ev,3,mode));
				hist_cc_reco_true_WithSpin->Fill(daughterWithSpin[0].E()+Ehad_true,cc_efficiency*vertex_weight*0.1783*prob_muflav(Ev,3,mode));

				hist_cc_ElCosTheta_true_WithSpin->Fill(daughterWithSpin[0].CosTheta(), cc_efficiency*vertex_weight*0.1783*prob_muflav(Ev,3,mode));
				if (daughterWithSpin[0].E() > 2.) {
					hist_cc_ElCosTheta_above2GeV_true_WithSpin->Fill(daughterWithSpin[0].CosTheta(), cc_efficiency*vertex_weight*0.1783*prob_muflav(Ev,3,mode));
				}
				hist_cc_ElPhi_true_WithSpin->Fill(daughterWithSpin[0].Phi(), cc_efficiency*vertex_weight*0.1783*prob_muflav(Ev,3,mode));
				hist_cc_El_vs_CosTheta_true_WithSpin->Fill(daughterWithSpin[0].E(), daughterWithSpin[0].CosTheta(), cc_efficiency*vertex_weight*0.1783*prob_muflav(Ev,3,mode));

				hist_PT_WithSpin->Fill(MissingEnergy.Et(),cc_efficiency*vertex_weight*0.1783*prob_muflav(Ev,3,mode));
				hist_twoneutrinos_PT_WithSpin->Fill(twoneutrinosWithSpin.Et(),cc_efficiency*vertex_weight*0.1783*prob_muflav(Ev,3,mode));



				Pxtot-=daughterWithSpin[0].Px();
				Pytot-=daughterWithSpin[0].Py();
				Pztot-=daughterWithSpin[0].Pz();

				//}

				//Fill Uniform
				El_smear = smear_energy_type(11, daughterUniform[0].E(), rangen);
				//if (El_smear > EM_thresh && vis_vertex && DUNE.is_fiducial(vertex_pos)) {

				Pxtot+=daughterUniform[0].Px();
				Pytot+=daughterUniform[0].Py();
				Pztot+=daughterUniform[0].Pz();

				MissingEnergy.SetXYZM(-Pxtot,-Pytot,-Pztot,0.);
				hist_lepton_PT_Uniform->Fill(daughterUniform[0].Pt(),prob_muflav(Ev,3,mode));
				hist_PT_vs_lepton_PT_Uniform->Fill(MissingEnergy.Et(),daughterUniform[0].Pt(),prob_muflav(Ev,3,mode));


				hist_cc_El_Uniform->Fill(El_smear,cc_efficiency*vertex_weight*0.1783*prob_muflav(Ev,3,mode));
				hist_cc_ehad_Uniform->Fill(Ehad,cc_efficiency*vertex_weight*0.1783*prob_muflav(Ev,3,mode));
				hist_cc_reco_Uniform->Fill(El_smear+Ehad,cc_efficiency*vertex_weight*0.1783*prob_muflav(Ev,3,mode));

				hist_cc_El_true_Uniform->Fill(daughterUniform[0].E(),cc_efficiency*vertex_weight*0.1783*prob_muflav(Ev,3,mode));
				hist_cc_ehad_true_Uniform->Fill(Ehad_true,cc_efficiency*vertex_weight*0.1783*prob_muflav(Ev,3,mode));
				hist_cc_reco_true_Uniform->Fill(daughterUniform[0].E()+Ehad_true,cc_efficiency*vertex_weight*0.1783*prob_muflav(Ev,3,mode));

				hist_cc_ElCosTheta_true_Uniform->Fill(daughterUniform[0].CosTheta(), cc_efficiency*vertex_weight*0.1783*prob_muflav(Ev,3,mode));
				hist_cc_ElPhi_true_Uniform->Fill(daughterUniform[0].Phi(), cc_efficiency*vertex_weight*0.1783*prob_muflav(Ev,3,mode));
				hist_cc_El_vs_CosTheta_true_Uniform->Fill(daughterUniform[0].E(), daughterUniform[0].CosTheta(), cc_efficiency*vertex_weight*0.1783*prob_muflav(Ev,3,mode));


				hist_PT_Uniform->Fill(MissingEnergy.Et(),cc_efficiency*vertex_weight*0.1783*prob_muflav(Ev,3,mode));
				hist_twoneutrinos_PT_Uniform->Fill(twoneutrinosUniform.Et(),cc_efficiency*vertex_weight*0.1783*prob_muflav(Ev,3,mode));



				Pxtot-=daughterUniform[0].Px();
				Pytot-=daughterUniform[0].Py();
				Pztot-=daughterUniform[0].Pz();

				//}




			}// end of nutau cc.



			// NC
			if (nc) {

				// cout << "# of gammas in evt : " << final_gamma_count+pion_gamma_count << endl;

				//hist_nofgamma->Fill(final_gamma_count+pion_gamma_count,vertex_weight);
				hist_sizeofgammas->Fill(gammas.size(),vertex_weight);

				hist_sizeofisphoton->Fill(photon_count,vertex_weight);

				// NC single gamma background
				if (background_photons.size()==1 && photon_count ==0){
					hist_nc_gamma->Fill(background_photons.at(0).lorentz.E(),egamma_misidrate*vertex_weight);
					hist_nc_reco->Fill(background_photons.at(0).lorentz.E()+Ehad,egamma_misidrate*vertex_weight);
					hist_nc_ehad->Fill(Ehad,egamma_misidrate*vertex_weight);

					hist_nc_gamma_true->Fill(background_photons_true.at(0).lorentz.E(),egamma_misidrate*vertex_weight);
					hist_nc_reco_true->Fill(background_photons_true.at(0).lorentz.E()+Ehad_true,egamma_misidrate*vertex_weight);
					hist_nc_ehad_true->Fill(Ehad_true,egamma_misidrate*vertex_weight);
					//std::cout<<"#######"<<background_photons[1].lorentz.E()<<" "<<background_photons[0].lorentz.E()<<std::endl;
					//std::cout<<"#######"<<background_photons.at(0).lorentz.E()<<std::endl;
					//if(background_photons.at(0).isPion) hist_backgroundsinglegamma_test->Fill(background_photons.at(0).lorentz.E(),egamma_misidrate*vertex_weight);
				}
			}//NC end.

		}//nrevertex end.


	}//event loop end.

	fnew->Write();

}


void genie_NC(TString filename){



	/*"The signal for nue apperance is an excess of charged-current(CC) nue and nuebar interactions over the expected background in the far detector. The background to nue appearance is composed of : (1) CC interactions of nue and nuebar intrinsic to the beam; (2) misidentified numu and numubar CC events; (3) neutral current (NC) backgrounds and (4) nutau and nutaubar CC events in which the tau's decay leptonically into electrons/positrons. NC and nutau backgrounds are due to interactions of higher-energy neutrinos but they contribute to backgrounds mainly at low energy, which is important for the sensitivity to CP violation."*/

	TRandom3 *rangen = new TRandom3(0);


//	TString filename = "gntp.0.numubar10k_gst";

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

	TFile *fnew = new TFile(TString("../gst0to40/out/")+filename+TString("_NC_study_out.root"),"recreate");

	TFile * fntuple = new TFile("DUNE_ntuple.root","UPDATE");



	std::vector<std::string> subchannels = {"nu_dune_elike_ncmisid"};
	std::vector<TTree*> list_o_trees;

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
	double Ev;
	double pxl;
	double pyl;
	double pzl;
	int nfpi0;
	int Np;
	int Npip;
	int Npim;
	int Nph;
	int No;
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
	double cc_efficiency=0.8;
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

	TH1F *hist_sizeofgammas = new TH1F("hist_sizeofgammas", "sizeofgammas", 10, 0., 10.);
	hist_sizeofgammas->Sumw2();

	TH1F *hist_PT = new TH1F("hist_PT", "PT", 80, 0., 2.);
	hist_PT->Sumw2();

	TH1F *hist_lepton_PT = new TH1F("hist_lepton_PT", "lepton_PT", 80, 0., 2.);
	hist_lepton_PT->Sumw2();

	TH2F *hist_PT_vs_lepton_PT = new TH2F("hist_PT_vs_lepton_PT", "PT_vs_lepton_EP", 80, 0., 2.,80.,0.,2.);
	hist_PT_vs_lepton_PT->Sumw2();

	TH1F *hist_twoneutrinos_PT = new TH1F("hist_twoneutrinos_PT", "twoneutrinos_PT", 80, 0., 2.);
	hist_twoneutrinos_PT->Sumw2();

	TH1F *hist_sizeofbackgroundgammas = new TH1F("hist_sizeofbackgroundgammas", "sizeofbackgroundgammas", 5, 0., 5.);
	hist_sizeofbackgroundgammas->Sumw2();
	TH1F *hist_sizeofisphoton = new TH1F("hist_sizeofisphoton", "sizeofisphoton", 5, 0., 5.);
	hist_sizeofisphoton->Sumw2();


	//TH1F *hist_singlegamma = new TH1F("hist_singlegamma", "singlegamma", 80, 0., 20.);
	//hist_singlegamma->Sumw2();



	//TH1F *hist_backgroundsinglegamma_test = new TH1F("hist_backgroundsinglegamma_test", "backgroundsinglegamma_test", 80, 0., 20.);
	//hist_backgroundsinglegamma_test->Sumw2();

	TH1F *hist_gamma1phi = new TH1F("hist_gamma1phi", "gamma1phi", 20., -5., 5.);
	hist_gamma1phi->Sumw2();
	TH1F *hist_gamma1costheta = new TH1F("hist_gamma1costheta", "gamma1costheta", 20., -1., 1.);
	hist_gamma1costheta->Sumw2();
	TH1F *hist_gamma1sintheta = new TH1F("hist_gamma1sintheta", "gamma1sintheta", 20., -1., 1.);
	hist_gamma1sintheta->Sumw2();

	TH1F *hist_angle_separation = new TH1F("hist_angle_separation", "angle_separation", 100., 0., 3.14);
	hist_angle_separation->Sumw2();

	TH1F *hist_angle_separation2 = new TH1F("hist_angle_separation2", "angle_separation2", 100., 0., 3.14);
	hist_angle_separation2->Sumw2();

	TH1F *hist_2gammainvmass = new TH1F("hist_2gammainvmass", "2gammainvmass", 20., 0., 0.2);
	hist_2gammainvmass->Sumw2();

	TH1F *hist_vis_vertex = new TH1F("hist_vis_vertex", "hist_vis_vertex", 2., 0., 2.);
	hist_vis_vertex->Sumw2();

	TH2F *hist_convlength_vs_E = new TH2F("hist_convlength_vs_E", "convlength_vs_E",20,0.,10.,20,0.,100.);
	hist_convlength_vs_E->Sumw2();

	TH2F *hist_dalitz = new TH2F("hist_dalitz", "dalitz",200,0.,5.,200,0.,5.);
	hist_dalitz->Sumw2();


	TH1F *hist_convlength = new TH1F("hist_convlength", "convlength",20,0.,100.);
	hist_convlength->Sumw2();

	TH1F *hist_cc_ehad = new TH1F("hist_cc_ehad", "cc_ehad", 80, 0., 20.);
	hist_cc_ehad->Sumw2();
	TH1F *hist_cc_reco = new TH1F("hist_cc_reco", "cc_reco", 80, 0., 20.);
	hist_cc_reco->Sumw2();
	TH1F *hist_cc_El = new TH1F("hist_cc_El", "cc_El", 80, 0., 20.);
	hist_cc_El->Sumw2();

	TH1F *hist_cc_ehad_true = new TH1F("hist_cc_ehad_true", "cc_ehad_true", 80, 0., 20.);
	hist_cc_ehad_true->Sumw2();
	TH1F *hist_cc_reco_true = new TH1F("hist_cc_reco_true", "cc_reco_true", 80, 0., 20.);
	hist_cc_reco_true->Sumw2();
	TH1F *hist_cc_El_true = new TH1F("hist_cc_El_true", "cc_El_true", 80, 0., 20.);
	hist_cc_El_true->Sumw2();



	TH1F *hist_nc_ehad = new TH1F("hist_nc_ehad", "nc_ehad", 80, 0., 20.);
	hist_nc_ehad->Sumw2();
	TH1F *hist_nc_reco = new TH1F("hist_nc_reco", "nc_reco", 80, 0., 20.);
	hist_nc_reco->Sumw2();
	TH1F *hist_nc_gamma = new TH1F("hist_nc_gamma", "nc_gamma", 80, 0., 20.);
	hist_nc_gamma->Sumw2();


	TH1F *hist_nc_ehad_true = new TH1F("hist_nc_ehad_true", "nc_ehad_true", 80, 0., 20.);
	hist_nc_ehad_true->Sumw2();
	TH1F *hist_nc_reco_true = new TH1F("hist_nc_reco_true", "nc_reco_true", 80, 0., 20.);
	hist_nc_reco_true->Sumw2();
	TH1F *hist_nc_gamma_true = new TH1F("hist_nc_gamma_true", "nc_gamma_true", 80, 0., 20.);
	hist_nc_gamma_true->Sumw2();



	SBN_detector DUNE_NC(5,false);

	for(int i=0; i < gstree->GetEntries(); i++){
		gstree->GetEntry(i);

		bool mode = true;
		if (neu<0){
			mode = false;
		}

		if (cc) {
			continue;
		}


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
					//cout << " proton mass : " << sqrt(energyf[j]*energyf[j] - pxf[j]*pxf[j]-pyf[j]*pyf[j]-pzf[j]*pzf[j])<< endl;
					kin_true=energyf[j]-MPROTON;
					//kin_smeared=smear_energy(kin_true, psmear, rangen);
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
					//cout << " pion mass : " << sqrt(energyf[j]*energyf[j] - pxf[j]*pxf[j]-pyf[j]*pyf[j]-pzf[j]*pzf[j])<< endl;
					kin_true=energyf[j]-MPION;
					//  kin_smeared=smear_energy(kin_true, pismear, rangen);
					kin_smeared=smear_energy_type(pdgf[j], kin_true, 1, rangen);
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
		hist_vis_vertex->Fill(vis_vertex);


		//Nrevertex. for stats.
		int Nrevertex = 25;
		for(int k=0;k<Nrevertex; k++){

			double vertex_weight = 1.0/((double) Nrevertex);


			//double vertex_pos[3];

			double vertex_pos_NC[3];

			//DUNE.random_pos(rangen, vertex_pos);
			DUNE_NC.random_pos(rangen, vertex_pos_NC);



			//photon flow
			std::vector<myphoton> gammas;
			std::vector<myphoton> gammas_true;

			//       std::vector<myphoton> gammas_smeared;

			// TLorentzVector final_gamma[10];
			//radiation gammas
			if (nf!=0){
				for(int j=0; j<nf;j++ ){

					if(pdgf[j]==22){ //energyf[j]>EM_thresh

						myphoton temp_vec;
						temp_vec.lorentz.SetPxPyPzE(pxf[j],pyf[j],pzf[j],energyf[j]);


						double temp_energy_smeared = smear_energy(energyf[j], EMsmear,rangen);

						if (temp_energy_smeared > EM_thresh) {

							gammas_true.push_back(temp_vec);

							temp_vec.lorentz.SetE(temp_energy_smeared);
							gammas.push_back(temp_vec);

							gammas.back().isPion=0;

						}
					}
				}
			}

			//pi0 gammas
			if (nf!=0){

				for(int k=0;k<nf;k++){

					if(pdgf[k]==111){

						double boost[3] = {pxf[k]/energyf[k],pyf[k]/energyf[k],pzf[k]/energyf[k]};//boost vector for pi0

						double x=0;
						double y=0;
						double z=0;

						rangen->Sphere(x,y,z,1.);//generate random direction on a unit sphere.

						myphoton gamma1, gamma2;


						gamma1.lorentz.SetPxPyPzE(x*pion_mass/2,y*pion_mass/2,z*pion_mass/2,pion_mass/2);// assigning gamma1 along the random direction
						gamma2.lorentz.SetPxPyPzE(-x*pion_mass/2,-y*pion_mass/2,-z*pion_mass/2,pion_mass/2);// assigning gamma2 to the opposite directon of gamma1, each of these has a half of pion rest energy.

						hist_gamma1phi->Fill(gamma1.lorentz.Phi());
						hist_gamma1costheta->Fill(cos(gamma1.lorentz.Theta()));
						hist_gamma1sintheta->Fill(sin(gamma1.lorentz.Theta()));
						//cos(gamma1.Theta());
						//cout << "gamma1 phi : " << gamma1.lorentz.Phi() << endl;

						gamma1.lorentz.Boost(pxf[k]/energyf[k],pyf[k]/energyf[k],pzf[k]/energyf[k]);//boost gamma1 and gamma2 to original pi0 frame.
						gamma2.lorentz.Boost(pxf[k]/energyf[k],pyf[k]/energyf[k],pzf[k]/energyf[k]);
						// check invariant amss of gamma1 gamm2 pair ==mpion
						hist_2gammainvmass->Fill((gamma1.lorentz+gamma2.lorentz).M());
						double angle_separation = gamma1.lorentz.Angle(gamma2.lorentz.BoostVector());
						double angle_separation2 = gamma1.lorentz.BoostVector().Angle(gamma2.lorentz.BoostVector());

						hist_angle_separation->Fill(angle_separation);
						hist_angle_separation2->Fill(angle_separation2);

						double temp_energy_smeared1= smear_energy(gamma1.lorentz.E(), EMsmear,rangen);
						double temp_energy_smeared2= smear_energy(gamma2.lorentz.E(), EMsmear,rangen);

						if (temp_energy_smeared1 > EM_thresh) {
							gamma1.isPion=1;
							gammas_true.push_back(gamma1);

							gamma1.lorentz.SetE(temp_energy_smeared1);
							gammas.push_back(gamma1);
							//pion_gamma_count++;

						}
						if (temp_energy_smeared2 > EM_thresh) {
							gamma2.isPion=1;

							gammas_true.push_back(gamma2);

							gamma2.lorentz.SetE(temp_energy_smeared2);
							gammas.push_back(gamma2);
							//pion_gamma_count++;

						}

					}


				}


			}

			if (gammas.size()>0){
				//cout << "# of gammas over EM_thresh : " << gammas.size() << ", first gamma energy : "<< gammas.at(0).lorentz.E()<< endl;
			}// we have everyting in gammas.

			double convlength[gammas.size()];
			double convpoint[gammas.size()][3];

			// cout << "vertex point x: " <<vertex_pos[0] << ", y: "<<vertex_pos[1]<<", z:"<< vertex_pos[2]<<  endl;


			for (int j=0; j< gammas.size(); j++){
				//convlength[j] = photon_conversion_length(gammas.at(j).lorentz.E(), rangen);
				convlength[j] = photon_conversion_length(gammas.at(j).lorentz.E(), rangen);
				gammas.at(j).convlength=convlength[j];

				double temp_3vec[3] = {gammas.at(j).lorentz.Px(),gammas.at(j).lorentz.Py(),gammas.at(j).lorentz.Pz()};


				int k = get_endpoint(vertex_pos_NC, convlength[j], temp_3vec, convpoint[j]);
				gammas.at(j).convpoint=convpoint[j];

				//  cout<<"conversion length: " << convlength[j]<< ", point x: " <<convpoint[j][0] << ", y: "<<convpoint[j][1]<<", z:"<< convpoint[j][2]<<  endl;

				hist_convlength_vs_E->Fill(gammas.at(j).lorentz.E(),gammas.at(j).convlength,vertex_weight);
				hist_convlength->Fill(gammas.at(j).convlength,vertex_weight);

			}



			std::vector<myphoton> background_photons;
			std::vector<myphoton> background_photons_true;

			int photon_count = 0;
			int background_count = 0;

			for (int m=0; m<gammas.size(); m++){
				if (vis_vertex){
					if (gammas.at(m).convlength > convlength_thresh){
						if(DUNE_NC.is_fiducial(gammas.at(m).convpoint)) photon_count++;
					}
					else{
						if(DUNE_NC.is_fiducial(gammas.at(m).convpoint)){

							background_count++;
							background_photons.push_back(gammas.at(m));
							background_photons_true.push_back(gammas_true.at(m));
						}
					}

				}
				else {
					if(DUNE_NC.is_fiducial(gammas.at(m).convpoint)){
						background_count++;
						background_photons.push_back(gammas.at(m));
						background_photons_true.push_back(gammas_true.at(m));
					}
				}
			}





			hist_sizeofbackgroundgammas->Fill(background_photons.size(),vertex_weight);

			double Pxtot = Pxhad;
			double Pytot = Pyhad;
			double Pztot = Pzhad;

			TLorentzVector MissingEnergy;





			// NC
			if (nc) {

				// cout << "# of gammas in evt : " << final_gamma_count+pion_gamma_count << endl;

				//hist_nofgamma->Fill(final_gamma_count+pion_gamma_count,vertex_weight);
				hist_sizeofgammas->Fill(gammas.size(),vertex_weight);

				hist_sizeofisphoton->Fill(photon_count,vertex_weight);

				// NC single gamma background
				if (background_photons.size()==1 && photon_count ==0){
					hist_nc_gamma->Fill(background_photons.at(0).lorentz.E(),egamma_misidrate*vertex_weight);
					hist_nc_reco->Fill(background_photons.at(0).lorentz.E()+Ehad,egamma_misidrate*vertex_weight);
					hist_nc_ehad->Fill(Ehad,egamma_misidrate*vertex_weight);

					hist_nc_gamma_true->Fill(background_photons_true.at(0).lorentz.E(),egamma_misidrate*vertex_weight);
					hist_nc_reco_true->Fill(background_photons_true.at(0).lorentz.E()+Ehad_true,egamma_misidrate*vertex_weight);
					hist_nc_ehad_true->Fill(Ehad_true,egamma_misidrate*vertex_weight);
					//std::cout<<"#######"<<background_photons[1].lorentz.E()<<" "<<background_photons[0].lorentz.E()<<std::endl;
					//std::cout<<"#######"<<background_photons.at(0).lorentz.E()<<std::endl;
					//if(background_photons.at(0).isPion) hist_backgroundsinglegamma_test->Fill(background_photons.at(0).lorentz.E(),egamma_misidrate*vertex_weight);
				
					m_Ereco = background_photons.at(0).lorentz.E()+Ehad;
					m_Etrue = Ev;
					m_l = 1300 ;	
					m_weight = egamma_misidrate*vertex_weight;
					m_nutype = neu;
					//std::cout<<"NCC: "<<m_Ereco<<" "<<m_Etrue<<" "<<m_l<<" "<<m_weight<<" "<<m_nutype<<std::endl;
					list_o_trees.at(0)->Fill();



				}
			}//NC end.

		}//nrevertex end.


	}//event loop end.

	fnew->cd();
	fnew->Write();

	fntuple->cd();
	for(auto &t: list_o_trees){
		t->Write();
	}
	fntuple->Purge();
	fntuple->Close();


}




void run_all_genie_study(){

	TString nutaubar = "gntp.0.nutaubar20k_gst";
	TString nutau = "gntp.0.nutau20k_gst";
	TString nue  = "gntp.0.nue20k_gst";
	TString nuebar = "gntp.0.nuebar20k_gst";
	TString numu  = "gntp.0.numu50k_gst";
	TString numubar = "gntp.0.numubar10k_gst";

	TString numu_nuebeam = "gntp.0.numuflux_nuebeam50k_gst";
	TString numubar_nuebarbeam = "gntp.0.numubarflux_nuebarbeam10k_gst";



	std::cout<<"Starting CC nue."<<std::endl;
	genie_study(nue);
	std::cout<<"Starting CC numu."<<std::endl;
	genie_study(numu);
	std::cout<<"Starting CC nutau."<<std::endl;
	genie_study(nutau);

	std::cout<<"Starting CC nuebar."<<std::endl;
	genie_study(nuebar);
	std::cout<<"Starting CC numubar."<<std::endl;
	genie_study(numubar);
	std::cout<<"Starting CC nutaubar."<<std::endl;
	genie_study(nutaubar);

	std::cout<<"Starting wierd CC numu_nuebeam."<<std::endl;
	genie_study(numu_nuebeam);
	std::cout<<"Starting wierd CC numubear_nuebarbeam."<<std::endl;
	genie_study(numubar_nuebarbeam);



	std::cout<<"Starting NC nue."<<std::endl;
	genie_NC(nue);
	std::cout<<"Starting NC numu."<<std::endl;
	genie_NC(numu);
	std::cout<<"Starting NC nutau."<<std::endl;
	genie_NC(nutau);

	std::cout<<"Starting NC nuebar."<<std::endl;
	genie_NC(nuebar);
	std::cout<<"Starting NC numubar."<<std::endl;
	genie_NC(numubar);
	std::cout<<"Starting NC nutaubar."<<std::endl;
	genie_NC(nutaubar);

	//std::cout<<"Starting wierd NC numu_nuebeam."<<std::endl;
	//genie_NC(numu_nuebeam);
	//std::cout<<"Starting wierd NC numubear_nuebarbeam."<<std::endl;
	//genie_NC(numubar_nuebarbeam);



}




