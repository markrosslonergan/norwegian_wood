#ifndef DETECTOR_H_
#define DETECTOR_H_

//#include "TRandom.h"
#include "TRandom3.h"
#include <vector>
#include <string>
#include "TLorentzVector.h"//yj


#include "params.h"

/*************************************************************
 *************************************************************
 *	TODO:
 *	    (4) overload detectors so I can just pass identifiers DONE
 ************************************************************
 ************************************************************/
	double smear_energy(double En, double Percen, TRandom3 * rangen);

    double smear_energy_type(int pdgid, double En, TRandom3 * rangen);//yj
    double smear_energy_type(int pdgid, double En, bool fully_contained, TRandom3 * rangen);//yj
	double smear_angle(double the, double ang, TRandom3 *rangen);
	double massive_smear_energy(double En, double Percen, TRandom3 * rangen,double mass);
	
	double muon_track_length(double El);
	double pion_track_length(double El);
	double photon_conversion_length(double ep, TRandom3 * r);
	double bethe(double beta);

	double CSDA_integrand(double Emu);
	double pion_containment(double posX, double posY, double posZ, TRandom3 * r);

	int get_endpoint(double *vertex,double track_L,double * pl,double *  endpoint);

bool isconvpointinfiducial(TLorentzVector gamma,double *vertex_pos,TRandom3 *rangen);//yj



class SBN_detector {
	double dh, dw, dl;


	public:
	char const * name;
	bool mumode;
	double height, width, length;
	double f_height, f_width, f_length;
        double volume;
	double f_volume;
	double mass;
	double f_mass;	
	double baseline;
	char const * fname; // location of root ntuple containing data file		
	char const * foscname;	//location of full oscillated
	char const * fbarname; // location of root ntuple containing data file		
	char const * fbaroscname;	//location of full oscillated 
	double potmodifier;
	int identifier;
	double proposal_modifier;

	//Constructors, if blank corresponds to background only. 
	SBN_detector (double h, double w, double l, double fh, double fw, double fl,double base);
	SBN_detector (int identifier, bool ismue = false);

	double osc_length(TRandom3 * rangen);

	bool is_active(double * pos);
	bool is_fiducial(double * pos);
	void random_pos(TRandom3 * rangen, double * vec);

	double track_length_escape(double * inside, double * outside);
	bool is_fully_contained(double *vertex,double * endpoint);
	


};



#endif
