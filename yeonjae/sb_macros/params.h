#ifndef params_H_
#define params_H_


/**************************
 * Bin and Spectra Parameters
 **************************/
#define N_m_bins 19
#define N_e_bins 11

#define N_e_spectra 7
#define N_m_spectra 2

#define N_dets 3
#define N_anti 2


/**************************
 * Physical Masses
 **************************/

#define MPROTON  0.93827//yj
#define MNEUTRON  0.93957//yj
#define MPION   0.13957
#define MKAON321    0.4936//yj
#define MKAON311    0.49767//yj
#define MKAON   0.493
#define MSIGMA3222  1.18937//yj
#define MSIGMA3122  1.11568//yj
#define MSIGMA3112  1.19744//yj
#define MSIGMA  1.189


/**************************
 *   Fit specific params
 ***************************/

//////////////////////
#define psmear  0.3//yj
#define nsmear 0.4//yj
#define pismear  0.30//yj
#define EMsmear  0.15//yj
#define MUsmear  0.30//yj
#define MUsmear_track  0.05//yj
#define osmear  0.3//yj

#define EMflat  0.02//yj
#define nflat   0.//yj
#define pflat   0.05//yj
#define oflat   0.05//yj

#define p_percent   0.1//yj

//////////////////////



#define n_thresh 0.05//yj
#define p_thresh  0.05//yj
#define pip_thresh  0.1//yj
#define pim_thresh  0.1//yj
#define other_thresh 0.05//yj
#define vertex_thresh  0.1//yj
#define EM_thresh  0.03//yj
#define MU_thresh  0.03//yj

/**************************
 *    Mode selection tools
 ***************************/

#define APP_ONLY 0
#define DIS_ONLY 1
#define BOTH_ONLY 2
#define WIERD_ONLY 3 	// This is combined, but with e_disapearance off.

#define NU_MODE 0
#define NU_NUBAR_MODE 1

/**************************
 *  	Detector Parameters
 **************************/
#define DET_SBND 0
#define DET_UBOONE 1
#define DET_ICARUS 2
#define DET_DUNE 4
#define DET_DUNE_NC 5  //yj

//Metric to imperial ton, tonne, tonnes
#define MET2IMP 1 //0.9071847



#endif
