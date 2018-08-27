
#include "TMatrix.h"
#include "TMatrixT.h"

static Int_t dune_nbins_nue = 44;//11
static Int_t dune_nbins_numu = 44;//19


static Int_t sbn_nbins_nue = 11;//11
static Int_t sbn_nbins_numu = 19;//19

//total 5+3-block histogram length (single detector) is
// fosc_nue, nue intrinsic, mu misID, tau misID, NC mis ID, instinsic numu, taumisID, NC misID_numu
// 11, 11, 11, 11, 11, 11, 11, 19, 19 = 115 bins

/*<subchannel name="fullosc" use="1" osc="21" />	
	  	<subchannel name="fulloscbar" use="1" osc="-21" />	
		<subchannel name="intrinsic" use="1" osc="11"/>
		<subchannel name="intrinsicbar" use="1" osc="-11"/>
		<subchannel name="mumisid" use="1" osc="22" />
		<subchannel name="mumisidbar" use="1" osc="-22" />
		<subchannel name="taumisid" use="1" osc="23"/>
		<subchannel name="taumisidbar" use="1" osc="-23"/>
		<subchannel name="ncmisid" use="1" osc="0"/>
		<subchannel name="ncmisidbar" use="1" osc="0"/>
</channel>
	
<channel name="mulike" use="1" numbins="22">
		<bins 
			edges="0.0 0.5 1.0 1.5 2.0 2.5 3.0 3.5 4.0 4.5 5.0 5.5 6.0 6.5 7.0 7.5 8.0 8.5 9.0 9.5 10.0  15 20"
		/>
		<subchannel name="intrinsic" use="1" osc="22"/>
		<subchannel name="intrinsicbar" use="1" osc="-22"/>
		<subchannel name="taumisid" use="1" osc="23"/>
		<subchannel name="taumisidbar" use="1" osc="-23"/>
		<subchannel name="ncmisid" use="1" osc="0"/>
		<subchannel name="ncmisidbar" use="1" osc="0"/>
</channel>
*/
static Int_t dune_nsubchannels_nue = 5;
static Int_t dune_nsubchannels_numu = 3;
static Int_t dune_ndetectors = 2;


static Int_t sbn_nsubchannels_nue = 7;
static Int_t sbn_nsubchannels_numu = 2;
static Int_t sbn_ndetectors = 3;


static Int_t nmodes =2;

static Int_t dune_nsingledet = dune_nsubchannels_nue*dune_nbins_nue + dune_nsubchannels_numu*dune_nbins_numu; //only neutrino mode
static Int_t dune_ntot = dune_ndetectors*(dune_nsubchannels_nue*dune_nbins_nue + dune_nsubchannels_numu*dune_nbins_numu)*nmodes;

static Int_t sbn_nsingledet = sbn_nsubchannels_nue*sbn_nbins_nue + sbn_nsubchannels_numu*sbn_nbins_numu; //only neutrino mode
static Int_t dune_ntot = sbn_ndetectors*(sbn_nsubchannels_nue*sbn_nbins_nue + sbn_nsubchannels_numu*sbn_nbins_numu)*nmodes;

//three detectors, neutrino and antineutrino mode
//cout << "Generating " << ntot << " by " << ntot << " fractional covariance matrix...\n";

static Int_t dune_pos_fosc = 0;
static Int_t dune_pos_instrinsic_nue = 1;
static Int_t dune_pos_mumisID = 2;
static Int_t dune_pos_taumisID = 3;
static Int_t dune_pos_NCmisID = 4;
static Int_t dune_pos_instrinsic_numu = 5;
static Int_t dune_pos_taumisID_numu = 6;
static Int_t dune_pos_NCmisID_numu = 7;



static Int_t sbn_pos_fosc = 0;
static Int_t sbn_pos_fosc_nubar = 1;
static Int_t sbn_pos_intrinsic_nue = 2;
static Int_t sbn_pos_mumisID = 3;
static Int_t sbn_pos_NCmisID = 4;
static Int_t sbn_pos_dirt = 5;
static Int_t sbn_pos_cosm = 6;
static Int_t sbn_pos_intrinsic_numu = 7;
static Int_t sbn_pos_NCmisID_numu = 8;



int size_increment = 1;

static Bool_t dune_is_nue_fosc(Int_t ii){
    return (((ii%dune_nsingledet)>=(dune_pos_fosc*dune_nbins_nue)) && ((ii%dune_nsingledet)<((dune_pos_fosc+size_increment)*dune_nbins_nue)));
}

static Bool_t dune_is_nue_intrinsic(Int_t ii){
    return (((ii%dune_nsingledet)>=(dune_pos_instrinsic_nue*dune_nbins_nue)) && ((ii%dune_nsingledet)<((dune_pos_instrinsic_nue+size_increment)*dune_nbins_nue)));
}

static Bool_t dune_is_nue_mumisID(Int_t ii){
    return (((ii%dune_nsingledet)>=(dune_pos_mumisID*dune_nbins_nue)) && ((ii%dune_nsingledet)<((dune_pos_mumisID+size_increment)*dune_nbins_nue)));
}

static Bool_t dune_is_nue_taumisID(Int_t ii){
    return (((ii%dune_nsingledet)>=(dune_pos_taumisID*dune_nbins_nue)) && ((ii%dune_nsingledet)<((dune_pos_taumisID+size_increment)*dune_nbins_nue)));
}

static Bool_t dune_is_nue_NCmisID(Int_t ii){
    return (((ii%dune_nsingledet)>=(dune_pos_NCmisID*dune_nbins_nue)) && ((ii%dune_nsingledet)<((dune_pos_NCmisID+size_increment)*dune_nbins_nue)));
}

static Bool_t dune_is_numu_instrinsic(Int_t ii){
    return (((ii%dune_nsingledet)>=(dune_nsubchannels_nue*dune_nbins_nue+(dune_pos_instrinsic_numu-dune_nsubchannels_nue)*dune_nbins_numu )) && ((ii%dune_nsingledet)<((dune_nsubchannels_nue*dune_nbins_nue+(dune_pos_instrinsic_numu-dune_nsubchannels_nue+size_increment)*dune_nbins_numu ))) );
}

static Bool_t dune_is_numu_taumisID(Int_t ii){
    return (((ii%dune_nsingledet)>=(dune_nsubchannels_nue*dune_nbins_nue+(dune_pos_taumisID_numu-dune_nsubchannels_nue)*dune_nbins_numu )) && ((ii%dune_nsingledet)<((dune_nsubchannels_nue*dune_nbins_nue+(dune_pos_taumisID_numu-dune_nsubchannels_nue+size_increment)*dune_nbins_numu ))) );
}

static Bool_t dune_is_numu_NCmisID(Int_t ii){
    return (((ii%dune_nsingledet)>=(dune_nsubchannels_nue*dune_nbins_nue+(dune_pos_NCmisID_numu-dune_nsubchannels_nue)*dune_nbins_numu )) && ((ii%dune_nsingledet)<((dune_nsubchannels_nue*dune_nbins_nue+(dune_pos_NCmisID_numu-dune_nsubchannels_nue+size_increment)*dune_nbins_numu ))) );
}

static Bool_t dune_is_antineutrino(Int_t ii){
	return (ii%dune_nbins_nue>21);
}

void covmx_generation_nu_nubar_xcheck_yeonjae_dune(){
    
    /*Int_t nbins_nue = 11;
     Int_t nbins_numu = 19;
     
     //total 9-block histogram length (single detector) is
     // fosc_nu, fosc_nubar, nueb_intrins, nueb_numu, nueb_NC, nueb_dirt, nueb_cosm, numu_intrins, numu_pipm
     // 11, 11, 11, 11, 11, 11, 11, 19, 19 = 115 bins
     Int_t nsubchannels_nue = 7;
     Int_t nsubchannels_numu = 2;
     Int_t ndetectors = 3;
     Int_t nmodes =2;
     
     Int_t nsingledet = nsubchannels_nue*nbins_nue + nsubchannels_numu*nbins_numu; //only neutrino mode
     Int_t ntot = ndetectors*(nsubchannels_nue*nbins_nue + nsubchannels_numu*nbins_numu)*nmodes; //three detectors, neutrino and antineutrino mode*/
    cout << "Generating " << dune_ntot << " by " << dune_ntot << " fractional covariance matrix...\n";
    
    //Sources of uncertainty (assumed independent/uncorrelated)
    
    
    //TMatrix dirtbkgdrate(dune_ntot,dune_ntot); TH2D *dirtmx = new TH2D("dirtmx","dirtmx",dune_ntot,0,dune_ntot,dune_ntot,0,dune_ntot); TH2D *dirtmxrho = new TH2D("dirtmxrho","dirtmxrho",dune_ntot,0,dune_ntot,dune_ntot,0,dune_ntot);
    //TMatrix cosmbkgdrate(dune_ntot,dune_ntot); TH2D *cosmmx = new TH2D("cosmmx","cosmmx",dune_ntot,0,dune_ntot,dune_ntot,0,dune_ntot); TH2D *cosmmxrho = new TH2D("cosmmxrho","cosmmxrho",dune_ntot,0,dune_ntot,dune_ntot,0,dune_ntot);
    TMatrixT<double> NCbkgdrate(dune_ntot,dune_ntot); TH2D *nuebNCmx = new TH2D("nuebNCmx","nuebNCmx",dune_ntot,0,dune_ntot,dune_ntot,0,dune_ntot); TH2D *nuebNCmxrho = new TH2D("nuebNCmxrho","nuebNCmxrho",dune_ntot,0,dune_ntot,dune_ntot,0,dune_ntot);
    
    TMatrixT<double> detsys(dune_ntot,dune_ntot); TH2D *detmx = new TH2D("detmx","detmx",dune_ntot,0,dune_ntot,dune_ntot,0,dune_ntot); TH2D *detmxrho = new TH2D("detmxrho","detmxrho",dune_ntot,0,dune_ntot,dune_ntot,0,dune_ntot);
    
    TMatrixT<double> xsec(dune_ntot,dune_ntot); TH2D *xsecmx = new TH2D("xsecmx","xsecmx",dune_ntot,0,dune_ntot,dune_ntot,0,dune_ntot); TH2D *xsecmxrho = new TH2D("xsecmxrho","xsecmxrho",dune_ntot,0,dune_ntot,dune_ntot,0,dune_ntot);
    TMatrixT<double> flux(dune_ntot,dune_ntot); TH2D *fluxmx = new TH2D("fluxmx","fluxmx",dune_ntot,0,dune_ntot,dune_ntot,0,dune_ntot); TH2D *fluxmxrho = new TH2D("fluxmxrho","fluxmxrho",dune_ntot,0,dune_ntot,dune_ntot,0,dune_ntot);
    
    TMatrixT<double> total(dune_ntot,dune_ntot); TH2D *totmx = new TH2D("totmx","totmx",dune_ntot,0,dune_ntot,dune_ntot,0,dune_ntot); TH2D *totmxrho = new TH2D("totmxrho","totmxrho",dune_ntot,0,dune_ntot,dune_ntot,0,dune_ntot);
    TMatrixT<double> total_nocorr(dune_ntot,dune_ntot); TH2D *totmx_nocorr = new TH2D("totmx_nocorr","totmx_nocorr",dune_ntot,0,dune_ntot,dune_ntot,0,dune_ntot);
    
    //Assume the following fractional sys. uncertainties and correlations for either neutrino mode or antineutrino mode
    
    //Dirt uncertainties: constrained through in-situ dirt-enhanced sample rate measurement
    //Double_t sigma_dirt = 0.15; //from SBN proposal
    //Double_t rho_dirt = 1.; //among only dirt background block (single detector)
    //assume neutrino and antineutrino are uncorrelated; antineutrino uncertainty is assumed to be the same
    
    //Cosmogenic background uncertainties: constrained through in-situ off-beam high-stats rate measurement
    Double_t sigma_cosm = 0.01; //from SBN proposal; argued as negligible but not given explicitly
    Double_t rho_cosm = 1.; //among only cosmogenic background block (single detector)
    //assume neutrino and antineutrino mode cosmogenic background is fully correlated;
    
    //NC backround uncertainties: constrained through in-situ NCpi0 rate measurement; from table X assuming 50% efficiency for SBND; then scaled by 1/R2 and fid. tonnage for uB and ICARUS
    Double_t sigma_NC_sbnd = 0.0024; //setting statistical uncertainty on NCpi0 sample in each detector as precision; from SBN proposal
    Double_t sigma_NC_ub = 0.0123; //roughly from scaling according to fiducial volume ratio (table xxxv, for nue analysis) and flux ratio
    Double_t sigma_NC_icarus = 0.0061; //roughly from scaling according to fiducial volume ratio (table xxxv, for nue analysis) and flux ratio
    Double_t rho_NC = 1.; //among only NC block (single detector)
    
    Double_t sigma_NC_ND = 0.0123;//yeonjae, arbitrary value
    Double_t sigma_NC_FD = 0.0024;//yeonjae, arbitrary value
    
    //assumed neutrino and atnineutrino are uncorrelated
    //antineutrino uncertainties are scaled according to statistics
    Double_t sigma_NC_sbnd_nubar = 0.0044;
    Double_t sigma_NC_ub_nubar = 0.023;
    Double_t sigma_NC_icarus_nubar = 0.011;
    
    Double_t sigma_NC_ND_nubar = 0.0123;//yeonjae, arbitrary value
    Double_t sigma_NC_FD_nubar = 0.0024;//yeonjae, arbitrary value
    
    //Detector uncertainties: assumed 2-3%, single detector
    Double_t sigma_det = 0.025; //from SBN proposal
    Double_t rho_det = 0.; //from SBN proposal
    //assumed neutrino and antineutrino fully correlated
    
    //Flux uncertainties:
    Double_t sigma_flux_nue = 0.153; //table III of SBN proposal, contributions added in quadrature
    Double_t sigma_flux_numu = 0.152;//0.152; //table III of SBN proposal, contributions added in quadrature

 
    Double_t rho_flux_numu_numu = 1.; //between fullosc, numu, numu_pipm, and nueb_numu
    Double_t rho_flux_nue_numu = 0.6; //between nue_intrins and the above
    Double_t rho_flux_nue_nue = 1.; //
    //correlations also hold from detector to detector
    Double_t rho_flux_nu_nubar = 0.4;
    //assume neutrino and antineutrino mode fluxes 40% correlated (because there's 40% wrong sign in nubar mode)
  
    //https://arxiv.org/pdf/1606.09550.pdf 
    //The νe and νe signal modes have independent normalization uncertainties of 2% each, while the νµ and νµ signal modes have independent normalization uncertainties of 5%
 
    //Cross-section uncertainties: assumed 20%
    Double_t sigma_CCxsec = 0.2;
    Double_t rho_CCxsec = 1.;

    Double_t sigma_NCxsec = 0.3;
    Double_t rho_NCxsec = 1.;
    Double_t rho_CC_NC = 0.5;
    
    //consider neutrino and antineutrino correlations
    //only for cross section; fully correlated CC and NC
    Double_t rho_CC_nu_nubar = 1.;
    Double_t rho_NC_nu_nubar = 1.;
    
    int i,j;
    
    Int_t onemode = dune_ntot/nmodes;
    
    Int_t flag_nu_nubar_correl = 1;
    
    
    
    
    for (i=0; i<dune_ntot; i++){
        for (j=0; j<dune_ntot; j++){
            
            if ((i<onemode && j<onemode) || (i>=onemode && j>=onemode) ){ //neutrino mode or antineutrino mode blocks
                
                //dirtbkgdrate(i,j) = 0.;
                //cosmbkgdrate(i,j) = 0.;
                NCbkgdrate(i,j) = 0.;
                detsys(i,j) = 0.;
                xsec(i,j) = 0.;
                flux(i,j) = 0.;
                total(i,j) = 0.;
                total_nocorr(i,j) = 0.;
                                //nueb_NC uncertainties -- applied to NC misID nue bkgd histogram only
                //neutrino mode:
              
                if (i<onemode && j<onemode){
                    if ( is_nue_NCmisID(i) && is_nue_NCmisID(j) && ((i/dune_nsingledet)==(j/dune_nsingledet)) ){
                        
                        if (((i%onemode)/dune_nsingledet)==0){//ND
                            if (i==j) NCbkgdrate(i,j) = sigma_NC_ND*sigma_NC_ND;
                            else NCbkgdrate(i,j) = rho_NC*sigma_NC_ND*sigma_NC_ND;
                        }
                        if (((i%onemode)/dune_nsingledet)==1){//FD
                            if (i==j) NCbkgdrate(i,j) = sigma_NC_FD*sigma_NC_FD;
                            else NCbkgdrate(i,j) = rho_NC*sigma_NC_FD*sigma_NC_FD;
                        }
                        
                    }
                }
                //antineutrino mode
                else if (i>=onemode && j>=onemode){
                    if ( is_nue_NCmisID(i) && is_nue_NCmisID(j) && ((i/dune_nsingledet)==(j/dune_nsingledet)) ){
                        if (((i%onemode)/dune_nsingledet)==0){//ND
                            if (i==j) NCbkgdrate(i,j) = sigma_NC_ND_nubar*sigma_NC_ND_nubar;
                            else NCbkgdrate(i,j) = rho_NC*sigma_NC_ND_nubar*sigma_NC_ND_nubar;
                        }
                        if (((i%onemode)/dune_nsingledet)==1){//FD
                            if (i==j) NCbkgdrate(i,j) = sigma_NC_FD_nubar*sigma_NC_FD_nubar;
                            else NCbkgdrate(i,j) = rho_NC*sigma_NC_FD_nubar*sigma_NC_FD_nubar;
                        }
                    }
                }
                
                //det uncertainties -- applied to all detectors and samples, but not the in-situ-constrained distributions, no correlations
                if ( i==j ){
                    detsys(i,j) = sigma_det*sigma_det;
                }
                
                
                
                //flux uncertainties -- applied to all detectors and samples, except the in-situ constrained distributions
                
                //yeonjae in-situ: measure the background directily to constraint the uncertatinty
                
                
                if ( !is_nue_intrinsic(i)  && !is_nue_intrinsic(j)  ) {
                    // numu = all except intrinsic nue and in-situ constrained mis-ID bins
		     

                    if (i==j){
			flux(i,j) = sigma_flux_numu*sigma_flux_numu;
		    }
                    else{
			flux(i,j) = rho_flux_numu_numu*sigma_flux_numu*sigma_flux_numu;
			}
                }
                
                if ( is_nue_intrinsic(i) && is_nue_intrinsic(j) ){ // intrinsic nue only
                    if (i==j) flux(i,j) = sigma_flux_nue*sigma_flux_nue;
                    else flux(i,j) = rho_flux_nue_nue*sigma_flux_nue*sigma_flux_nue;
                }
                
                else if ( is_nue_intrinsic(i) && !is_nue_intrinsic(j) ){ //intrinsic nue correlations with numu
                    flux(i,j) = rho_flux_nue_numu*sigma_flux_nue*sigma_flux_numu;
                }
                
                else if ( is_nue_intrinsic(j)  ){ // itrinsic nue correlations with numu //yeonjae
                    flux(i,j) = rho_flux_nue_numu*sigma_flux_numu*sigma_flux_nue;
                }
                
                //xsec uncertainties -- applied to all detectors and samples, except the in-situ constrained distributions
                if ( true ){ // exclude in-situ constrained backgrounds
                    if ( !is_numu_NCmisID(i) && !is_numu_NCmisID(j) ){ // all charged current here
                        if (i==j) xsec(i,j) = sigma_CCxsec*sigma_CCxsec;
                        else xsec(i,j) = rho_CCxsec*sigma_CCxsec*sigma_CCxsec;
                    }
                    if ( is_numu_NCmisID(i) && is_numu_NCmisID(j) ){ // all neutral current here
                        if (i==j) xsec(i,j) = sigma_NCxsec*sigma_NCxsec;
                        else xsec(i,j) = rho_NCxsec*sigma_NCxsec*sigma_NCxsec;
                    }
                    if ( is_numu_NCmisID(i) && !is_numu_NCmisID(j) ){
                        xsec(i,j) = rho_CC_NC*sigma_CCxsec*sigma_NCxsec;
                    }
                    if ( !is_numu_NCmisID(i) && is_numu_NCmisID(j) ){
                        xsec(i,j) = rho_CC_NC*sigma_NCxsec*sigma_CCxsec;
                    }
                }
                
            }//end if i && j in neutrino mode block or antineutrino mode block
            
            if (flag_nu_nubar_correl==1){
                
                if ((i<onemode && j>=onemode) || (i>=onemode && j<onemode) ){ //neutrino mode or antineutrino mode blocks
                    
                    //dirtbkgdrate(i,j) = 0.;
                    //cosmbkgdrate(i,j) = 0.;
                    NCbkgdrate(i,j) = 0.;
                    detsys(i,j) = 0.;
                    xsec(i,j) = 0.;
                    flux(i,j) = 0.;
                    total(i,j) = 0.;
                    total_nocorr(i,j) = 0.;
                    
                    // fosc_nu, fosc_nubar, nueb_intrins, nueb_numu, nueb_NC, nueb_dirt, nueb_cosm, numu_intrins, numu_pipm
                    // 11, 11, 11, 11, 11, 11, 11, 19, 19 = 115 bins
                    
                    //dirt uncertainties -- applied to dirt histogram only
                    //assumed uncorrelated in neutrino and antineutrino
                    
                    //cosmogenic uncertainties -- applied to cosmogenic bkgd histogram only
                    //assumed fully correlated in neutrino and antineutrino
                    
                    /*
                    if ( is_nue_cosmogenic(i%onemode) && is_nue_cosmogenic(j%onemode) && (((i%onemode)/dune_nsingledet)==((j%onemode)/dune_nsingledet)) ){
                        //yeonjae cosmogenic and inthe same detector
                        if (i==(j+onemode) || j==(i+onemode)) cosmbkgdrate(i,j) = sigma_cosm*sigma_cosm;
                        else cosmbkgdrate(i,j) = rho_cosm*sigma_cosm*sigma_cosm;
                    }
                    */
                     
                    //nueb_NC uncertainties -- applied to NC misID nue bkgd histogram only
                    //assumed uncorrelated in neutrino and antineutrino
                    
                    //det uncertainties -- applied to all detectors and samples, but not the in-situ-constrained distributions, no correlations
                    if ( i==(j+onemode) || j==(i+onemode) ){
                        detsys(i,j) = sigma_det*sigma_det;
                    }
                    
                    //flux uncertainties -- applied to all detectors and samples, except the in-situ constrained distributions
                    //assumed fully correlated only between fullosc nu in nu and nubar mode, and between fullosc nubar in nu and nubar mode
                    //      if ( (((i%dune_nsingledet)<(2*dune_nbins_nue)) && ((j%dune_nsingledet)<(2*dune_nbins_nue))) ){//remove this if statement if we want fully correlated in nu and nubar
                    if ( !is_nue_intrinsic(i) &&
                        !is_nue_intrinsic(j) ){ // numu = all except intrinsic nue and in-situ constrained mis-ID bins
                        if (i==j) flux(i,j) = sigma_flux_numu*sigma_flux_numu;
                        else flux(i,j) = rho_flux_numu_numu*rho_flux_nu_nubar*sigma_flux_numu*sigma_flux_numu;
                    }
                    if ( is_nue_intrinsic(i) && is_nue_intrinsic(j) ){ // intrinsic nue only
                        if (i==j) flux(i,j) = sigma_flux_nue*sigma_flux_nue;
                        else flux(i,j) = rho_flux_nue_nue*rho_flux_nu_nubar*sigma_flux_nue*sigma_flux_nue;
                    }
                    else if ( is_nue_intrinsic(i) && !is_nue_intrinsic(j) ){ //intrinsic nue correlations with numu
                        flux(i,j) = rho_flux_nue_numu*rho_flux_nu_nubar*sigma_flux_nue*sigma_flux_numu;
                    }
                    else if ( is_nue_intrinsic(j) ){ // itrinsic nue correlations with numu
                        flux(i,j) = rho_flux_nue_numu*rho_flux_nu_nubar*sigma_flux_numu*sigma_flux_nue;
                        //      }
                    }
                    
                    //xsec uncertainties -- applied to all detectors and samples, except the in-situ constrained distributions
                    if ( true ){ // exclude in-situ constrained backgrounds
                        if (!is_numu_NCmisID(i) && !is_numu_NCmisID(j)){ // all charged current here
                            if (i==j) xsec(i,j) = sigma_CCxsec*sigma_CCxsec;
                            else xsec(i,j) = rho_CCxsec*sigma_CCxsec*sigma_CCxsec;
                        }
                        if (is_numu_NCmisID(i) && is_numu_NCmisID(j)){ // all neutral current here
                            if (i==j) xsec(i,j) = sigma_NCxsec*sigma_NCxsec;
                            else xsec(i,j) = rho_NCxsec*sigma_NCxsec*sigma_NCxsec;
                        }
                        if (is_numu_NCmisID(i) && !is_numu_NCmisID(j)){
                            xsec(i,j) = rho_CC_NC*sigma_CCxsec*sigma_NCxsec;
                        }
                        if (!is_numu_NCmisID(i) && is_numu_NCmisID(j)){
                            xsec(i,j) = rho_CC_NC*sigma_NCxsec*sigma_CCxsec;
                        }
                    }
                 
                    
                
                 }//end if i && j in neutrino-antineutrino correlation block
        
            }//end if correlation flag on
            
        }//end loop over j
    }//end loop over i

    //Here we will fix the CS due to fact the upper half of all bins is antineutrinos
// these numbers from page 68 on https://arxiv.org/pdf/0806.1449.pdf
   double flux_corr = 0.5;
   double flux_sigma = 0.25/0.15;

   double xsec_corr = 0.5;
   double xsec_sigma = 1.5;


    for(int i=0;i<dune_ntot; i++){
    for(int j=0;j<dune_ntot; j++){
		if(i==j && is_nue_fosc(i)){ flux(i,j) = flux(i,j)*1.19;
		}else if(i==j){
			flux(i,j) += 0.012;
		}
		
		if( is_antineutrino(i) && is_antineutrino(j)){
			flux(i,j) = flux(i,j)*flux_sigma*flux_sigma; 
			xsec(i,j) = xsec(i,j)*xsec_sigma*xsec_sigma; 
		}
		if( (is_antineutrino(i) && !is_antineutrino(j)) || (is_antineutrino(j) && !is_antineutrino(i))){
			flux(i,j) = flux(i,j)*flux_sigma*flux_corr;
			xsec(i,j) = xsec(i,j)*xsec_sigma*xsec_corr;
		}

    }
    }




    total =  NCbkgdrate + detsys + flux + xsec; //dirtbkgdrate + cosmbkgdrate +
   


 
    for (i=0; i<dune_ntot; i++){
        for (j=0; j<dune_ntot; j++){
            //dirtmx ->SetBinContent(i+1,j+1,dirtbkgdrate(i,j)); dirtmxrho->SetBinContent(i+1,j+1,dirtbkgdrate(i,j)/sqrt(dirtbkgdrate(i,i)*dirtbkgdrate(j,j)));
            //cosmmx ->SetBinContent(i+1,j+1,cosmbkgdrate(i,j)); cosmmxrho->SetBinContent(i+1,j+1,cosmbkgdrate(i,j)/sqrt(cosmbkgdrate(i,i)*cosmbkgdrate(j,j)));
            nuebNCmx ->SetBinContent(i+1,j+1,NCbkgdrate(i,j)); nuebNCmxrho->SetBinContent(i+1,j+1,NCbkgdrate(i,j)/sqrt(NCbkgdrate(i,i)*NCbkgdrate(j,j)));
            detmx ->SetBinContent(i+1,j+1,detsys(i,j));  detmxrho ->SetBinContent(i+1,j+1,detsys(i,j)/sqrt(detsys(i,i)*detsys(j,j)));
            fluxmx ->SetBinContent(i+1,j+1,flux(i,j)); fluxmxrho ->SetBinContent(i+1,j+1,flux(i,j)/sqrt(flux(i,i)*flux(j,j)));
            xsecmx ->SetBinContent(i+1,j+1,xsec(i,j));    xsecmxrho ->SetBinContent(i+1,j+1,xsec(i,j)/sqrt(xsec(i,i)*xsec(j,j)));
            totmx ->SetBinContent(i+1,j+1,NCbkgdrate(i,j)+detsys(i,j)+xsec(i,j)+flux(i,j));
            totmxrho ->SetBinContent(i+1,j+1,(NCbkgdrate(i,j)+detsys(i,j)+xsec(i,j)+flux(i,j))/sqrt((NCbkgdrate(i,i)+detsys(i,i)+xsec(i,i)+flux(i,i))*(NCbkgdrate(j,j)+detsys(j,j)+xsec(j,j)+flux(j,j))));
            total(i,j) = totmx->GetBinContent(i+1,j+1);
            if (i==j){
                totmx_nocorr ->SetBinContent(i+1,j+1,NCbkgdrate(i,j)+detsys(i,j)+xsec(i,j)+flux(i,j));
                total_nocorr(i,j) = totmx_nocorr->GetBinContent(i+1,j+1);
            }
        }
    }
    
    cout << "total determinant: " << total.Determinant() << " , " << "total_nocorr determinant: "<< total_nocorr.Determinant() << endl;
    
    gStyle->SetPadRightMargin(0.2);
    gStyle->SetOptStat(0);
    gStyle->SetGridStyle(1);
    
    TCanvas *c = new TCanvas("c","c",800,800);
    c->Divide(3,3);
    c->cd(1); gPad->SetGrid();
    //dirtmx->Draw("colz");
    c->cd(2); gPad->SetGrid();
    //cosmmx->Draw("colz");
    c->cd(3);
    nuebNCmx->Draw("colz");
    c->cd(4);
    detmx->Draw("colz");
    c->cd(5);
    xsecmx->Draw("colz");
    c->cd(6);
    fluxmx->Draw("colz");
    c->cd(8);
    totmx_nocorr->Draw("colz");
    c->cd(9);
    totmx->Draw("colz");
    
    
    TCanvas *d = new TCanvas("d","d",800,800);
    d->Divide(3,3);
    d->cd(1); gPad->SetGrid();
    //dirtmxrho->Draw("colz");
    d->cd(2); gPad->SetGrid();
    //cosmmxrho->Draw("colz");
    d->cd(3);
    nuebNCmxrho->Draw("colz");
    d->cd(4);
    detmxrho->Draw("colz");
    d->cd(5);
    xsecmxrho->Draw("colz");
    d->cd(6);
    fluxmxrho->Draw("colz");
    d->cd(9);
    totmxrho->Draw("colz");
    
    
    if (dune_nbins_nue<3 && dune_nbins_numu<4){
        
        /*printf("\n========================\n");
        printf("       Dirt Covmx");
        printf("\n========================\n");
        printf("\t");
        for (i=0; i<dune_ntot; i++){
            printf("%i\t",i%dune_nsingledet);
        }
        printf("\n");
        for (i=0; i<dune_ntot; i++){
            printf("%i\t",i%dune_nsingledet);
            for (j=0; j<dune_ntot; j++){
                printf("%5.1e\t",dirtbkgdrate(i,j));
            }
            printf("\n");
        }*/
        
        
        /*printf("\n========================\n");
        printf("       Cosm Covmx");
        printf("\n========================\n");
        printf("\t");
        for (i=0; i<dune_ntot; i++){
            printf("%i\t",i%dune_nsingledet);
        }
        printf("\n");
        for (i=0; i<dune_ntot; i++){
            printf("%i\t",i%dune_nsingledet);
            for (j=0; j<dune_ntot; j++){
                printf("%5.1e\t",cosmbkgdrate(i,j));
            }
            printf("\n");
        }*/
        
        printf("\n========================\n");
        printf("       nueb_NC Covmx");
        printf("\n========================\n");
        printf("\t");
        for (i=0; i<dune_ntot; i++){
            printf("%i\t",i%dune_nsingledet);
        }
        printf("\n");
        for (i=0; i<dune_ntot; i++){
            printf("%i\t",i%dune_nsingledet);
            for (j=0; j<dune_ntot; j++){
                printf("%5.1e\t",NCbkgdrate(i,j));
            }
            printf("\n");
        }
        
        printf("\n========================\n");
        printf("       Det Covmx");
        printf("\n========================\n");
        printf("\t");
        for (i=0; i<dune_ntot; i++){
            printf("%i\t",i%dune_ntot);
        }
        printf("\n");
        for (i=0; i<dune_ntot; i++){
            printf("%i\t",i%dune_nsingledet);
            for (j=0; j<dune_ntot; j++){
                printf("%5.1e\t",detsys(i,j));
            }
            printf("\n");
        }
        
        printf("\n========================\n");
        printf("       Flux Covmx");
        printf("\n========================\n");
        printf("\t");
        for (i=0; i<dune_ntot; i++){
            printf("%i\t",i%dune_ntot);
        }
        printf("\n");
        for (i=0; i<dune_ntot; i++){
            printf("%i\t",i%dune_nsingledet);
            for (j=0; j<dune_ntot; j++){
                printf("%5.1e\t",flux(i,j));
            }
            printf("\n");
        }
        
        
        printf("\n========================\n");
        printf("       Xsec Covmx");
        printf("\n========================\n");
        printf("\t");
        for (i=0; i<dune_ntot; i++){
            printf("%i\t",i%dune_ntot);
        }
        printf("\n");
        for (i=0; i<dune_ntot; i++){
            printf("%i\t",i%dune_nsingledet);
            for (j=0; j<dune_ntot; j++){
                printf("%5.1e\t",xsec(i,j));
            }
            printf("\n");
        }
        
    }
    
    char name[500];
    sprintf(name,"covariance_matrices_xcheck_%ix%i.root",dune_ntot,dune_ntot);
    TFile *f = new TFile(name,"RECREATE");
    
    //dirtmx->Write();
    //cosmmx->Write();
    nuebNCmx->Write();
    detmx->Write();
    fluxmx->Write();
    xsecmx->Write();
    totmx->Write();
    totmx_nocorr->Write();
    
    //dirtmxrho->Write();
    //cosmmxrho->Write();
    nuebNCmxrho->Write();
    detmxrho->Write();
    fluxmxrho->Write();
    xsecmxrho->Write();
    totmxrho->Write();
    
    //dirtbkgdrate.Write();
    //cosmbkgdrate.Write();
//    NCbkgdrate.Write();
//    detsys.Write();
//    flux.Write();
//    xsec.Write();
    total.Write();
//    total_nocorr.Write();
    
    f->Close();
     
    
    
}
