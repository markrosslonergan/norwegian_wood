//definitions of functions here.







int DecayIt(TLorentzVector parent, TRandom3 *rangen, TLorentzVector daughters[3]){
    // this version neglects muon polarization, and electron mass
    // assumes the pure V-A coupling
    // the Neutrinos are correclty V-A.
    
    // get parent mass
    double parentmass = parent.M();
    int N_DAUGHTER = 3;
    
    // sum of daughter particles' mass
    double daughtermass[N_DAUGHTER];
    double sumofdaughtermass = 0.0;
    for (int index=0; index<N_DAUGHTER; index++){
        daughtermass[index] = daughters[index].M();
        sumofdaughtermass+=daughtermass[index];
    }
    
    // parent particle is at rest
    // calcualte daughter momemtum
    double daughtermomentum[N_DAUGHTER];
    // calcualte electron energy
    double xmax = (1.0+daughtermass[0]*daughtermass[0]/parentmass/parentmass);
    double x;
    
    double Ee, Ene;
    
    double gam;
    double EMax=parentmass/2-daughtermass[0];
    
    int MAX_LOOP =1000;
    // generating random energy
    for (int loop1=0; loop1<MAX_LOOP; loop1++){
        Ee=rangen->Uniform(1.);// check this line
        for(int loop2=0; loop2<MAX_LOOP;loop2++){
            x=xmax*rangen->Uniform(1.);
            gam = rangen->Uniform(1.);
            if (gam <= x*(1.-x)) break;
            x = xmax;
        }
        Ene = x;
        if (Ene >= (1.-Ee)) break;
        Ene = 1.-Ee;
    }
    double Enm=(2.-Ee-Ene);
    
    //initialisation of rotation parameters
    
    double costheta, sintheta, rphi, rtheta, rpsi;
    costheta = 1.-2./Ee-2./Ene+2./Ene/Ee;
    sintheta= sqrt(1.-costheta*costheta);
    
    rphi=TMath::TwoPi()*rangen->Uniform(1.);
    rtheta=TMath::ACos(2.*rangen->Uniform(1.)-1.);
    rpsi=TMath::TwoPi()*rangen->Uniform(1.);
    
    TRotation rot;
    rot.SetXEulerAngles(rphi,rtheta,rpsi);
    
    //electron 0
    daughtermomentum[0]=sqrt(Ee*Ee*EMax*EMax+2.0*Ee*EMax*daughtermass[0]);
    TVector3 direction0(0.0,0.0,1.0);
    
    direction0*=daughtermomentum[0];
    
    direction0 *= rot;
    daughters[0].SetPxPyPzE(direction0.X(),direction0.Y(),direction0.Z(),sqrt(daughtermomentum[0]*daughtermomentum[0]+daughtermass[0]*daughtermass[0]));
    
    //electron neutrino
    daughtermomentum[1]=sqrt(Ene*Ene*EMax*EMax+2.0*Ene*EMax * daughtermass[1]);
    TVector3 direction1(sintheta,0.0,costheta);
    
    direction1*=daughtermomentum[1];
    
    direction1 *= rot;
    
    daughters[1].SetPxPyPzE(direction1.X(),direction1.Y(),direction1.Z(),sqrt(daughtermomentum[1]*daughtermomentum[1]+daughtermass[1]*daughtermass[1]));
    
    //tau neutrino
    daughtermomentum[2] = sqrt(Enm*Enm*EMax*EMax +2.0*Enm*EMax*daughtermass[2]);
    TVector3 direction2(-Ene/Enm*sintheta,0,-Ee/Enm-Ene/Enm*costheta);
    
    direction2*=daughtermomentum[2];
    
    direction2 *= rot;
    
    daughters[2].SetPxPyPzE(direction2.X(),direction2.Y(),direction2.Z(),sqrt(daughtermomentum[2]*daughtermomentum[2]+daughtermass[2]*daughtermass[2]));
    
    return 0;
    //double EMASS =
    
    
}




double R_c(double x,double omega){
    
    int n_max = (int)(100.*x);
    
    if(n_max<10)n_max=10;
    
    double L2 = 0.0;
    
    for(int n=1; n<=n_max; n++){
        L2 += pow(x,n)/(n*n);
    }
    
    double r_c;
    
    r_c = 2.*L2-(TMath::Pi()*TMath::Pi()/3.)-2.;
    r_c = r_c + omega * (1.5+2.*log((1.-x)/x));
    r_c = r_c - log(x)*(2.*log(x)-1.);
    r_c = r_c + (3.*log(x)-1.-1./x)*log(1.-x);
    
    return r_c;
}




double F_c(double x, double x0 , double omega)
{
    
    double f_c;
    
    f_c = (5.+17.*x-34.*x*x)*(omega+log(x))-22.*x+34.*x*x;
    f_c = (1.-x)/(3.*x*x)*f_c;
    f_c = (6.-4.*x)*R_c(x,omega)+(6.-6.*x)*std::log(x) + f_c;
    f_c = ((1./137.)/TMath::TwoPi()) * (x*x-x0*x0) * f_c;
    
    return f_c;
}

double F_theta(double x, double x0,double omega)
{
    double f_theta;
    
    f_theta = (1.+x+34*x*x)*(omega+TMath::Log(x))+3.-7.*x-32.*x*x;
    f_theta = f_theta + ((4.*(1.-x)*(1.-x))/x)*TMath::Log(1.-x);
    f_theta = (1.-x)/(3.*x*x) * f_theta;
    f_theta = (2.-4.*x)*R_c(x,omega)+(2.-6.*x)*TMath::Log(x)-f_theta;
    
    
    /*static constexpr double Avogadro = 6.02214179e+23/mole;
     
     //
     // c   = 299.792458 mm/ns
     // c^2 = 898.7404 (mm/ns)^2
     //
     static constexpr double c_light   = 2.99792458e+8 * m/s;
     static constexpr double c_squared = c_light * c_light;
     
     //
     // h     = 4.13566e-12 MeV*ns
     // hbar  = 6.58212e-13 MeV*ns
     // hbarc = 197.32705e-12 MeV*mm
     //
     static constexpr double c_light   = 2.99792458e+8 * m/s;
     static constexpr double c_squared = c_light * c_light;
     
     static constexpr double h_Planck      = 6.62606896e-34 * joule*s;
     static constexpr double hbar_Planck   = h_Planck/twopi;
     static constexpr double hbarc         = hbar_Planck * c_light;
     static constexpr double hbarc_squared = hbarc * hbarc;
     
     //
     //
     //
     static constexpr double eplus = 1. ;// positron charge
     static constexpr double electron_charge = - eplus; // see SystemOfUnits.h
     static constexpr double e_squared = eplus * eplus;
     static constexpr double epsilon0 = 1./(c_squared*mu0);
     static constexpr double elm_coupling           = e_squared/(4*pi*epsilon0);
     static constexpr double fine_structure_const   = elm_coupling/hbarc;
     */
    double fine_structure_const = 1./137.;
    
    f_theta = (fine_structure_const/TMath::TwoPi()) * (x*x-x0*x0) * f_theta;
    
    return f_theta;
}



int DecayItWithSpin(TLorentzVector parent, TRandom3 *rangen, TLorentzVector daughters[3], TVector3 parent_polarization){
    // this version neglects muon polarization, and electron mass
    // assumes the pure V-A coupling
    // the Neutrinos are correclty V-A.
    
    // get parent mass
    double parentmass = parent.M();
    
    int N_DAUGHTER = 3;
    
    // sum of daughter particles' mass
    double daughtermass[N_DAUGHTER];
    double sumofdaughtermass = 0.0;
    for (int index=0; index<N_DAUGHTER; index++){
        daughtermass[index] = daughters[index].M();
        sumofdaughtermass+=daughtermass[index];
    }
    
    double EMASS = 0.00051;
    
    double michel_rho = 0.75;
    double michel_delta = 0.75;
    double michel_xsi = 1.00;
    double michel_eta = 0.00;
    
    double rndm, x, ctheta;
    
    double FG;
    double FG_max = 2.00;
    
    double W_mue = (parentmass*parentmass+EMASS*EMASS)/(2.*parentmass);
    double x0 = EMASS/W_mue;
    
    double x0_squared = x0*x0;
    
    // ***********************************************
    //    x0 <= x <= 1. and -1 <= y <= 1
    //
    //    F(x,y) = f(x)*g(x,y); g(x,y) = 1.+g(x)*y
    // ***********************************************
    
    // ***** sampling F(x,y) directly (brute force) *****
    
    int MAX_LOOP=10000;
    for (int loop_count =0; loop_count<MAX_LOOP; loop_count++){
        rndm = rangen->Uniform(1.);
        x = x0 + rndm*(1.-x0);
        double x_squared = x*x;
        
        double F_IS, F_AS, G_IS, G_AS;
        
        F_IS = 1./6.*(-2.*x_squared+3.*x-x0_squared);
        F_AS = 1./6.*sqrt(x_squared-x0_squared)*(2.*x-2.+sqrt(1.-x0_squared));
        
        G_IS = 2./9.*(michel_rho-0.75)*(4.*x_squared-3.*x-x0_squared);
        G_IS = G_IS + michel_eta*(1.-x)*x0;
        
        G_AS = 3.*(michel_xsi-1.)*(1.-x);
        G_AS = G_AS+2.*(michel_xsi*michel_delta-0.75)*(4.*x-4.+sqrt(1.-x0_squared));
        G_AS = 1./9.*sqrt(x_squared-x0_squared)*G_AS;
        
        F_IS = F_IS + G_IS;
        F_AS = F_AS + G_AS;
        
        
        // *** Radiative Corrections ***
        double omega = log(parentmass/EMASS);
        double R_IS = F_c(x,x0,omega);
        
        double F = 6.*F_IS + R_IS/sqrt(x_squared-x0_squared);
        
        // *** Radiative Corrections ***
        double R_AS = F_theta(x,x0,omega);
        rndm = rangen->Uniform(1.);
        ctheta = 2.*rndm-1.;
        double G = 6.*F_AS - R_AS/sqrt(x_squared-x0_squared);
        
        FG = sqrt(x_squared-x0_squared)*F*(1.+(G/F)*ctheta);
        
        if(FG>FG_max){
            cout << "***Problem in Tau Decay *** : FG > FG_max"<<endl;
            FG_max =FG;
        }
        
        rndm = rangen->Uniform(1.);
        
        if (FG >= rndm*FG_max) break;
    }
    
    double energy = x*W_mue;
    rndm = rangen->Uniform(1.);
    
    
    double phi = TMath::TwoPi()*rndm;
    
    if(energy < EMASS) energy = EMASS;
    
    double daughtermomentum[3];
    
    daughtermomentum[0] = sqrt(energy*energy - EMASS*EMASS);
    
    double stheta = sqrt(1.-ctheta*ctheta);
    double cphi = TMath::Cos(phi);
    double sphi = TMath::Sin(phi);
    
    //Coordinates of the decay position with respect to the muon spin
    
    double px = stheta*cphi;
    double py = stheta*sphi;
    double pz = ctheta;
    
    TVector3 direction0(px, py, pz);
    
    direction0.RotateUz(parent_polarization);
    
    direction0 *= daughtermomentum[0];
    
    daughters[0].SetPxPyPzE(direction0.X(),direction0.Y(),direction0.Z(),sqrt(daughtermomentum[0]*daughtermomentum[0]+daughtermass[0]*daughtermass[0]));
    
    // daughter 1, 2 (neutrinos)
    // create neutrinos in the CM frame of two neutrinos
    double energy2 = parentmass*(1.0 - x/2.0);
    double vmass   = sqrt((energy2-daughtermomentum[0])*(energy2+daughtermomentum[0]));
    double beta = -1.0*daughtermomentum[0]/energy2;
    double costhetan = 2.*rangen->Uniform(1.)-1.0;
    double sinthetan = sqrt((1.0-costhetan)*(1.0+costhetan));
    double phin  = TMath::TwoPi()*rangen->Uniform(1.);
    double sinphin = TMath::Sin(phin);
    double cosphin = TMath::Cos(phin);
    
    TVector3 direction1 (sinthetan*cosphin,sinthetan*sinphin,costhetan);
    
    direction1 *= vmass/2.;
    
    daughters[1].SetPxPyPzE(direction1.X(),direction1.Y(),direction1.Z(),vmass/2.);
    daughters[2].SetPxPyPzE(-1.0*direction1.X(),-1.0*direction1.Y(),-1.0*direction1.Z(),vmass/2.);
    
    // boost to the muon rest frame
    daughters[1].Boost( direction0.x()*beta, direction0.y()*beta, direction0.z()*beta);
    daughters[2].Boost( direction0.x()*beta, direction0.y()*beta, direction0.z()*beta);
    
    return 0;
}




int threebodydecay(TRandom3 *rangen, double mP, double m0, double m1, double m2, double p0[4], double p1[4], double p2[4]){
    
    //mP: mass of parent particle
    
    double mass_sum_daughters =m0+m1+m2;
    
    double mommax=0.0;
    double momsum=0.0;
    double rd1=0;
    double rd = 0;
    double rd2 = 0;
    double energy = 0;
    
    double absP0 = 0;
    double absP1 = 0;
    double absP2 = 0;
    
    do{
        // TRandom3->Uniform(0,1);
        rd1 = rangen->Uniform(0,1);
        rd2 = rangen->Uniform(0,1);
        if(rd2>rd1){
            rd=rd1;rd1=rd2;rd2=rd;
        };
        mommax=0.0;
        momsum=0.0;
        
        //Daughter 0
        energy = rd2*(mP - mass_sum_daughters);
        absP0 = sqrt(energy*energy + 2.0*energy*m0);
        if(absP0 > mommax) mommax = absP0;
        momsum = momsum + absP0;
        
        // Daughter 1
        energy = (1.0 - rd1)*(mP - mass_sum_daughters);
        absP1 = sqrt(energy*energy + 2.0*energy*m1);
        if(absP1 > mommax) mommax = absP1;
        momsum = momsum + absP1;
        
        // Daughter 2
        energy = (rd1 - rd2)*(mP - mass_sum_daughters);
        absP2 = sqrt(energy*energy + 2.0*energy*m2);
        if(absP2 > mommax) mommax = absP2;
        momsum = momsum + absP2;
        
        
        
    } while (mommax > momsum-mommax);
    
    double costheta, sintheta, phi,sinphi,cosphi;
    double costhetan, sinthetan,phin,sinphin,cosphin;
    
    costheta=2.0*rangen->Uniform(0,1)-1.0;
    sintheta= sqrt((1-costheta)*(1+costheta));
    phi=2*3.14159*rangen->Uniform(0,1);
    sinphi=sin(phi);
    cosphi=cos(phi);
    
    std::vector<double > tp0 = {sqrt(m0*m0+absP0*absP0),absP0*sintheta*cosphi, absP0*sintheta*sinphi, absP0*costheta};
    
    p0[1]= absP0*sintheta*cosphi;
    p0[2]= absP0*sintheta*sinphi;
    p0[3]= absP0*costheta;
    
    
    costhetan = (absP1*absP1-absP2*absP2-absP0*absP0)/(2.0*absP2*absP0);
    sinthetan = sqrt((1-costhetan)*(1+costhetan));
    phin = 2*3.14159*rangen->Uniform(0,1);
    sinphin=sin(phin);
    cosphin=cos(phin);
    
    std::vector<double > d2 = {sinthetan*cosphin*costheta*cosphi - sinthetan*sinphin*sinphi +   costhetan*sintheta*cosphi, sinthetan*cosphin*costheta*sinphi + sinthetan*sinphin*cosphi + costhetan*sintheta*sinphi, -sinthetan*cosphin*sintheta + costhetan*costheta};
    
    p0[0] = sqrt(m0*m0+absP0*absP0);
    p1[0] = sqrt(m1*m1+absP1*absP1);
    double E2 = sqrt(m2*m2+absP2*absP2);
    
    std::vector<double > vp2 = {E2,absP2*d2[0],absP2*d2[1],absP2*d2[2]};
    
    p2[0]=vp2[0];
    p2[1]=vp2[1];
    p2[2]=vp2[2];
    p2[3]=vp2[3];
    
    
    p1[1] = -(tp0[1]+vp2[1]);
    p1[2] = -(tp0[2]+vp2[2]);
    p1[3] = -(tp0[3]+vp2[3]);
    return 0;
    
    
    
}



double smear_energy(double En, double Percen, TRandom3 * rangen){
    double ans = 0;
    while(ans <=0){
        ans = rangen->Gaus(En,Percen*En/sqrt(En));
    }
    return ans;
}

double smear_energy_type(int pdgid, double En, TRandom3 *rangen){
    double ans = -1;
    if (abs(pdgid) == 13){
        while(ans <0){
            ans = rangen->Gaus(En, MUsmear*En);
        }
        return ans;
    }
    else if (abs(pdgid) == 211){
        while(ans <0){
            ans = rangen->Gaus(En, pismear*En);
        }
        return ans;
    }
    else if (abs(pdgid) == 11 || pdgid ==22){
        double temp_sigma = sqrt( pow(EMflat*En,2) + pow(EMsmear*sqrt(En),2) );
        while(ans <0){
            ans = rangen->Gaus(En, temp_sigma);
        }
        return ans;
    }
    else if(pdgid == 2212){
        
        if (En < 0.4) {
            while(ans <0){
                ans = rangen -> Gaus(En, p_percent*En);
                if(ans<=0){
                    //cout << "case : slow proton" << endl;
                }
            }
            return ans;
        }
        else {
            double temp_sigma1 = sqrt( pow(pflat*En,2)+ pow(psmear*sqrt(En),2) );
            
            while(ans <0){
                ans = rangen->Gaus(En, temp_sigma1);
            }
            return ans;
        }
        
    }
    else if (pdgid == 2112){
        while(ans <0){
            ans = rangen->Gaus(En, nsmear*sqrt(En));
        }
        return ans;
    }
    else if (abs(pdgid) == 321 || pdgid == 311 || pdgid == 3222 || pdgid == 3112 || pdgid == 3122 ){
        
        double temp_sigma2 = sqrt( pow(oflat*En,2)+ pow(osmear*sqrt(En),2) );
        while(ans <0){
            ans = rangen->Gaus(En, temp_sigma2);
        }
        //cout << "case : other" << endl;
        return ans;
        
    }
    else {
        cout << "case : wrong call" << endl;
        //abort();
        return 0;
    }
    
}


double smear_energy_type(int pdgid, double En ,bool fully_contained, TRandom3 *rangen){
    double ans = -1;
    if (abs(pdgid) == 13){
        if(fully_contained){
            while(ans <0){
                ans = rangen->Gaus(En, MUsmear_track*sqrt(En));//MUsmear_track = 0.05;
            }
        }
        else{
            while(ans <0){
                ans = rangen->Gaus(En, MUsmear*En);//MUsmear = 0.3;
            }
        }
        return ans;
    }
    else if (abs(pdgid) == 211){
        
        if(fully_contained){
            while(ans <0){
                ans = rangen->Gaus(En, 0.15*sqrt(En));
            }
        }
        else{
            while(ans <0){
                ans = rangen->Gaus(En, pismear*En);
            }
        }
        return ans;
    }
    else {
        cout << "case : wrong call" << endl;
        //abort();
        return 0;
    }
    
}


double smear_angle(double th, double an, TRandom3 * rangen){
    double ans = 0;
    ans = rangen->Gaus(th,an);
    return ans;
}
double massive_smear_energy(double En, double Percen, TRandom3 * rangen, double mass){
    double ans = 0;
    while(ans <=mass){
        ans = rangen->Gaus(En,Percen*En/sqrt(En));
    }
    return ans;
}


double muon_track_length(double El){
    double	ans = 0.0;
    std::vector<double> list ={0.000925764,0.022103,0.0648654,0.122844,0.192052,0.269917,0.354701,0.445186,0.540483,0.639929,0.743014,0.849338,0.958581,1.07048,1.18483,1.30145,1.42017,1.54089,1.66347,1.78783,1.91388,2.04155,2.17076,2.30146,2.43359,2.5671,2.70195,2.83811,2.97552,3.11416,3.25399,3.39499,3.53713,3.68038,3.82471,3.97012,4.11657,4.26405,4.41253,4.56201,4.71246,4.86386,5.01621,5.16949,5.32368,5.47878,5.63476,5.79162,5.94935,6.10794,6.26737,6.42763,6.58873,6.75063,6.91335,7.07687,7.24118,7.40627,7.57213,7.73877,7.90617,8.07432,8.24321,8.41285,8.58323,8.75433,8.92615,9.09869,9.27194,9.4459,9.62056,9.79591,9.97195,10.1487,10.3261,10.5042,10.6829,10.8623,11.0424,11.2232,11.4046,11.5866,11.7693,11.9526,12.1366,12.3212,12.5065,12.6923,12.8788,13.0659,13.2537,13.442,13.631,13.8205,14.0107,14.2015,14.3929,14.5848,14.7774,14.9705,15.1643,15.3586,15.5535,15.749,15.9451,16.1417,16.339,16.5368,16.7351,16.934,17.1335,17.3336,17.5342,17.7354,17.9371,18.1394,18.3422,18.5456,18.7495,18.954,19.159,19.3646,19.5707,19.7773,19.9845,20.1922,20.4004,20.6092,20.8185,21.0283,21.2387,21.4495,21.6609,21.8729,22.0853,22.2983,22.5117,22.7257,22.9402,23.1553,23.3708,23.5868,23.8034,24.0205,24.238,24.4561,24.6747,24.8937,25.1133,25.3334,25.5539,25.775,25.9966,26.2186,26.4412,26.6642,26.8878,27.1118,27.3363,27.5613,27.7868,28.0128,28.2392,28.4662,28.6936,28.9215,29.1499,29.3787,29.6081,29.8379,30.0682,30.299,30.5302,30.7619,30.9941,31.2268,31.4599,31.6936,31.9276,32.1622,32.3972,32.6327,32.8686,33.105,33.3419,33.5793,33.8171,34.0553,34.2941,34.5332,34.7729,35.013,35.2536,35.4946,35.7361};
    
    /*	TF1 f("CSDA", "CSDA_integrand(x)", 0.105, El);
   		ROOT::Math::WrappedTF1 wf1(f);
     
   		// Create the Integrator
     ROOT::Math::Integrator ig(ROOT::Math::IntegrationOneDim::kNONADAPTIVE);
     
     
     // Set parameters of the integration
     ig.SetFunction(wf1, false);
     ig.SetRelTolerance(0.001);
     ans = ig.Integral(0.105,El);
     //	std::cout << "integral result is " << ans <<std::endl;
     
     return ans*100; //for meters
     */
    if(floor((El-0.11)/0.02)+2 >= list.size()){
        
        ans = list.back()*100;
        
        return ans;
    }
    
    double y0 = 100*list[(int)floor((El-0.11)/0.02)];
    double y1 = 100*list[(int)floor((El-0.11)/0.02)+1];
    double x0 = 0.11+floor((El-0.11)/0.02)*0.02;
    double x1 = 0.11+(floor((El-0.11)/0.02)+1)*0.02;
    
    ans = y0+(y1-y0)*(El-x0)/(x1-x0);
    return ans;
    
}

double pion_track_length(double El){
    double	ans = 0.0;
    std::vector<double> list ={0.000925764,0.022103,0.0648654,0.122844,0.192052,0.269917,0.354701,0.445186,0.540483,0.639929,0.743014,0.849338,0.958581,1.07048,1.18483,1.30145,1.42017,1.54089,1.66347,1.78783,1.91388,2.04155,2.17076,2.30146,2.43359,2.5671,2.70195,2.83811,2.97552,3.11416,3.25399,3.39499,3.53713,3.68038,3.82471,3.97012,4.11657,4.26405,4.41253,4.56201,4.71246,4.86386,5.01621,5.16949,5.32368,5.47878,5.63476,5.79162,5.94935,6.10794,6.26737,6.42763,6.58873,6.75063,6.91335,7.07687,7.24118,7.40627,7.57213,7.73877,7.90617,8.07432,8.24321,8.41285,8.58323,8.75433,8.92615,9.09869,9.27194,9.4459,9.62056,9.79591,9.97195,10.1487,10.3261,10.5042,10.6829,10.8623,11.0424,11.2232,11.4046,11.5866,11.7693,11.9526,12.1366,12.3212,12.5065,12.6923,12.8788,13.0659,13.2537,13.442,13.631,13.8205,14.0107,14.2015,14.3929,14.5848,14.7774,14.9705,15.1643,15.3586,15.5535,15.749,15.9451,16.1417,16.339,16.5368,16.7351,16.934,17.1335,17.3336,17.5342,17.7354,17.9371,18.1394,18.3422,18.5456,18.7495,18.954,19.159,19.3646,19.5707,19.7773,19.9845,20.1922,20.4004,20.6092,20.8185,21.0283,21.2387,21.4495,21.6609,21.8729,22.0853,22.2983,22.5117,22.7257,22.9402,23.1553,23.3708,23.5868,23.8034,24.0205,24.238,24.4561,24.6747,24.8937,25.1133,25.3334,25.5539,25.775,25.9966,26.2186,26.4412,26.6642,26.8878,27.1118,27.3363,27.5613,27.7868,28.0128,28.2392,28.4662,28.6936,28.9215,29.1499,29.3787,29.6081,29.8379,30.0682,30.299,30.5302,30.7619,30.9941,31.2268,31.4599,31.6936,31.9276,32.1622,32.3972,32.6327,32.8686,33.105,33.3419,33.5793,33.8171,34.0553,34.2941,34.5332,34.7729,35.013,35.2536,35.4946,35.7361};
    
   /* 
    TF1 f("CSDA", "CSDA_integrand(x)", 0.105, El);
   		ROOT::Math::WrappedTF1 wf1(f);
     
   		// Create the Integrator
     ROOT::Math::Integrator ig(ROOT::Math::IntegrationOneDim::kNONADAPTIVE);
     
     
     // Set parameters of the integration
     ig.SetFunction(wf1, false);
     ig.SetRelTolerance(0.001);
     ans = ig.Integral(0.105,El);
     //std::cout << "integral result is " << ans <<std::endl;
     */
    
    
    
    if(floor((El-0.11)/0.02)+2 >= list.size()){
        
        ans = list.back()*100;
        
        return ans;
    }
    
    
    double y0 = 100*list[(int)floor((El-0.11)/0.02)];
    double y1 = 100*list[(int)floor((El-0.11)/0.02)+1];
    double x0 = 0.11+floor((El-0.11)/0.02)*0.02;
    double x1 = 0.11+(floor((El-0.11)/0.02)+1)*0.02;
    
    ans = y0+(y1-y0)*(El-x0)/(x1-x0);
   // return ans;
    
    //std::cout << "interpolation result is " << ans<<std::endl;

    return ans; //for meters
    
}



double photon_conversion_length(double ep, TRandom3 * r){
    
    double ans= 0.0;
    double att_len = 0;
    
    if(ep < 0.01)
    {
        att_len = 25.0;
    }
    else if (ep <= 1.0)
    {
        att_len = exp(0.321799*ep)*(-391.977+402.232/(pow(ep,0.00797561)));
    }
    else {
        att_len = 15.0;
    }
    
    
    ans = r->Exp(att_len);
    //std::cout<<ep<<" "<<att_len<<" "<<ans/100.0<<std::endl;
    
    return ans;
    
}



double bethe(double beta){
    
    if(beta>1 || beta <0){std::cout<<"ERROR: beta must be 0<b<1 in bethe!"<<std::endl;}
    
    double me = 0.51;//in MeV
    //double c = 3.0E10;//in cm
    double mu = 105.;// in MeV
    double K  = 0.307075; // MeV mol^-1 cm^2
    double Z = 18;
    //double Z = 29;
    double A = 40;
    //double A = 63.5;
    double I = 10*Z;//mean excitation energy
    double d = 0;
    double z = -1;
    double gam = 1/sqrt(1-beta*beta);
    double wmax = 2*me*beta*beta*gam*gam/(1+2*gam*me/mu+pow(me/mu,2));
    
    return 0.001*fabs(K*z*z*Z/(A*beta*beta)*(0.5*log(2*me*pow(beta*gam,2)*wmax/(I*I))-beta*beta-d/2));//in GeV
    
    
}


double CSDA_integrand(double Emu){
    double beta = sqrt(Emu*Emu-0.105*0.105)/Emu;
    double rho = 1.3954;
    return 1.0/(bethe(beta)*rho);
    
}

double pion_containment(double posX, double posY, double posZ, TRandom3 * r)
{
    
    
    
    return 0;
}

int get_endpoint(double *vertex,double track_L, double * pl,double *  endpoint){
    double pnorm = sqrt(pl[0]*pl[0]+pl[1]*pl[1]+pl[2]*pl[2]);
    
    
    endpoint[0] = vertex[0]+track_L*pl[0]/pnorm;
    endpoint[1] = vertex[1]+track_L*pl[1]/pnorm;
    endpoint[2] = vertex[2]+track_L*pl[2]/pnorm;
    
    
    return 1;
}


SBN_detector::SBN_detector(double h, double w, double l, double fh, double fw, double fl, double base){
    height = h;
    length = l;
    width = w;
    volume = h*l*w;
    
    f_height = fh;
    f_length = fl;
    f_width = fw;
    f_volume = fh*fl*fw;
    
    baseline=base;
    
    dh = (h-fh)/2.0;
    dl = (l-fl)/2.0;
    dw = (w-fw)/2.0;
    
    proposal_modifier =1;
    
    fname ="../rootfiles/ntuple.ICARUS.root";
    potmodifier = 1.0;
    identifier = DET_ICARUS;
    
}

SBN_detector::SBN_detector(int ident, bool ismu ){
    
    // add case(DET_DUNE) ismu says is muon or is electron..
    
    //lengh units are centi-meters
    switch (ident)
    {
        case(DET_DUNE):
            height = 1200;
            length = 5750;
            width = 1440;
            volume = height*length*width;
            
            //f_height = 1140;
            //f_length = 5500;
            //f_width =1330 ;
            //f_volume = f_height*f_width*f_length;
            
            //if(ismu){
            
            
            
            //}
            
            
            //dh = (height-f_height)/2.0;
            //dl = (length-f_length)/2.0;
            //dw = (width-f_width)/2.0;
            
            //dh = 82.4;//94.8725
            //dl = 82.4;
            //dw = 82.4;
            dh = 90.5619;
            dl = 90.5619;
            dw = 90.5619;
            
            
            name = "DUNE";
            fname ="../rootfiles/ntuple.DUNE.root";
            //foscname ="../rootfiles/ntuple.SBND_fullosc.root";
            fbarname ="../rootfiles/NUBAR_MODE/ntuple.DUNE.root";
            //fbaroscname ="../rootfiles/NUBAR_MODE/ntuple.SBND_fullosc.root";
            
            potmodifier = 1.0;
            identifier = DET_DUNE;
            
            baseline = 1300*1000;
            
            break;
            
        case(DET_DUNE_NC):
            height = 1200;
            length = 230;
            width = 1440;
            volume = height*length*width;
            
            //f_height = 1140;
            //f_length = 5500;
            //f_width =1330 ;
            //f_volume = f_height*f_width*f_length;
            
            //if(ismu){
            // fiducial volume / fiducial surface area = (1200-94.8725)*(1450-94.8725)*(5800-94.8725) / 2( (1200-94.8725)*(1450-94.8725)+ (1450-94.8725)*(5800-94.8725)+(1200-94.8725)*(5800-94.8725)) ~= 275 cm
            
            // 150*6*3.5*2.3/(2*(600*350+350*230+230*600))
            
            
            // surface area 1 : 31067314 cm^2 = 3.1067*10^7
            // surface area 2 : 108120000 cm^2 = 1.0812*10^8
            
            
            //}
            
            
            //dh = (height-f_height)/2.0;
            //dl = (length-f_length)/2.0;
            //dw = (width-f_width)/2.0;
            
            //dh = 82.4;//94.8725
            //dl = 82.4;
            //dw = 82.4;
            //dh = 30.;//13.5206
            //dl = 30.;
            //dw = 30.;
            dh = 25.4909;
            dl = 25.4909;
            dw = 25.4909;
            
            
            name = "DUNE_NC";
            fname ="../rootfiles/ntuple.DUNE.root";
            //foscname ="../rootfiles/ntuple.SBND_fullosc.root";
            fbarname ="../rootfiles/NUBAR_MODE/ntuple.DUNE.root";
            //fbaroscname ="../rootfiles/NUBAR_MODE/ntuple.SBND_fullosc.root";
            
            potmodifier = 1.0*25.;//yj
            identifier = DET_DUNE_NC;
            
            baseline = 1300*1000;
            
            break;
            
            
            
        case(DET_SBND):
            height = 400;
            length = 500;
            width = 2*200;
            volume = height*length*width;
            mass =112;
            
            if(ismu)
            {
                f_height = 370;
                f_length = 405;
                f_width =2*183.5 ;
                f_mass = 77.0;
                
                proposal_modifier = 0.8*5212690.0/(6.6*887966.0*MET2IMP) ;
                
            }
            else
            {
                
                f_height =350;
                f_length = 420;
                f_width = 2*173.5 ;
                f_mass = 71.4;
                
                proposal_modifier = 36798.0/(6.6*6931*MET2IMP) ;
            }
            
            
            f_volume = f_height*f_length*f_width;
            
            
            baseline=110;
            
            
            
            dh = (height-f_height)/2.0;
            dl = (length-f_length)/2.0;
            dw = (width-f_width)/2.0;
            
            name = "SBND";
            fname ="../rootfiles/ntuple.SBND.root";
            foscname ="../rootfiles/ntuple.SBND_fullosc.root";
            fbarname ="../rootfiles/NUBAR_MODE/ntuple.SBND.root";
            fbaroscname ="../rootfiles/NUBAR_MODE/ntuple.SBND_fullosc.root";
            
            potmodifier = 1.0;
            identifier = DET_SBND;
            break;
        case(DET_UBOONE):
            height = 233.0;
            length = 1037.0;
            width = 256.0;
            mass = 86.6;
            volume = height*length*width;
            
            if(ismu)
            {
                
                f_height = 203;
                f_length = 942;
                f_width = 226;
                f_mass = 60.5;
                proposal_modifier =1.1 *  0.907*mass/89.0*173302/(6.6*26182);
                //proposal_modifier =  5212690.0/(6.6*887966.0*MET2IMP) ;
            }
            else
            {
                f_height = 183.0;
                f_length = 957.0;
                f_width = 206.0;
                f_mass = 47.9;
                
                //proposal_modifier =36798.0/(6.6*6931*MET2IMP) ;
                proposal_modifier =0.9*  0.90*86.6/89.0*1469/(6.6*258.7);
                //proposal_modifier   =1469/1417.9;
            }
            
            f_volume = f_height*f_length*f_width;
            
            baseline=470.0;
            
            
            dh = (height-f_height)/2.0;
            dl = (length-f_length)/2.0;
            dw = (width-f_width)/2.0;
            
            name = "uBooNE";
            fname ="../rootfiles/ntuple.uBooNE.root";
            foscname ="../rootfiles/ntuple.uBooNE_fullosc.root";
            fbarname ="../rootfiles/NUBAR_MODE/ntuple.uBooNE.root";
            fbaroscname ="../rootfiles/NUBAR_MODE/ntuple.uBooNE_fullosc.root";
            potmodifier = 2.0;
            identifier = DET_UBOONE;
            break;
        case(DET_ICARUS):
            
            height = 316.0;
            length = 1795;
            width = 4*150;
            volume = height*length*width;
            mass = 476;
            
            if(ismu)
            {
                
                f_height = 286;
                f_length = 1700;
                f_width = 4*133;
                f_mass = 363;
                //proposal_modifier =0.9*mass/476*173302*MET2IMP/(6.6*28182);
                proposal_modifier =1.1 *0.907*mass/476*173302/(6.6*26182);
                //proposal_modifier =  5212690.0/(6.6*887966.0*MET2IMP) ;
                //proposal_modifier = 173302*87*MET2IMP/(6.6*28182*89);
                
                
                
            }
            else
            {
                
                f_height = 266;
                f_length = 1715;
                f_width = 4*123.5;
                f_mass = 315;
                proposal_modifier =0.9*  0.90*86.6/89.0*1469/(6.6*258.7);
                //proposal_modifier = 1469*87*MET2IMP/(6.6*238.7*89);
            }
            
            
            
            
            f_volume = f_height*f_length*f_width;
            
            baseline=600;
            
            dh = (height-f_height)/2.0;
            dl = (length-f_length)/2.0;
            dw = (width-f_width)/2.0;
            
            name = "ICARUS";
            fname ="../rootfiles/ntuple.ICARUS.root";
            foscname ="../rootfiles/ntuple.ICARUS_fullosc.root";
            fbarname ="../rootfiles/NUBAR_MODE/ntuple.ICARUS.root";
            fbaroscname ="../rootfiles/NUBAR_MODE/ntuple.ICARUS_fullosc.root";
            
            potmodifier = 1.0;
            identifier = DET_ICARUS;
            break;
            
        default:
            std::cout<<"#ERROR: SBN_detector::SBN_detector(int ident) basefun.c, ident not one of 1,2,3. "<<std::endl;
            exit(EXIT_FAILURE);
    }
    
}

bool SBN_detector::is_active(double * pos){
    bool ans = false;
    
    if(0.0 <= pos[0] && pos[0]<= height && 0.0 <= pos[1] && pos[1]<= width && 0.0 <= pos[2] && pos[2]<= length  )
    {
        ans = true;
    }
    
    
    return ans;
}

bool SBN_detector::is_fiducial(double * pos){
    bool ans = false;
    
    if(dh <= pos[0] && pos[0]<= height-dh && dw <= pos[1] && pos[1]<= width-dw && dl <= pos[2] && pos[2]<= length-dl  )
    {
        ans = true;
    }
    
    
    return ans;
}

/*bool SBN_detector::is_fiducial(double * pos){
 //bool ans = false;
 
 //if(dh <= pos[0] && pos[0]<= height-dh && dw <= pos[1] && pos[1]<= width-dw && dl <= pos[2] && pos[2]<= length-dl  )
 //{
 //ans = true;
 //}
 bool x = false;
 bool y = false;
 bool z = false;
 
 if ( dh <= pos[0] && pos[0]<= height-dh ){
 x= true;
 }
 if ( (30<=pos[1]&&pos[1]<=695) || (755<=pos[1]&&pos[1]<=1420) ){
 y= true;
 }
 if ( dl <= pos[2] && pos[2]<= length-dl ){
 z = true;
 }
 
 //cout << (x&&y&&z) << endl;
 
 return (x && y && z);
 
 
 
 }*/

void SBN_detector::random_pos(TRandom3 * rangen, double * vec){
    
    vec[0]=rangen->Uniform(height);
    vec[1]=rangen->Uniform(width);
    vec[2]=rangen->Uniform(length);
    
}


bool SBN_detector::is_fully_contained(double *vertex,double * endpoint){
    bool ans = false;
    if( SBN_detector::is_active(vertex)&& SBN_detector::is_active(endpoint)){
        
        ans = true;
    }
    
    
    return ans;
}


double SBN_detector::track_length_escape(double * in, double *out){
    // Returns in CM
    
    if(is_active(out) || !is_active(in) )
    {
        std::cout<<"#ERROR: SBN_detector::track_length_escape. in array is outside OR out array is inside!"<<std::endl;
        return 0;
    }
    
    
    double a[3] = {in[0],in[1],in[2]};
    double m[3] = {out[0]-in[0],out[1]-in[1],out[2]-in[2]};
    double minlen= 1e12;
    int npos = 0;
    
    double t[6] = {-a[0]/m[0], (-a[0]+height)/m[0],-a[1]/m[1],(-a[1]+width)/m[1],-a[2]/m[2],(-a[2]+length)/m[2]};
    std::vector<double > tpos;
    
    for(int i=0;i<6;i++)
    {
        if(t[i]>=0)
        {
            npos++;
            tpos.push_back(t[i]);
        }
    }
    
    for(int j=0; j<npos; j++){
        
        double pt[3] = {a[0]+m[0]*tpos[j],a[1]+m[1]*tpos[j],a[2]+m[2]*tpos[j]};
        double len = sqrt(pow(pt[0]-in[0],2)+pow(pt[1]-in[1],2)+pow(pt[2]-in[2],2));
        
        if(len<=minlen)
        {
            minlen=len;
        }
        
    }
    
    
    return minlen;
    
    
}


double SBN_detector::osc_length(TRandom3 * rangen){
    
    return baseline - rangen->Uniform(50.0);
    
    
}

