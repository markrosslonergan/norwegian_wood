#include "TCanvas.h"

#include "TF1.h"
#include "TMath.h"

#define TAU 3
#define MUON 2
#define EL 1



double prob_eflav(double E, int flav, bool mode){
    
    //probability reference : 1505.01826.pdf
    if (!mode) {
        E *= -1.;
    }
    
    
    double L  = 1300.;
    
    double theta23 = 0.738;
    double theta12 = 0.5843;
    double theta13 = 0.148;
    
    double delta = 0.;//
    
    double m31sq = 2.457e-3;
    
    //double m31sq = -2.52e-3; //for inverted mass hierarchy
    double m21sq = 7.5e-5;
    
    double s12 = TMath::Sin(theta12);
    double s13 = TMath::Sin(theta13);
    double s23 = TMath::Sin(theta23);
    double c12 = TMath::Cos(theta12);
    double c13 = TMath::Cos(theta13);
    double c23 = TMath::Cos(theta23);
    
    double mrensq = m31sq-s12*s12*m21sq;
    double eps = m21sq/mrensq;
    double a = 1.52e-4*0.4957*2.8*E;
    //a*=2.;
    //a=0.;
    double lambdaP = 0.5*((mrensq+a)+TMath::Sign(1.0,mrensq)*TMath::Sqrt((mrensq-a)*(mrensq-a)+4*s13*s13*a*mrensq))+eps*mrensq*s12*s12;
    double lambdaM = 0.5*((mrensq+a)-TMath::Sign(1.0,mrensq)*TMath::Sqrt((mrensq-a)*(mrensq-a)+4*s13*s13*a*mrensq))+eps*mrensq*s12*s12;
    double lambdaN = c12*c12*eps*mrensq;
    
    
    double c2phi = (mrensq*TMath::Cos(2.0*theta13)-a)/(lambdaP-lambdaM);
    
    double sphisq = (1-c2phi)/2;
    double cphisq = (1+c2phi)/2;
    
    double S = 1.;
    double Jr = c12*s12*c23*s23*c13*c13*s13;
    
    double croneckerdelta, APM, APN, AMN;
    
    double BPM, C, BPN, BMN;
    
    switch (flav) {
            
        case 1:
            croneckerdelta = 1.;
            APM = -1.;
            APN = 0.;
            AMN = 0.;
            
            BPM = 0.;
            C = 0.;
            BPN = 0.;
            BMN = 0.;
            break;
        case 2:
            croneckerdelta = 0.;
            APM = s23*s23;
            APN = 0.;
            AMN = 0.;
            
            BPM = 1.;
            C = 1.;
            BPN = 0.;
            BMN = 0.;
            break;
        case 3:
            croneckerdelta = 0.;
            APM = c23*c23;
            APN = 0.;
            AMN = 0.;
            
            BPM = -1.;
            C = -1.;
            BPN = 0.;
            BMN = 0.;
            break;
        default:
            break;
    }
    
    
    
    double term0 = croneckerdelta;
    
    double term1 = 4*(APM*sphisq*cphisq+eps*BPM*Jr*mrensq*mrensq*((lambdaP-lambdaM)-(mrensq-a))/((lambdaP-lambdaM)*(lambdaP-lambdaM)*(lambdaP-lambdaN)));
    double freq1 = TMath::Sin((lambdaP-lambdaM)*1.27*L/E)*TMath::Sin((lambdaP-lambdaM)*1.27*L/E);
    
    double term2 = 4*(APN*cphisq+eps*BPN*(Jr/(c13*c13))*mrensq*((lambdaP-lambdaM)-(mrensq+a))/((lambdaP-lambdaM)*(lambdaP-lambdaN)));
    double freq2 = TMath::Sin((lambdaP-lambdaN)*1.27*L/E)*TMath::Sin((lambdaP-lambdaN)*1.27*L/E);
    
    double term3 = 4*(AMN*sphisq+eps*BMN*(Jr/(c13*c13))*mrensq*((lambdaP-lambdaM)+(mrensq+a))/((lambdaP-lambdaM)*(lambdaM-lambdaN)));
    double freq3 = TMath::Sin((lambdaM-lambdaN)*1.27*L/E)*TMath::Sin((lambdaM-lambdaN)*1.27*L/E);
    double term4_1 = 8*eps*Jr*mrensq*mrensq*mrensq/((lambdaP-lambdaM)*(lambdaP-lambdaN)*(lambdaM-lambdaN));
    double freq4_1 = TMath::Sin((lambdaP-lambdaM)*1.27*L/E)*TMath::Sin((lambdaM-lambdaN)*1.27*L/E)*C*TMath::Cos((lambdaP-lambdaN)*1.27*L/E);
    
    double prob = term0 + term1*freq1 + term2*freq2 + term3*freq3 + term4_1*freq4_1;
    //pro
    
    
    return prob;
}

double prob_muflav(double E, int flav, bool mode){
    
    //probability reference : 1505.01826.pdf
    
    if (!mode) {
        E *= -1.;
    }
    
    double L  = 1300.;
    
    double theta23 = 41.6*TMath::Pi()/180.;
    theta23 = 45.0*TMath::Pi()/180.;
    double theta12 = 33.5*TMath::Pi()/180.;
    double theta13 = 8.5*TMath::Pi()/180.;
    
    double delta = 0.;//
    
    double m31sq = 2.52e-3;
    //double m31sq = -2.52e-3; //for inverted mass hierarchy
    double m21sq = 7.5e-5;
    
    double s12 = TMath::Sin(theta12);
    double s13 = TMath::Sin(theta13);
    double s23 = TMath::Sin(theta23);
    double c12 = TMath::Cos(theta12);
    double c13 = TMath::Cos(theta13);
    double c23 = TMath::Cos(theta23);
    
    double mrensq = m31sq-s12*s12*m21sq;
    double eps = m21sq/mrensq;
    double a = 1.52e-4*0.4957*2.8*E;
    
    //a*=2.;
    //a=0.;
    double lambdaP = 0.5*((mrensq+a)+TMath::Sign(1.0,mrensq)*TMath::Sqrt((mrensq-a)*(mrensq-a)+4*s13*s13*a*mrensq))+eps*mrensq*s12*s12;
    double lambdaM = 0.5*((mrensq+a)-TMath::Sign(1.0,mrensq)*TMath::Sqrt((mrensq-a)*(mrensq-a)+4*s13*s13*a*mrensq))+eps*mrensq*s12*s12;
    double lambdaN = c12*c12*eps*mrensq;
    
    
    //double phi = 0.5*TMath::ASin(mrensq*TMath::Sin(2.0*theta13)/(lambdaP-lambdaM));
    //double cphi = TMath::Cos(phi);
    //double sphi = TMath::Sin(phi);
    double c2phi = (mrensq*TMath::Cos(2.0*theta13)-a)/(lambdaP-lambdaM);
    
    double sphisq = (1-c2phi)/2;
    double cphisq = (1+c2phi)/2;
    
    double S = 1.;
    double Jr = c12*s12*c23*s23*c13*c13*s13;
    
    double croneckerdelta, APM, APN, AMN;
    
    double BPM, C, BPN, BMN;
    
    switch (flav) {
            
        case 1:
            croneckerdelta = 0.;
            APM = s23*s23;
            APN = 0.;
            AMN = 0.;
            
            BPM = 1.;
            C = 1.;
            BPN = 0.;
            BMN = 0.;
            break;
        case 2:
            croneckerdelta = 1.;
            APM = -s23*s23*s23*s23;
            APN = -s23*s23*c23*c23;
            AMN = -s23*s23*c23*c23;
            
            BPM = TMath::Cos(2.0*theta23)-1;
            C = TMath::Cos(2.0*theta23)-1;
            BPN = TMath::Cos(2.0*theta23);
            BMN = TMath::Cos(2.0*theta23);
            break;
        case 3:
            croneckerdelta = 0.;
            APM = -s23*s23*c23*c23;
            APN = s23*s23*c23*c23;
            AMN = s23*s23*c23*c23;
            
            BPM = -TMath::Cos(2.0*theta23);
            C = -TMath::Cos(2.0*theta23);
            BPN = -TMath::Cos(2.0*theta23);
            BMN = -TMath::Cos(2.0*theta23);
            break;
        default:
            break;
    }
    
    
    
    double term0 = croneckerdelta;
    //double term1 = 4*(APM*sphi*sphi*cphi*cphi + eps*BPM*Jr*mrensq*mrensq*((lambdaP-lambdaM)-(mrensq-a))/(lambdaP-lambdaM)/(lambdaP-lambdaM)/(lambdaP-lambdaN));
    //double freq1 = TMath::Sin(1.27*(lambdaP-lambdaM)*L/E)*TMath::Sin(1.27*(lambdaP-lambdaM)*L/E);
    //double term2 = 4*(APN*cphi*cphi + eps*BPN*(Jr/c13/c13)*mrensq*((lambdaP-lambdaM)-(mrensq+a))/(lambdaP-lambdaM)/(lambdaP-lambdaN));
    //double freq2 = TMath::Sin(1.27*(lambdaP-lambdaN)*L/E)*TMath::Sin(1.27*(lambdaP-lambdaN)*L/E);
    //double term3 = 4*(AMN*sphi*sphi+eps*BMN*(Jr/c13/c13)*mrensq*((lambdaP-lambdaM)+(mrensq+a))/(lambdaP-lambdaM)/(lambdaM-lambdaN));
    
    //double freq3 = TMath::Sin(1.27*(lambdaM-lambdaN)*L/E)*TMath::Sin(1.27*(lambdaM-lambdaN)*L/E);
    //double term4 = 8*eps*Jr*mrensq*mrensq*mrensq/(lambdaP-lambdaM)/(lambdaP-lambdaN)/(lambdaM-lambdaN)*TMath::Sin(1.27*(lambdaP-lambdaM)*L/E)*TMath::Sin(1.27*(lambdaM-lambdaN)*L/E);
    //double freq4 = C*TMath::Cos(1.27*(lambdaP-lambdaN)*L/E);
    
    //double prob = term0 + term1*freq1 + term2*freq2 + term3*freq3 + term4*freq4;
    
    //double term0 = croneckerdelta;
    double term1 = 4*(APM*sphisq*cphisq+eps*BPM*Jr*mrensq*mrensq*((lambdaP-lambdaM)-(mrensq-a))/((lambdaP-lambdaM)*(lambdaP-lambdaM)*(lambdaP-lambdaN)));
    double freq1 = TMath::Sin((lambdaP-lambdaM)*1.27*L/E)*TMath::Sin((lambdaP-lambdaM)*1.27*L/E);
    
    double term2 = 4*(APN*cphisq+eps*BPN*(Jr/(c13*c13))*mrensq*((lambdaP-lambdaM)-(mrensq+a))/((lambdaP-lambdaM)*(lambdaP-lambdaN)));
    double freq2 = TMath::Sin((lambdaP-lambdaN)*1.27*L/E)*TMath::Sin((lambdaP-lambdaN)*1.27*L/E);
    
    double term3 = 4*(AMN*sphisq+eps*BMN*(Jr/(c13*c13))*mrensq*((lambdaP-lambdaM)+(mrensq+a))/((lambdaP-lambdaM)*(lambdaM-lambdaN)));
    double freq3 = TMath::Sin((lambdaM-lambdaN)*1.27*L/E)*TMath::Sin((lambdaM-lambdaN)*1.27*L/E);
    double term4_1 = 8*eps*Jr*mrensq*mrensq*mrensq/((lambdaP-lambdaM)*(lambdaP-lambdaN)*(lambdaM-lambdaN));
    double freq4_1 = TMath::Sin((lambdaP-lambdaM)*1.27*L/E)*TMath::Sin((lambdaM-lambdaN)*1.27*L/E)*C*TMath::Cos((lambdaP-lambdaN)*1.27*L/E);
    
    double prob = term0 + term1*freq1 + term2*freq2 + term3*freq3 + term4_1*freq4_1;
    //pro
    
    
    return prob;
}





void genie_tau_test(){
    
    bool mode = true;

    TF1 *mue_prob = new TF1("mue_prob", "prob_muflav(x,EL,1)", 0.1, 100.);
    mue_prob->SetNpx(100000);
    TCanvas *c1 = new TCanvas("c1");
    c1->cd();
    mue_prob->Draw();
    gPad->SetLogx();
    c1->SaveAs("c1_1.pdf");
    
    TF1 *mumu_prob = new TF1("mumu_prob", "prob_muflav(x,MUON,1)", 0.1, 100.);
    mumu_prob->SetNpx(100000);
    TCanvas *c2 = new TCanvas("c2");
    c2->cd();
    mumu_prob->Draw();
    gPad->SetLogx();
    c2->SaveAs("c2_1.pdf");
    
    TF1 *mutau_prob = new TF1("mutau_prob", "prob_muflav(x,TAU,1)", 0.1, 100.);
    mutau_prob->SetNpx(100000);
    TCanvas *c3 = new TCanvas("c3");
    c3->cd();
    mutau_prob->Draw();
    gPad->SetLogx();
    c3->SaveAs("c3_1.pdf");
    
    TF1 *ee_prob = new TF1("ee_prob", "prob_eflav(x,EL,1)", 0.1, 100.);
    ee_prob->SetNpx(100000);
    TCanvas *c4 = new TCanvas("c4");
    c4->cd();
    ee_prob->Draw();
    gPad->SetLogx();
    c4->SaveAs("c4_1.pdf");
    
    TF1 *emu_prob = new TF1("emu_prob", "prob_eflav(x,MUON,1)", 0.1, 100.);
    emu_prob->SetNpx(100000);
    TCanvas *c5 = new TCanvas("c5");
    c5->cd();
    emu_prob->Draw();
    gPad->SetLogx();
    c5->SaveAs("c5_1.pdf");
    
    TF1 *etau_prob = new TF1("etau_prob", "prob_eflav(x,TAU,1)", 0.1, 100.);
    etau_prob->SetNpx(100000);
    TCanvas *c6 = new TCanvas("c6");
    c6->cd();
    etau_prob->Draw();
    gPad->SetLogx();
    c6->SaveAs("c6_1.pdf");
    
    
    
    
    TF1 *esum_prob = new TF1("esum_prob", "prob_eflav(x,EL,1)+prob_eflav(x,MUON,1)+prob_eflav(x,TAU,1)", 0.1, 100.);
    esum_prob->SetNpx(100000);
    TCanvas *c7 = new TCanvas("c7");
    c7->cd();
    esum_prob->Draw();
    gPad->SetLogx();
    c7->SaveAs("c7_1.pdf");
    
    TF1 *musum_prob = new TF1("musum_prob", "prob_muflav(x,EL,1)+prob_muflav(x,MUON,1)+prob_muflav(x,TAU,1)", 0.1, 100.);
    musum_prob->SetNpx(100000);
    TCanvas *c8 = new TCanvas("c8");
    c8->cd();
    musum_prob->Draw();
    gPad->SetLogx();
    c8->SaveAs("c8_1.pdf");
    
    

    
    
    
}
