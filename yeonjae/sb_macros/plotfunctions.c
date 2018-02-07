
#include "TFile.h"
#include "TCanvas.h"
#include "TH1F.h"

#include "THStack.h"
#include "TLegend.h"


void plot_cut_test_weak(){
    
    TFile *fnue = new TFile("test_nue_CC_slimmed_af.root");
    TFile *fnutau = new TFile("test_nutau_CC_slimmed_af.root");
    
    TH1F *hist_nuecc_PT = (TH1F*) fnue->Get("hist_PT_weak");
    TH1F *hist_nuecc_PTsq = (TH1F*) fnue->Get("hist_PTsq_weak");
    TH1F *hist_nuecc_PTcube = (TH1F*) fnue->Get("hist_PTcube_weak");
    TH1F *hist_nuecc_PTquartic = (TH1F*) fnue->Get("hist_PTquartic_weak");
    TH1F *hist_nuecc_PTsqOverPLsq = (TH1F*) fnue->Get("hist_PTsqOverPLsq_weak");
    TH1F *hist_nuecc_PTsqOverPL = (TH1F*) fnue->Get("hist_PTsqOverPL_weak");
    
    TH1F *hist_nutaucc_PT = (TH1F*) fnutau->Get("hist_PT_weak");
    TH1F *hist_nutaucc_PTsq = (TH1F*) fnutau->Get("hist_PTsq_weak");
    TH1F *hist_nutaucc_PTcube = (TH1F*) fnutau->Get("hist_PTcube_weak");
    TH1F *hist_nutaucc_PTquartic = (TH1F*) fnutau->Get("hist_PTquartic_weak");
    TH1F *hist_nutaucc_PTsqOverPLsq = (TH1F*) fnutau->Get("hist_PTsqOverPLsq_weak");
    TH1F *hist_nutaucc_PTsqOverPL = (TH1F*) fnutau->Get("hist_PTsqOverPL_weak");
    
    
    
    hist_nuecc_PT->SetLineColor(kRed);
    hist_nuecc_PTsq->SetLineColor(kRed);
    hist_nuecc_PTcube->SetLineColor(kRed);
    hist_nuecc_PTquartic->SetLineColor(kRed);
    hist_nuecc_PTsqOverPLsq->SetLineColor(kRed);
    hist_nuecc_PTsqOverPL->SetLineColor(kRed);
    
    hist_nuecc_PT->SetLineWidth(2);
    hist_nuecc_PTsq->SetLineWidth(2);
    hist_nuecc_PTcube->SetLineWidth(2);
    hist_nuecc_PTquartic->SetLineWidth(2);
    hist_nuecc_PTsqOverPLsq->SetLineWidth(2);
    hist_nuecc_PTsqOverPL->SetLineWidth(2);
    
    hist_nutaucc_PT->SetLineWidth(2);
    
    hist_nutaucc_PT->SetLineStyle(1);
    hist_nutaucc_PTsq->SetLineWidth(2);
    hist_nutaucc_PTcube->SetLineWidth(2);
    hist_nutaucc_PTquartic->SetLineWidth(2);
    hist_nutaucc_PTsqOverPLsq->SetLineWidth(2);
    hist_nutaucc_PTsqOverPL->SetLineWidth(2);
    
    
    
    TCanvas *c1 = new TCanvas ("c1","c1",800,600);
    //c1->SetLogy();
    c1->cd();
    hist_nuecc_PT->DrawNormalized();
    hist_nutaucc_PT->DrawNormalized("histsame");
    TLegend *leg1 = new TLegend(0.6,0.6,0.9,0.9);
    leg1->AddEntry(hist_nuecc_PT,"nuecc_PT","f");
    leg1->AddEntry(hist_nutaucc_PT,"nutaucc_PT","f");
    leg1->Draw();
    c1->SaveAs("norm_PT_compare_norm.pdf");
    
   
    TCanvas *c2 = new TCanvas ("c2","c2",800,600);
    //c1->SetLogy();
    c2->cd();
    hist_nuecc_PTsq->DrawNormalized();
    hist_nutaucc_PTsq->DrawNormalized("histsame");
    TLegend *leg2 = new TLegend(0.6,0.6,0.9,0.9);
    leg2->AddEntry(hist_nuecc_PTsq,"nuecc_PTsq","f");
    leg2->AddEntry(hist_nutaucc_PTsq,"nutaucc_PTsq","f");
    leg2->Draw();
    c2->SaveAs("norm_PTsq_compare_norm.pdf");
    
    TCanvas *c3 = new TCanvas ("c3","c3",800,600);
    //c1->SetLogy();
    c3->cd();
    hist_nuecc_PTcube->DrawNormalized();
    hist_nutaucc_PTcube->DrawNormalized("histsame");
    TLegend *leg3 = new TLegend(0.6,0.6,0.9,0.9);
    leg3->AddEntry(hist_nuecc_PTcube,"nuecc_PTcube","f");
    leg3->AddEntry(hist_nutaucc_PTcube,"nutaucc_PTcube","f");
    leg3->Draw();
    c3->SaveAs("norm_PTcube_compare_norm.pdf");
    
    
    TCanvas *c4 = new TCanvas ("c4","c4",800,600);
    //c1->SetLogy();
    c4->cd();
    hist_nuecc_PTquartic->DrawNormalized();
    hist_nutaucc_PTquartic->DrawNormalized("histsame");
    TLegend *leg4 = new TLegend(0.6,0.6,0.9,0.9);
    leg4->AddEntry(hist_nuecc_PTquartic,"nuecc_PTquartic","f");
    leg4->AddEntry(hist_nutaucc_PTquartic,"nutaucc_PTquartic","f");
    leg4->Draw();
    c4->SaveAs("norm_PTquartic_compare_norm.pdf");
    
    
    TCanvas *c5 = new TCanvas ("c5", "c5", 800, 600);
    c5->cd();
    hist_nuecc_PTsqOverPLsq->DrawNormalized();
    hist_nutaucc_PTsqOverPLsq->DrawNormalized("histsame");
    TLegend *leg5 = new TLegend(0.6,0.6,0.9,0.9);
    leg5->AddEntry(hist_nuecc_PTquartic,"nuecc_PTsqOverPLsq","f");
    leg5->AddEntry(hist_nutaucc_PTquartic,"nutaucc_PTsqOverPLsq","f");
    leg5->Draw();
    c5->SaveAs("norm_PTsqOverPLsq_compare_norm.pdf");
    
    TCanvas *c6 = new TCanvas ("c6", "c6", 800, 600);
    c6->cd();
    hist_nuecc_PTsqOverPL->DrawNormalized();
    hist_nutaucc_PTsqOverPL->DrawNormalized("histsame");
    TLegend *leg6 = new TLegend(0.6,0.6,0.9,0.9);
    leg6->AddEntry(hist_nuecc_PTquartic,"nuecc_PTsqOverPL","f");
    leg6->AddEntry(hist_nutaucc_PTquartic,"nutaucc_PTsqOverPL","f");
    leg6->Draw();
    c6->SaveAs("norm_PTsqOverPL_compare_norm.pdf");
    
    //bf
    //int nbinsPT = 9;
    //int nbinsPTsq = 6;
    //int nbinsPTcube = 5;
    //int nbinsPTquartic = 4;
    //int nbinsPTsqOverPLsq = 3;
    //int nbinsPTsqOverPL = 18;
    
    //af
    int nbinsPT = 8;
    int nbinsPTsq = 6;
    int nbinsPTcube = 4;
    int nbinsPTquartic = 2;
    int nbinsPTsqOverPLsq = 4;
    int nbinsPTsqOverPL = 8;
    
    
    //cout << hist_nuecc_PT->GetBinContent(1) << endl;
    
    hist_nuecc_PT->Scale(1./hist_nuecc_PT->GetEntries());
    hist_nuecc_PTsq->Scale(1./hist_nuecc_PTsq->GetEntries());
    hist_nuecc_PTcube->Scale(1./hist_nuecc_PTcube->GetEntries());
    hist_nuecc_PTquartic->Scale(1./hist_nuecc_PTquartic->GetEntries());
    hist_nuecc_PTsqOverPLsq->Scale(1./hist_nuecc_PTsqOverPLsq->GetEntries());
    hist_nuecc_PTsqOverPL->Scale(1./hist_nuecc_PTsqOverPL->GetEntries());
    
    hist_nutaucc_PT->Scale(1./hist_nutaucc_PT->Integral());
    hist_nutaucc_PTsq->Scale(1./hist_nutaucc_PTsq->Integral());
    hist_nutaucc_PTcube->Scale(1./hist_nutaucc_PTcube->Integral());
    hist_nutaucc_PTquartic->Scale(1./hist_nutaucc_PTquartic->Integral());
    hist_nutaucc_PTsqOverPLsq->Scale(1./hist_nutaucc_PTsqOverPLsq->Integral());
    hist_nutaucc_PTsqOverPL->Scale(1./hist_nutaucc_PTsqOverPL->Integral());
    
    
    double nuecc_PT_integral = 0.;
    double nuecc_PTsq_integral = 0.;
    double nuecc_PTcube_integral = 0.;
    double nuecc_PTquartic_integral = 0.;
    double nuecc_PTsqOverPLsq_integral = 0.;
    double nuecc_PTsqOverPL_integral = 0.;
    
    double nutaucc_PT_integral = 0.;
    double nutaucc_PTsq_integral = 0.;
    double nutaucc_PTcube_integral = 0.;
    double nutaucc_PTquartic_integral = 0.;
    double nutaucc_PTsqOverPLsq_integral = 0.;
    double nutaucc_PTsqOverPL_integral = 0.;
    
    
    for (int i=1; i<=nbinsPT; i++){
        nuecc_PT_integral+=hist_nuecc_PT->GetBinContent(i);
        nutaucc_PT_integral+=hist_nutaucc_PT->GetBinContent(i);
    }
    for (int i=1; i<=nbinsPTsq; i++){
        nuecc_PTsq_integral+=hist_nuecc_PTsq->GetBinContent(i);
        nutaucc_PTsq_integral+=hist_nutaucc_PTsq->GetBinContent(i);
    }
    for (int i=1; i<=nbinsPTcube; i++){
        nuecc_PTcube_integral+=hist_nuecc_PTcube->GetBinContent(i);
        nutaucc_PTcube_integral+=hist_nutaucc_PTcube->GetBinContent(i);
    }
    for (int i=1; i<=nbinsPTquartic; i++){
        nuecc_PTquartic_integral+=hist_nuecc_PTquartic->GetBinContent(i);
        nutaucc_PTquartic_integral+=hist_nutaucc_PTquartic->GetBinContent(i);
    }
    for (int i=1; i<=nbinsPTsqOverPLsq; i++){
        nuecc_PTsqOverPLsq_integral+=hist_nuecc_PTsqOverPLsq->GetBinContent(i);
        nutaucc_PTsqOverPLsq_integral+=hist_nutaucc_PTsqOverPLsq->GetBinContent(i);
    }
    for (int i=1; i<=nbinsPTsqOverPL; i++){
        nuecc_PTsqOverPL_integral+=hist_nuecc_PTsqOverPL->GetBinContent(i);
        nutaucc_PTsqOverPL_integral+=hist_nutaucc_PTsqOverPL->GetBinContent(i);
    }
    
    cout << "pT nue: " << nuecc_PT_integral << " , nutau: " << nutaucc_PT_integral << endl;
    cout << "pTsq nue: " << nuecc_PTsq_integral << " , nutau: " << nutaucc_PTsq_integral << endl;
    cout << "pTcube nue: " << nuecc_PTcube_integral << " , nutau: " << nutaucc_PTcube_integral << endl;
    cout << "pTquartic nue: " << nuecc_PTquartic_integral << " , nutau: " << nutaucc_PTquartic_integral << endl;
    cout << "pTsqOverPLsq nue: " << nuecc_PTsqOverPLsq_integral << " , nutau: " << nutaucc_PTsqOverPLsq_integral << endl;
    cout << "pTsqOverPL nue: " << nuecc_PTsqOverPL_integral << " , nutau: " << nutaucc_PTsqOverPL_integral << endl;
    
    
    
    //cout << hist_nuecc_PT->GetBinContent(1)+hist_nuecc_PT->GetBinContent(2)+hist_nuecc_PT->GetBinContent(3)+hist_nuecc_PT->GetBinContent(4)+hist_nuecc_PT->GetBinContent(5) << endl;
    //cout << hist_nutaucc_PT->GetBinContent(1)+hist_nutaucc_PT->GetBinContent(2)+hist_nutaucc_PT->GetBinContent(3)+hist_nutaucc_PT->GetBinContent(4)+hist_nutaucc_PT->GetBinContent(5) << endl;
    
    
    //cout << hist_nuecc_PTsq->GetBinContent(1)+hist_nuecc_PTsq->GetBinContent(2)+hist_nuecc_PTsq->GetBinContent(3)+hist_nuecc_PTsq->GetBinContent(4)+hist_nuecc_PTsq->GetBinContent(5)+hist_nuecc_PTsq->GetBinContent(6) << endl;
    //cout << hist_nutaucc_PTsq->GetBinContent(1)+hist_nutaucc_PTsq->GetBinContent(2)+hist_nutaucc_PTsq->GetBinContent(3)+hist_nutaucc_PTsq->GetBinContent(4)+hist_nutaucc_PTsq->GetBinContent(5)+hist_nutaucc_PTsq->GetBinContent(6) << endl;
    
    //cout << hist_nuecc_PTcube->GetBinContent(1)+hist_nuecc_PTcube->GetBinContent(2)+hist_nuecc_PTcube->GetBinContent(3)+hist_nuecc_PTcube->GetBinContent(4) << endl;
    //cout << hist_nutaucc_PTcube->GetBinContent(1)+hist_nutaucc_PTcube->GetBinContent(2)+hist_nutaucc_PTcube->GetBinContent(3)+hist_nutaucc_PTcube->GetBinContent(4) << endl;
    
    //cout << hist_nuecc_PTquartic->GetBinContent(1)+hist_nuecc_PTquartic->GetBinContent(2) << endl;
    //cout << hist_nutaucc_PTquartic->GetBinContent(1)+hist_nutaucc_PTquartic->GetBinContent(2) << endl;
    
    
    
}



void plot_cut_test(){

    TFile *fnue = new TFile("test_nue_CC_slimmed_af.root");
    TFile *fnutau = new TFile("test_nutau_CC_slimmed_af.root");
    
    TH1F *hist_nuecc_PT = (TH1F*) fnue->Get("hist_PT");
    TH1F *hist_nuecc_PTsq = (TH1F*) fnue->Get("hist_PTsq");
    TH1F *hist_nuecc_PTcube = (TH1F*) fnue->Get("hist_PTcube");
    TH1F *hist_nuecc_PTquartic = (TH1F*) fnue->Get("hist_PTquartic");
    TH1F *hist_nuecc_PTsqOverPLsq = (TH1F*) fnue->Get("hist_PTsqOverPLsq");
    TH1F *hist_nuecc_PTsqOverPL = (TH1F*) fnue->Get("hist_PTsqOverPL");
    
    TH1F *hist_nutaucc_PT = (TH1F*) fnutau->Get("hist_PT");
    TH1F *hist_nutaucc_PTsq = (TH1F*) fnutau->Get("hist_PTsq");
    TH1F *hist_nutaucc_PTcube = (TH1F*) fnutau->Get("hist_PTcube");
    TH1F *hist_nutaucc_PTquartic = (TH1F*) fnutau->Get("hist_PTquartic");
    TH1F *hist_nutaucc_PTsqOverPLsq = (TH1F*) fnutau->Get("hist_PTsqOverPLsq");
    TH1F *hist_nutaucc_PTsqOverPL = (TH1F*) fnutau->Get("hist_PTsqOverPL");
    
    
    
    hist_nuecc_PT->SetLineColor(kRed);
    hist_nuecc_PTsq->SetLineColor(kRed);
    hist_nuecc_PTcube->SetLineColor(kRed);
    hist_nuecc_PTquartic->SetLineColor(kRed);
    hist_nuecc_PTsqOverPLsq->SetLineColor(kRed);
    hist_nuecc_PTsqOverPL->SetLineColor(kRed);
    
    hist_nuecc_PT->SetLineWidth(2);
    hist_nuecc_PTsq->SetLineWidth(2);
    hist_nuecc_PTcube->SetLineWidth(2);
    hist_nuecc_PTquartic->SetLineWidth(2);
    hist_nuecc_PTsqOverPLsq->SetLineWidth(2);
    hist_nuecc_PTsqOverPL->SetLineWidth(2);
    
    hist_nutaucc_PT->SetLineWidth(2);
    hist_nutaucc_PTsq->SetLineWidth(2);
    hist_nutaucc_PTcube->SetLineWidth(2);
    hist_nutaucc_PTquartic->SetLineWidth(2);
    hist_nutaucc_PTsqOverPLsq->SetLineWidth(2);
    hist_nutaucc_PTsqOverPL->SetLineWidth(2);

    
    
    TCanvas *c1 = new TCanvas ("c1","c1",800,600);
    //c1->SetLogy();
    c1->cd();
    hist_nuecc_PT->DrawNormalized();
    hist_nutaucc_PT->DrawNormalized("same");
    TLegend *leg1 = new TLegend(0.6,0.6,0.9,0.9);
    leg1->AddEntry(hist_nuecc_PT,"nuecc_PT","f");
    leg1->AddEntry(hist_nutaucc_PT,"nutaucc_PT","f");
    leg1->Draw();
    c1->SaveAs("norm_PT_compare_norm.pdf");
    
    
    TCanvas *c2 = new TCanvas ("c2","c2",800,600);
    //c1->SetLogy();
    c2->cd();
    hist_nuecc_PTsq->DrawNormalized();
    hist_nutaucc_PTsq->DrawNormalized("same");
    TLegend *leg2 = new TLegend(0.6,0.6,0.9,0.9);
    leg2->AddEntry(hist_nuecc_PTsq,"nuecc_PTsq","f");
    leg2->AddEntry(hist_nutaucc_PTsq,"nutaucc_PTsq","f");
    leg2->Draw();
    c2->SaveAs("norm_PTsq_compare_norm.pdf");

    TCanvas *c3 = new TCanvas ("c3","c3",800,600);
    //c1->SetLogy();
    c3->cd();
    hist_nuecc_PTcube->DrawNormalized();
    hist_nutaucc_PTcube->DrawNormalized("same");
    TLegend *leg3 = new TLegend(0.6,0.6,0.9,0.9);
    leg3->AddEntry(hist_nuecc_PTcube,"nuecc_PTcube","f");
    leg3->AddEntry(hist_nutaucc_PTcube,"nutaucc_PTcube","f");
    leg3->Draw();
    c3->SaveAs("norm_PTcube_compare_norm.pdf");
    
    
    TCanvas *c4 = new TCanvas ("c4","c4",800,600);
    //c1->SetLogy();
    c4->cd();
    hist_nuecc_PTquartic->DrawNormalized();
    hist_nutaucc_PTquartic->DrawNormalized("same");
    TLegend *leg4 = new TLegend(0.6,0.6,0.9,0.9);
    leg4->AddEntry(hist_nuecc_PTquartic,"nuecc_PTquartic","f");
    leg4->AddEntry(hist_nutaucc_PTquartic,"nutaucc_PTquartic","f");
    leg4->Draw();
    c4->SaveAs("norm_PTquartic_compare_norm.pdf");
    
    
    TCanvas *c5 = new TCanvas ("c5", "c5", 800, 600);
    c5->cd();
    hist_nuecc_PTsqOverPLsq->DrawNormalized();
    hist_nutaucc_PTsqOverPLsq->DrawNormalized("same");
    TLegend *leg5 = new TLegend(0.6,0.6,0.9,0.9);
    leg5->AddEntry(hist_nuecc_PTquartic,"nuecc_PTsqOverPLsq","f");
    leg5->AddEntry(hist_nutaucc_PTquartic,"nutaucc_PTsqOverPLsq","f");
    leg5->Draw();
    c5->SaveAs("norm_PTsqOverPLsq_compare_norm.pdf");
    
    TCanvas *c6 = new TCanvas ("c6", "c6", 800, 600);
    c6->cd();
    hist_nuecc_PTsqOverPL->DrawNormalized();
    hist_nutaucc_PTsqOverPL->DrawNormalized("same");
    TLegend *leg6 = new TLegend(0.6,0.6,0.9,0.9);
    leg6->AddEntry(hist_nuecc_PTquartic,"nuecc_PTsqOverPL","f");
    leg6->AddEntry(hist_nutaucc_PTquartic,"nutaucc_PTsqOverPL","f");
    leg6->Draw();
    c6->SaveAs("norm_PTsqOverPL_compare_norm.pdf");

    //bf
    //int nbinsPT = 9;
    //int nbinsPTsq = 6;
    //int nbinsPTcube = 5;
    //int nbinsPTquartic = 4;
    //int nbinsPTsqOverPLsq = 3;
    //int nbinsPTsqOverPL = 18;
    
    //af
    int nbinsPT = 9;
    int nbinsPTsq = 7;
    int nbinsPTcube = 4;
    int nbinsPTquartic = 2;
    int nbinsPTsqOverPLsq = 3;
    int nbinsPTsqOverPL = 19;
    
    
    //cout << hist_nuecc_PT->GetBinContent(1) << endl;
    
    hist_nuecc_PT->Scale(1./hist_nuecc_PT->GetEntries());
    hist_nuecc_PTsq->Scale(1./hist_nuecc_PTsq->GetEntries());
    hist_nuecc_PTcube->Scale(1./hist_nuecc_PTcube->GetEntries());
    hist_nuecc_PTquartic->Scale(1./hist_nuecc_PTquartic->GetEntries());
    hist_nuecc_PTsqOverPLsq->Scale(1./hist_nuecc_PTsqOverPLsq->GetEntries());
    hist_nuecc_PTsqOverPL->Scale(1./hist_nuecc_PTsqOverPL->GetEntries());
    
    hist_nutaucc_PT->Scale(1./hist_nutaucc_PT->GetEntries());
    hist_nutaucc_PTsq->Scale(1./hist_nutaucc_PTsq->GetEntries());
    hist_nutaucc_PTcube->Scale(1./hist_nutaucc_PTcube->GetEntries());
    hist_nutaucc_PTquartic->Scale(1./hist_nutaucc_PTquartic->GetEntries());
    hist_nutaucc_PTsqOverPLsq->Scale(1./hist_nutaucc_PTsqOverPLsq->GetEntries());
    hist_nutaucc_PTsqOverPL->Scale(1./hist_nutaucc_PTsqOverPL->GetEntries());

    
    double nuecc_PT_integral = 0.;
    double nuecc_PTsq_integral = 0.;
    double nuecc_PTcube_integral = 0.;
    double nuecc_PTquartic_integral = 0.;
    double nuecc_PTsqOverPLsq_integral = 0.;
    double nuecc_PTsqOverPL_integral = 0.;
    
    double nutaucc_PT_integral = 0.;
    double nutaucc_PTsq_integral = 0.;
    double nutaucc_PTcube_integral = 0.;
    double nutaucc_PTquartic_integral = 0.;
    double nutaucc_PTsqOverPLsq_integral = 0.;
    double nutaucc_PTsqOverPL_integral = 0.;

    
    for (int i=1; i<=nbinsPT; i++){
        nuecc_PT_integral+=hist_nuecc_PT->GetBinContent(i);
        nutaucc_PT_integral+=hist_nutaucc_PT->GetBinContent(i);
    }
    for (int i=1; i<=nbinsPTsq; i++){
        nuecc_PTsq_integral+=hist_nuecc_PTsq->GetBinContent(i);
        nutaucc_PTsq_integral+=hist_nutaucc_PTsq->GetBinContent(i);
    }
    for (int i=1; i<=nbinsPTcube; i++){
        nuecc_PTcube_integral+=hist_nuecc_PTcube->GetBinContent(i);
        nutaucc_PTcube_integral+=hist_nutaucc_PTcube->GetBinContent(i);
    }
    for (int i=1; i<=nbinsPTquartic; i++){
        nuecc_PTquartic_integral+=hist_nuecc_PTquartic->GetBinContent(i);
        nutaucc_PTquartic_integral+=hist_nutaucc_PTquartic->GetBinContent(i);
    }
    for (int i=1; i<=nbinsPTsqOverPLsq; i++){
        nuecc_PTsqOverPLsq_integral+=hist_nuecc_PTsqOverPLsq->GetBinContent(i);
        nutaucc_PTsqOverPLsq_integral+=hist_nutaucc_PTsqOverPLsq->GetBinContent(i);
    }
    for (int i=1; i<=nbinsPTsqOverPL; i++){
        nuecc_PTsqOverPL_integral+=hist_nuecc_PTsqOverPL->GetBinContent(i);
        nutaucc_PTsqOverPL_integral+=hist_nutaucc_PTsqOverPL->GetBinContent(i);
    }
    
    cout << "pT nue: " << nuecc_PT_integral << " , nutau: " << nutaucc_PT_integral << endl;
    cout << "pTsq nue: " << nuecc_PTsq_integral << " , nutau: " << nutaucc_PTsq_integral << endl;
    cout << "pTcube nue: " << nuecc_PTcube_integral << " , nutau: " << nutaucc_PTcube_integral << endl;
    cout << "pTquartic nue: " << nuecc_PTquartic_integral << " , nutau: " << nutaucc_PTquartic_integral << endl;
    cout << "pTsqOverPLsq nue: " << nuecc_PTsqOverPLsq_integral << " , nutau: " << nutaucc_PTsqOverPLsq_integral << endl;
    cout << "pTsqOverPL nue: " << nuecc_PTsqOverPL_integral << " , nutau: " << nutaucc_PTsqOverPL_integral << endl;
    

    
    //cout << hist_nuecc_PT->GetBinContent(1)+hist_nuecc_PT->GetBinContent(2)+hist_nuecc_PT->GetBinContent(3)+hist_nuecc_PT->GetBinContent(4)+hist_nuecc_PT->GetBinContent(5) << endl;
    //cout << hist_nutaucc_PT->GetBinContent(1)+hist_nutaucc_PT->GetBinContent(2)+hist_nutaucc_PT->GetBinContent(3)+hist_nutaucc_PT->GetBinContent(4)+hist_nutaucc_PT->GetBinContent(5) << endl;
    
    
    //cout << hist_nuecc_PTsq->GetBinContent(1)+hist_nuecc_PTsq->GetBinContent(2)+hist_nuecc_PTsq->GetBinContent(3)+hist_nuecc_PTsq->GetBinContent(4)+hist_nuecc_PTsq->GetBinContent(5)+hist_nuecc_PTsq->GetBinContent(6) << endl;
    //cout << hist_nutaucc_PTsq->GetBinContent(1)+hist_nutaucc_PTsq->GetBinContent(2)+hist_nutaucc_PTsq->GetBinContent(3)+hist_nutaucc_PTsq->GetBinContent(4)+hist_nutaucc_PTsq->GetBinContent(5)+hist_nutaucc_PTsq->GetBinContent(6) << endl;
    
    //cout << hist_nuecc_PTcube->GetBinContent(1)+hist_nuecc_PTcube->GetBinContent(2)+hist_nuecc_PTcube->GetBinContent(3)+hist_nuecc_PTcube->GetBinContent(4) << endl;
    //cout << hist_nutaucc_PTcube->GetBinContent(1)+hist_nutaucc_PTcube->GetBinContent(2)+hist_nutaucc_PTcube->GetBinContent(3)+hist_nutaucc_PTcube->GetBinContent(4) << endl;
    
    //cout << hist_nuecc_PTquartic->GetBinContent(1)+hist_nuecc_PTquartic->GetBinContent(2) << endl;
    //cout << hist_nutaucc_PTquartic->GetBinContent(1)+hist_nutaucc_PTquartic->GetBinContent(2) << endl;

    
    
}




void plotfunctions(){

}


void plot_polarcompare(){
    TFile *forward = new TFile("test_nutau_CC_tauproduction_test_polarization_forward.root");
    TFile *backward = new TFile("test_nutau_CC_tauproduction_test_polarization_backward.root");
    TFile *transverse = new TFile("test_nutau_CC_tauproduction_test_polarization_transverse.root");
    
    TH1F *hist_cc_El_forward = (TH1F*)forward->Get("hist_cc_El_true_WithSpin");
    TH1F *hist_cc_El_backward = (TH1F*)backward->Get("hist_cc_El_true_WithSpin");
    TH1F *hist_cc_El_transverse = (TH1F*)transverse->Get("hist_cc_El_true_WithSpin");
    
    TH1F *hist_cc_ElCosTheta_forward = (TH1F*)forward->Get("hist_cc_ElCosTheta_true_WithSpin");
    TH1F *hist_cc_ElCosTheta_backward = (TH1F*)backward->Get("hist_cc_ElCosTheta_true_WithSpin");
    TH1F *hist_cc_ElCosTheta_transverse = (TH1F*)transverse->Get("hist_cc_ElCosTheta_true_WithSpin");
    
    TH1F *hist_cc_ElCosTheta_above2GeV_forward = (TH1F*)forward->Get("hist_cc_ElCosTheta_above2GeV_true_WithSpin");
    TH1F *hist_cc_ElCosTheta_above2GeV_backward = (TH1F*)backward->Get("hist_cc_ElCosTheta_above2GeV_true_WithSpin");
    TH1F *hist_cc_ElCosTheta_above2GeV_transverse = (TH1F*)transverse->Get("hist_cc_ElCosTheta_above2GeV_true_WithSpin");
    
    
    hist_cc_El_forward->SetTitle("forward spin");
    hist_cc_El_backward->SetTitle("backward spin");
    hist_cc_El_backward->SetTitle("transverse spin");
    
    
    hist_cc_El_forward->SetFillColor(kBlue);
    hist_cc_El_forward->SetLineColor(kBlue);
    hist_cc_ElCosTheta_forward->SetFillColor(kBlue);
    hist_cc_ElCosTheta_forward->SetLineColor(kBlue);
    hist_cc_ElCosTheta_above2GeV_forward->SetFillColor(kBlue);
    hist_cc_ElCosTheta_above2GeV_forward->SetLineColor(kBlue);
    
    
    
    hist_cc_El_backward->SetFillColor(kGreen+3);
    hist_cc_El_backward->SetLineColor(kGreen+3);
    hist_cc_ElCosTheta_backward->SetFillColor(kGreen+3);
    hist_cc_ElCosTheta_backward->SetLineColor(kGreen+3);
    hist_cc_ElCosTheta_above2GeV_backward->SetFillColor(kGreen+3);
    hist_cc_ElCosTheta_above2GeV_backward->SetLineColor(kGreen+3);
    
    
    hist_cc_El_transverse->SetFillColor(kRed);
    hist_cc_El_transverse->SetLineColor(kRed);
    hist_cc_ElCosTheta_transverse->SetFillColor(kRed);
    hist_cc_ElCosTheta_transverse->SetLineColor(kRed);
    hist_cc_ElCosTheta_above2GeV_transverse->SetFillColor(kRed);
    hist_cc_ElCosTheta_above2GeV_transverse->SetLineColor(kRed);
    
    TCanvas *c1 = new TCanvas ("c1","c1",800,600);
    //c1->SetLogy();
    c1->cd();
    
    hist_cc_El_transverse->Draw();
    hist_cc_El_backward->Draw("same");
    hist_cc_El_forward->Draw("same");
    
    TLegend *leg = new TLegend(0.6,0.6,0.9,0.9);
    leg->AddEntry(hist_cc_El_transverse,"transverse","f");
    leg->AddEntry(hist_cc_El_backward,"backward","f");
    leg->AddEntry(hist_cc_El_forward,"forward","f");
    
    leg->Draw();
    
    TCanvas *c2 = new TCanvas ("c2","c2",800,600);
    //c1->SetLogy();
    c2->cd();
    
    hist_cc_ElCosTheta_forward->Draw();
    hist_cc_ElCosTheta_backward->Draw("same");
    hist_cc_ElCosTheta_transverse->Draw("same");
    
    TLegend *leg2 = new TLegend(0.6,0.6,0.9,0.9);
    leg2->AddEntry(hist_cc_ElCosTheta_transverse,"transverse","f");
    leg2->AddEntry(hist_cc_ElCosTheta_backward,"backward","f");
    leg2->AddEntry(hist_cc_ElCosTheta_forward,"forward","f");
    
    leg2->Draw();
    
    TCanvas *c3 = new TCanvas ("c3","c3",800,600);
    //c1->SetLogy();
    c3->cd();
    
    hist_cc_ElCosTheta_above2GeV_forward->Draw();
    hist_cc_ElCosTheta_above2GeV_backward->Draw("same");
    hist_cc_ElCosTheta_above2GeV_transverse->Draw("same");
    
    TLegend *leg3 = new TLegend(0.6,0.6,0.9,0.9);
    leg3->AddEntry(hist_cc_ElCosTheta_above2GeV_transverse,"transverse","f");
    leg3->AddEntry(hist_cc_ElCosTheta_above2GeV_backward,"backward","f");
    leg3->AddEntry(hist_cc_ElCosTheta_above2GeV_forward,"forward","f");
    
    leg3->Draw();
    
    
    
}



void stackaa(){
    TFile *fnue = new TFile("nuehistograms.root");
    TFile *fnumu = new TFile("numuhistograms.root");
    
    TFile *fnuenc = new TFile("test_nue_NC.root");
    TFile *fnumunc = new TFile("test_numu_NC.root");
    TFile *fnuebarnc = new TFile("test_nuebar_NC.root");
    TFile *fnumubarnc = new TFile("test_numubar_NC.root");
    
    TFile *fnutau = new TFile("test_nutau.root");
    
    TFile *fnuebar = new TFile("nuebarhistograms.root");
    TFile *fnumubar = new TFile("numubarhistograms.root");
    
    TH1F *hist_nuecc = (TH1F*)fnue->Get("hist_cc_reco");
    TH1F *hist_nuebarcc = (TH1F*)fnuebar->Get("hist_cc_reco");
    TH1F *hist_numucc = (TH1F*)fnumu->Get("hist_cc_reco");
    TH1F *hist_numubarcc = (TH1F*)fnumubar->Get("hist_cc_reco");
    
    TH1F *hist_nutaucc = (TH1F*)fnutau->Get("hist_cc_reco");
    
    
    //hist_mu->SetTitle("hist_mu_reco");
    //TH1F *hist_gamma = (TH1F*)f->Get("hist_gamma");
    TH1F *hist_nuenc = (TH1F*)fnuenc->Get("hist_nc_reco");
    TH1F *hist_nuebarnc = (TH1F*)fnuebarnc->Get("hist_nc_reco");
    TH1F *hist_numunc = (TH1F*)fnumunc->Get("hist_nc_reco");
    TH1F *hist_numubarnc = (TH1F*)fnumubarnc->Get("hist_nc_reco");
    
    
    
    hist_nutaucc->SetFillColor(kBlue+9);
    hist_nutaucc->SetLineColor(kBlue+9);
    
    hist_nuecc->SetFillColor(kAzure+9);
    hist_nuecc->SetLineColor(kAzure+9);
    hist_nuebarcc->SetFillColor(kAzure+9);
    hist_nuebarcc->SetLineColor(kAzure+9);
    
    hist_numucc->SetFillColor(kGreen+3);
    hist_numucc->SetLineColor(kGreen+3);
    hist_numubarcc->SetFillColor(kGreen+3);
    hist_numubarcc->SetLineColor(kGreen+3);
    
    hist_nuenc->SetFillColor(kBlue);
    hist_nuenc->SetLineColor(kBlue);
    hist_nuebarnc->SetFillColor(kBlue);
    hist_nuebarnc->SetLineColor(kBlue);
    
    hist_numunc->SetFillColor(kRed);
    hist_numunc->SetLineColor(kRed);
    hist_numubarnc->SetFillColor(kRed);
    hist_numubarnc->SetLineColor(kRed);
    
    
    // for 40 kt
    //double numu_scale = 3.9767;
    //double nue_scale = 0.0465198;
    //double numubar_scale = 0.130666;
    //double nuebar_scale = 0.0058103;
    
    //for 54.096 kt
    double numu_scale = 5.378;
    double nue_scale = 0.06291;
    double numubar_scale = 0.1767;
    double nuebar_scale = 0.007858;
    
    hist_nuecc->Scale(nue_scale);
    hist_nuenc->Scale(nue_scale);
    hist_nuebarcc->Scale(nuebar_scale);
    hist_nuebarnc->Scale(nuebar_scale);
    hist_numucc->Scale(numu_scale);
    hist_numunc->Scale(numu_scale);
    hist_numubarcc->Scale(numubar_scale);
    hist_numubarnc->Scale(numubar_scale);
    
    
    hist_nutaucc->Scale(numu_scale*.5);
    
    
    double nuecc_integral = 0.;
    double nuenc_integral = 0.;
    double numucc_integral = 0.;
    
    double nutaucc_integral = 0.;
    
    double numunc_integral = 0.;
    double nuebarcc_integral = 0.;
    double nuebarnc_integral = 0.;
    double numubarcc_integral = 0.;
    double numubarnc_integral = 0.;
    
    for (int i=3; i< 33; i++){
        //cout << hist_nuebackgroundsinglegamma->GetBinContent(i) << endl;
        //nuenccontent0580+=hist_nuebackgroundsinglegamma->GetBinContent(i);
        nuecc_integral+=hist_nuecc->GetBinContent(i);
        nuenc_integral+=hist_nuenc->GetBinContent(i);
        nuebarcc_integral+=hist_nuebarcc->GetBinContent(i);
        nuebarnc_integral+=hist_nuebarnc->GetBinContent(i);
        numucc_integral+=hist_numucc->GetBinContent(i);
        
        nutaucc_integral+=hist_nutaucc->GetBinContent(i);
        
        
        numunc_integral+=hist_numunc->GetBinContent(i);
        numubarcc_integral+=hist_numubarcc->GetBinContent(i);
        numubarnc_integral+=hist_numubarnc->GetBinContent(i);
        
        
    }
    cout << "integral nuecc: "<< nuecc_integral<<", nuenc: "<<nuenc_integral<<", nuebarcc: "<<nuebarcc_integral<<", nuebarnc: "<<nuebarnc_integral << endl;
    cout << "integral numucc: "<< numucc_integral<<", numunc: "<<numunc_integral<<", numubarcc: "<<numubarcc_integral<<", numubarnc: "<<numubarnc_integral << endl;
    cout << "integral nue+nubar cc: " << nuecc_integral+nuebarcc_integral << ", NC: "<<nuenc_integral+nuebarnc_integral+numunc_integral+numubarnc_integral << ", integral numu+numubar cc: "<< numucc_integral+numubarcc_integral << endl;
    
    cout << "integral nutaucc: "<< nutaucc_integral<< endl;
    
    THStack *hs = new THStack("hs","nue beam CC + NC + numu muon misid");
    
    hs->Add(hist_numucc,"histo");
    hs->Add(hist_nutaucc,"histo");
    
    hs->Add(hist_numubarcc,"histo");
    
    hs->Add(hist_nuenc,"histo");
    hs->Add(hist_nuebarnc,"histo");
    hs->Add(hist_numunc,"histo");
    hs->Add(hist_numubarnc,"histo");
    
    hs->Add(hist_nuecc,"histo");
    hs->Add(hist_nuebarcc,"histo");
    
    
    TCanvas *c1 = new TCanvas ("c1","c1",800,600);
    //c1->SetLogy();
    c1->cd();
    
    hs->Draw();
    hs->GetXaxis()->SetRangeUser(0.5,8.);
    
    
    TLegend *leg = new TLegend(0.6,0.6,0.9,0.9);
    leg->AddEntry(hist_nuecc,"hist_nuecc","f");
    leg->AddEntry(hist_nuebarcc,"hist_nuebarcc","f");
    leg->AddEntry(hist_numucc,"hist_numucc","f");
    
    leg->AddEntry(hist_nutaucc,"hist_nutaucc","f");
    
    leg->AddEntry(hist_numubarcc,"hist_numubarcc","f");
    leg->AddEntry(hist_nuenc,"hist_nuenc","f");
    leg->AddEntry(hist_nuebarnc,"hist_nuebarnc","f");
    leg->AddEntry(hist_numunc,"hist_numunc","f");
    leg->AddEntry(hist_numubarnc,"hist_numubarnc","f");
    leg->Draw();
    
    c1->SaveAs("stack_active_normalized_tau.pdf");
    
    
    
    
    /*hist_intrinsic->Scale(nue_scale);
     
     hist_nuebackgroundsinglegamma->Scale(nue_scale);
     
     hist_numubackgroundsinglegamma->Scale(numu_scale);
     hist_mu->Scale(numu_scale);
     
     
     //hist_pi0->SetFillColor(kGreen);
     //hist_pi0->SetLineColor(kGreen);
     
     TCanvas *c1 = new TCanvas ("c1","c1",800,600);
     //c1->SetLogy();
     c1->cd();
     
     //hist_backgroundsinglegamma->Scale();
     
     THStack *hs = new THStack("hs","nue beam CC + NC + numu muon misid");
     
     hs->Add(hist_mu,"histo");
     hs->Add(hist_nuebackgroundsinglegamma,"histo");
     hs->Add(hist_numubackgroundsinglegamma,"histo");
     hs->Add(hist_intrinsic,"histo");
     //hs->Add(hist_pi0,"histo");
     
     //hs->Setlogy();
     hs->Draw();
     
     TLegend *leg = new TLegend(0.6,0.6,0.9,0.9);
     leg->AddEntry(hist_intrinsic,"nue CC intrinsic","f");   // h1 and h2 are histogram pointers
     leg->AddEntry(hist_nuebackgroundsinglegamma,"nue NC","f");
     leg->AddEntry(hist_numubackgroundsinglegamma,"numu NC","f");
     leg->AddEntry(hist_mu,"numu CC muon misID","f");
     //leg->AddEntry(hist_pi0,"Pi0*0.05","f");
     leg->Draw();
     
     c1->SaveAs("stack_nue1.pdf");
     
     
     TCanvas *c2 = new TCanvas ("c2","c2",800,600);
     //c1->SetLogy();
     c2->cd();
     
     hs->GetXaxis()->SetRangeUser(0.5,8.);
     
     hs->GetXaxis()->SetTitle("E_{reco}[GeV]");
     hs->GetYaxis()->SetTitle("Events/0.25[GeV]");
     hs->Draw();
     
     //TLegend *leg = new TLegend(0.6,0.6,0.9,0.9);
     //leg->AddEntry(hist_intrinsic,"CC intrinsic","f");   // h1 and h2 are histogram pointers
     //leg->AddEntry(hist_nuebackgroundsinglegamma,"nue NC single gamma","f");
     //leg->AddEntry(hist_numubackgroundsinglegamma,"numu NC single gamma","f");
     //leg->AddEntry(hist_pi0,"Pi0*0.05","f");
     leg->Draw();
     
     
     
     c2->SaveAs("stack_nue2.pdf");*/
    
    
}
void stackbb(){
    TFile *fnue = new TFile("test_nue_CC_prob.root");
    TFile *fnumu = new TFile("test_numu_CC.root");
    
    TFile *fnuenc = new TFile("test_nue_NC.root");
    TFile *fnumunc = new TFile("test_numu_NC.root");
    TFile *fnuebarnc = new TFile("test_nuebar_NC.root");
    TFile *fnumubarnc = new TFile("test_numubar_NC.root");
    
    TFile *fnutau = new TFile("test_nutau_CC.root");
    TFile *fnutaubar = new TFile("test_nutaubar_CC.root");
    
    TFile *fnutaunc = new TFile("test_nutau_NC.root");
    TFile *fnutaubarnc = new TFile("test_nutaubar_NC.root");
    
    TFile *fnuebar = new TFile("test_nuebar_CC_prob.root");
    TFile *fnumubar = new TFile("test_numubar_CC.root");
    
    
    
    
    TH1F *hist_nuecc = (TH1F*)fnue->Get("hist_cc_reco");
    TH1F *hist_nuebarcc = (TH1F*)fnuebar->Get("hist_cc_reco");
    TH1F *hist_nuenc = (TH1F*)fnuenc->Get("hist_nc_reco");
    TH1F *hist_nuebarnc = (TH1F*)fnuebarnc->Get("hist_nc_reco");
    
    TH1F *hist_numucc = (TH1F*)fnumu->Get("hist_cc_reco");
    TH1F *hist_numubarcc = (TH1F*)fnumubar->Get("hist_cc_reco");
    TH1F *hist_numunc = (TH1F*)fnumunc->Get("hist_nc_reco");
    TH1F *hist_numubarnc = (TH1F*)fnumubarnc->Get("hist_nc_reco");
    
    
    TH1F *hist_nutaucc = (TH1F*)fnutau->Get("hist_cc_reco");
    TH1F *hist_nutaubarcc = (TH1F*)fnutaubar->Get("hist_cc_reco");
    TH1F *hist_nutaunc = (TH1F*)fnutaunc->Get("hist_nc_reco");
    TH1F *hist_nutaubarnc = (TH1F*)fnutaubarnc->Get("hist_nc_reco");
    
    
    
    
    hist_nutaucc->SetFillColor(kBlue+9);
    hist_nutaucc->SetLineColor(kBlue+9);
    hist_nutaubarcc->SetFillColor(kBlue+9);
    hist_nutaubarcc->SetLineColor(kBlue+9);
    
    hist_nuecc->SetFillColor(kAzure+9);
    hist_nuecc->SetLineColor(kAzure+9);
    hist_nuebarcc->SetFillColor(kAzure+9);
    hist_nuebarcc->SetLineColor(kAzure+9);
    
    hist_numucc->SetFillColor(kGreen+3);
    hist_numucc->SetLineColor(kGreen+3);
    hist_numubarcc->SetFillColor(kGreen+3);
    hist_numubarcc->SetLineColor(kGreen+3);
    
    hist_nuenc->SetFillColor(kBlue);
    hist_nuenc->SetLineColor(kBlue);
    hist_nuebarnc->SetFillColor(kBlue);
    hist_nuebarnc->SetLineColor(kBlue);
    
    hist_numunc->SetFillColor(kRed);
    hist_numunc->SetLineColor(kRed);
    hist_numubarnc->SetFillColor(kRed);
    hist_numubarnc->SetLineColor(kRed);
    
    
    hist_nutaunc->SetFillColor(kRed+9);
    hist_nutaunc->SetLineColor(kRed+9);
    hist_nutaubarnc->SetFillColor(kRed+9);
    hist_nutaubarnc->SetLineColor(kRed+9);
    
    
    
    //for 54.096 kt
    // double numu_scale = 5.378;
    // double nue_scale = 0.06291;
    // double numubar_scale = 0.1767;
    // double nuebar_scale = 0.007858;
    
    //for 55.641 kt
    double numu_scale = 5.53168;
    double nue_scale = 0.06471;
    double numubar_scale = 0.18176;
    double nuebar_scale = 0.008082;
    
    hist_nuecc->Scale(nue_scale);
    hist_nuenc->Scale(nue_scale);
    hist_nuebarcc->Scale(nuebar_scale);
    hist_nuebarnc->Scale(nuebar_scale);
    
    hist_numucc->Scale(numu_scale);
    hist_numunc->Scale(numu_scale);
    hist_numubarcc->Scale(numubar_scale);
    hist_numubarnc->Scale(numubar_scale);
    
    
    hist_nutaucc->Scale(numu_scale*.5);
    hist_nutaunc->Scale(numu_scale*.5);
    hist_nutaubarcc->Scale(numubar_scale*.5);
    hist_nutaubarnc->Scale(numubar_scale*.5);
    
    
    double nuecc_integral = 0.;
    double nuenc_integral = 0.;
    double nuebarcc_integral = 0.;
    double nuebarnc_integral = 0.;
    
    double numucc_integral = 0.;
    double numunc_integral = 0.;
    double numubarcc_integral = 0.;
    double numubarnc_integral = 0.;
    
    double nutaucc_integral = 0.;
    double nutaubarcc_integral = 0.;
    double nutaunc_integral = 0.;
    double nutaubarnc_integral = 0.;
    
    for (int i=3; i< 33; i++){
        //cout << hist_nuebackgroundsinglegamma->GetBinContent(i) << endl;
        //nuenccontent0580+=hist_nuebackgroundsinglegamma->GetBinContent(i);
        nuecc_integral+=hist_nuecc->GetBinContent(i);
        nuenc_integral+=hist_nuenc->GetBinContent(i);
        nuebarcc_integral+=hist_nuebarcc->GetBinContent(i);
        nuebarnc_integral+=hist_nuebarnc->GetBinContent(i);
        
        numucc_integral+=hist_numucc->GetBinContent(i);
        numunc_integral+=hist_numunc->GetBinContent(i);
        numubarcc_integral+=hist_numubarcc->GetBinContent(i);
        numubarnc_integral+=hist_numubarnc->GetBinContent(i);
        
        nutaucc_integral+=hist_nutaucc->GetBinContent(i);
        nutaunc_integral+=hist_nutaunc->GetBinContent(i);
        nutaubarcc_integral+=hist_nutaubarcc->GetBinContent(i);
        nutaubarnc_integral+=hist_nutaubarnc->GetBinContent(i);
        
        
        
        
        
    }
    cout << "integral nuecc: "<< nuecc_integral<<", nuenc: "<<nuenc_integral<<", nuebarcc: "<<nuebarcc_integral<<", nuebarnc: "<<nuebarnc_integral << endl;
    cout << "integral numucc: "<< numucc_integral<<", numunc: "<<numunc_integral<<", numubarcc: "<<numubarcc_integral<<", numubarnc: "<<numubarnc_integral << endl;
    cout << "integral nue+nubar cc: " << nuecc_integral+nuebarcc_integral << ", NC: "<<nuenc_integral+nuebarnc_integral+numunc_integral+numubarnc_integral << ", integral numu+numubar cc: "<< numucc_integral+numubarcc_integral << endl;
    
    cout << "integral nutaucc: "<< nutaucc_integral<< ", nutaunc : " << nutaunc_integral<< ", nutaubarcc: "<< nutaubarcc_integral<< ", nutaubarnc : " << nutaubarnc_integral<< endl;
    
    THStack *hs = new THStack("hs","nue beam CC + NC + numu muon misid");
    
    hs->Add(hist_numucc,"histo");
    //hs->Add(hist_nutaucc,"histo");
    
    hs->Add(hist_numubarcc,"histo");
    
    hs->Add(hist_nuenc,"histo");
    hs->Add(hist_nuebarnc,"histo");
    hs->Add(hist_numunc,"histo");
    hs->Add(hist_numubarnc,"histo");
    
    hs->Add(hist_nuecc,"histo");
    hs->Add(hist_nuebarcc,"histo");
    
    
    TCanvas *c1 = new TCanvas ("c1","c1",800,600);
    //c1->SetLogy();
    c1->cd();
    
    hs->Draw();
    hs->GetXaxis()->SetRangeUser(0.5,8.);
    
    
    TLegend *leg = new TLegend(0.6,0.6,0.9,0.9);
    leg->AddEntry(hist_nuecc,"hist_nuecc","f");
    leg->AddEntry(hist_nuebarcc,"hist_nuebarcc","f");
    leg->AddEntry(hist_numucc,"hist_numucc","f");
    
    leg->AddEntry(hist_nutaucc,"hist_nutaucc","f");
    
    leg->AddEntry(hist_numubarcc,"hist_numubarcc","f");
    leg->AddEntry(hist_nuenc,"hist_nuenc","f");
    leg->AddEntry(hist_nuebarnc,"hist_nuebarnc","f");
    leg->AddEntry(hist_numunc,"hist_numunc","f");
    leg->AddEntry(hist_numubarnc,"hist_numubarnc","f");
    leg->Draw();
    
    c1->SaveAs("stack_active_normalized_tau_bb.pdf");
    
    
    
}





void stack_study(){
    
    TFile *fnue = new TFile("../gst0to40/out/gntp.0.nue20k_gst_study_out.root");
    TFile *fnuebar = new TFile("../gst0to40/out/gntp.0.nuebar20k_gst_study_out.root");
    TFile *fnumu = new TFile("../gst0to40/out/gntp.0.numu50k_gst_study_out.root");
    TFile *fnumubar = new TFile("../gst0to40/out/gntp.0.numubar10k_gst_study_out.root");
    TFile *fnutau = new TFile("../gst0to40/out/gntp.0.nutau20k_gst_study_out.root");
    TFile *fnutaubar = new TFile("../gst0to40/out/gntp.0.nutaubar20k_gst_study_out.root");
    TFile *fnuesig = new TFile("../gst0to40/out/gntp.0.numuflux_nuebeam50k_gst_study_out.root");
    TFile *fnuebarsig = new TFile("../gst0to40/out/gntp.0.numubarflux_nuebarbeam10k_gst_study_out.root");
    
    //TFile *fnue = new TFile("../gst0to40/out/gntp.0.nue20k_gst_NC_study_out.root");
    //TFile *fnuebar = new TFile("../gst0to40/out/gntp.0.nuebar20k_gst_NC_study_out.root");
    //TFile *fnumu = new TFile("../gst0to40/out/gntp.0.numu50k_gst_NC_study_out.root");
    //TFile *fnumubar = new TFile("../gst0to40/out/gntp.0.numubar10k_NC_gst_study_out.root");
    //TFile *fnutau = new TFile("../gst0to40/out/gntp.0.nutau20k_gst_NC_study_out.root");
    //TFile *fnutaubar = new TFile("../gst0to40/out/gntp.0.nutaubar20k_NC_gst_study_out.root");

    
    
    TFile *fnuenc = new TFile("../gst0to40/out/gntp.0.nue20k_gst_NC_study_out.root");
    TFile *fnuebarnc = new TFile("../gst0to40/out/gntp.0.nuebar20k_gst_NC_study_out.root");
    TFile *fnumunc = new TFile("../gst0to40/out/gntp.0.numu50k_gst_NC_study_out.root");
    TFile *fnumubarnc = new TFile("../gst0to40/out/gntp.0.numubar10k_gst_NC_study_out.root");
    TFile *fnutaunc = new TFile("../gst0to40/out/gntp.0.nutau20k_gst_NC_study_out.root");
    TFile *fnutaubarnc = new TFile("../gst0to40/out/gntp.0.nutaubar20k_gst_NC_study_out.root");
   
    
    
    
    
    TH1F *hist_nuecc = (TH1F*)fnue->Get("hist_cc_reco");
    TH1F *hist_nuebarcc = (TH1F*)fnuebar->Get("hist_cc_reco");
    TH1F *hist_nuenc = (TH1F*)fnuenc->Get("hist_nc_reco");
    TH1F *hist_nuebarnc = (TH1F*)fnuebarnc->Get("hist_nc_reco");
    
    TH1F *hist_numucc = (TH1F*)fnumu->Get("hist_cc_reco");
    TH1F *hist_numubarcc = (TH1F*)fnumubar->Get("hist_cc_reco");
    TH1F *hist_numunc = (TH1F*)fnumunc->Get("hist_nc_reco");
    TH1F *hist_numubarnc = (TH1F*)fnumubarnc->Get("hist_nc_reco");
    
    
    TH1F *hist_nutaucc = (TH1F*)fnutau->Get("hist_cc_reco");
    TH1F *hist_nutaubarcc = (TH1F*)fnutaubar->Get("hist_cc_reco");
    
    //TH1F *hist_nutaucc = (TH1F*)fnutau->Get("hist_cc_reco_pt05cut");
    //TH1F *hist_nutaubarcc = (TH1F*)fnutaubar->Get("hist_cc_reco_pt05cut");
    
    TH1F *hist_nutaunc = (TH1F*)fnutaunc->Get("hist_nc_reco");
    TH1F *hist_nutaubarnc = (TH1F*)fnutaubarnc->Get("hist_nc_reco");
    
    TH1F *hist_nuesig = (TH1F*)fnuesig->Get("hist_cc_reco");
    TH1F *hist_nuebarsig = (TH1F*)fnuebarsig->Get("hist_cc_reco");
    
    
    //hist_nuesig->SetFillColor(kBlack);
    hist_nuesig->SetLineStyle(2);
    hist_nuesig->SetLineWidth(2);
    hist_nuesig->SetLineColor(kBlack);
    
    hist_nuebarsig->SetFillColor(kBlack);
    hist_nuebarsig->SetLineColor(kBlack);

    
    
    //hist_nuecc->SetFillColor(kBlue);
    
    
    hist_nuecc->SetLineColor(kBlue);
    hist_nuecc->SetLineStyle(2);
    hist_nuecc->SetLineWidth(2);
    
    hist_nuebarcc->SetFillColor(kAzure+9);
    hist_nuebarcc->SetLineColor(kAzure+9);
    
    //hist_numucc->SetFillColor(kGreen+3);
    hist_numucc->SetLineColor(kGreen+3);
    hist_numucc->SetLineStyle(2);
    hist_numucc->SetLineWidth(2);
    //hist_numubarcc->SetFillColor(kGreen+3);
    //hist_numubarcc->SetLineColor(kGreen+3);
    
    
    //hist_nuenc->SetFillColor(kRed);
    
    hist_nuenc->SetLineColor(kRed);
    hist_nuenc->SetLineStyle(2);
    hist_nuenc->SetLineWidth(2);
    
    //hist_nuebarnc->SetFillColor(kRed);
    //hist_nuebarnc->SetLineColor(kRed);
    
    hist_numunc->SetFillColor(kRed);
    hist_numunc->SetLineColor(kRed);
    hist_numubarnc->SetFillColor(kRed);
    hist_numubarnc->SetLineColor(kRed);
    
    
    hist_nutaunc->SetFillColor(kRed);
    hist_nutaunc->SetLineColor(kRed);
    
    hist_nutaubarnc->SetFillColor(kRed);
    hist_nutaubarnc->SetLineColor(kRed);
    
   // hist_nutaucc->SetFillColor(kGray+2);
    hist_nutaucc->SetLineColor(kGray+2);
    hist_nutaucc->SetLineStyle(2);
    hist_nutaucc->SetLineWidth(2);
   // hist_nutaubarcc->SetFillColor(kGray+2);
   // hist_nutaubarcc->SetLineColor(kGray+2);

    
    
    
    
    //for 54.096 kt
    // double numu_scale = 5.378;
    // double nue_scale = 0.06291;
    // double numubar_scale = 0.1767;
    // double nuebar_scale = 0.007858;
    
    //for 55.641 kt
    double numu_scale = 5.53168;
    double nue_scale = 0.06471;
    double numubar_scale = 0.18176;
    double nuebar_scale = 0.008082;
    
    double sig_scale = numu_scale/5.;
    
    
    
    
    hist_nuecc->Scale(nue_scale*0.5);
    hist_nuenc->Scale(nue_scale*0.5);
    
    hist_nuebarcc->Scale(nuebar_scale*0.5);
    hist_nuebarnc->Scale(nuebar_scale*0.5);
    
    hist_numucc->Scale(numu_scale*0.2);
    hist_numunc->Scale(numu_scale*0.2);
    
    hist_numubarcc->Scale(numubar_scale);
    hist_numubarnc->Scale(numubar_scale);
    
    
    hist_nutaucc->Scale(numu_scale*.5);
    hist_nutaunc->Scale(numu_scale*.5);
    
    hist_nutaubarcc->Scale(numubar_scale*.5);
    hist_nutaubarnc->Scale(numubar_scale*.5);
    
    hist_nuesig->Scale(sig_scale);
    hist_nuebarsig->Scale(numubar_scale);
    
    
    double nuecc_integral = 0.;
    double nuenc_integral = 0.;
    double nuebarcc_integral = 0.;
    double nuebarnc_integral = 0.;
    
    double numucc_integral = 0.;
    double numunc_integral = 0.;
    double numubarcc_integral = 0.;
    double numubarnc_integral = 0.;
    
    double nutaucc_integral = 0.;
    double nutaubarcc_integral = 0.;
    double nutaunc_integral = 0.;
    double nutaubarnc_integral = 0.;
    
    double sig_integral = 0.;
    double sigbar_integral = 0.;
    
    
    double bkg_nutau_effi_ratio[40] = {
        0.185177101265301,
        0.0569097745791012,
        0.142302952010498,
        0.271548091251609,
        0.365215909373156,
        0.486457828608005,
        0.463130167051557,
        0.49821981228389,
        0.508886165054105,
        0.523416189301326,
        0.560997511737675,
        0.518695827424634,
        0.545142155498725,
        0.530004233403559,
        0.562003158877504,
        0.509813049542501,
        0.509139483299998,
        0.486598707120586,
        0.48577769000742,
        0.466012371584107,
        0.442657261327112,
        0.444240085377149,
        0.483963356356026,
        0.42227539095115,
        0.43934426885497,
        0.375047996507843,
        0.407865995461402,
        0.363710583961687,
        0.390076613505933,
        0.357954552368202,
        0.396268301879026,
        0.387585758062179,
        0.385305674134454,
        0.306145466796258,
        0.380393505599226,
        0.374360284108998,
        0.413174986327187,
        0.372564212374476,
        0.349898958801367,
        0.334581841047651};
    
    double sig_nue_effi_ratio[40] = {
        0.521657779382294,
        0.942923207519666,
        0.889683700153627,
        0.911539031741566,
        0.966729129379551,
        0.978183764359334,
        0.98112755900117,
        0.986379904525923,
        0.988931224173579,
        0.986405364710914,
        0.98894157766952,
        0.991419120944965,
        0.986495120514319,
        0.983898143443015,
        0.983880964903094,
        0.976496791989119,
        0.974032690840062,
        0.973973405815723,
        0.968892013833642,
        0.966360898273846,
        0.956858588056045,
        0.949320337116643,
        0.954163860342054,
        0.949450889612176,
        0.947072966044237,
        0.936950797360869,
        0.922385043903225,
        0.921928862645288,
        0.912912829355829,
        0.912636981103151,
        0.910803981014623,
        0.90542620165338,
        0.912877628360244,
        0.90726553947602,
        0.892808253437137,
        0.905400646664422,
        0.900718823866867,
        0.920413369689987,
        0.898669628730067,
        0.906748296678402};
    
    double sig_nuebar_effi_ratio[40] = {
        0.467598209176704,
        0.905333140336751,
        0.94675699636117,
        0.975334171104503,
        0.986367401602673,
        0.989121705406261,
        0.991702315457223,
        0.991740641269727,
        0.989318036530991,
        0.989350893606755,
        0.989258242835556,
        0.989291376349116,
        0.984408364788412,
        0.984596049045815,
        0.984490472875525,
        0.97483447993455,
        0.972367149824261,
        0.972380639343171,
        0.965231268271163,
        0.965173266414924,
        0.962709979925019,
        0.960461471728456,
        0.958224330051069,
        0.962427877418409,
        0.960574039401596,
        0.953653817437637,
        0.948517286326605,
        0.948318943897705,
        0.948257198898437,
        0.94609077251976,
        0.939081541068118,
        0.93872434222039,
        0.945899496457523,
        0.941256104576436,
        0.941200830480307,
        0.945837893896248,
        0.94329295786115,
        0.948085175550173,
        0.943439865156784,
        0.934131393019274};
    
    
    
    double bkg_nue_effi_ratio[40] = {
        0.57637154779581,
        0.940354238004384,
        0.904928491193115,
        0.929295362547548,
        0.975467243856832,
        0.986160849280521,
        0.987226465066376,
        0.99084924663688,
        0.98828589943228,
        0.98932976182368,
        0.991886710089514,
        0.986304888140026,
        0.987851027174437,
        0.980294153288184,
        0.976846530015149,
        0.965178931563289,
        0.956401597513423,
        0.96099507727344,
        0.956426650123147,
        0.95104734649484,
        0.942147585689838,
        0.932628551161413,
        0.93264084542889,
        0.930816391108137,
        0.929545470219414,
        0.922777892086839,
        0.910427377732449,
        0.905238785635381,
        0.897634373665486,
        0.912382350939481,
        0.904059130071588,
        0.899374361362336,
        0.902402387962717,
        0.909700773718198,
        0.896579607415514,
        0.906522916845458,
        0.909459581049356,
        0.916525085592965,
        0.893978041480055,
        0.89196689153647
    };
    
    double bkg_NC_effi_ratio[40] = {
        0.18047465517046,
        0.101151907725579,
        0.232938787430535,
        0.337010457724242,
        0.648029150535554,
        0.501602340213943,
        0.645180823288289,
        0.931730203198994,
        0.761877907924052,
        0.776838758789639,
        0.372765819839769,
        0.84718572080561,
        0.65233019873389,
        0.593016077460186,
        0.507954119678036,
        0.6212656967605,
        0.57345478472973,
        0.620662846465288,
        0.507703967080246,
        0.559117027037734,
        0.348976890933601,
        0.349160950475233,
        0.372580987352924,
        0.232986701371153,
        0.279438769960226,
        0.207053277186504,
        0.155284049293359,
        0.310483034741583,
        0.20695622854475,
        0.232985115308552,
        0.310737887138773,
        0.132925008629442,
        0.0932994577334251,
        0.186445714386841,
        0.349048037607936,
        0.399259898770444,
        0.18651271829387,
        0.310834997254636,
        0.103361131972826,
        0.116469686089845};
        
        
    
    
    bool apply_effi = true;

    
    if(apply_effi){
        for (int k = 0; k<40; k++){
            hist_nutaucc->SetBinContent(k+1,hist_nutaucc->GetBinContent(k+1)*bkg_nutau_effi_ratio[k]);
            hist_nutaubarcc->SetBinContent(k+1,hist_nutaubarcc->GetBinContent(k+1)*bkg_nutau_effi_ratio[k]);
            
            hist_nuesig->SetBinContent(k+1,hist_nuesig->GetBinContent(k+1)*sig_nue_effi_ratio[k]);
            hist_nuebarsig->SetBinContent(k+1,hist_nuebarsig->GetBinContent(k+1)*sig_nuebar_effi_ratio[k]);
            
            hist_nuecc->SetBinContent(k+1,hist_nuecc->GetBinContent(k+1)*bkg_nue_effi_ratio[k]);
            hist_nuebarcc->SetBinContent(k+1,hist_nuebarcc->GetBinContent(k+1)*bkg_nue_effi_ratio[k]);
            
            hist_nuenc->SetBinContent(k+1,hist_nuenc->GetBinContent(k+1)*bkg_NC_effi_ratio[k]);
            hist_nuebarnc->SetBinContent(k+1,hist_nuebarnc->GetBinContent(k+1)*bkg_NC_effi_ratio[k]);
            hist_numunc->SetBinContent(k+1,hist_numunc->GetBinContent(k+1)*bkg_NC_effi_ratio[k]);
            hist_numubarnc->SetBinContent(k+1,hist_numubarnc->GetBinContent(k+1)*bkg_NC_effi_ratio[k]);
            hist_nutaunc->SetBinContent(k+1,hist_nutaunc->GetBinContent(k+1)*bkg_NC_effi_ratio[k]);
            hist_nutaubarnc->SetBinContent(k+1,hist_nutaubarnc->GetBinContent(k+1)*bkg_NC_effi_ratio[k]);
            
            
        }
        
    }
    
    for (int i=3; i< 33; i++){
        sig_integral+=hist_nuesig->GetBinContent(i);
        sigbar_integral+=hist_nuebarsig->GetBinContent(i);
        //cout << hist_nuebackgroundsinglegamma->GetBinContent(i) << endl;
        //nuenccontent0580+=hist_nuebackgroundsinglegamma->GetBinContent(i);
        nuecc_integral+=hist_nuecc->GetBinContent(i);
        nuenc_integral+=hist_nuenc->GetBinContent(i);
        nuebarcc_integral+=hist_nuebarcc->GetBinContent(i);
        nuebarnc_integral+=hist_nuebarnc->GetBinContent(i);
        
        numucc_integral+=hist_numucc->GetBinContent(i);
        numunc_integral+=hist_numunc->GetBinContent(i);
        numubarcc_integral+=hist_numubarcc->GetBinContent(i);
        numubarnc_integral+=hist_numubarnc->GetBinContent(i);
        
        nutaucc_integral+=hist_nutaucc->GetBinContent(i);
        nutaunc_integral+=hist_nutaunc->GetBinContent(i);
        nutaubarcc_integral+=hist_nutaubarcc->GetBinContent(i);
        nutaubarnc_integral+=hist_nutaubarnc->GetBinContent(i);
        
    }

    cout << "integral sig nue : " << sig_integral << ", nuebar : "<< sigbar_integral<< endl;
    cout << "integral nuecc: "<< nuecc_integral<<",  nuebarcc: "<<nuebarcc_integral<< endl;
    cout << "integral numucc: "<< numucc_integral<<", numubarcc: "<<numubarcc_integral<< endl;
    cout << "integral nue+nubar cc: " << nuecc_integral+nuebarcc_integral << ", integral numu+numubar cc: "<< numucc_integral+numubarcc_integral << endl;
    
    cout << "integral nutaucc: "<< nutaucc_integral<< ", nutaubarcc: "<< nutaubarcc_integral<< endl;
    cout << "integral nuenc : " << nuenc_integral << ", nuebarnc: "<<nuebarnc_integral<<endl;
    cout << "integral numunc : " << numunc_integral << ", numubarnc: "<<numubarnc_integral<<endl;
    cout << "integral nutaunc : " << nutaunc_integral << ", nutaubarnc: "<<nutaubarnc_integral<<endl;
    
    
    hist_nuesig->Add(hist_nuebarsig);
    hist_nuecc->Add(hist_nuebarcc);
    
    hist_nuenc->Add(hist_nuebarnc);
    hist_nuenc->Add(hist_numunc);
    hist_nuenc->Add(hist_numubarnc);
    hist_nuenc->Add(hist_nutaunc);
    hist_nuenc->Add(hist_nutaubarnc);
    
    hist_nutaucc->Add(hist_nutaubarcc);
    
    hist_numucc->Add(hist_numubarcc);
    
    
    THStack *hs = new THStack("hs","DUNE nue appearance");
    
    
    hs->Add(hist_numucc,"histo");
    //hs->Add(hist_numubarcc,"histo");
    hs->Add(hist_nutaucc,"histo");
    //hs->Add(hist_nutaubarcc,"histo");
    //hs->Add(hist_nutaunc,"histo");
    //hs->Add(hist_nutaubarnc,"histo");
    hs->Add(hist_nuenc,"histo");
    //hs->Add(hist_nuebarnc,"histo");
    //hs->Add(hist_numunc,"histo");
    //hs->Add(hist_numubarnc,"histo");
    hs->Add(hist_nuecc,"hist");
    //hs->Add(hist_nuebarcc,"histo");
    hs -> Add(hist_nuesig,"hist");
    //hs -> Add(hist_nuebarsig,"hist");
   
 
    /*****************mark adds to TFile***************/
    TFile *fm = new TFile("DUNE_bf.root","RECREATE");
    fm->cd();
    TH1D * signal = (TH1D*)hist_nuesig->Clone("nu_dune_elike_signal");
    TH1D * intrinsic = (TH1D*)hist_nuecc->Clone("nu_dune_elike_intrinsic");
    TH1D * ncmisid = (TH1D*)hist_nuenc->Clone("nu_dune_elike_ncmisid");
    TH1D * mumisid = (TH1D*)hist_numucc->Clone("nu_dune_elike_mumisid");
    TH1D * taumisid = (TH1D*)hist_nutaucc->Clone("nu_dune_elike_taumisid");
    signal->Write();
    intrinsic->Write();
    ncmisid->Write();
    mumisid->Write();
    taumisid->Write();
    fm->Close();



 
    TCanvas *c1 = new TCanvas ("c1","c1",600,600);
    //c1->SetLogy();
    c1->cd();
    
    hs->Draw();
    hs->GetXaxis()->SetRangeUser(0.5,8.);
    
    
    TLegend *leg = new TLegend(0.6,0.6,0.9,0.9);
    leg->AddEntry(hist_nuesig,"nue signal","l");
    leg->AddEntry(hist_nuecc,"nue+nuebar CC","l");
    leg->AddEntry(hist_nuenc,"NC","l");
    leg->AddEntry(hist_nutaucc,"nutau+nutaubar CC","l");
    

    //leg->AddEntry(hist_nutaubarcc,"hist_nutaubarcc","l");

    //    leg->AddEntry(hist_nuebarcc,"hist_nuebarcc","l");
    
    
       // leg->AddEntry(hist_nutaubarnc,"hist_nutaubarnc","l");
    //leg->AddEntry(hist_nutaunc,"hist_nutaunc","l");
    leg->AddEntry(hist_numucc,"numu+numubar CC","l");

    
    //leg->AddEntry(hist_numubarcc,"hist_numubarcc","l");
        //leg->AddEntry(hist_nuebarnc,"hist_nuebarnc","l");
    //leg->AddEntry(hist_numunc,"hist_numunc","l");
    //leg->AddEntry(hist_numubarnc,"hist_numubarnc","l");
    leg->Draw();
    
    c1->SaveAs("stack_active_normalized_tau_cc.pdf");
    
    
    
}



void stack_mu(){
    
    //TFile *fnumusig = new TFile("../gst0to40/out/gntp.0.numu50k_gst_genie_mu_out.root");
    //TFile *fnumubarcc = new TFile("../gst0to40/out/gntp.0.numubar10k_gst_genie_mu_out.root");
    
    TFile *fnumusig = new TFile("../gst0to40/out_mu_200cm/gntp.0.numu50k_gst_genie_mu_out.root");
    TFile *fnumubarcc = new TFile("../gst0to40/out_mu_200cm/gntp.0.numubar10k_gst_genie_mu_out.root");
    
    TFile *fnutaucc = new TFile("../gst0to40/out/gntp.0.nutau20k_gst_genie_mu_out.root");
    TFile *fnutaubarcc = new TFile("../gst0to40/out/gntp.0.nutaubar20k_gst_genie_mu_out.root");
    
    TH1F *hist_numusig = (TH1F*)fnumusig->Get("hist_cc_reco");
    TH1F *hist_numunc = (TH1F*)fnumusig->Get("hist_nc_reco");
    TH1F *hist_numusigtrue = (TH1F*)fnumusig->Get("hist_cc_reco_true");
    TH1F *hist_nutaucc = (TH1F*)fnutaucc->Get("hist_cc_reco");
    TH1F *hist_nutaubarcc = (TH1F*)fnutaubarcc->Get("hist_cc_reco");
    
    TH1F *hist_numubar = (TH1F*)fnumubarcc->Get("hist_cc_reco");
   // TH1F *hist_numusigtrue = (TH1F*)fnumusig->Get("hist_cc_reco_true");

    
    //hist_numusig->SetFillColor(kPink);
    hist_numusig->SetLineColor(kBlack);
    hist_numusig->SetLineStyle(2);
    hist_numusig->SetLineWidth(2);
    
    hist_numubar->SetLineColor(kGreen+1);
    hist_numubar->SetLineStyle(2);
    hist_numubar->SetLineWidth(2);
    
    hist_numunc->SetLineColor(kRed);
    hist_numunc->SetLineStyle(2);
    hist_numunc->SetLineWidth(2);

    
    hist_numusigtrue->SetLineColor(kBlack);
    hist_numusigtrue->SetLineStyle(2);
    hist_numusigtrue->SetLineWidth(2);


    
    //hist_nutaucc->SetFillColor(kGray+2);
    hist_nutaucc->SetLineColor(kGray+2);
    hist_nutaucc->SetLineStyle(2);
    hist_nutaucc->SetLineWidth(2);
    //hist_nutaubarcc->SetFillColor(kGreen+3);
    hist_nutaubarcc->SetLineColor(kGray+2);
    hist_nutaubarcc->SetLineStyle(2);
    hist_nutaubarcc->SetLineWidth(2);

    
    //for 54.096 kt
    // double numu_scale = 5.378;
    // double nue_scale = 0.06291;
    // double numubar_scale = 0.1767;
    // double nuebar_scale = 0.007858;
    
    //for 55.641 kt
    double numu_scale = 5.53168;
    double nue_scale = 0.06471;
    double numubar_scale = 0.18176;
    double nuebar_scale = 0.008082;
    
    double sig_scale = numu_scale/5.;
    

    
    hist_numusig->Scale(sig_scale);
    hist_numunc->Scale(sig_scale);
    hist_numusigtrue->Scale(sig_scale);

    hist_nutaucc->Scale(numu_scale*0.5);
    hist_nutaubarcc->Scale(numubar_scale*0.5);
    hist_numubar->Scale(numubar_scale);
    
    double sig_integral = 0.;
    double sig_integral_true = 0.;
    double nutaucc_integral = 0.;
    double nutaubarcc_integral = 0.;
    double numubarcc_integral = 0.;
    double numunc_integral = 0.;
    
    
    double bkg_nutau_effi_ratio[40] = {1.,
        1.,
        0.325581395348837,
        0.397050482132729,
        0.425581395348838,
        0.455324357405141,
        0.48062015503876,
        0.475631951466128,
        0.53359173126615,
        0.516015796401931,
        0.552501761804088,
        0.544920440636476,
        0.52093023255814,
        0.551573187414501,
        0.55098389982111,
        0.625116279069768,
        0.605306256141501,
        0.622851365015167,
        0.674091057975762,
        0.721288014311271,
        0.824127906976744,
        0.787226657410621,
        0.837209302325581,
        0.901610017889087,
        0.900436046511628,
        0.887949260042283,
        0.945736434108528,
        0.913728432108027,
        0.850712678169541,
        0.877414268821443,
        0.976744186046511,
        0.913728432108027,
        0.960990247561891,
        0.850712678169541,
        0.801804928844152,
        0.871556350626118,
        0.96124031007752,
        0.993023255813955,
        1.,
        0.837209302325581};
    
    //double sig_nue_effi_ratio[40] = {
    //    };
    
    //double sig_nuebar_effi_ratio[40] = {
    //    };
    
    
    
    //double bkg_nue_effi_ratio[40] = {
    //        };
    
    double bkg_NC_effi_ratio[40] = {
    1.,
    1.,
    0.48837209302324,
    0.488372093023265,
    0.651162790697676,
    0.651162790697676,
    0.488372093023265,
    0.488372093023265,
    0.488372093023265,
    0.390697674418608,
    0.325581395348838,
    0.781395348837217,
    0.488372093023257,
    0.390697674418608,
    0.390697674418608,
    0.488372093023257,
    0.488372093023257,
    0.418604651162796,
    0.418604651162796,
    0.418604651162796,
    0.732558139534889,
    0.488372093023257,
    0.418604651162796,
    0.418604651162796,
    0.488372093023259,
    0.813953488372078,
    0.651162790697683,
    0.434108527131789,
    0.293023255813956,
    0.390697674418608,
    0.48837209302325,
    0.558139534883728,
    0.558139534883728,
    0.355179704016919,
    0.325581395348842,
    0.325581395348842,
    0.177589852008459,
    0.355179704016919,
    0.390697674418608,
    0.217054263565894
        };
    
    
    
    
    bool apply_effi = true;
    
    
    if(apply_effi){
        for (int k = 0; k<40; k++){
            hist_nutaucc->SetBinContent(k+1,hist_nutaucc->GetBinContent(k+1)*bkg_nutau_effi_ratio[k]);
            //hist_nutaubarcc->SetBinContent(k+1,hist_nutaubarcc->GetBinContent(k+1)*bkg_nutau_effi_ratio[k]);
            
            //hist_nuesig->SetBinContent(k+1,hist_nuesig->GetBinContent(k+1)*sig_nue_effi_ratio[k]);
            //hist_nuebarsig->SetBinContent(k+1,hist_nuebarsig->GetBinContent(k+1)*sig_nuebar_effi_ratio[k]);
            
            //hist_nuecc->SetBinContent(k+1,hist_nuecc->GetBinContent(k+1)*bkg_nue_effi_ratio[k]);
            //hist_nuebarcc->SetBinContent(k+1,hist_nuebarcc->GetBinContent(k+1)*bkg_nue_effi_ratio[k]);
            
            hist_numunc->SetBinContent(k+1,hist_numunc->GetBinContent(k+1)*bkg_NC_effi_ratio[k]);
            //hist_nuebarnc->SetBinContent(k+1,hist_nuebarnc->GetBinContent(k+1)*bkg_NC_effi_ratio[k]);
            //hist_numunc->SetBinContent(k+1,hist_numunc->GetBinContent(k+1)*bkg_NC_effi_ratio[k]);
            //hist_numubarnc->SetBinContent(k+1,hist_numubarnc->GetBinContent(k+1)*bkg_NC_effi_ratio[k]);
            //hist_nutaunc->SetBinContent(k+1,hist_nutaunc->GetBinContent(k+1)*bkg_NC_effi_ratio[k]);
            //hist_nutaubarnc->SetBinContent(k+1,hist_nutaubarnc->GetBinContent(k+1)*bkg_NC_effi_ratio[k]);
            
            
            
            
            
        }
        
    }

    
    hist_nutaucc->Add(hist_nutaubarcc);
    
    for (int i=3; i< 81; i++){
        numunc_integral+=hist_numunc->GetBinContent(i);
        sig_integral+=hist_numusig->GetBinContent(i);
        sig_integral_true+=hist_numusigtrue->GetBinContent(i);
        nutaucc_integral+=hist_nutaucc->GetBinContent(i);
        //nutaubarcc_integral+=hist_nutaubarcc->GetBinContent(i);
        numubarcc_integral+=hist_numubar->GetBinContent(i);
    }
    cout << "integral sig numu : " << sig_integral<< ", true : "<<sig_integral_true << endl;
    cout << "integral bkg nutaucc : " << nutaucc_integral << endl;
    //cout << "integral bkg nutaubarcc : " << nutaubarcc_integral << endl;
    cout << "integral bkg numubar : " << numubarcc_integral << endl;
    cout << "integral NC numu : "<< numunc_integral<<endl;



 
    /*****************mark adds to TFile***************/
    TFile *fm = new TFile("DUNE_bf.root","UPDATE");
    fm->cd();
    TH1D * signal = (TH1D*)hist_numusig->Clone("nu_dune_mulike_signal");
    TH1D * ncmisid = (TH1D*)hist_numunc->Clone("nu_dune_mulike_ncmisid");
    TH1D * mubar = (TH1D*)hist_numubar->Clone("nu_dune_mulike_barsignal");
    TH1D * taumisid = (TH1D*)hist_nutaucc->Clone("nu_dune_mulike_taumisid");
    signal->Write();
    ncmisid->Write();
    mubar->Write();
    taumisid->Write();
    fm->Close();






    THStack *hs = new THStack("hs","mu dis. stack");
    TCanvas *c1 = new TCanvas ("c1","c1",600,600);
    c1->cd();
    
    //hs->Add(hist_nutaucc,"histo");
    
    //hist_numusig->GetXaxis()->SetRangeUser(0.5,8.0);
    //hist_numubar->GetXaxis()->SetRangeUser(0.5,8.0);

    hs->Add(hist_numubar,"histo");
    hs->Add(hist_nutaucc,"histo");

    hs->Add(hist_numunc,"histo");
    hs->Add(hist_numusig,"histo");
    
    
    //hs->SetRangeUser(0.5,8.0);

    
    hs->Draw();
    
    TLegend *leg = new TLegend(0.6,0.6,0.9,0.9);
    leg->AddEntry(hist_numusig,"signal numu CC","l");
    leg->AddEntry(hist_numunc,"NC","l");
    leg->AddEntry(hist_nutaucc,"nutau+nutaubar CC","l");

    leg->AddEntry(hist_numubar,"bkg numubar CC","l");
    
    leg->Draw();
    //leg->AddEntry(hist_nuenc,"NC","l");
    //leg->AddEntry(hist_nutaucc,"nutau+nutaubar CC","l");
    
    c1->SaveAs("stack_mu.pdf");
    
    //TCanvas *c2 = new TCanvas ("c2","c2",600,600);
    //c2->cd();
    
    //hist_numusigtrue->Draw("histo");
    
    
    
}



void effi_bf(){

    TFile *fnutau = new TFile("gntp.0.nutau20k_gst_study_out.root");
   
    TH1F *hist_nutaucc = (TH1F*)fnutau->Get("hist_cc_reco");
    TH1F *hist_effi_nutaucc = (TH1F*)fnutau->Get("hist_cc_reco_effi");
    
    for (int i=1; i< 81; i++){
        hist_nutaucc->SetBinContent(i, hist_nutaucc->GetBinContent(i)/hist_effi_nutaucc->GetBinContent(i));
        cout << hist_nutaucc->GetBinContent(i) << endl;
    }
    hist_nutaucc->Rebin(2);
    hist_nutaucc->Draw("histo");
    
    
}




void flux_normalization(){
    
    
    TFile *fnue = new TFile("g4lbne_FHC_FD.root");
    TH1F *numucc_flux = (TH1F*)fnue->Get("numu_cceventrate");
    TH1F *numubarcc_flux = (TH1F*)fnue->Get("numubar_cceventrate");
    TH1F *nuecc_flux = (TH1F*)fnue->Get("nue_cceventrate");
    TH1F *nuebarcc_flux = (TH1F*)fnue->Get("nuebar_cceventrate");
    
    TH1F *numunc_flux = (TH1F*)fnue->Get("numu_nceventrate");
    TH1F *numubarnc_flux = (TH1F*)fnue->Get("numubar_nceventrate");
    TH1F *nuenc_flux = (TH1F*)fnue->Get("nue_nceventrate");
    TH1F *nuebarnc_flux = (TH1F*)fnue->Get("nuebar_nceventrate");
    
    
    double numucc_flux_integral = 0.;
    double numubarcc_flux_integral = 0.;
    double nuecc_flux_integral = 0.;
    double nuebarcc_flux_integral = 0.;
    double numunc_flux_integral = 0.;
    double numubarnc_flux_integral = 0.;
    double nuenc_flux_integral = 0.;
    double nuebarnc_flux_integral = 0.;
    
    for(int i = 1; i < 89; i++){
        
        numucc_flux_integral += numucc_flux->GetBinContent(i);
        numubarcc_flux_integral += numubarcc_flux->GetBinContent(i);
        nuecc_flux_integral += nuecc_flux->GetBinContent(i);
        nuebarcc_flux_integral += nuebarcc_flux->GetBinContent(i);
        
        numunc_flux_integral += numunc_flux->GetBinContent(i);
        numubarnc_flux_integral += numubarnc_flux->GetBinContent(i);
        
        nuenc_flux_integral += nuenc_flux->GetBinContent(i);
        nuebarnc_flux_integral += nuebarnc_flux->GetBinContent(i);
        
        
        cout << numucc_flux->GetBinContent(i) << ", " << numucc_flux_integral << endl;
        
    }
    
    
    
    cout << "muoncc : " << numucc_flux_integral << ", electroncc : " << nuecc_flux_integral << endl;
    cout << "numubarcc : " << numubarcc_flux_integral << ", enubarcc : " << nuebarcc_flux_integral << endl;
    
    cout << "muonnc : " << numunc_flux_integral << ", electronnc : " << nuenc_flux_integral << endl;
    cout << "numubarnc : " << numubarnc_flux_integral << ", enubarnc : " << nuebarnc_flux_integral << endl;
    
    //double KT = 40.;
    double KT = 55.641;
    double POT = 1.47E21*3.5;
    
    numucc_flux_integral *= KT*POT;
    nuecc_flux_integral *= KT*POT;
    numubarcc_flux_integral *= KT*POT;
    nuebarcc_flux_integral *= KT*POT;
    numunc_flux_integral *= KT*POT;
    nuenc_flux_integral *= KT*POT;
    numubarnc_flux_integral *= KT*POT;
    nuebarnc_flux_integral *= KT*POT;
    
    
    cout << "cc: muon evt # : " << numucc_flux_integral << ", electron evt # : " << nuecc_flux_integral << "mu nu bar evt # : " << numubarcc_flux_integral << ", e nu bar evt # : " << nuebarcc_flux_integral<< endl;
    
    cout << "nc : muon evt # : " << numunc_flux_integral << ", electron evt # : " << nuenc_flux_integral << "mu nu bar evt # : " << numubarnc_flux_integral << ", e nu bar evt # : " << nuebarnc_flux_integral << endl;
    
    cout << "muon evt # : " << numucc_flux_integral+numunc_flux_integral << ", electron evt # : " << nuecc_flux_integral+nuenc_flux_integral << "mu nu bar evt # : " << numubarcc_flux_integral+numubarnc_flux_integral << ", e nu bar evt # : " << nuebarcc_flux_integral+nuebarnc_flux_integral<< endl;
    
    
    
}// load hists and plot here
