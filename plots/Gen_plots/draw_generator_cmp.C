
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <iostream>
#include <cstring>
#include <algorithm>
#include <vector>

#include "TROOT.h"
#include "TFile.h"
#include "TObject.h"
#include "TRatioPlot.h"
#include "TH1F.h"
#include "TProfile.h"
#include "TStyle.h"
#include "TCanvas.h"
#include "TPostScript.h"
#include "TLorentzVector.h"
#include "TVector3.h"
#include "TFitter.h"
#include "TSystem.h"



void make_hists_from_tree(TTree *t1, TH1F *h_m, TH1F *h_cost){
    //read event data
    h_m->Sumw2();
    h_cost->Sumw2();
    Long64_t size  =  t1->GetEntries();
    int nSelected=0;
    double root2 = sqrt(2);
    
    Int_t lep1_id, lep2_id;
    TLorentzVector *lep_pls = 0;
    TLorentzVector *lep_mns = 0;
    TLorentzVector cm;
    t1->SetBranchAddress("lep_pls", &lep_pls);
    t1->SetBranchAddress("lep_mns", &lep_mns);
    t1->SetBranchAddress("lep1_id", &lep1_id);
    t1->SetBranchAddress("lep2_id", &lep2_id);
    for (int i=0; i<size; i++) {
        t1->GetEntry(i);
        cm = *lep_pls + *lep_mns;
        //printf("M= %.2f \n", cm.M());
        if(cm.M() >= 150. && (lep1_id < 14) && (lep2_id < 14) && cm.M() <=200.){
            nSelected++;
            h_m->Fill(cm.M());

            double mu_p_pls = (lep_pls->E()+lep_pls->Pz())/root2;
            double mu_p_min = (lep_pls->E()-lep_pls->Pz())/root2;
            double mu_m_pls = (lep_mns->E()+lep_mns->Pz())/root2;
            double mu_m_min = (lep_mns->E()-lep_mns->Pz())/root2;
            double qt2 = cm.Px()*cm.Px()+cm.Py()*cm.Py();
            double cm_m2 = cm.M2();
            double gen_cost = 2*(mu_m_pls*mu_p_min - mu_m_min*mu_p_pls)/sqrt(cm_m2*(cm_m2 + qt2));
            //pls and mns leptons are backwards in ttree, so flip sign
            gen_cost = -gen_cost;

            //flip sign for reconstruction
            if(cm.Pz() < 0.) gen_cost = -gen_cost;
            h_cost->Fill(gen_cost);
        }
        



    }
    printf("Selected %i events \n", nSelected);
    h_cost->Scale(1./h_cost->Integral());
    h_m->Scale(1./h_m->Integral());


    t1->ResetBranchAddresses();
    return;
}

void make_ratio_plot(char title[80], TH1F* h1, char h1_label[80], TH1F* h2, char h2_label[80], char ratio_label[80], 
        char axis_label[80], bool logy=false){

    TCanvas *c = new TCanvas("c_m", "Histograms", 200, 10, 900, 700);
    TPad *pad1 = new TPad("pad1", "pad1", 0.,0.3,0.98,1.);
    pad1->SetBottomMargin(0);
    pad1->Draw();
    pad1->cd();
    if(logy) pad1->SetLogy();
    //h1->SetMinimum(0.1);
    h1->Draw("hist E");
    gStyle->SetEndErrorSize(4);
    h2->Draw("hist E same");


    gStyle->SetLegendBorderSize(0);
    TLegend *leg1 = new TLegend(0.15, 0.7, 0.35, 0.85);
    leg1->AddEntry(h1, h1_label, "l");
    leg1->AddEntry(h2, h2_label, "l");
    leg1->Draw();

    //gPad->BuildLegend();
    c->cd();
    TPad *pad2 = new TPad("pad2", "pad2", 0.,0,.98,0.3);
    //pad2->SetTopMargin(0);
    pad2->SetBottomMargin(0.2);
    pad2->SetGridy();
    pad2->Draw();
    pad2->cd();

    auto ratio = (TH1F *) h1->Clone("h_ratio");
    ratio->SetMinimum(0.01);
    ratio->SetMaximum(2.0);
    ratio->Sumw2();
    ratio->SetStats(0);
    ratio->Divide(h2);
    ratio->SetMarkerStyle(21);
    ratio->SetLineColor(kBlack);
    ratio->Draw("ep");
    c->cd();

    ratio->SetTitle("");
    // Y axis ratio plot settings
    ratio->GetYaxis()->SetTitle(ratio_label);
    ratio->GetYaxis()->SetNdivisions(505);
    ratio->GetYaxis()->SetTitleSize(20);
    ratio->GetYaxis()->SetTitleFont(43);
    ratio->GetYaxis()->SetTitleOffset(1.2);
    ratio->GetYaxis()->SetLabelFont(43); // Absolute font size in pixel (precision 3)
    ratio->GetYaxis()->SetLabelSize(15);
    // X axis ratio plot settings
    ratio->GetXaxis()->SetTitle(axis_label);
    ratio->GetXaxis()->SetTitleSize(20);
    ratio->GetXaxis()->SetTitleFont(43);
    ratio->GetXaxis()->SetTitleOffset(3.);
    ratio->GetXaxis()->SetLabelFont(43); // Absolute font size in pixel (precision 3)
    ratio->GetXaxis()->SetLabelSize(20);

    //lumi_sqrtS = "";       // used with iPeriod = 0, e.g. for simulation-only plots (default is an empty string)
    //int iPeriod = 4; 
    //CMS_lumi(pad1, iPeriod, 33 );
    c->Print(title);
    return;
}

void draw_generator_cmp(){
    gStyle->SetOptStat(0);
    TFile *f_mad = TFile::Open("../generator_stuff/mass_binned_200k.root");
    TTree *t_mad = (TTree *)f_mad->Get("T_lhe");

    TFile *f_pwg = TFile::Open("../generator_stuff/powheg_m150_evts.root");
    TTree *t_pwg = (TTree *)f_pwg->Get("T_lhe");

    TH1F *h_pwg_m = new TH1F("h_pwg_m", "POWHEG vs. aMC@NLO DY Samples (150 < M < 200)", 20, 150, 200);
    TH1F *h_mad_m = new TH1F("h_mad_m", "POWHEG vs. aMC@NLO DY Samples (150 < M < 200)", 20, 150, 200);
    TH1F *h_pwg_cost = new TH1F("h_pwg_cost", "POWHEG vs. aMC@NLO DY Samples (150 < M < 200)", 20, -1., 1.);
    TH1F *h_mad_cost = new TH1F("h_mad_cost", "POWHEG vs. aMC@NLO DY Samples (150 < M < 200)", 20, -1., 1.);

    make_hists_from_tree(t_mad, h_mad_m, h_mad_cost);
    make_hists_from_tree(t_pwg, h_pwg_m, h_pwg_cost);

    h_mad_m->SetLineColor(kRed);
    h_mad_cost->SetLineColor(kRed);
    h_pwg_m->SetLineColor(kBlue);
    h_pwg_cost->SetLineColor(kBlue);

    h_mad_m->SetLineWidth(3);
    h_mad_cost->SetLineWidth(3);
    h_pwg_m->SetLineWidth(3);
    h_pwg_cost->SetLineWidth(3);

    make_ratio_plot("pwg_vs_mad_m150_m_cmp.pdf", h_mad_m, "amc@NLO",h_pwg_m, "POWHEG", "aMC/POWHEG", "M (GeV)", true);
    make_ratio_plot("pwg_vs_mad_m150_cost_cmp.pdf", h_mad_cost, "amc@NLO",h_pwg_cost, "POWHEG", "aMC/POWHEG", "cos(#theta_{*})", false);

    return;
}

