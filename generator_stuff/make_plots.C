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



void make_hist_from_tree(TTree *t1, TH1F *h, Double_t xsec){
    //read event data
    h->Sumw2();
    Long64_t size  =  t1->GetEntries();
    Double_t lumi = 1000* 35.9;
    Double_t evt_weight = lumi*xsec / size;
    printf("%.3e \n", evt_weight);
    int nSelected=0;
    
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
        if(cm.M() >= 100. && cm.M() < 200.){
            nSelected++;
            h->Fill(cm.Pt(), evt_weight);
        }
        



    }
    printf("Selected %i events \n", nSelected);


    t1->ResetBranchAddresses();
    return;
}


void make_plots(){
    gStyle->SetOptStat(0);
    TFile *f_unbinned = TFile::Open("mass_unbinned_2M.root");
    TTree *t_unbinned = (TTree *)f_unbinned->Get("T_lhe");

    TFile *f_binned = TFile::Open("mass_binned_100k.root");
    TTree *t_binned = (TTree *)f_binned->Get("T_lhe");

    TH1F *h_pt_un = new TH1F("h_pt_un", "Binned vs. UnBinned DY Samples (100<M<200)", 20, 100, 200);
    TH1F *h_pt = new TH1F("h_pt", "", 20, 100, 200);

    make_hist_from_tree(t_unbinned, h_pt_un, 5941.0);
    make_hist_from_tree(t_binned, h_pt, 226.6);

    h_pt_un->Scale(1./h_pt_un->Integral());
    h_pt->Scale(1./h_pt->Integral());

    TCanvas * c1 = new TCanvas("c1", "", 800, 800);
    TPad *pad1 = new TPad("pad1c", "pad1", 0.,0.3,0.98,1.);
    pad1->SetBottomMargin(0.01);
    pad1->Draw();
    pad1->SetLogy();
    pad1->cd();
    h_pt->SetLineColor(kRed);
    h_pt_un->SetLineColor(kBlue);
    h_pt->SetLineWidth(3);
    h_pt_un->SetLineWidth(3);
    h_pt_un->Draw("hist E");
    h_pt->Draw("hist E same");

    TLegend *leg3 = new TLegend(0.5, 0.65, 0.75, 0.8);
    leg3->AddEntry(h_pt, "DY Mass binned (M-100to200)", "l");
    leg3->AddEntry(h_pt_un, "DY Unbinned (M-50)", "l");
    leg3->Draw();

    c1->Update();

    c1->cd();
    TPad *pad2 = new TPad("cost_pad2", "pad2", 0.,0,.98,0.3);
    //pad2->SetTopMargin(0);
    pad2->SetBottomMargin(0.2);
    pad2->SetGridy();
    pad2->Draw();
    pad2->cd();
    auto h_ratio = (TH1F *) h_pt->Clone("h_ratio");
    h_ratio->Divide(h_pt_un);
    h_ratio->SetMarkerStyle(21);
    h_ratio->SetLineColor(kBlack);
    h_ratio->Draw("ep");
    h_ratio->SetMinimum(0.);
    h_ratio->SetMaximum(2.0);

    h_ratio->GetYaxis()->SetTitle("Binned/Unbinned");
    h_ratio->GetYaxis()->SetNdivisions(505);
    h_ratio->GetYaxis()->SetTitleSize(20);
    h_ratio->GetYaxis()->SetTitleFont(43);
    h_ratio->GetYaxis()->SetTitleOffset(1.2);
    h_ratio->GetYaxis()->SetLabelFont(43); // Absolute font size in pixel (precision 3)
    h_ratio->GetYaxis()->SetLabelSize(15);
    //h axis cost_ratio plot settings
    h_ratio->GetXaxis()->SetTitle("dilepton Mass (ee, #mu#mu, #tau#tau) (GeV)");
    h_ratio->GetXaxis()->SetTitleSize(20);
    h_ratio->GetXaxis()->SetTitleFont(43);
    h_ratio->GetXaxis()->SetTitleOffset(3.);
    h_ratio->GetXaxis()->SetLabelFont(43); // Absolute font size in pixel (precision 3)
    h_ratio->GetXaxis()->SetLabelSize(20);

    return;
}

