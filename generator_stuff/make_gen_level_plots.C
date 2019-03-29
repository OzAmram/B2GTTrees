
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <iostream>
#include <cstring>
#include <algorithm>
#include <vector>

#include "TROOT.h"
#include "TFile.h"
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
#include "TClonesArray.h"
#include "make_gen_level_plots.h"




void fill_gen_hist(TTree *t, TH1F *h_m, TH1F *h_pt){
    Long64_t size  =  t->GetEntries();
    float weight, gen_pt[100], gen_eta[100], gen_phi[100], gen_mass[100];
    int gen_status[100], gen_id[100], gen_parent[100], gen_size;

    TClonesArray *gen_parts = new TClonesArray("TGenParticle");
    //TGenEventInfo *gen_evt;

    int OUTGOING = 23;
    //Particle ID's
    int ELECTRON = 11; 
    int MUON = 13;
    int PHOTON = 22;
    int Z=23;
    int GLUON = 21;
    int TAU = 15;
    int PROTON = 2212;

    TLorentzVector lep1, lep2;
    bool got_lep1, got_lep2;

    //t->SetBranchAddress("GenEvtInfo", &gen_evt);
    t->SetBranchAddress("weight", &weight);
    t->SetBranchAddress("GenParticle", &gen_parts);
    //t->SetBranchAddress("GenParticle.eta", &gen_eta);
    //t->SetBranchAddress("GenParticle.phi", &gen_phi);
    //t->SetBranchAddress("GenParticle.mass", &gen_mass);
    //t->SetBranchAddress("GenParticle@.GetEntries()", &gen_size);
    for(int i=0; i<size; i++){
        got_lep1 = got_lep2 = false;
        t->GetEntry(i);
        //weight = gen_evt->weight;
        gen_size = gen_parts->GetEntries();
        for(int j=0; j<gen_size; j++){
            TGenParticle *genp = (TGenParticle *)((*gen_parts)[j]);
            if((abs(genp->pdgId) == MUON || abs(genp->pdgId) == ELECTRON) &&
                    genp->status == OUTGOING){
                if(!got_lep1){
                    got_lep1= true;
                }
                else if(!got_lep2){
                    lep2.SetPtEtaPhiM(genp->pt, genp->eta, genp->phi, genp->mass);
                    got_lep1= true;
                }
                else{printf("Extra lepton! \n");}
            }
        }
        TLorentzVector cm = lep1 + lep2;
        if(cm.M() > 100.){
            h_m->Fill(cm.M(), weight);
            h_pt->Fill(cm.Pt(), weight);
        }
    }
    h_m->Scale(1./h_m->Integral());
    h_pt->Scale(1./h_pt->Integral());

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
    h1->Draw("hist E");
    //m_stack->SetMaximum(65000);
    h1->SetMinimum(0.1);
    gStyle->SetEndErrorSize(4);
    h2->Draw("hist E same");


    gStyle->SetLegendBorderSize(0);
    TLegend *leg1 = new TLegend(0.25, 0.05, 0.4, 0.2);
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
    c->Print(title);
    return;
}


void make_gen_level_plots(){
    gStyle->SetOptStat(0);
    TFile *f_unbinned = TFile::Open("gen_level/DY_M50_gen.root");
    TTree *t_unbinned = (TTree *)f_unbinned->Get("Events");

    TFile *f_binned = TFile::Open("gen_level/DY_M100_gen.root");
    TTree *t_binned = (TTree *)f_binned->Get("Events");

    TH1F *h_unbin_m = new TH1F("h_unbin_m", "", 20, 100, 200);
    TH1F *h_unbin_pt = new TH1F("h_unbin_pt", "", 40, 0, 400);

    TH1F *h_bin_m = new TH1F("h_bin_m", "", 20, 100, 200);
    TH1F *h_bin_pt = new TH1F("h_bin_pt", "", 40, 0, 400);
    fill_gen_hist(t_unbinned, h_unbin_m, h_unbin_pt);
    fill_gen_hist(t_binned, h_bin_m, h_bin_pt);



    make_ratio_plot("gen_m_cmp.pdf", h_bin_m, "Binned MC (DY_M-100to200)",h_unbin_m, "Unbinned MC (DY_M50)", "Binned/Unbinned", "M (GeV)", true);
    make_ratio_plot("gen_pt_cmp.pdf", h_bin_pt, "Binned MC (DY_M-100to200)",h_unbin_pt, "Unbinned MC (DY_M50)", "Binned/Unbinned", "Pt (GeV)", true);
    return;
}
