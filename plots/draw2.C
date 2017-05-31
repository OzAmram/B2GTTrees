//perform fits to Reconstructed MuMu data to extract Asym

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
#include "TH1F.h"
#include "TProfile.h"
#include "TStyle.h"
#include "TCanvas.h"
#include "TPostScript.h"
#include "TLorentzVector.h"
#include "TVector3.h"
#include "TFitter.h"
#include "TSystem.h"
#include "Math/Functor.h"

TH1F* h_integral(TH1F *h, char *title){
    int nbins = h->GetSize() - 2; //remove over and underflow bins
    Double_t min = h->GetBinLowEdge(1); 
    Double_t max = h->GetBinLowEdge(nbins+1);
    Double_t entries = h->GetEntries();
    TH1F *h2 = new TH1F(title,"",nbins, min, max);
    //printf("nbins %i min %f max %f \n",nbins, min, max);
    Double_t sum =0;
    for(int i=1; i<nbins+2; i++){
        sum += h->GetBinContent(i)/entries;
        //printf("sum %f \n", sum);
        h2->SetBinContent(i,sum);
    }
    Float_t rightmax = 1.1*h2->GetMaximum();
    Float_t scale = gPad->GetUymax()/rightmax;
    h2->SetLineWidth(4);
    h2->SetLineColor(kBlue);
    h2->Scale(scale);
    h2->Draw("same");
    TGaxis *axis = new TGaxis(gPad->GetUxmax(),gPad->GetUymin(),
                      gPad->GetUxmax(), gPad->GetUymax(), 0,rightmax,510,"+L");
    axis->SetLineColor(kBlue);
    axis->SetLabelColor(kBlue);
    axis->Draw();


    return h2;
}



void draw2(){
    TFile *f_data = TFile::Open("../analyze/output_files/ttbar_background_mar2.root");
    TH1F *h_m = new TH1F("h_m", "mass above 150 gev, combination of samples; mass (GeV)", 30, 150, 300);

    TH1F *h_cost1 = new TH1F("cost1", "Cos(#theta) for M in [150, 300] (Gev);cos(#theta)", 25, -1,1);


    TH1F *h_cost_cut = new TH1F("cost_cut", "Cos(#theta) for M in [150, 180] (Gev): met_pt <80, mu1_pt > 40, mu2_pt>20 ;cos(#theta)", 
            25, -1,1);
    TH1F *h_mu1_pt = new TH1F("mu1_pt", "Pt of leading muon for M in [150,180]", 50, 0,200);
    TH1F *h_mu2_pt = new TH1F("mu2_pt", "Pt of 2nd muon for M in [150,180]", 50, 0,200);

    TH1F *h_met_pt = new TH1F("met_pt", "DY MC Missing pt for M in [150,300]", 50, 0,150);
    
    TH1F *h_jet1_pt = new TH1F("jet1_pt", "ttbar Pt of leading jet for M in [150,300], met<50", 50, 0,200);
    TH1F *h_jet2_pt = new TH1F("jet2_pt", "ttbar Pt of 2nd jet for M in [150,300], met<50", 50, 0,200);
    //read event data
    TTree *t1 = (TTree *)f_data->Get("T_data"); //Tree of data
    Long64_t size  =  t1->GetEntries();
    Double_t m, xF, cost, mu1_pt, mu2_pt, jet1_pt, jet2_pt, gen_weight;
    Float_t met_pt;
    t1->SetBranchAddress("m", &m);
    t1->SetBranchAddress("xF", &xF);
    t1->SetBranchAddress("cost", &cost);
    t1->SetBranchAddress("mu1_pt", &mu1_pt);
    t1->SetBranchAddress("mu2_pt", &mu2_pt);
    t1->SetBranchAddress("jet1_pt", &jet1_pt);
    t1->SetBranchAddress("jet2_pt", &jet2_pt);
    t1->SetBranchAddress("met_pt", &met_pt);
    t1->SetBranchAddress("gen_weight", &gen_weight);

    int nBack = 0;
    int nEvents = 0;
    int nForward = 0;
    int nSelected = 0;
    int nF = 0;
    int nB=0;

    for (int i=0; i<size; i++) {
        t1->GetEntry(i);
        if(m >= 150. && m < 300.){
            h_cost1->Fill(cost, gen_weight);
            h_m->Fill(m,gen_weight);
            h_mu1_pt->Fill(mu1_pt,gen_weight);
            h_mu2_pt->Fill(mu2_pt,gen_weight);
            h_jet1_pt->Fill(jet1_pt, gen_weight);
            h_jet2_pt->Fill(jet2_pt,gen_weight);
            h_met_pt->Fill(met_pt, gen_weight);

            nEvents++;
            if(met_pt < 80. && mu1_pt > 40. && mu2_pt > 20.){
                if(cost > 0) nForward++;
                if(cost < 0) nBack++;
                if(cost > 0 && xF > 0.1) nF++;
                if(cost < 0 && xF>0.1) nB++;
                nSelected++;
                h_cost_cut->Fill(cost);
            }
        }
    }
    float AFB = float(nForward - nBack)/float(nForward + nBack);
    float AFB_highx = float(nF - nB)/float(nF + nB);

    printf("high x count %i, AFB %0.3f \n", nF+nB, AFB_highx);

    printf("Running through %i events. Counting AFB is %0.3f. There were %i selected. Selection efficiency was %0.3f\n", 
            nEvents, AFB, nSelected, float(nSelected)/float(nEvents));

    TCanvas *c_cost1 = new TCanvas("c_cost1", "Histograms", 200, 10, 900, 700);
    c_cost1->cd();
    h_cost1->Draw();
    h_cost1->SetFillColor(35);


    TCanvas *c_met_pt = new TCanvas("c_met_pt", "Histograms", 200, 10, 900, 700);
    c_met_pt->cd();
    h_met_pt->SetFillColor(35);
    h_met_pt->Draw();
    c_met_pt->Update();
    TH1F *met_pt_int = h_integral(h_met_pt, "met_pt_int");

    TCanvas *c_m = new TCanvas("c_m", "Histograms", 200, 10, 900, 700);
    c_m->cd();
    h_m->Draw();

    TCanvas *c_jet1_pt = new TCanvas("c_jet1_pt", "Histograms", 200, 10, 900, 700);
    c_jet1_pt->cd();
    h_jet1_pt->Draw();

    TCanvas *c_jet2_pt = new TCanvas("c_jet2_pt", "Histograms", 200, 10, 900, 700);
    c_jet2_pt->cd();
    h_jet2_pt->Draw();
    /*
    h_cost_cut->Draw("same");
    h_cost_cut->SetLineWidth(3);
    h_cost_cut->SetLineColor(kRed);
    c_cost1->Update();

    leg = new TLegend(0.1,0.6,0.5,0.9);
    leg->AddEntry(h_cost1, "No extra cuts", "f");
    leg->AddEntry(h_cost_cut, "met_pt < 80, mu1_pt > 40, mu2_pt > 20", "f");
    leg->Draw();
    */

    //TCanvas *c_cost_cut = new TCanvas("c_cost_cut", "Histograms", 200, 10, 900, 700);
    //c_cost_cut->cd();

/*

    TCanvas *c_mu1_pt = new TCanvas("c_mu1_pt", "Histograms", 200, 10, 900, 700);
    c_mu1_pt->cd();
    h_mu1_pt->SetFillColor(35);
    h_mu1_pt->Draw();
    c_mu1_pt->Update();
    TH1F *mu1_pt_int = h_integral(h_mu1_pt, "mu1_pt_int");

    TCanvas *c_mu2_pt = new TCanvas("c_mu2_pt", "Histograms", 200, 10, 900, 700);
    c_mu2_pt->cd();
    h_mu2_pt->SetFillColor(35);
    h_mu2_pt->Draw();
    c_mu2_pt->Update();
    TH1F *mu2_pt_int = h_integral(h_mu2_pt, "mu2_pt_int");


    //TCanvas *c_cost2 = new TCanvas("c_cost2", "Histograms", 200, 10, 900, 700);
    TCanvas *c_jet1_pt = new TCanvas("c_jet1_pt", "Histograms", 200, 10, 900, 700);
    c_jet1_pt->cd();
    h_jet1_pt->Draw();


*/


 
}

    
    
