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


void read_tree(TTree *t1, TH1F *h_m, TH1F *h_cost){
    //read event data
    Long64_t size  =  t1->GetEntries();
    Double_t m, xF, cost, mu1_pt, mu2_pt, jet1_cmva, jet2_cmva, gen_weight;
    Float_t cost_pt;
    t1->SetBranchAddress("m", &m);
    t1->SetBranchAddress("xF", &xF);
    t1->SetBranchAddress("cost", &cost);
    t1->SetBranchAddress("jet2_CMVA", &jet2_cmva);
    t1->SetBranchAddress("jet1_CMVA", &jet1_cmva);
    t1->SetBranchAddress("gen_weight", &gen_weight);

    for (int i=0; i<size; i++) {
        t1->GetEntry(i);
        if(m >= 150. && jet1_cmva < 0. && jet2_cmva < 0){
            h_m->Fill(m,gen_weight);
            h_cost->Fill(cost, gen_weight);


        }
    }
}


void draw_stack(){
    TFile *f_mc = TFile::Open("../analyze/output_files/DYToLL_mc_2016_may9.root");
    TTree *t_mc = (TTree *)f_mc->Get("T_data");



    TH1F *mc_m = new TH1F("mc_m", "MC Signal (qqbar, qglu, qbarglu)", 30, 150, 1000);
    mc_m->SetFillColor(kRed);
    mc_m->SetMarkerStyle(21);
    mc_m->SetMarkerColor(kRed);

    TH1F *mc_cost = new TH1F("mc_cost", "MC Signal (qqbar, qglu, qbarglu)", 40, -1,1);
    mc_cost->SetFillColor(kRed);
    mc_cost->SetMarkerStyle(21);
    mc_cost->SetMarkerColor(kRed);


    read_tree(t_mc, mc_m, mc_cost);

    TCanvas *c1 = new TCanvas("c1", "Histograms", 200, 10, 900, 700);
    mc_m->Draw("h");




}



