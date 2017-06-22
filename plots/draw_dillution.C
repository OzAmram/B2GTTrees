
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





void draw_dillution(){
    const int nBins = 10;
    TH1F *h_xF = new TH1F("h_xf", "xF distribution, M>150", nBins, 0, 0.4);
    TH1F *h_Nc = new TH1F("h_Nc", "Number Correct; xF", nBins, 0, 0.4);
    TH1F *h_Ni = new TH1F("h_Ni", "Number Incorrect", nBins, 0, 0.4);


    //read event data
    TFile* f_mc = (TFile*) TFile::Open("../analyze/output_files/DYToLL_mc_2016_jun13.root");
    TTree *t1 = (TTree *) f_mc ->Get("T_data");
    Long64_t size  =  t1->GetEntries();
    Double_t m, xF, cost, cost_st, mu1_pt, mu2_pt, jet1_pt, jet2_pt, gen_weight;
    Float_t met_pt;
    t1->SetBranchAddress("m", &m);
    t1->SetBranchAddress("xF", &xF);
    t1->SetBranchAddress("cost", &cost);
    t1->SetBranchAddress("cost_st", &cost_st);
    t1->SetBranchAddress("mu1_pt", &mu1_pt);
    t1->SetBranchAddress("mu2_pt", &mu2_pt);
    t1->SetBranchAddress("jet1_pt", &jet1_pt);
    t1->SetBranchAddress("jet2_pt", &jet2_pt);
    t1->SetBranchAddress("met_pt", &met_pt);
    t1->SetBranchAddress("gen_weight", &gen_weight);


    for (int i=0; i<size; i++) {
        t1->GetEntry(i);
        if(m >= 150. ){
            //printf("%0.2f \n", xF);
            h_xF->Fill(xF, gen_weight);
            Double_t ratio = cost_st/cost;
            if(ratio > 0) h_Nc->Fill(xF);
            if(ratio < 0) h_Ni->Fill(xF);

            }
        }
    
    TCanvas *c1 = new TCanvas("c1", "canvas", 200,10, 900,700);
    h_xF->SetFillColor(kBlue);
    h_xF->Draw();
    c1->Update();

    TCanvas *c2 = new TCanvas("c2", "canva", 100,100, 700,700);
    h_Nc ->SetLineColor(kBlue);
    h_Nc->SetStats(kFALSE);
    h_Nc ->SetLineWidth(2);
    h_Ni ->SetLineColor(kRed);
    h_Ni ->SetLineWidth(2);
    h_Nc ->Draw();
    h_Ni ->Draw("same");
    c2->Update();
    c2->BuildLegend();

    Double_t Nc, Ni, dillu[nBins], dillu_error[nBins], xF_center[nBins]; 

    for (int i=1; i < nBins; i++){
        Nc = h_Nc->GetBinContent(i);
        Ni = h_Ni->GetBinContent(i);
        xF_center[i] = h_Nc->GetBinCenter(i);
        dillu[i] = (Nc - Ni)/(Nc + Ni);
        dillu_error[i] = sqrt(pow(1./(Nc+Ni) - Nc/(Nc+Ni)/(Nc+Ni), 2)*Nc + 
                              pow(1./(Nc+Ni) - Ni/(Nc+Ni)/(Nc+Ni), 2)*Ni);
    }

    TGraph *g_dillu = new TGraph(nBins, xF_center, dillu);

    TCanvas *c3 = new TCanvas("c3", "canvas", 200,10, 900,700);
    g_dillu->Draw("AP");
    g_dillu->SetMarkerStyle(20);
    g_dillu->SetTitle("Dillution Factor: (N_c - N_i)/(N_c + N_i)");
    g_dillu->GetXaxis()->SetTitle("xF");

    c3->Update();



    }
