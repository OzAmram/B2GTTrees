
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
#include "TH2F.h"
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

int draw_mar23(){
    TFile *f_ttbar = TFile::Open("../analyze/output_files/ttbar_background_mar23.root");
    TFile *f_mc = TFile::Open("../analyze/output_files/DYToLL_mc_mar23.root");

    TTree *t_ttbar = (TTree *)f_ttbar->Get("T_data"); //Tree of data
    Long64_t size  =  t_ttbar->GetEntries();

    Double_t jet1_csv, jet2_csv, jet1_cmva, jet2_cmva, gen_weight;
    t_ttbar->SetBranchAddress("jet2_CMVA", &jet2_cmva);
    t_ttbar->SetBranchAddress("jet1_CMVA", &jet1_cmva);
    t_ttbar->SetBranchAddress("jet2_CSV", &jet2_csv);
    t_ttbar->SetBranchAddress("jet1_CSV", &jet1_csv);
    t_ttbar->SetBranchAddress("gen_weight", &gen_weight);


    Int_t nbins = 110;
    Double_t min = -10.;
    Double_t max = 1.;
    TH1F *h_ttbar_csv1 = new TH1F("h_ttbar_csv1", "CSV of leading jet from TTbar sample", nbins, min,max);
    TH1F *h_ttbar_csv2 = new TH1F("h_ttbar_csv2", "CSV of 2nd jet from TTbar sample", nbins, min,max);
    TH1F *h_mc_csv1 = new TH1F("h_mc_csv1", "CSV of leading jet from DY sample", nbins, min,max);
    TH1F *h_mc_csv2 = new TH1F("h_mc_csv2", "CSV of 2nd jet from DY sample", nbins, min,max);


    nbins = 200;
    min = -1.;
    max = -.5;
    TH1F *h_ttbar_cmva1 = new TH1F("h_ttbar_cmva1", "CMVA of leading jet from TTbar sample", nbins, min,max);
    TH1F *h_ttbar_cmva2 = new TH1F("h_ttbar_cmva2", "CMVA of 2nd jet from TTbar sample", nbins, min,max);
    TH1F *h_mc_cmva1 = new TH1F("h_mc_cmva1", "CMVA of leading jet from DY sample", nbins, min,max);
    TH1F *h_mc_cmva2 = new TH1F("h_mc_cmva2", "CMVA of 2nd jet from DY sample", nbins, min,max);

    TH2F *h_ttbar_cmva = new TH2F("h_ttbar_cmva", "CVMA of 1st & 2nd jet",
            nbins,min,max, nbins,min,max);


    for (int i=0; i<size; i++) {
        t_ttbar->GetEntry(i);

        h_ttbar_csv1->Fill(jet1_csv, gen_weight);
        h_ttbar_csv2->Fill(jet2_csv, gen_weight);
        h_ttbar_cmva1->Fill(jet1_cmva, gen_weight);
        h_ttbar_cmva2->Fill(jet2_cmva, gen_weight);

        h_ttbar_cmva->Fill(jet1_cmva, jet2_cmva, gen_weight);

    }

    TTree *t_mc = (TTree *) f_mc->Get("T_data");
    size = t_mc->GetEntries();
    
    t_mc->SetBranchAddress("jet2_CMVA", &jet2_cmva);
    t_mc->SetBranchAddress("jet1_CMVA", &jet1_cmva);
    t_mc->SetBranchAddress("jet2_CSV", &jet2_csv);
    t_mc->SetBranchAddress("jet1_CSV", &jet1_csv);
    t_mc->SetBranchAddress("gen_weight", &gen_weight);
    for (int i=0; i<size; i++) {
        t_mc->GetEntry(i);

        h_mc_csv1->Fill(jet1_csv, gen_weight);
        h_mc_csv2->Fill(jet2_csv, gen_weight);
        h_mc_cmva1->Fill(jet1_cmva, gen_weight);
        h_mc_cmva2->Fill(jet2_cmva, gen_weight);

    }
    h_mc_csv1->Scale(1/h_mc_csv1->Integral());
    h_mc_csv2->Scale(1/h_mc_csv2->Integral());
    h_mc_cmva1->Scale(1/h_mc_cmva1->Integral());
    h_mc_cmva2->Scale(1/h_mc_cmva2->Integral());

    h_ttbar_csv1->Scale(1/h_ttbar_csv1->Integral());
    h_ttbar_csv2->Scale(1/h_ttbar_csv2->Integral());
    h_ttbar_cmva1->Scale(1/h_ttbar_cmva1->Integral());
    h_ttbar_cmva2->Scale(1/h_ttbar_cmva2->Integral());

    TCanvas *c1 = new TCanvas("c1", "Histograms", 200, 10, 900, 700);
    h_mc_csv1->Draw();
    h_mc_csv1->SetFillColor(kRed);
    h_ttbar_csv1->Draw("same");
    h_ttbar_csv1->SetFillColor(kBlue);
    c1->Update();
    gPad->BuildLegend();

    TCanvas *c2 = new TCanvas("c2", "Histograms", 200, 10, 900, 700);
    h_mc_csv2->Draw();
    h_mc_csv2->SetFillColor(kRed);
    h_ttbar_csv2->Draw("same");
    h_ttbar_csv2->SetFillColor(kBlue);
    c2->Update();
    gPad->BuildLegend();


    TCanvas *c3 = new TCanvas("c3", "Histograms", 200, 10, 900, 700);
    h_mc_cmva1->Draw();
    h_mc_cmva1->SetFillColor(kRed);
    TCanvas *c3p = new TCanvas("c3p", "Histograms", 200, 10, 900, 700);
    h_ttbar_cmva1->Draw();
    h_ttbar_cmva1->SetFillColor(kBlue);
    c3->Update();
    c3p->Update();
    //gPad->BuildLegend();

    TCanvas *c4 = new TCanvas("c4", "Histograms", 200, 10, 900, 700);
    h_mc_cmva2->Draw();
    h_mc_cmva2->SetLineColor(kRed);
    h_ttbar_cmva2->Draw("same");
    h_ttbar_cmva2->SetFillColor(kBlue);
    c4->Update();
    gPad->BuildLegend();

    TCanvas *c5 = new TCanvas("c5", "Histograms", 200, 10, 900, 700);
    h_ttbar_cmva->Draw();
    return;
}


