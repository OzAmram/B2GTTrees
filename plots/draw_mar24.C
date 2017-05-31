
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

int draw_mar24(){
    TFile *f_data = TFile::Open("../analyze/output_files/DYToLL_data_2015_mar23.root");

    TTree *t_data = (TTree *)f_data->Get("T_data"); //Tree of data
    Long64_t size  =  t_data->GetEntries();

    Double_t jet1_csv, jet2_csv, jet1_cmva, jet2_cmva, gen_weight;
    t_data->SetBranchAddress("jet2_CMVA", &jet2_cmva);
    t_data->SetBranchAddress("jet1_CMVA", &jet1_cmva);
    t_data->SetBranchAddress("jet2_CSV", &jet2_csv);
    t_data->SetBranchAddress("jet1_CSV", &jet1_csv);


    Int_t nbins = 110;
    Double_t min = -10.;
    Double_t max = 1.;
    TH1F *h_data_csv1 = new TH1F("h_data_csv1", "CSV of leading jet from data sample", nbins, min,max);
    TH1F *h_data_csv2 = new TH1F("h_data_csv2", "CSV of 2nd jet from data sample", nbins, min,max);


    nbins = 40;
    min = -1.;
    max = 1.;
    TH1F *h_data_cmva1 = new TH1F("h_data_cmva1", "CMVA of leading jet from data sample", nbins, min,max);
    TH1F *h_data_cmva2 = new TH1F("h_data_cmva2", "CMVA of 2nd jet from data sample", nbins, min,max);


    for (int i=0; i<size; i++) {
        t_data->GetEntry(i);

        h_data_csv1->Fill(jet1_csv);
        h_data_csv2->Fill(jet2_csv);
        h_data_cmva1->Fill(jet1_cmva);
        h_data_cmva2->Fill(jet2_cmva);


    }

    h_data_csv2->Scale(1/h_data_csv2->Integral());
    h_data_cmva1->Scale(1/h_data_cmva1->Integral());
    h_data_cmva2->Scale(1/h_data_cmva2->Integral());

    TCanvas *c1 = new TCanvas("c1", "Histograms", 200, 10, 900, 700);
    h_data_csv1->Draw();
    h_data_csv1->SetFillColor(kBlue);
    c1->Update();

    TCanvas *c2 = new TCanvas("c2", "Histograms", 200, 10, 900, 700);
    h_data_csv2->Draw();
    h_data_csv2->SetFillColor(kBlue);
    c2->Update();


    TCanvas *c3 = new TCanvas("c3", "Histograms", 200, 10, 900, 700);
    h_data_cmva1->Draw();
    h_data_cmva1->SetFillColor(kBlue);
    c3->Update();

    TCanvas *c4 = new TCanvas("c4", "Histograms", 200, 10, 900, 700);
    h_data_cmva2->Draw();
    h_data_cmva2->SetFillColor(kBlue);
    c4->Update();

    return;
}


