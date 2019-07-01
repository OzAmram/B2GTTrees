
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
#include "../../Utils/HistMaker.C"
//

void SetErrors(TH2D *h_rate, TH2D *h_total){
    int nBins_x = h_rate->GetXaxis()->GetNbins();
    int nBins_y = h_rate->GetYaxis()->GetNbins();
    for (int i=1; i <= nBins_x; i++){
        for (int j=1; j <= nBins_y; j++){
            //binomial distribution divided by n
            Double_t n = h_total->GetBinContent(i,j);
            Double_t r = h_rate->GetBinContent(i,j);
            Double_t err = sqrt(r*(1-r)/n);
            if(r <=0){
                if (j != 0){
                    printf("Substituting 0 rate with 1 lower pt bin \n");
                    r = h_rate->GetBinContent(i, j-1);
                    h_rate->SetBinContent(i,j,r);
                    Double_t prev_err = h_rate->GetBinError(i, j-1);
                    err = 1.4 * prev_err;
                    printf(" i j prev error %i %i %.2f \n", i,j, prev_err );
                }

            }
            printf("n r Err is %.0f %.2f %.2f \n", n,r,err);
            if(r == 0) err = 0.2;
            h_rate->SetBinError(i,j, err);
        }
    }
    return;
}

void SetErrors(TH1D *h_rate, TH1D *h_total){
    int nBins = h_rate->GetXaxis()->GetNbins();
    for (int i=0; i <= nBins; i++){
        //binomial distribution divided by n
        Double_t n = h_total->GetBinContent(i);
        Double_t p = h_rate->GetBinContent(i);
        Double_t err = sqrt(p*(1-p)/n);
        printf("Err is %.2f \n", err);
        if(p == 0) err = 0.2;
        h_rate->SetBinError(i, err);
    }
    return;
}

void SetCorrectedRate(TH2D *h_data_rate, TH2D *h_data_total, TH2D *h_contam_rate, TH2D* h_contam_total){
    int nBins_x = h_data_rate->GetXaxis()->GetNbins();
    int nBins_y = h_data_rate->GetYaxis()->GetNbins();
    printf("nx,ny = %i %i \n", nBins_x, nBins_y);
    //printf("Get size %i \n", nBins);
    for (int i=1; i <= nBins_x; i++){
        for (int j=1; j <= nBins_y; j++){
            Double_t r_old = h_data_rate->GetBinContent(i,j);
            Double_t n_old = h_data_total->GetBinContent(i,j);
            Double_t p_old = n_old * r_old;
            if(r_old == 0 || n_old == 0) continue;

            Double_t r_contam = h_contam_rate->GetBinContent(i,j);
            Double_t n_contam = h_contam_total->GetBinContent(i,j);
            Double_t p_contam = r_contam * n_contam;

            Double_t eta_center = h_data_total->GetXaxis()->GetBinCenter(i);
            Double_t pt_center = h_data_total->GetYaxis()->GetBinCenter(j);

            Double_t p_new = p_old - p_contam;
            if(p_new <= 0) p_new = 0;
            Double_t n_new = n_old - n_contam;
            Double_t r_new = p_new/n_new;
            printf("Bin (%.3f, %.0f): (r, n,p) Old (%.2f, %.0f, %.0f) New (%.2f, %.0f, %.0f) \n", eta_center, pt_center, r_old, n_old, p_old, r_new, n_new, p_new);

            h_data_rate->SetBinContent(i,j, r_new);
            h_data_total->SetBinContent(i,j, n_new);
        }
    }
        

}

void construct_fakerate_template(TH2D *h_rate, TH2D *h_total, TTree *t, int flag = FLAG_MUONS, bool isData = true){
    Double_t pt, eta, gen_weight;
    Bool_t pass;
    Double_t lumi;
    if(flag == FLAG_MUONS){
        t->SetBranchAddress("mu_pt", &pt);
        t->SetBranchAddress("mu_eta", &eta);
        lumi = mu_lumi16;
    }
    else{
        t->SetBranchAddress("el_pt", &pt);
        t->SetBranchAddress("el_eta", &eta);
        lumi = el_lumi16;
    }
    t->SetBranchAddress("pass", &pass);
    if(!isData) t->SetBranchAddress("gen_weight", &gen_weight);
    Long64_t size  =  t->GetEntries();
    for (int i=0; i<size; i++) {
        t->GetEntry(i);
        Double_t fill_weight = 1;
        if(!isData) fill_weight = gen_weight;
        h_total->Fill(eta,pt, fill_weight);
        if(pass) h_rate->Fill(eta,pt, fill_weight);
    }
    h_rate->Divide(h_total);
    if(!isData) h_total->Scale(1000*lumi);
    return;
}
    
        
    


void draw_fakerate(){
    /*

    TFile *f = TFile::Open("../analyze/FakeRate/root_files/2016/SingleElectron_data_fakerate_june25.root");
    TFile *f_mc = TFile::Open("../analyze/FakeRate/root_files/2016/SingleElectron_mc_fakerate_contam_june25.root");
    int FLAG = FLAG_ELECTRONS;
    Float_t pt_bins[] = {10,20,35,60,100, 1000};
    int n_pt_bins = 5;
    */

    ///*
    TFile *f = TFile::Open("../analyze/FakeRate/root_files/2016/SingleMuon_data_fakerate_june25.root");
    TFile *f_mc = TFile::Open("../analyze/FakeRate/root_files/2016/SingleMuon_mc_fakerate_contam_june25.root");
    int FLAG = FLAG_MUONS;
    Float_t pt_bins[] = {10,20,26,35,1000};
    int n_pt_bins = 4;
    //*/

    //
    Float_t eta_bins[] = {0,1.479,2.4};
    int n_eta_bins = 2;

    bool write_out = true;

    TH2D *h_rate = new TH2D("h_rate", "passing iso cuts", n_eta_bins, eta_bins, n_pt_bins, pt_bins);
    TH2D *h_total = new TH2D("h_total", "passing iso cuts", n_eta_bins, eta_bins, n_pt_bins, pt_bins);

    TH2D *h_contam_rate = new TH2D("h_contam_rate", "passing iso cuts", n_eta_bins, eta_bins, n_pt_bins, pt_bins);
    TH2D *h_contam_total = new TH2D("h_contam_total", "passing iso cuts", n_eta_bins, eta_bins, n_pt_bins, pt_bins);

    TTree *t_data = (TTree *) f->Get("T_data");
    TTree *t_mc = (TTree *) f_mc->Get("T_data");



    construct_fakerate_template(h_rate, h_total, t_data, FLAG, true);
    construct_fakerate_template(h_contam_rate, h_contam_total, t_mc, FLAG,  false);

    //TH2D* h_rate = (TH2D *)f->Get("h_rate"); 
    //TH2D* h_total = (TH2D *)f->Get("h_total"); 

    TH2D* h_rate_new = (TH2D *)h_rate->Clone("h_rate_new");
    TH2D* h_total_new = (TH2D *)h_total->Clone("h_total_new");


    /*
    TH2D* h_contam_pass = (TH2D *)f_mc->Get("h_pass"); 
    TH2D* h_contam_total = (TH2D *)f_mc->Get("h_total"); 
    h_contam_pass->Scale(1000*tot_lumi);
    h_contam_total->Scale(1000*tot_lumi);
    TH2D* h_contam_rate = h_contam_pass->Clone("h_conam_rate");
    h_contam_rate->Divide(h_contam_total);
    */

    SetCorrectedRate(h_rate_new, h_total_new, h_contam_rate, h_contam_total);

    SetErrors(h_rate_new, h_total_new);
    if(write_out){
        TFile *f_new = TFile::Open("../analyze/FakeRate/root_files/2016/SingleMuon_data_fakerate_corrected_june25.root", "RECREATE");
        //TFile *f_new = TFile::Open("../analyze/FakeRate/root_files/2016/SingleElectron_data_fakerate_corrected_june25.root", "RECREATE");
        f_new->cd();
        h_rate_new->Write();
        h_total_new->Write();
        printf("Wrote out corrected rate \n");
    }

    TH1D *rate_barrel = h_rate_new->ProjectionY("rate_barrel", 1,1, "e");
    TH1D *rate_endcap = h_rate_new->ProjectionY("rate_endcap", 2,2, "e");

    TH1D *total_barrel = h_total_new->ProjectionY("total_bar", 1,1, "e");
    TH1D *total_endcap = h_total_new->ProjectionY("total_endcap", 2,2, "e");


    

    //SetErrors(rate_barrel, total_barrel);
    //SetErrors(rate_endcap, total_endcap);

    

    TCanvas *c1 = new TCanvas("c1", "Histograms", 200, 10, 900, 700);
    c1->SetLogx();
    if(FLAG == FLAG_MUONS) rate_barrel->SetTitle("Muons fakerate");
    else rate_barrel->SetTitle("Electrons fakerate");
    rate_barrel->SetStats(0);
    rate_barrel->SetLineWidth(3);
    rate_barrel ->Draw("E1");
    rate_barrel->GetXaxis()->SetTitle("p_t (Gev)");
    rate_barrel->GetXaxis()->SetRangeUser(10,500);
    rate_barrel->GetYaxis()->SetRangeUser(0, 0.7);
    c1->Update();

    //TCanvas *c2 = new TCanvas("c2", "Histograms", 200, 10, 900, 700);
    //rate_endcap->SetTitle("Fake-Rate Muons with Trigger: End Cap");
    rate_endcap->SetStats(0);
    rate_endcap->SetLineWidth(3);
    rate_endcap->Draw("E1 same");
    rate_endcap->SetLineColor(kRed);
    //rate_endcap->GetXaxis()->SetRangeUser(10,100);
    rate_endcap->GetYaxis()->SetRangeUser(0, 0.5);
    TLegend *leg1 = new TLegend(0.5, 0.65, 0.75, 0.8);
    leg1->AddEntry(rate_barrel, "Barrel", "f");
    leg1->AddEntry(rate_endcap, "Endcap",  "f");
    leg1->Draw();
    



    /*
    TCanvas *c3 = new TCanvas("c3", "Histograms", 200, 10, 900, 700);
    total_barrel->SetTitle("Number of events");
    total_barrel->SetStats(0);
    total_barrel->SetLineWidth(3);
    total_barrel ->Draw("hist");
    total_barrel->GetXaxis()->SetTitle("p_t (Gev)");
    total_barrel->GetXaxis()->SetRangeUser(10,100);
    c3->SetLogy();
    c3->SetLogx();

    total_endcap->SetTitle("Number of events");
    total_endcap->SetStats(0);
    total_endcap->SetLineWidth(3);
    total_endcap->SetLineColor(kRed);
    total_endcap->Draw("hist same");
    total_endcap->GetXaxis()->SetTitle("p_t (Gev)");
    total_endcap->GetXaxis()->SetRangeUser(10,100);
    printf("Totals: Bar %.0f End %.0f \n", total_barrel->Integral(), total_endcap->Integral());
    TLegend *leg2 = new TLegend(0.5, 0.65, 0.75, 0.8);
    leg2->AddEntry(total_barrel, "Barrel",  "f");
    leg2->AddEntry(total_endcap, "Endcap" , "f");
    leg2->Draw();
    */
}

    
    
