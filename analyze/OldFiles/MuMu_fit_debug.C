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
#include "TPostScript.h"
#include "TLorentzVector.h"
#include "TVector3.h"
#include "TFitter.h"
#include "TSystem.h"
//#include"Minuit2/Minuit2Minimizer.h"
#include "Math/Functor.h"



int n_xf_bins = 3;
float xf_max = 0.5;
Float_t xf_bins[] = {0., 0.1, 0.2, 0.3, 0.4};
int n_m_bins = 1;
float m_max = 1000.;
int n_cost_bins = 10;
Float_t cost_bins[] = {-1.0, -.8, -.6, -.4, -.2, 0., 0.2, 0.4, 0.6, 0.8, 1.0};

float m_low =  150;
float m_high = 700;
Double_t alpha = 0.05;

TH2F *h_asym, *h_sym, *h_nosig, *h_ttbar, *h_data;
TH2F *h_mc_count, *h_nosig_count, *h_ttbar_count, *h_data_count;
//void *x_axp, *y_axp, *z_axp;

//define axis globally for convinience 
TAxis *x_ax, *y_ax, *z_ax;

//MC templates
TFile* f_mc = (TFile*) TFile::Open("output_files/DYToLL_mc_2016_apr10.root");
TTree *t_mc = (TTree *) f_mc ->Get("T_data");
TTree *t_nosig = (TTree *) f_mc ->Get("T_back");
TFile* f_ttbar = (TFile*) TFile::Open("output_files/ttbar_background_apr10.root");
TTree *t_ttbar = (TTree *) f_ttbar ->Get("T_data");

TFile *f_data = TFile::Open("output_files/DYToLL_data_2016_apr4.root");
TTree *t_data = (TTree *)f_data->Get("T_data"); 


vector<double> v_xF;
vector<double> v_m;
vector<double> v_cost;
unsigned int nDataEvents;
bool first_run = false;

Double_t get_prob(Double_t xF, Double_t cost, TH2F *h){
    TAxis* x_ax =  h->GetXaxis();
    TAxis *y_ax =  h->GetYaxis();
    //binning the same in all templates
    int xbin = x_ax->FindBin(xF);
    int ybin = y_ax->FindBin(cost);

    return h->GetBinContent(xbin, ybin);
}

Double_t get_bin_width(Double_t xF, Double_t cost, TH2F *h){
    TAxis *x_ax =  h->GetXaxis();
    TAxis *y_ax =  h->GetYaxis();
    int xbin = x_ax->FindBin(xF);
    int ybin = y_ax->FindBin(cost);
    Double_t xbin_width = x_ax->GetBinWidth(xbin);
    Double_t ybin_width = y_ax->GetBinWidth(ybin);

    //return xbin_width * ybin_width;
    return 1.0;
}

// fcn passes back f = - 2*ln(L), (where L is the likelihood of the event)
void fcn(int& npar, double* deriv, double& f, double par[], int flag){
    double lnL = 0.0;

    for (int i=0; i<nDataEvents; i++){
        Double_t p_sym = get_prob(v_xF[i], v_cost[i], h_sym);
        Double_t p_asym = get_prob(v_xF[i],  v_cost[i], h_asym);
        Double_t p_nosig = get_prob(v_xF[i], v_cost[i], h_nosig);
        Double_t p_ttbar = get_prob(v_xF[i], v_cost[i], h_ttbar);

        double AFB = par[0];
        double r_nosig = par[1];
        double r_ttbar = par[2];
        double prob = r_nosig*p_nosig + r_ttbar*p_ttbar + (1 - r_nosig - r_ttbar) * (p_sym + AFB*p_asym);
        if(prob > 1) printf("Warning prob is too big \n");
        prob = max(prob, 1e-20);
        if(prob >0.0) lnL += log(prob);
        else printf("Warning, prob is negative \n");
    }
    f = -2.0 * lnL;

}

Double_t  template_fcn(Double_t *x, Double_t *par){
    Double_t xx = x[0];
    Double_t yy = x[1];

    Double_t p_sym = get_prob(xx, yy,  h_sym);
    Double_t p_asym = get_prob(xx, yy, h_asym);
    Double_t p_nosig = get_prob(xx, yy,  h_nosig);
    Double_t p_ttbar = get_prob(xx, yy,  h_ttbar);

    Double_t AFB = par[0];
    Double_t r_nosig = par[1];
    Double_t r_ttbar = par[2];
    Double_t prob = r_nosig*p_nosig + r_ttbar*p_ttbar + (1 - r_nosig - r_ttbar) * (p_sym + AFB*p_asym);
    return prob;
}


int gen_template(TTree *t1, TH2F* h, TH2F* h_count,  bool is_data=false){
    Long64_t nEntries  =  t1->GetEntries();
    //printf("size is %i \n", nEntries);
    Double_t m, xF, cost, gen_weight, jet1_cmva, jet2_cmva;
    Float_t met_pt;
    t1->SetBranchAddress("m", &m);
    t1->SetBranchAddress("xF", &xF);
    t1->SetBranchAddress("cost", &cost);
    t1->SetBranchAddress("jet1_CMVA", &jet1_cmva);
    t1->SetBranchAddress("jet2_CMVA", &jet2_cmva);
    t1->SetBranchAddress("met_pt", &met_pt);
    if(!is_data) t1->SetBranchAddress("gen_weight", &gen_weight);
    else{
        gen_weight = 1.0;
        printf("total sample size is %i \n", (int) nEntries);
    }
    int nEvents = 0;
    int n=0;
    for (int i=0; i<nEntries; i++) {
        t1->GetEntry(i);
        if(m >= m_low && m <= m_high && met_pt < 50. && jet1_cmva < 0. && jet2_cmva < 0.){
            n++;
            Double_t bin_width = get_bin_width(xF, cost, h);
            h->Fill(xF, cost, gen_weight/bin_width); 
            h_count ->Fill(xF, cost, 1);
            if(is_data){
                nDataEvents++;
                v_m.push_back(m);
                v_xF.push_back(xF);
                v_cost.push_back(cost);
            }
            nEvents++;
        }
    }
    //printf("N MC is %i \n", n);
    if(is_data) printf("Fit will be run on %i events \n", nEvents);
    if(is_data) h->Sumw2();
    h->Scale(1./h->Integral());
    t1->ResetBranchAddresses();
    return 0;
}


int gen_mc_template(TTree *t1, TH2F* h_sym, TH2F *h_asym, TH2F *h_count, bool print=false){
    Long64_t nEntries  =  t1->GetEntries();
    //printf("size is %i \n", nEntries);
    Double_t m, xF, cost, gen_weight, reweight, jet1_cmva, jet2_cmva, cost_st;
    Float_t met_pt;
    t1->SetBranchAddress("m", &m);
    t1->SetBranchAddress("xF", &xF);
    t1->SetBranchAddress("cost", &cost);
    t1->SetBranchAddress("cost_st", &cost_st);
    t1->SetBranchAddress("jet1_CMVA", &jet1_cmva);
    t1->SetBranchAddress("jet2_CMVA", &jet2_cmva);
    t1->SetBranchAddress("met_pt", &met_pt);
    t1->SetBranchAddress("gen_weight", &gen_weight);
    //t1->SetBranchAddress("reweight", &reweight);
    int n = 0;
    for (int i=0; i<nEntries; i++) {
        t1->GetEntry(i);
        if(m >= m_low && m <= m_high && met_pt < 50. && jet1_cmva < 0. && jet2_cmva < 0. ){
            reweight = (4./3.)*cost_st*(2. + alpha)/
                       (1. + cost_st*cost_st + alpha*(1.- cost_st*cost_st));
            n++;
            Double_t bin_width = get_bin_width(xF, cost, h_sym);
            h_sym->Fill(xF, cost, gen_weight/bin_width); 
            h_sym->Fill(xF, -cost, gen_weight/bin_width); 
            h_asym->Fill(xF, cost, reweight * gen_weight/bin_width);
            h_asym->Fill(xF, -cost, -reweight * gen_weight/bin_width);
            h_count->Fill(xF, cost, 1);
        }
    }
    printf("N sym is %i \n", n);
    float norm = h_sym -> Integral();
    h_sym->Scale(1./norm);
    h_asym->Scale(1./norm);
    t1->ResetBranchAddresses();
    return 0;
}


int gen_ttbar_template(TTree *t1, TH2F* h, TH2F* h_count,  bool is_data=false){
    Long64_t nEntries  =  t1->GetEntries();
    Double_t m, xF, cost, gen_weight, jet1_cmva, jet2_cmva;
    Float_t met_pt;
    t1->SetBranchAddress("m", &m);
    t1->SetBranchAddress("xF", &xF);
    t1->SetBranchAddress("cost", &cost);
    t1->SetBranchAddress("jet1_CMVA", &jet1_cmva);
    t1->SetBranchAddress("jet2_CMVA", &jet2_cmva);
    t1->SetBranchAddress("met_pt", &met_pt);
    if(!is_data) t1->SetBranchAddress("gen_weight", &gen_weight);
    else{
        gen_weight = 1.0;
        printf("total sample size is %i \n", (int) nEntries);
    }
    int nEvents = 0;
    for (int i=0; i<nEntries; i++) {
        t1->GetEntry(i);
        if(m >= m_low && m <= m_high && met_pt < 50.){
            Double_t bin_width = get_bin_width(xF, cost, h);
            h->Fill(xF, cost, gen_weight/bin_width); 
            h_count ->Fill(xF, cost, 1);
            nEvents++;
        }
    }
    h->Scale(1./h->Integral());
    t1->ResetBranchAddresses();
    return 0;
}

void setup(){
    //setup global variables
    //TH1::SetDefaultSumw2(kTRUE);
    h_mc_count = new TH2F("h_mc_count", "Events in bins for MC templates",
            n_xf_bins, xf_bins, n_cost_bins, cost_bins);
    h_sym_count = new TH2F("h_sym_count", "Events in bins for MC templates",
            n_xf_bins, xf_bins, n_cost_bins, cost_bins);
    h_mc = new TH2F("h_mc", "MC DY distribution",
            n_xf_bins, xf_bins, n_cost_bins, cost_bins);
    h_sym = new TH2F("h_sym", "Symmetric template of mc (xF, cost_r) ",
            n_xf_bins, xf_bins, n_cost_bins, cost_bins);
    h_asym = new TH2F("h_asym", "Asymmetric template of mc (xF cost_r) ",
            n_xf_bins, xf_bins, n_cost_bins, cost_bins);
    h_ttbar = new TH2F("h_ttbar", "TTBar with met < 50 and CMVA < 0",
            n_xf_bins, xf_bins, n_cost_bins, cost_bins);
    h_ttbar_count = new TH2F("h_ttbar_count", "Events in bins for ttbar template",
            n_xf_bins, xf_bins, n_cost_bins, cost_bins);
    h_nosig_count = new TH2F("h_nosig_count", "Events in bins for nosig template",
            n_xf_bins, xf_bins, n_cost_bins, cost_bins);
    h_nosig = new TH2F("h_nosig", "Template of no asymetry drell-yan events from mc (xF, cost_r)",
            n_xf_bins, xf_bins, n_cost_bins, cost_bins);

    h_data = new TH2F("h_data", "Data template of (x_f, cost_r) ",
            n_xf_bins, xf_bins, n_cost_bins, cost_bins);
    h_data_count = new TH2F("h_data_count", "Data template of (x_f,  c_r)",
            n_xf_bins, xf_bins, n_cost_bins, cost_bins);
    gen_mc_template(t_mc, h_sym, h_asym, h_sym_count);
    gen_template(t_mc, h_mc, h_mc_count);
    gen_template(t_nosig, h_nosig, h_nosig_count);
    gen_ttbar_template(t_ttbar, h_ttbar, h_ttbar_count);

    gen_template(t_data, h_data, h_data_count, true);
    printf("Finishing setup \n");
    return;
}

void plot_things(Double_t AFB, Double_t r_ttbar, Double_t r_nosig){
    
    TH1D *h_data_cost = h_data->ProjectionY("h_data_cost");
    TH1D *h_ttbar_cost = h_ttbar->ProjectionY("h_ttbar_cost");
    TH1D *h_sym_cost = h_sym->ProjectionY("h_sym_cost");
    TH1D *h_asym_cost = h_asym->ProjectionY("h_asym_cost");
    TH1D *h_nosig_cost = h_nosig->ProjectionY("h_nosig_cost");

    TString title = TString::Format("Fitted cost vs data, AFB=%0.3f, M in [%1.0f, %1.0f]", 
                                    AFB, m_low, m_high);


    TH1D *h_fit = new TH1D("h_fit", title, n_cost_bins, cost_bins);

    
    /*
    TCanvas *c1 = new TCanvas("c1", "canvas", 200,10, 900,700);
    h_ttbar_nocut_count->SetBit(TH2::kNoStats);
    h_ttbar_nocut_count->Draw("colz");
    h_data->SetMarkerStyle(20);
    h_data->SetMarkerSize(0.5);
    h_data->Draw("scat same");
    c1->Update();


    TCanvas *c2 = new TCanvas("c2", "canva", 100,100, 700,700);
    TH1D *h_ttbar_xf_nocut = h_ttbar_nocut->ProjectionX("h_ttbar_nocut_xf");
    h_ttbar_xf_nocut->SetLineColor(kBlue);
    h_ttbar_xf_nocut->SetLineWidth(0.2);
    h_ttbar_xf_nocut->Draw();
    TH1D *h_ttbar_xf = h_ttbar->ProjectionX("h_ttbar_xf");
    h_ttbar_xf->SetLineColor(kRed);
    h_ttbar_xf->SetLineWidth(0.2);
    h_ttbar_xf->Draw("same");
    c2->Update();
    c2->BuildLegend();

    TCanvas *c3 = new TCanvas("c3", "canva", 100,100, 700,700);
    TH1D *h_ttbar_cost_nocut = h_ttbar_nocut->ProjectionY("h_ttbar_cost_nocut");
    h_ttbar_cost_nocut->SetBit(TH1::kNoStats);
    h_ttbar_cost_nocut->SetLineColor(kBlue);
    h_ttbar_cost_nocut->SetLineWidth(0.2);
    h_ttbar_cost_nocut->Draw("same");
    TH1D *h_ttbar_cost = h_ttbar->ProjectionY("h_ttbar_cost");
    h_ttbar_cost->SetLineColor(kRed);
    h_ttbar_cost->SetLineWidth(0.2);
    h_ttbar_cost->Draw("same");
    c3->Update();
    c3->BuildLegend();
    */



    h_ttbar_cost ->Scale(r_ttbar);
    h_nosig_cost ->Scale(r_nosig);
    h_sym_cost ->Scale((1. - r_nosig - r_ttbar));
    h_asym_cost ->Scale((1. - r_nosig - r_ttbar)*AFB);

    h_fit->Add(h_ttbar_cost);
    h_fit->Add(h_sym_cost);
    h_fit->Add(h_asym_cost);
    
    h_ttbar_cost ->SetFillColor(kGreen);
    h_nosig_cost ->SetFillColor(kYellow);
    h_sym_cost ->SetFillColor(kBlue);
    h_asym_cost ->SetFillColor(kRed);

    THStack *cost_fit = new THStack("cost_fit", title);
    cost_fit->Add(h_nosig_cost);
    cost_fit->Add(h_ttbar_cost);
    cost_fit->Add(h_sym_cost);
    cost_fit->Add(h_asym_cost);

    
    TCanvas *c5 = new TCanvas("c5", "data vs fit", 200, 10, 900,700);
    c5->cd();
    h_data_cost->SetMarkerStyle(kFullCircle);
    h_data_cost->SetMarkerColor(1);
    cost_fit->SetMaximum(0.25);
    cost_fit->Draw();
    h_data_cost->Draw("E1 same");
    c5->Update();
    gPad->BuildLegend();
    

    TCanvas *c6 = new TCanvas("c6", "data vs fit", 200, 10, 900,700);
    c6->cd();
    h_fit->SetMaximum(0.25);
    h_fit->SetFillColor(kBlue);
    h_fit->Draw();
    h_data_cost->Draw("P same");
    c6->Update();
    gPad->BuildLegend();


    /*
    TCanvas *c6 = new TCanvas("c6", "data vs fit", 200, 10, 900,700);
    h_fit->SetFillColor(kBlue);
    h_fit->SetMaximum(0.25);
    h_data_cost->SetMarkerStyle(kFullCircle);
    h_data_cost->SetMarkerColor(1);
    h_fit->Draw();
    h_data_cost->Draw("P same");
    c6->Update();
    */





    TCanvas *c32 = new TCanvas("c32", "canva", 100,100, 700,700);
    //h_data_cost->SetBit(TH1::kNoStats);
    h_data_cost->Draw();



    TCanvas *c33 = new TCanvas("c33", "canva", 100,100, 700,700);
    //h_asym_cost->SetBit(TH1::kNoStats);
    h_asym_cost->Draw();

    TCanvas *c34 = new TCanvas("c34", "canva", 100,100, 700,700);
    //h_asym_cost->SetBit(TH1::kNoStats);
    h_sym_cost->Draw();


}

void MuMu_fit_debug(){
    setup();
    Double_t AFB_fit, r_ttbar_fit, r_nosig_fit;
    AFB_fit = 0.618;
    r_ttbar_fit = 0.1238;
    r_nosig_fit = 2.644e-10;
    plot_things(AFB_fit, r_ttbar_fit, r_nosig_fit);

    printf("Integrals are %f %f %f %f %f \n", h_data->Integral(), h_sym->Integral(), 
                                           h_asym->Integral(), h_ttbar->Integral(), 
                                           h_nosig->Integral());


    /*
    float AFB_start = 0.4;
    float AFB_start_error = 0.1;
    float AFB_max = 0.75;
    float r_ttbar_start = 0.05;
    float r_ttbar_start_error = 0.05;
    float r_ttbar_max = 0.2;
    float r_nosig_start = 0.005;
    float r_nosig_start_error = 0.005;
    float r_nosig_max = 0.02;

    TVirtualFitter * minuit = TVirtualFitter::Fitter(0,3);
    minuit->SetFCN(fcn);
    minuit->SetParameter(0,"AFB", AFB_start, AFB_start_error, -AFB_max, AFB_max);
    minuit->SetParameter(1,"r_nosig", r_nosig_start, r_nosig_start_error, 0, r_nosig_max);
    minuit->SetParameter(2,"r_ttbar", r_ttbar_start, r_ttbar_start_error, 0, r_ttbar_max);
    Double_t arglist[100];
    arglist[0] = 10000.;
    minuit->ExecuteCommand("MIGRAD", arglist,0);


    printf("Trying template fit \n\n\n");
    Int_t npar = 3;
    Int_t ndim = 2;
    TF2 *fit_fcn = new TF2("F1", &template_fcn, 0., 1., -1.,1., npar, ndim);
    fit_fcn->SetParNames("AFB", "r_nosig", "r_ttbar");
    fit_fcn->SetParameter(0, AFB_start); //AFB
    fit_fcn->SetParLimits(0, -AFB_max, AFB_max);
    fit_fcn->SetParError(0, AFB_start_error);
    fit_fcn->SetParameter(1, r_nosig_start); //r_nosig
    fit_fcn->SetParLimits(1, 0, r_nosig_max);
    fit_fcn->SetParError(1, r_nosig_start_error);
    fit_fcn->SetParameter(2, r_ttbar_start); //r_ttbar
    fit_fcn->SetParLimits(2, 0, r_ttbar_max);
    fit_fcn->SetParError(2, r_ttbar_start_error);
    h_data->Sumw2();
    h_data->Fit(fit_fcn, "WL M");
    */




    /*
    f_data->Close();
    f_ttbar->Close();
    f_mc->Close();
    */

}



