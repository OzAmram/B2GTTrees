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



int n_xf_bins = 5;
float xf_max = 1.0;
Float_t xf_bins[] = {0., 0.05, 0.1, 0.15, 0.25, 1.0};
int n_m_bins = 1;
float m_max = 1000.;
int n_cost_bins = 10;
Float_t cost_bins[] = {-1.0, -.8, -.6, -.4, -.2, 0., 0.2, 0.4, 0.6, 0.8, 1.0};



float m_low = 150;
float m_high = 100000;

Double_t alpha = 0.05;

Double_t med_btag = 0.4432;

TH2F *h_asym, *h_sym, *h_nosig, *h_ttbar, *h_data, *h_mc;
TH2F *h_mc_count, *h_nosig_count, *h_ttbar_count, *h_data_count;
//void *x_axp, *y_axp, *z_axp;

//define axis globally for convinience 
TAxis *x_ax, *y_ax, *z_ax;

//MC templates
TFile* f_mc = (TFile*) TFile::Open("output_files/DYToLL_mc_2016_may9.root");
TTree *t_mc = (TTree *) f_mc ->Get("T_data");
TTree *t_nosig = (TTree *) f_mc ->Get("T_back");
TFile* f_ttbar = (TFile*) TFile::Open("output_files/ttbar_background_may9.root");
TTree *t_ttbar = (TTree *) f_ttbar ->Get("T_data");

TFile *f_data = TFile::Open("output_files/DYToLL_data_samesign_may11.root");
TTree *t_data = (TTree *)f_data->Get("T_data"); 


vector<double> v_xF;
vector<double> v_m;
vector<double> v_cost;
unsigned int nDataEvents;

Double_t get_prob(Double_t xF, Double_t cost, TH2F *h){
    TAxis* x_ax =  h->GetXaxis();
    TAxis *y_ax =  h->GetYaxis();
    //binning the same in all templates
    int xbin = x_ax->FindBin(xF);
    int ybin = y_ax->FindBin(cost);

    return h->GetBinContent(xbin, ybin);
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
        if(m >= m_low && m <= m_high && met_pt < 50. && jet1_cmva < med_btag 
                && jet2_cmva < med_btag){
            n++;
            h->Fill(xF, cost, gen_weight); 
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
    if(is_data) printf("Fit will be run on %i events in mass range [%1.0f, %1.0f] \n", 
            nEvents, m_low, m_high);
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
    int n = 0;
    for (int i=0; i<nEntries; i++) {
        t1->GetEntry(i);
        if(m >= m_low && m <= m_high && met_pt < 50.
                && jet1_cmva < med_btag && jet2_cmva < med_btag){
            reweight = (4./3.)*cost_st*(2. + alpha)/
                (1. + cost_st*cost_st + alpha*(1.- cost_st*cost_st));
            n++;
            h_sym->Fill(xF, cost, gen_weight); 
            h_sym->Fill(xF, -cost, gen_weight); 
            h_asym->Fill(xF, cost, reweight * gen_weight);
            h_asym->Fill(xF, -cost, -reweight * gen_weight);
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
        if(m >= m_low && m <= m_high && met_pt < 50. 
                && jet1_cmva < med_btag && jet2_cmva < med_btag){
            h->Fill(xF, cost, gen_weight); 
            h_count ->Fill(xF, cost, 1);
            nEvents++;
        }
    }
    printf("N ttbar events %i \n", nEvents);
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
    h_mc = new TH2F("h_mc", "MC DY distribution xF > 0.15",
            n_xf_bins, xf_bins, n_cost_bins, cost_bins);
    h_sym = new TH2F("h_sym", "Symmetric template of mc (xF, cost_r) xF>0.15",
            n_xf_bins, xf_bins, n_cost_bins, cost_bins);
    h_asym = new TH2F("h_asym", "Asymmetric template of mc (xF cost_r) xF > 0.15",
            n_xf_bins, xf_bins, n_cost_bins, cost_bins);
    h_ttbar = new TH2F("h_ttbar", "TTBar with met < 50 and CMVA < 0",
            n_xf_bins, xf_bins, n_cost_bins, cost_bins);
    h_ttbar_count = new TH2F("h_ttbar_count", "Events in bins for ttbar template",
            n_xf_bins, xf_bins, n_cost_bins, cost_bins);
    h_nosig_count = new TH2F("h_nosig_count", "Events in bins for nosig template",
            n_xf_bins, xf_bins, n_cost_bins, cost_bins);
    h_nosig = new TH2F("h_nosig", "Template of no asymetry drell-yan events from mc (xF, cost_r)",
            n_xf_bins, xf_bins, n_cost_bins, cost_bins);

    h_data = new TH2F("h_data", "Data template of (x_f, cost_r) xF > 0.15",
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

void MuMu_fit_new(){
    setup();

    printf("Integrals are %f %f %f %f %f \n", h_data->Integral(), h_sym->Integral(), 
                                           h_asym->Integral(), h_ttbar->Integral(), 
                                           h_nosig->Integral());



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
    h_data->Fit(fit_fcn, "WL M N");


    Double_t AFB_fit, r_ttbar_fit, r_nosig_fit;
    Double_t AFB_fit_error, r_ttbar_fit_error, r_nosig_fit_error;
    AFB_fit = minuit->GetParameter(0); 
    AFB_fit_error = minuit->GetParError(0);
    r_ttbar_fit = minuit->GetParameter(2); 
    r_ttbar_fit_error = minuit->GetParError(2);
    r_nosig_fit = minuit->GetParameter(1); 
    r_nosig_fit_error = minuit->GetParError(1);



    printf("\n\n\n Fit on M=[%.0f, %.0f], %i Events: AFB = %0.3f +/- %0.3f r_ttbar = %0.3f +/- %0.3f r_nosig = %0.4f +/- %0.4f \n", 
                m_low, m_high, nDataEvents, AFB_fit, AFB_fit_error, r_ttbar_fit, r_ttbar_fit_error, r_nosig_fit, r_nosig_fit_error);


    /*
    f_data->Close();
    f_ttbar->Close();
    f_mc->Close();
    */

}



