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
int n_cost_bins = 10;
Float_t cost_bins[] = {-1.0, -.8, -.6, -.4, -.2, 0., 0.2, 0.4, 0.6, 0.8, 1.0};

Float_t m_bins[] = {150, 200, 300, 400, 600, 1000};
const int n_m_bins = int(sizeof(m_bins)/sizeof(m_bins[0])) - 1; 
int n_fitted_events[n_m_bins];

Float_t m_min = m_bins[0];
Float_t m_max = m_bins[n_m_bins];
Float_t m_low, m_high;




Double_t alpha = 0.05;

TH3F *h_asym, *h_sym, *h_nosig, *h_ttbar, *h_data, *h_mc;

//define axis globally for convinience 
TAxis *x_ax, *y_ax, *z_ax;

//MC templates
TFile* f_mc = (TFile*) TFile::Open("output_files/DYToLL_mc_2016_apr11.root");
TTree *t_mc = (TTree *) f_mc ->Get("T_data");
TTree *t_nosig = (TTree *) f_mc ->Get("T_back");
TFile* f_ttbar = (TFile*) TFile::Open("output_files/ttbar_background_apr11.root");
TTree *t_ttbar = (TTree *) f_ttbar ->Get("T_data");

TFile *f_data = TFile::Open("output_files/DYToLL_data_2016_apr11.root");
TTree *t_data = (TTree *)f_data->Get("T_data"); 


vector<double> v_xF;
vector<double> v_m;
vector<double> v_cost;

Double_t AFB_fit[n_m_bins], AFB_err[n_m_bins], r_ttbar_fit[n_m_bins], r_ttbar_err[n_m_bins];

unsigned int nDataEvents;

Double_t get_prob(Double_t xF, Double_t cost, Double_t m, TH3F *h){
    TAxis* x_ax =  h->GetXaxis();
    TAxis *y_ax =  h->GetYaxis();
    TAxis *z_ax =  h->GetZaxis();

    int xbin = x_ax->FindBin(xF);
    int ybin = y_ax->FindBin(cost);
    int zbin = z_ax->FindBin(m);

    return h->GetBinContent(xbin, ybin, zbin);
}


// fcn passes back f = - 2*ln(L), (where L is the likelihood of the event)
void fcn(int& npar, double* deriv, double& f, double par[], int flag){
    double lnL = 0.0;

    int nFittedEvents = 0;
    for (int i=0; i<nDataEvents; i++){
        if(v_m[i] >= m_low && v_m[i] <= m_high){
            nFittedEvents++;
            Double_t p_sym = get_prob(v_xF[i], v_cost[i], v_m[i], h_sym);
            Double_t p_asym = get_prob(v_xF[i],  v_cost[i], v_m[i], h_asym);
            Double_t p_nosig = get_prob(v_xF[i], v_cost[i], v_m[i], h_nosig);
            Double_t p_ttbar = get_prob(v_xF[i], v_cost[i], v_m[i], h_ttbar);

            double AFB = par[0];
            double r_nosig = par[1];
            double r_ttbar = par[2];
            double prob = r_nosig*p_nosig + r_ttbar*p_ttbar + (1 - r_nosig - r_ttbar) * (p_sym + AFB*p_asym);
            if(prob > 1) printf("Warning prob is too big \n");
            prob = max(prob, 1e-20);
            if(prob >0.0) lnL += log(prob);
            else printf("Warning, prob is negative \n");
        }
    }
    int j=0;
    while(m_bins[j] != m_low) j++;
    n_fitted_events[j] = nFittedEvents;


    f = -2.0 * lnL;

}


int gen_template(TTree *t1, TH3F* h, bool is_data=false){
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
        if(m >= m_min && m <= m_max && met_pt < 50. && jet1_cmva < 0. && jet2_cmva < 0.){
            n++;
            h->Fill(xF, cost, m, gen_weight); 
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
            nEvents, m_low, m_max);
    h->Scale(1./h->Integral());
    t1->ResetBranchAddresses();
    return 0;
}


int gen_mc_template(TTree *t1, TH3F* h_sym, TH3F *h_asym, bool print=false){
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
        if(m >= m_min && m <= m_max && met_pt < 50. && jet1_cmva < 0. && jet2_cmva < 0.){
            reweight = (4./3.)*cost_st*(2. + alpha)/
                (1. + cost_st*cost_st + alpha*(1.- cost_st*cost_st));
            n++;
            h_sym->Fill(xF, cost, m, gen_weight); 
            h_sym->Fill(xF, -cost, m, gen_weight); 
            h_asym->Fill(xF, cost, m, reweight * gen_weight);
            h_asym->Fill(xF, -cost, m, -reweight * gen_weight);
        }
    }
    printf("N sym is %i \n", n);
    float norm = h_sym -> Integral();
    h_sym->Scale(1./norm);
    h_asym->Scale(1./norm);
    t1->ResetBranchAddresses();
    return 0;
}


int gen_ttbar_template(TTree *t1, TH3F* h, bool is_data=false){
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
        if(m >= m_min && m <= m_max && met_pt < 50.){
            h->Fill(xF, cost, m, gen_weight); 
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
    h_mc = new TH3F("h_mc", "MC DY distribution xF > 0.15",
            n_xf_bins, xf_bins, n_cost_bins, cost_bins, n_m_bins, m_bins);
    h_sym = new TH3F("h_sym", "Symmetric template of mc (xF, cost_r, m) xF>0.15",
            n_xf_bins, xf_bins, n_cost_bins, cost_bins, n_m_bins, m_bins);
    h_asym = new TH3F("h_asym", "Asymmetric template of mc (xF cost_r, m) xF > 0.15",
            n_xf_bins, xf_bins, n_cost_bins, cost_bins, n_m_bins, m_bins);
    h_ttbar = new TH3F("h_ttbar", "TTBar with met < 50 and CMVA < 0",
            n_xf_bins, xf_bins, n_cost_bins, cost_bins, n_m_bins, m_bins);
    h_nosig = new TH3F("h_nosig", "Template of no asymetry drell-yan events from mc (xF, cost_r,m )",
            n_xf_bins, xf_bins, n_cost_bins, cost_bins, n_m_bins, m_bins);

    h_data = new TH3F("h_data", "Data template of (x_f, cost_r,m) xF > 0.15",
            n_xf_bins, xf_bins, n_cost_bins, cost_bins, n_m_bins, m_bins);
    gen_mc_template(t_mc, h_sym, h_asym);
    gen_template(t_mc, h_mc);
    gen_template(t_nosig, h_nosig);
    gen_ttbar_template(t_ttbar, h_ttbar);

    gen_template(t_data, h_data, true);
    h_data->Sumw2();
    printf("Finishing setup \n");
    return;
}

void MuMu_fit_all(){
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

    TVirtualFitter *minuit;
    for(int i = 0; i < n_m_bins; i++){

        m_low = m_bins[i];
        m_high = m_bins[i+1];

        printf("I=%i. Fitting in range %.0f %.0f \n", i, m_low, m_high);


        minuit = TVirtualFitter::Fitter(0,3);
        minuit->SetFCN(fcn);
        minuit->SetParameter(0,"AFB", AFB_start, AFB_start_error, -AFB_max, AFB_max);
        minuit->SetParameter(1,"r_nosig", r_nosig_start, r_nosig_start_error, 0, r_nosig_max);
        minuit->SetParameter(2,"r_ttbar", r_ttbar_start, r_ttbar_start_error, 0, r_ttbar_max);
        Double_t arglist[100];
        arglist[0] = 10000.;
        minuit->ExecuteCommand("MIGRAD", arglist,0);



        
        AFB_fit[i] = minuit->GetParameter(0);
        AFB_err[i] = minuit->GetParError(0);
        r_ttbar_fit[i] = minuit->GetParameter(2); 
        r_ttbar_err[i] = minuit->GetParError(2); 
        minuit->Clear();
    }


    for(int i=0; i<n_m_bins; i++){
        printf("I=%i M=[%.0f, %.0f] AFB = %0.2f +/- %0.2f r_ttbar = %0.2f +/- %0.2f \n", 
                i, m_bins[i], m_bins[i+1], AFB_fit[i], AFB_err[i], r_ttbar_fit[i], r_ttbar_err[i]);

    }

    /*
    f_data->Close();
    f_ttbar->Close();
    f_mc->Close();
    */

}



