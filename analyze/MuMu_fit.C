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
#include"Minuit2/Minuit2Minimizer.h"
#include "Math/Functor.h"



bool rebin = true;

int n_xf_bins = 5;
float xf_max = 0.5;
int n_m_bins = 1;
float m_max = 1000.;
int n_cost_bins = 10;

float m_low = 200;
float m_high = 300;

TH3F *h_asym, *h_sym, *h_nosig, *h_ttbar;
//void *x_axp, *y_axp, *z_axp;

//define axis globally for convinience 
TAxis *x_ax, *y_ax, *z_ax;

//MC templates
TFile* f_mc = (TFile*) TFile::Open("output_files/DYToLL_mc_2016_apr10.root");
TTree *t_mc = (TTree *) f_mc ->Get("T_data");
TTree *t_nosig = (TTree *) f_mc ->Get("T_back");
TFile* f_ttbar = (TFile*) TFile::Open("output_files/ttbar_background_apr10.root");
TTree *t_ttbar = (TTree *) f_ttbar ->Get("T_data");

int gen_template(TTree *t1, TH3F* h, bool print=false);
int gen_mc_template(TTree *t1, TH3F* h_sym, TH3F* h_asym, bool print=false);

vector<double> v_xF;
vector<double> v_m;
vector<double> v_cost;
unsigned int nDataEvents;
bool first_run = false;

Double_t get_prob(Double_t xF, Double_t m, Double_t cost, TH3F *h){
    TAxis* x_ax =  h_sym->GetXaxis();
    TAxis *y_ax =  h_sym->GetYaxis();
    TAxis *z_ax =  h_sym->GetZaxis();
    //binning the same in all templates
    int xbin = x_ax->FindBin(xF);
    int ybin = y_ax->FindBin(m);
    int zbin = z_ax->FindBin(cost);

    return h->GetBinContent(xbin, ybin,zbin);
}

// fcn passes back f = - 2*ln(L), (where L is the likelihood of the event)
void fcn(int& npar, double* deriv, double& f, double par[], int flag){
    double lnL = 0.0;

    for (int i=0; i<nDataEvents; i++){
        Double_t p_sym = get_prob(v_xF[i], v_m[i], v_cost[i], h_sym);
        Double_t p_asym = get_prob(v_xF[i], v_m[i], v_cost[i], h_asym);
        Double_t p_nosig = get_prob(v_xF[i], v_m[i], v_cost[i], h_nosig);
        Double_t p_ttbar = get_prob(v_xF[i], v_m[i], v_cost[i], h_ttbar);

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

double template_fcn(Double_t *x, Double_t *par){
    Double_t xx = x[0];
    Double_t yy = x[1];
    Double_t zz = x[2];

    Double_t p_sym = get_prob(xx, yy, zz, h_sym);
    Double_t p_asym = get_prob(xx, yy, zz, h_asym);
    Double_t p_nosig = get_prob(xx, yy, zz, h_nosig);
    Double_t p_ttbar = get_prob(xx, yy, zz, h_ttbar);

    double AFB = par[0];
    double r_nosig = par[1];
    double r_ttbar = par[2];
    double prob = r_nosig*p_nosig + r_ttbar*p_ttbar + (1 - r_nosig - r_ttbar) * (p_sym + AFB*p_asym);
    return prob;
}

int gen_template(TTree *t1, TH3F* h, bool print=false){
    Long64_t size  =  t1->GetEntries();
    Double_t m, xF, cost, gen_weight, jet1_cmva, jet2_cmva;
    Float_t met_pt;
    t1->SetBranchAddress("m", &m);
    t1->SetBranchAddress("xF", &xF);
    t1->SetBranchAddress("cost", &cost);
    t1->SetBranchAddress("gen_weight", &gen_weight);
    t1->SetBranchAddress("jet1_CMVA", &jet1_cmva);
    t1->SetBranchAddress("jet2_CMVA", &jet2_cmva);
    t1->SetBranchAddress("met_pt", &met_pt);
    int nBack = 0;
    int nForward = 0;
    int nB=0;
    int nF=0;
    int nEmpties=0;
    for (int i=0; i<size; i++) {
        t1->GetEntry(i);
        if(m >= m_low && m <= m_high && met_pt < 50. && jet1_cmva < 0. && jet2_cmva < 0.){
            h->Fill(xF, m, cost, gen_weight); 
        }
    }
    h->Scale(1./h->Integral());
    return 0;
}

int gen_mc_template(TTree *t1, TH3F* h_sym, TH3F *h_asym, bool print=false){
    Long64_t size  =  t1->GetEntries();
    Double_t m, xF, cost, gen_weight, reweight, jet1_cmva, jet2_cmva;
    Float_t met_pt;
    t1->SetBranchAddress("m", &m);
    t1->SetBranchAddress("xF", &xF);
    t1->SetBranchAddress("cost", &cost);
    t1->SetBranchAddress("gen_weight", &gen_weight);
    t1->SetBranchAddress("reweight", &reweight);
    t1->SetBranchAddress("jet1_CMVA", &jet1_cmva);
    t1->SetBranchAddress("jet2_CMVA", &jet2_cmva);
    t1->SetBranchAddress("met_pt", &met_pt);
    for (int i=0; i<size; i++) {
        t1->GetEntry(i);
        if(m >= m_low && m <= m_high && met_pt < 50. && jet1_cmva < 0. && jet2_cmva < 0.){
            h_sym->Fill(xF, m, cost, gen_weight); 
            h_sym->Fill(xF, m, -cost, gen_weight); 
            h_asym->Fill(xF, m, cost, reweight * gen_weight);
            h_asym->Fill(xF, m, -cost, -reweight * gen_weight);
        }
    }
    float norm = h_sym -> Integral();
    h_sym->Scale(1/norm);
    h_asym->Scale(1/norm);
    return 0;
}
void setup(){
    //setup global variables
    if(!rebin){
        //symmetric and anti-symmetric templates
        h_asym = (TH3F *) f_mc->Get("f_asym");
        h_sym = (TH3F *) f_mc->Get("f_sym");
        //DY events with no asymmetry (no signal) (gluglu, qq or qbarqbar)
        h_nosig = (TH3F *) f_mc->Get("f_back");
        //tbar background 
        h_ttbar = (TH3F *)f_ttbar->Get("ttbar_back"); 
        h_sym->GetXaxis();
        h_sym->GetYaxis();
        h_sym->GetZaxis();
    }
    else{
        h_sym = new TH3F("h_sym", "Symmetric template of mc (xF, m, cost_r)",
                n_xf_bins, 0, xf_max, 1,0,m_max, n_cost_bins, -1.,1.);
        h_asym = new TH3F("h_asym", "Asymmetric template of mc (xF, m, cost_r)",
                n_xf_bins, 0, xf_max, 1,0,m_max, n_cost_bins, -1.,1.);
        h_ttbar = new TH3F("h_ttbar", "Template of TTBar background (xF, m, cost_r)",
                n_xf_bins, 0, xf_max, 1,0,m_max, n_cost_bins, -1.,1.);
        h_nosig = new TH3F("h_nosig", "Template of no asymetry drell-yan events from mc (xF, m, cost_r)",
                n_xf_bins, 0, xf_max, 1,0,m_max, n_cost_bins, -1.,1.);
        gen_mc_template(t_mc, h_sym, h_asym);
        gen_template(t_ttbar, h_ttbar);
        gen_template(t_nosig, h_nosig);
    }
    printf("Finishing setup \n");
    return;
}

void MuMu_fit(){
    setup();
    TFile *f_data = TFile::Open("output_files/DYToLL_data_2016_apr4.root");

    //read event data
    TTree *t1 = (TTree *)f_data->Get("T_data"); 
    Long64_t size  =  t1->GetEntries();
    printf("Total sampe size is %i \n", (int) size);

    TH3F* h_data = new TH3F("h_data", "Data template of (x_f, M, c_r)",
                       n_xf_bins, 0, xf_max, 1,0,m_max, n_cost_bins, -1.,1.);
    Double_t m, xF, cost;
    Float_t met_pt;
    t1->SetBranchAddress("m", &m);
    t1->SetBranchAddress("xF", &xF);
    t1->SetBranchAddress("cost", &cost);
    int nBack = 0;
    int nForward = 0;
    int nB=0;
    int nF=0;
    int nEmpties=0;
    for (int i=0; i<size; i++) {
        t1->GetEntry(i);
        if(m >= m_low && m <= m_high){
            v_m.push_back(m);
            v_xF.push_back(xF);
            v_cost.push_back(cost);
            h_data->Fill(xF, m, cost); 
            nDataEvents++;
            Double_t prob = get_prob(xF,m,cost, h_sym);
            if(prob < 1e-20) nEmpties++;
            if(cost <0) nBack++;
            if(cost >0) nForward++;
            if(cost <0 && abs(xF) > 0.1) nB++;
            if(cost >0 && abs(xF) > 0.1) nF++;
        }
    }
    h_data->Scale(1./h_data->Integral());
    float AFB = float(nForward - nBack)/float(nForward + nBack);
    float AFB_highx = float(nF - nB)/float(nF + nB);

    printf("Integrals are %f %f %f %f %f \n", h_data->Integral(), h_sym->Integral(), 
                                           h_asym->Integral(), h_ttbar->Integral(), 
                                           h_nosig->Integral());


    printf("Running fit on %i events (%i hit empty bins).\n" 
            "Counting AFB was %0.3f. High x total was %i,  AFB was %0.3f\n", 
            nDataEvents, nEmpties, AFB, nF+nB, AFB_highx);

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
    int npar = 3;
    int ndim = 3;
    TF3 *fit_fcn = new TF3("F1", &template_fcn, 0., xf_max, 0.,m_max,-1.,1., npar, ndim);
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




    f_data->Close();
    f_ttbar->Close();
    f_mc->Close();

}



