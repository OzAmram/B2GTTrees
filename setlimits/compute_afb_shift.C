
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
#include "../analyze/TemplateMaker.C"
#include "../analyze/AFB_fit/root_files.h""



float m_low;
float m_high;

bool print = true;


Double_t med_btag = 0.4432;

TH2F *h_elel_asym, *h_elel_sym, *h_elel_back,  *h_elel_data, *h_elel_mc;
TH2F *h_elel_mc_count, *h_elel_sym_count;
TH2F *h_mumu_asym, *h_mumu_sym, *h_mumu_back,  *h_mumu_data, *h_mumu_mc;
TH2F *h_mumu_mc_count, *h_mumu_sym_count;
//void *x_axp, *y_axp, *z_axp;

//define axis globally for convinience 
TAxis *x_ax, *y_ax, *z_ax;



vector<double> v_elel_xF;
vector<double> v_elel_cost;
vector<double> v_mumu_xF;
vector<double> v_mumu_cost;
unsigned int nElEl_DataEvents;
unsigned int nMuMu_DataEvents;

Double_t get_prob(Double_t xF, Double_t cost, TH2F *h){
    TAxis* x_ax =  h->GetXaxis();
    TAxis *y_ax =  h->GetYaxis();
    //binning the same in all templates
    int xbin = x_ax->FindBin(xF);
    int ybin = y_ax->FindBin(cost);

    return h->GetBinContent(xbin, ybin);
}



void setup(){
    //setup global variables
    //TH1::SetDefaultSumw2(kTRUE);
    h_elel_mc_count = new TH2F("h_elel_mc_count", "Events in bins for MC templates",
            n_xf_bins, xf_bins, n_cost_bins, cost_bins);
    h_elel_mc_count->SetDirectory(0);
    h_elel_sym_count = new TH2F("h_elel_sym_count", "Events in bins for MC templates",
            n_xf_bins, xf_bins, n_cost_bins, cost_bins);
    h_elel_sym_count->SetDirectory(0);
    h_elel_sym = new TH2F("h_elel_sym", "Symmetric template of mc",
            n_xf_bins, xf_bins, n_cost_bins, cost_bins);
    h_elel_sym->SetDirectory(0);
    h_elel_asym = new TH2F("h_elel_asym", "Asymmetric template of mc",
            n_xf_bins, xf_bins, n_cost_bins, cost_bins);
    h_elel_asym->SetDirectory(0);
    h_elel_back = new TH2F("h_elel_back", "Combined background template",
            n_xf_bins, xf_bins, n_cost_bins, cost_bins);
    h_elel_back->SetDirectory(0);
    h_elel_data = new TH2F("h_elel_data", "Data template of (x_f, cost_r)",
            n_xf_bins, xf_bins, n_cost_bins, cost_bins);
    h_elel_data->SetDirectory(0);



    gen_mc_template(t_elel_mc, alpha, h_elel_sym, h_elel_asym, h_elel_sym_count, m_low, m_high, FLAG_ELECTRONS);

    h_mumu_mc_count = new TH2F("h_mumu_mc_count", "Events in bins for MC templates",
            n_xf_bins, xf_bins, n_cost_bins, cost_bins);
    h_mumu_mc_count->SetDirectory(0);
    h_mumu_sym_count = new TH2F("h_mumu_sym_count", "Events in bins for MC templates",
            n_xf_bins, xf_bins, n_cost_bins, cost_bins);
    h_mumu_sym_count->SetDirectory(0);
    h_mumu_sym = new TH2F("h_mumu_sym", "Symmetric template of mc",
            n_xf_bins, xf_bins, n_cost_bins, cost_bins);
    h_mumu_sym->SetDirectory(0);
    h_mumu_asym = new TH2F("h_mumu_asym", "Asymmetric template of mc",
            n_xf_bins, xf_bins, n_cost_bins, cost_bins);
    h_mumu_asym->SetDirectory(0);
    h_mumu_back = new TH2F("h_mumu_back", "Combined background template",
            n_xf_bins, xf_bins, n_cost_bins, cost_bins);
    h_mumu_back->SetDirectory(0);

    h_mumu_data = new TH2F("h_mumu_data", "Data template of (x_f, cost_r)",
            n_xf_bins, xf_bins, n_cost_bins, cost_bins);
    h_mumu_data->SetDirectory(0);


    gen_mc_template(t_mumu_mc, alpha, h_mumu_sym, h_mumu_asym, h_mumu_sym_count, m_low, m_high, FLAG_MUONS);

    printf("Finishing setup \n");
    return;
}


void compute_afb_shift(){
    Double_t AFB = 0.6;
    init();
    f_m100_dilu = (TFile*) TFile::Open("dilus/M150_dilus.root");
    TH1F *utype_dilu = (TH1F *) f_m100_dilu ->Get("utype_dilu");
    TH1F *dtype_dilu = (TH1F *) f_m100_dilu ->Get("dtype_dilu");
    TH1F *mix_dilu = (TH1F *) f_m100_dilu ->Get("mix_dilu");

    TH1F *h_dilu_ratio = (TH1F *)utype_dilu->Clone("h_dilu_ratio");
    h_dilu_ratio->Divide(mix_dilu);
    for(int i=0; i<1; i++){
        printf("Starting loop \n");
        m_low = m_bins[i];
        m_high = m_bins[i+1];
        alpha = alphas[i];

        setup();

        TH2F *f = (TH2F *) h_mumu_sym->Clone("f");
        TH2F *g = (TH2F *) h_mumu_sym->Clone("f");
        f->Add(h_mumu_asym, AFB);
        TH2F *h_asym_dilu = (TH2F *) h_mumu_asym->Clone("h_mumu_asym_dilu");
        for(int i=0; i<=n_xf_bins; i++){
            for(int j=0; j<n_cost_bins; j++){
                Double_t p = h_asym_dilu->GetBinContent(i,j);
                Double_t dilu_ratio = h_dilu_ratio->GetBinContent(i);
                h_asym_dilu->SetBinContent(i,j, p*dilu_ratio);
            }
        }
        g->Add(h_asym_dilu, AFB);
        TH2F* f_sqrd = (TH2F *) f->Clone("f2");
        f_sqrd->Multiply(f);

        TH2F *M = (TH2F *) g->Clone("M");
        M->Divide(f_sqrd);
        M->Multiply(h_mumu_asym);
        M->Multiply(h_mumu_asym);
        Double_t M_int = M->Integral();

        TH2F *A = (TH2F *) g->Clone("A");
        A->Multiply(h_mumu_asym);
        A->Divide(h_mumu_asym);
        Double_t delta_AFB = (1/M_int) * A->Integral();
        printf("For M_low %.0f, delta AFB %.2f \n", m_low, delta_AFB);
    }
}



