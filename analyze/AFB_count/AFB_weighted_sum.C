
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
#include "../TemplateMaker.C"



int n_xf_bins = 5;
float xf_max = 1.0;
Float_t xf_bins[] = {0., 0.05, 0.1, 0.15, 0.25, 1.0};
int n_cost_bins = 10;
Float_t cost_bins[] = {-1.0, -.8, -.6, -.4, -.2, 0., 0.2, 0.4, 0.6, 0.8, 1.0};
Float_t m_bins[] = {150,200,250,350,500,700,100000};
Double_t alphas[6] = {0.0981, 0.0703, 0.0480, 0.0386, 0.0148, 0.0180};
Double_t alpha;

int FLAG = FLAG_ELECTRONS;

//int n_xf_bins = 4;
//float xf_max = 1.0;
//Float_t xf_bins[] = {0., 0.04, 0.08, 0.14, 1.0};
//float m_max = 1000.;
//int n_cost_bins = 8;
//Float_t cost_bins[] = {-1.0, -.75, -.5, -.25, 0., 0.25, 0.5, 0.75, 1.0};


float m_low;
float m_high;
//alpha = 0.0981;

bool print = true;


Double_t med_btag = 0.4432;

TH2F *h_asym, *h_sym, *h_back,  *h_data, *h_mc;
TH2F *h_mc_count, *h_sym_count;
//void *x_axp, *y_axp, *z_axp;

//define axis globally for convinience 
TAxis *x_ax, *y_ax, *z_ax;

//MC templates
TFile* f_mc = (TFile*) TFile::Open("output_files/ElEl_DY_sep25.root");
TTree *t_mc = (TTree *) f_mc ->Get("T_data");
TTree *t_nosig = (TTree *) f_mc ->Get("T_back");
TFile* f_back = (TFile*) TFile::Open("output_files/ElEl_combined_back_sep25.root");
TTree *t_back = (TTree *) f_back ->Get("T_data");

TFile *f_data = TFile::Open("output_files/SingleElectron_data_sep22.root");
TTree *t_data = (TTree *)f_data->Get("T_data"); 
/*
TFile* f_mc = (TFile*) TFile::Open("output_files/MuMu_DY_sep8.root");
TTree *t_mc = (TTree *) f_mc ->Get("T_data");
TTree *t_nosig = (TTree *) f_mc ->Get("T_back");
TFile* f_back = (TFile*) TFile::Open("output_files/MuMu_combined_back_aug30.root");
TTree *t_back = (TTree *) f_back ->Get("T_data");

TFile *f_data = TFile::Open("output_files/SingleMuon_data_aug28.root");
TTree *t_data = (TTree *)f_data->Get("T_data"); 
*/


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
void count_AFB(double *AFB, double *AFB_unc){
    Double_t D_F=0;
    Double_t D_B=0;
    Double_t N_F=0;
    Double_t N_B=0;

    Double_t D_FF = 0;
    Double_t N_FF = 0;
    Double_t C_FF = 0;

    Double_t D_BB = 0;
    Double_t N_BB = 0;
    Double_t C_BB = 0;

    for (int i =0; i nDataEvents; i++){
        Double_t cost = v_cost[i];
        Double_t cost2 = cost*cost;
        Double_t h = alpha * (1 - cost2);
        Double_t w_D = 0.5 *cost2 / (pow((1 + cost2 + h),3));
        Double_t w_N = 0.5 *abs(cost) / (pow((1 + cost2 + h),2));

        if(cost >0){
            D_F += w_D;
            N_F += w_N;
            D_FF += w_D*w_D;
            N_FF += w_N*w_N;
            C_FF += w_N*w_D;
        }
        else{
            D_B += w_D;
            N_B += w_N;
            D_BB += w_D*w_D;
            N_BB += w_N*w_N;
            C_BB += w_N*w_D;
        }
    }
    AFB = (3.0/8.0) * (N_F - N_B)/ (D_F + D_B);
    Double_t term1 = (9.0/64.0)*(N_FF + N_BB);
    Double_t term2 = (AFB*AFB*(D_FF + D_BB));
    Double_t term3 = (-3.0/4.0)*(C_FF - C_BB);
    AFB_unc = (1./(D_F + D_B)) * sqrt(term1 + term2 + term3);
}





void setup(){
    //setup global variables
    //TH1::SetDefaultSumw2(kTRUE);
    printf("Starting setup \n");
    h_mc_count = new TH2F("h_mc_count", "Events in bins for MC templates",
            n_xf_bins, xf_bins, n_cost_bins, cost_bins);
    h_mc_count->SetDirectory(0);
    h_sym_count = new TH2F("h_sym_count", "Events in bins for MC templates",
            n_xf_bins, xf_bins, n_cost_bins, cost_bins);
    h_sym_count->SetDirectory(0);
    h_sym = new TH2F("h_sym", "Symmetric template of mc (xF, cost_r) xF>0.15",
            n_xf_bins, xf_bins, n_cost_bins, cost_bins);
    h_sym->SetDirectory(0);
    h_asym = new TH2F("h_asym", "Asymmetric template of mc (xF cost_r) xF > 0.15",
            n_xf_bins, xf_bins, n_cost_bins, cost_bins);
    h_asym->SetDirectory(0);
    h_back = new TH2F("h_back", "Combined background template",
            n_xf_bins, xf_bins, n_cost_bins, cost_bins);

    h_back->SetDirectory(0);
    h_data = new TH2F("h_data", "Data template of (x_f, cost_r) xF > 0.15",
            n_xf_bins, xf_bins, n_cost_bins, cost_bins);

    h_data->SetDirectory(0);
    printf("Generating templates \n");

    //gen_mc_template(t_mc, alpha, h_sym, h_asym, h_sym_count, m_low, m_high, FLAG);
    TTree *ts[2] = {t_back, t_nosig};
    //gen_combined_background_template(2, ts, h_back, m_low, m_high, FLAG);

    nDataEvents = gen_data_template(t_data, h_data, &v_m, &v_xF, &v_cost, m_low, m_high);
    //f_data->Close();
    //f_back->Close();
    //f_mc->Close();
    printf("Finishing setup \n");
    return;
}


void cleanup(){
    //delete h_mc_count;
    //delete h_sym_count;
    //delete h_sym;
    //delete h_asym;
    //delete h_back;
    //delete h_data;
    v_m.clear();
    v_cost.clear();
    v_xF.clear();
    printf("Finishing cleanup\n");
}

void single_fit_all(){
    int n_m_bins = 6;
    Double_t AFB_fit[n_m_bins], AFB_fit_err[n_m_bins], r_back_fit[n_m_bins], r_back_fit_err[n_m_bins];
    unsigned int nEvents[n_m_bins];

    for(int i=0; i<n_m_bins; i++){
        printf("Starting loop \n");
        m_low = m_bins[i];
        m_high = m_bins[i+1];
        alpha = alphas[i];
        Double_t AFB, AFB_unc;
        count_AFB(&AFB, &AFB_unc);


        setup();


        nEvents[i] = nDataEvents;

        AFB_fit[i] = AFB;
        AFB_fit_err[i] = AFB_unc;



        cleanup();
    }
    for(int i=0; i<n_m_bins; i++){
        printf("\n\n\n Fit on M=[%.0f, %.0f], %i Events: AFB = %0.3f +/- %0.3f \n", 
                    m_bins[i], m_bins[i+1], nEvents[i], AFB_fit[i], AFB_fit_err[i]);

    }

}



