
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
Float_t m_bins[] = {150,200,250,300,400, 500,600, 700, 850, 100000};
//                   150-200  200-250 250-350         350-500 500-700 700+
Double_t alphas[9] = {0.0981, 0.0703, 0.0480, 0.0480, 0.0386, 0.0148, 0.0148, 0.0180, 0.0180};
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

TH2F *h_elel_asym, *h_elel_sym, *h_elel_back,  *h_elel_data, *h_elel_mc;
TH2F *h_elel_mc_count, *h_elel_sym_count;
TH2F *h_mumu_asym, *h_mumu_sym, *h_mumu_back,  *h_mumu_data, *h_mumu_mc;
TH2F *h_mumu_mc_count, *h_mumu_sym_count;
//void *x_axp, *y_axp, *z_axp;

//define axis globally for convinience 
TAxis *x_ax, *y_ax, *z_ax;

//MC templates
TFile* f_elel_mc = (TFile*) TFile::Open("output_files/ElEl_DY_sep25.root");
TTree *t_elel_mc = (TTree *) f_elel_mc ->Get("T_data");
TTree *t_elel_nosig = (TTree *) f_elel_mc ->Get("T_back");
TFile* f_elel_back = (TFile*) TFile::Open("output_files/ElEl_combined_back_sep25.root");
TTree *t_elel_back = (TTree *) f_elel_back ->Get("T_data");

TFile *f_elel_data = TFile::Open("output_files/SingleElectron_data_sep22.root");
TTree *t_elel_data = (TTree *)f_elel_data->Get("T_data"); 

TFile *f_elel_QCD = TFile::Open("../analyze/output_files/ElEl_QCD_est_nov2.root");
TTree *t_elel_QCD = (TTree *)f_elel_QCD->Get("T_data");

TFile *f_elel_WJets = TFile::Open("../analyze/output_files/ElEl_WJets_est_nov2.root");
TTree *t_elel_WJets = (TTree *)f_elel_WJets->Get("T_data");

TFile *f_elel_WJets_contam = TFile::Open("../analyze/FakeRate/root_files/ElEl_fakerate_WJets_MC_dec4.root");
TTree *t_elel_WJets_contam = (TTree *)f_elel_WJets_contam->Get("T_data");

TFile *f_elel_QCD_contam = TFile::Open("../analyze/FakeRate/root_files/ElEl_fakerate_QCD_MC_dec4.root");
TTree *t_elel_QCD_contam = (TTree *)f_elel_QCD_contam->Get("T_data");
////////////////////////////////////////

TFile* f_mumu_mc = (TFile*) TFile::Open("output_files/MuMu_DY_sep8.root");
TTree *t_mumu_mc = (TTree *) f_mumu_mc ->Get("T_data");
TTree *t_mumu_nosig = (TTree *) f_mumu_mc ->Get("T_back");
TFile* f_mumu_back = (TFile*) TFile::Open("output_files/MuMu_combined_back_aug30.root");
TTree *t_mumu_back = (TTree *) f_mumu_back ->Get("T_data");

TFile *f_mumu_data = TFile::Open("output_files/SingleMuon_data_aug28.root");
TTree *t_mumu_data = (TTree *)f_mumu_data->Get("T_data"); 

TFile *f_mumu_QCD = TFile::Open("../analyze/output_files/MuMu_QCD_est_nov2.root");
TTree *t_mumu_QCD = (TTree *)f_mumu_QCD->Get("T_data");

TFile *f_mumu_WJets = TFile::Open("../analyze/output_files/MuMu_WJets_est_Nov2.root");
TTree *t_mumu_WJets = (TTree *)f_mumu_WJets->Get("T_data");

TFile *f_mumu_WJets_contam = TFile::Open("../analyze/FakeRate/root_files/MuMu_fakerate_Wjets_MC_dec4.root");
TTree *t_mumu_WJets_contam = (TTree *)f_mumu_WJets_contam->Get("T_data");

TFile *f_mumu_QCD_contam = TFile::Open("../analyze/FakeRate/root_files/MuMu_fakerate_QCD_MC_dec4.root");
TTree *t_mumu_QCD_contam = (TTree *)f_mumu_QCD_contam->Get("T_data");


vector<double> v_elel_xF;
vector<double> v_elel_m;
vector<double> v_elel_cost;
vector<double> v_mumu_xF;
vector<double> v_mumu_m;
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


// fcn passes back f = - 2*ln(L), (where L is the likelihood of the event)
void fcn(int& npar, double* deriv, double& f, double par[], int flag){
    double lnL = 0.0;
    int elel_misses = 0;
    int mumu_misses = 0;
    if(print) printf("\n \n \n ");

    double AFB = par[0];
    double r_elel_back = par[1];
    double r_mumu_back = par[2];
    for (int i=0; i<nElEl_DataEvents; i++){
        Double_t p_elel_sym = get_prob(v_elel_xF[i], v_elel_cost[i], h_elel_sym);
        Double_t p_elel_asym = get_prob(v_elel_xF[i],  v_elel_cost[i], h_elel_asym);
        Double_t p_elel_back = get_prob(v_elel_xF[i], v_elel_cost[i], h_elel_back);



        double elel_prob = r_elel_back*p_elel_back + (1 - r_elel_back) * (p_elel_sym + AFB*p_elel_asym);
        if(elel_prob > 1) printf("Warning prob is too big \n");
        if(print && p_elel_sym < 1e-20){
            elel_misses++;
            //printf(" Warning p_sym is 0 or negative! for elel bin xf: %0.2f cost: %1.2f \n", v_elel_xF[i], v_elel_cost[i]);
            if(p_elel_back < 1e-20) printf("p_back Is also 0 or negative! \n");
            if(elel_prob < 1e-20) printf("Warning prob is also 0 or negative! \n");
            p_elel_sym = 1e-20;
            elel_prob = r_elel_back*p_elel_back + (1 - r_elel_back) * (p_elel_sym + AFB*p_elel_asym);
        }
        //if(prob < 1e-20) printf(" Warning prob is 0 or negative! for bin xf: %0.2f cost: %1.2f \n", v_xF[i], v_cost[i]);
        elel_prob = max(elel_prob, 1e-20);
        if(elel_prob >0.0) lnL += log(elel_prob);
    }
    for (int i=0; i<nMuMu_DataEvents; i++){
        Double_t p_mumu_sym = get_prob(v_mumu_xF[i], v_mumu_cost[i], h_mumu_sym);
        Double_t p_mumu_asym = get_prob(v_mumu_xF[i],  v_mumu_cost[i], h_mumu_asym);
        Double_t p_mumu_back = get_prob(v_mumu_xF[i], v_mumu_cost[i], h_mumu_back);



        double mumu_prob = r_mumu_back*p_mumu_back + (1 - r_mumu_back) * (p_mumu_sym + AFB*p_mumu_asym);
        if(mumu_prob > 1) printf("Warning prob is too big \n");
        if(print && p_mumu_sym < 1e-20){
            mumu_misses++;
            printf(" Warning p_sym is 0 or negative! for mumu bin xf: %0.2f cost: %1.2f \n", v_mumu_xF[i], v_mumu_cost[i]);
            if(p_mumu_back < 1e-20) printf("p_back Is also 0 or negative! \n");
            if(mumu_prob < 1e-20) printf("Warning prob is also 0 or negative! \n");
            p_mumu_sym = 1e-20;
            mumu_prob = r_mumu_back*p_mumu_back + (1 - r_mumu_back) * (p_mumu_sym + AFB*p_mumu_asym);
        }
        //if(prob < 1e-20) printf(" Warning prob is 0 or negative! for bin xf: %0.2f cost: %1.2f \n", v_xF[i], v_cost[i]);
        mumu_prob = max(mumu_prob, 1e-20);
        if(mumu_prob >0.0) lnL += log(mumu_prob);
    }
    f = -2.0 * lnL;
    if(print) {
        printf("ElEl: %i misses out of %i events \n\n\n", elel_misses, nElEl_DataEvents);
        printf("MuMu: %i misses out of %i events \n\n\n", elel_misses, nMuMu_DataEvents);
        print = false;
    }

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
    TTree *elel_ts[2] = {t_elel_back, t_elel_nosig};

    gen_fakes_template(t_elel_WJets, t_elel_QCD, t_elel_WJets_contam, t_elel_QCD_contam, h_elel_back, m_low, m_high, FLAG_ELECTRONS);
    gen_combined_background_template(2, elel_ts, h_elel_back, m_low, m_high, FLAG_ELECTRONS);

    nElEl_DataEvents = gen_data_template(t_elel_data, h_elel_data, &v_elel_m, &v_elel_xF, &v_elel_cost, m_low, m_high);

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
    TTree *mumu_ts[2] = {t_mumu_back, t_mumu_nosig};
    gen_fakes_template(t_mumu_WJets, t_mumu_QCD, t_mumu_WJets_contam, t_mumu_QCD_contam, h_mumu_back, m_low, m_high, FLAG_MUONS);
    gen_combined_background_template(2, mumu_ts, h_mumu_back, m_low, m_high, FLAG_MUONS);

    nMuMu_DataEvents = gen_data_template(t_mumu_data, h_mumu_data, &v_mumu_m, &v_mumu_xF, &v_mumu_cost, m_low, m_high);
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
    v_elel_m.clear();
    v_elel_cost.clear();
    v_elel_xF.clear();
    v_mumu_m.clear();
    v_mumu_cost.clear();
    v_mumu_xF.clear();
    printf("Finishing cleanup\n");
}
void combined_fit_all(){
    int n_m_bins = 9;
    Double_t AFB_fit[n_m_bins], AFB_fit_err[n_m_bins], r_elel_back_fit[n_m_bins], r_elel_back_fit_err[n_m_bins], 
             r_mumu_back_fit[n_m_bins], r_mumu_back_fit_err[n_m_bins];

    unsigned int nElElEvents[n_m_bins];
    unsigned int nMuMuEvents[n_m_bins];

    Double_t r_mumu_back_starts[] = {0.124, 0.143, 0.16, 0.16, 0.17, 0.17, 0.17, 0.1, 0.1};
    Double_t r_elel_back_starts[] = {0.1, 0.16, 0.2, 0.2, 0.17, 0.16, 0.16, 0.13, 0.13};

    for(int i=0; i<n_m_bins; i++){
        printf("Starting loop \n");
        m_low = m_bins[i];
        m_high = m_bins[i+1];
        alpha = alphas[i];

        setup();

        /*
        printf("ElEL: Integrals are %f %f %f %f  \n", h_elel_data->Integral(), h_elel_sym->Integral(), 
                                               h_elel_asym->Integral(), h_elel_back->Integral() );
        printf("MuMu: Integrals are %f %f %f %f  \n", h_mumu_data->Integral(), h_mumu_sym->Integral(), 
                                               h_mumu_asym->Integral(), h_mumu_back->Integral() );
        h_elel_sym->Print();
        h_mumu_sym->Print();
        */



        float AFB_start = 0.6;
        float AFB_start_error = 0.1;
        float AFB_max = 0.75;
        float r_back_start = 0.12;
        float r_back_start_error = 0.02;
        float r_back_max = 0.6;

        TVirtualFitter * minuit = TVirtualFitter::Fitter(0,3);
        minuit->SetFCN(fcn);
        minuit->SetParameter(0,"AFB", AFB_start, AFB_start_error, -AFB_max, AFB_max);
        minuit->SetParameter(1,"r_elel_back", r_elel_back_starts[i], r_back_start_error, 0, r_back_max);
        minuit->SetParameter(2,"r_mumu_back", r_mumu_back_starts[i], r_back_start_error, 0, r_back_max);
        Double_t arglist[100];
        arglist[0] = 10000.;
        minuit->ExecuteCommand("MIGRAD", arglist,0);
        Double_t up = 1.0;
        minuit->SetErrorDef(up);
        arglist[0] = 0.;
        minuit->ExecuteCommand("MINOS", arglist, 0);




        AFB_fit[i] = minuit->GetParameter(0); 
        AFB_fit_err[i] = minuit->GetParError(0);
        r_elel_back_fit[i] = minuit->GetParameter(1); 
        r_elel_back_fit_err[i] = minuit->GetParError(1);
        r_mumu_back_fit[i] = minuit->GetParameter(2); 
        r_mumu_back_fit_err[i] = minuit->GetParError(2);

        nElElEvents[i] = nElEl_DataEvents;
        nMuMuEvents[i] = nMuMu_DataEvents;



        cleanup();
    }
    for(int i=0; i<n_m_bins; i++){
        printf("\n Fit on M=[%.0f, %.0f], %i ElEl Events, %i MuMu Events: AFB = %0.3f +/- %0.3f r_elel_back = %0.3f +/- %0.3f r_mumu_back =%0.3f +/- %0.3f \n", 
                    m_bins[i], m_bins[i+1], nElElEvents[i], nMuMuEvents[i], AFB_fit[i], AFB_fit_err[i], r_elel_back_fit[i], r_elel_back_fit_err[i], r_mumu_back_fit[i], r_mumu_back_fit_err[i]);

    }

}



