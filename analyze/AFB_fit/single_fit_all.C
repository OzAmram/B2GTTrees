
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



int FLAG = FLAG_ELECTRONS;
//int FLAG = FLAG_MUONS;

//MC templates
// /*
TFile* f_mc = (TFile*) TFile::Open("output_files/ElEl_DY_nov25.root");
TTree *t_mc = (TTree *) f_mc ->Get("T_data");
TTree *t_nosig = (TTree *) f_mc ->Get("T_back");
TFile* f_back = (TFile*) TFile::Open("output_files/ElEl_combined_back_sep25.root");
TTree *t_back = (TTree *) f_back ->Get("T_data");

TFile *f_data = TFile::Open("output_files/SingleElectron_data_sep22.root");
TTree *t_data = (TTree *)f_data->Get("T_data"); 

TFile *f_QCD = TFile::Open("../analyze/output_files/ElEl_QCD_est_nov2.root");
TTree *t_QCD = (TTree *)f_QCD->Get("T_data");

TFile *f_WJets = TFile::Open("../analyze/output_files/ElEl_WJets_est_nov2.root");
TTree *t_WJets = (TTree *)f_WJets->Get("T_data");

TFile *f_WJets_contam = TFile::Open("../analyze/output_files/ElEl_fakerate_WJets_MC_nov2.root");
TTree *t_WJets_contam = (TTree *)f_WJets_contam->Get("T_data");
// */
 /*
TFile* f_mc = (TFile*) TFile::Open("output_files/MuMu_DY_sep8.root");
TTree *t_mc = (TTree *) f_mc ->Get("T_data");
TTree *t_nosig = (TTree *) f_mc ->Get("T_back");
TFile* f_back = (TFile*) TFile::Open("output_files/MuMu_combined_back_aug30.root");
TTree *t_back = (TTree *) f_back ->Get("T_data");

TFile *f_data = TFile::Open("output_files/SingleMuon_data_aug28.root");
TTree *t_data = (TTree *)f_data->Get("T_data"); 

TFile *f_QCD = TFile::Open("../analyze/output_files/MuMu_QCD_est_nov2.root");
TTree *t_QCD = (TTree *)f_QCD->Get("T_data");

TFile *f_WJets = TFile::Open("../analyze/output_files/MuMu_WJets_est_Nov2.root");
TTree *t_WJets = (TTree *)f_WJets->Get("T_data");

TFile *f_WJets_contam = TFile::Open("../analyze/output_files/MuMu_fakerate_Wjets_MC_nov2.root");
TTree *t_WJets_contam = (TTree *)f_WJets_contam->Get("T_data");
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
void fcn(int& npar, double* deriv, double& f, double par[], int flag){
    double lnL = 0.0;
    int misses = 0;
    if(print) printf("\n \n \n ");

    for (int i=0; i<nDataEvents; i++){
        Double_t p_sym = get_prob(v_xF[i], v_cost[i], h_sym);
        Double_t p_asym = get_prob(v_xF[i],  v_cost[i], h_asym);
        Double_t p_back = get_prob(v_xF[i], v_cost[i], h_back);



        double AFB = par[0];
        double r_back = par[1];
        double prob = r_back*p_back + (1 - r_back) * (p_sym + AFB*p_asym);
        if(prob > 1) printf("Warning prob is too big \n");
        if(print && p_sym < 1e-20){
            misses++;
            printf(" Warning p_sym is 0 or negative! for bin xf: %0.2f cost: %1.2f \n", v_xF[i], v_cost[i]);
            if(p_back < 1e-20) printf("p_back Is also 0 or negative! \n");
            if(prob < 1e-20) printf("Warning prob is also 0 or negative! \n");
            p_sym = 1e-20;
            prob = r_back*p_back + (1 - r_back) * (p_sym + AFB*p_asym);
        }
        //if(prob < 1e-20) printf(" Warning prob is 0 or negative! for bin xf: %0.2f cost: %1.2f \n", v_xF[i], v_cost[i]);
        prob = max(prob, 1e-20);
        if(prob >0.0) lnL += log(prob);
    }
    f = -2.0 * lnL;
    if(print) {
        printf("%i misses out of %i events \n\n\n", misses, nDataEvents);
        print = false;
    }

}

Double_t  template_fcn(Double_t *x, Double_t *par){
    Double_t xx = x[0];
    Double_t yy = x[1];

    Double_t p_sym = get_prob(xx, yy,  h_sym);
    Double_t p_asym = get_prob(xx, yy, h_asym);
    Double_t p_back = get_prob(xx, yy,  h_back);

    Double_t AFB = par[0];
    Double_t r_back = par[1];
    Double_t prob = r_back*p_back + (1 - r_back) * (p_sym + AFB*p_asym);
    return prob;
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

    gen_mc_template(t_mc, alpha, h_sym, h_asym, h_sym_count, m_low, m_high, FLAG);
    TTree *ts[2] = {t_back, t_nosig};
    gen_combined_background_template(2, ts, h_back, t_WJets, t_QCD, t_WJets_contam, m_low, m_high, FLAG);

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

        setup();

        printf("Integrals are %f %f %f %f  \n", h_data->Integral(), h_sym->Integral(), 
                                               h_asym->Integral(), h_back->Integral() );
        h_sym->Print();



        float AFB_start = 0.6;
        float AFB_start_error = 0.1;
        float AFB_max = 0.75;
        float r_back_start = 0.12;
        float r_back_start_error = 0.04;
        float r_back_max = 0.6;

        TVirtualFitter * minuit = TVirtualFitter::Fitter(0,2);
        minuit->SetFCN(fcn);
        minuit->SetParameter(0,"AFB", AFB_start, AFB_start_error, -AFB_max, AFB_max);
        minuit->SetParameter(1,"r_back", r_back_start, r_back_start_error, 0, r_back_max);
        Double_t arglist[100];
        arglist[0] = 10000.;
        minuit->ExecuteCommand("MIGRAD", arglist,0);
        Double_t up = 1.0;
        minuit->SetErrorDef(up);
        arglist[0] = 0.;
        minuit->ExecuteCommand("MINOS", arglist, 0);


        /*
        printf("Trying template fit \n\n\n");
        Int_t npar = 2;
        Int_t ndim = 2;
        TF2 *fit_fcn = new TF2("F1", &template_fcn, 0., 1., -1.,1., npar, ndim);
        fit_fcn->SetParNames("AFB", "r_back");
        fit_fcn->SetParameter(0, AFB_start); //AFB
        fit_fcn->SetParLimits(0, -AFB_max, AFB_max);
        fit_fcn->SetParError(0, AFB_start_error);
        fit_fcn->SetParameter(1, r_back_start); //r_back
        fit_fcn->SetParLimits(1, 0, r_back_max);
        fit_fcn->SetParError(1, r_back_start_error);
        h_data->Sumw2();
        h_data->Fit(fit_fcn, "WL M N");
        */


        AFB_fit[i] = minuit->GetParameter(0); 
        AFB_fit_err[i] = minuit->GetParError(0);
        r_back_fit[i] = minuit->GetParameter(1); 
        r_back_fit_err[i] = minuit->GetParError(1);

        nEvents[i] = nDataEvents;



        cleanup();
    }
    for(int i=0; i<n_m_bins; i++){
        printf("\n\n\n Fit on M=[%.0f, %.0f], %i Events: AFB = %0.3f +/- %0.3f r_back = %0.3f +/- %0.3f \n", 
                    m_bins[i], m_bins[i+1], nEvents[i], AFB_fit[i], AFB_fit_err[i], r_back_fit[i], r_back_fit_err[i]);

    }

}



