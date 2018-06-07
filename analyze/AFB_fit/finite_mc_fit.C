

//perform fits to Reconstructed MuMu data to extract Asym

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <iostream>
#include <cstring>
#include <algorithm>
#include <vector>

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
#include "../TemplateMaker_systematics.C"
//#include "../TemplateMaker.C"
#include "root_files.h"




const TString fout_name("AFB_fit/fit_results/m_bins/MuMu_fit_finite_mc_stat_june07.root");



float m_low;
float m_high;

bool print = true;


Double_t med_btag = 0.4432;

TH2F *h_elel_asym, *h_elel_sym, *h_elel_back,  *h_elel_data, *h_elel_dilu, *h_elel_fitted_mc;
TH2F *h_elel_mc_count, *h_elel_sym_count, *h_elel_mc, *h_elel_mc_errs;
TH2F *h_mumu_asym, *h_mumu_sym, *h_mumu_back,  *h_mumu_data, *h_mumu_dilu, *h_mumu_fitted_mc;
TH2F *h_mumu_mc_count, *h_mumu_sym_count, *h_mumu_mc, *h_mumu_mc_errs;
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

void print_hist(TH2F *h){
    printf("\n");
    for(int i=1; i<= n_xf_bins; i++){
        for(int j=1; j<= n_cost_bins; j++){
            printf("%.2e ",    h->GetBinContent(i,j));
        }
        printf("\n");
    }

}

void gen_fit_templates_from_mc_and_dilu(TH2F *h_sym, TH2F *h_asym, TH2F *h_mc, TH2F* h_dilu){
    for(int i=1; i<=n_xf_bins; i++){    
        for(int j=1; j<=n_cost_bins; j++){
            Float_t bin_weight = h_mc->GetBinContent(i,j);
            Double_t xf_center = h_mc->GetXaxis()->GetBinCenter(i);
            Double_t cost_center = h_mc->GetYaxis()->GetBinCenter(j);
            Double_t dilu_factor = h_dilu->GetBinContent(i,j);

                
            Double_t reweight = (4./3.)*cost_center*(2. + alpha)/
                (1. + cost_center*cost_center + alpha*(1.- cost_center*cost_center));


            //printf("cost,  xf ,rw rw_avg: %.2f,, %.2f %.2f %.2f \n", cost_center, xf_center,reweight, reweight_avg);
           
            h_sym->Fill(xf_center, cost_center, bin_weight);
            h_sym->Fill(xf_center, -cost_center, bin_weight);

            h_asym->Fill(xf_center, cost_center, reweight *bin_weight* dilu_factor);
            h_asym->Fill(xf_center, -cost_center, -reweight *bin_weight * dilu_factor);

            h_asym->Fill(xf_center, -cost_center, reweight *bin_weight * (1-dilu_factor));
            h_asym->Fill(xf_center, cost_center, -reweight *bin_weight *(1-dilu_factor));
        }
    }
    Double_t h_int = h_sym->Integral();
    h_sym->Scale(1./h_int);
    h_asym->Scale(1./h_int);
    return;
}
void gen_fit_templates_from_finite_mc(TH2F *h_sym, TH2F *h_asym, TH2F *h_mc_corr, TH2F* h_mc_inc){
    for(int i=1; i<=n_xf_bins; i++){    
        for(int j=1; j<=n_cost_bins; j++){
            Float_t bin_weight_corr = h_mc_corr->GetBinContent(i,j);
            Float_t bin_weight_inc = h_mc_inc->GetBinContent(i,j);
            Double_t xf_center = h_mc_inc->GetXaxis()->GetBinCenter(i);
            Double_t cost_center = h_mc_inc->GetYaxis()->GetBinCenter(j);

                
            Double_t reweight = (4./3.)*cost_center*(2. + alpha)/
                (1. + cost_center*cost_center + alpha*(1.- cost_center*cost_center));


            //printf("cost,  xf ,rw rw_avg: %.2f,, %.2f %.2f %.2f \n", cost_center, xf_center,reweight, reweight_avg);
           
            h_sym->Fill(xf_center, cost_center, bin_weight_corr + bin_weight_inc);
            h_sym->Fill(xf_center, -cost_center, bin_weight_corr + bin_weight_inc);

            h_asym->Fill(xf_center, cost_center, reweight *bin_weight_corr);
            h_asym->Fill(xf_center, -cost_center, -reweight *bin_weight_corr);

            h_asym->Fill(xf_center, -cost_center, reweight *bin_weight_inc);
            h_asym->Fill(xf_center, cost_center, -reweight *bin_weight_inc);
        }
    }
    Double_t h_int = h_sym->Integral();
    h_sym->Scale(1./h_int);
    h_asym->Scale(1./h_int);
    return;
}
void setHistParams(TVirtualFitter *minuit, TH2F *h_mc, int start_idx){
    char base_str[20] = "hist_bin_%i_%i";
    char param_str[20];
    print_hist(h_mc);
    for(int i=0; i< n_xf_bins; i++){
        for(int j=0; j< n_cost_bins; j++){
            int param_idx = start_idx + i*n_cost_bins + j;
            sprintf(param_str, base_str, i, j);
            Float_t start = h_mc->GetBinContent(i+1, j+1);
            Float_t unc = h_mc->GetBinError(i+1, j+1);
            //if(param_idx > 40) printf("Start, unc = %.2e %.2e \n", start, unc);
            //printf("Unc ratio is %.2f \n", unc/start);
            minuit->SetParameter(param_idx,param_str, start, unc, 0.,1e4);
        }
    }
    return;
}
void update_templates(double *par, TH2F *h_mc, TH2F *h_sym, TH2F *h_asym, TH2F *h_dilu, int start_idx){
    for(int i=0; i< n_xf_bins; i++){
        for(int j=0; j< n_cost_bins; j++){
            int param_idx = start_idx + i*n_cost_bins + j;
            h_mc->SetBinContent(i+1, j+1, par[param_idx]);
        }
    }

    gen_fit_templates_from_mc_and_dilu(h_sym, h_asym, h_mc, h_dilu);
}

double get_hists_likelihood(TH2F *h_mc, TH2F *h_fitted_mc){
    double sum = 0;
    for(int i=0; i< n_xf_bins; i++){
        for(int j=0; j< n_cost_bins; j++){
            double obs = h_mc->GetBinContent(i+1, j+1);
            double unc = h_mc->GetBinError(i+1, j+1);
            double mean = h_fitted_mc->GetBinContent(i+1, j+1);
            sum += pow((obs-mean)/unc,2);
        }
    }
    return sum;
}


// fcn passes back f = - 2*ln(L), (where L is the likelihood of the event)
void fcn(int& npar, double* deriv, double& f, double par[], int flag){
    f=0;
    double lnL = 0.0;
    double lnL_hists = 0.0;
    int misses = 0;
    if(print) printf("\n \n \n ");

    double AFB = par[0];
    double r_back = par[1];


    /*
    update_templates(par, h_elel_fitted_mc, h_elel_sym, h_elel_asym, h_elel_dilu, 2);
    //print_hist(h_elel_fitted_mc);
    for (int i=0; i<nElEl_DataEvents; i++){
        Double_t p_sym = get_prob(v_elel_xF[i], v_elel_cost[i], h_elel_sym);
        Double_t p_asym = get_prob(v_elel_xF[i],  v_elel_cost[i], h_elel_asym);
        Double_t p_back = get_prob(v_elel_xF[i], v_elel_cost[i], h_elel_back);



        if(p_back < 0){ 
            //Printf("P_back < 0! %.2e %.2f %.2f \n", p_back, v_xF[i], v_cost[i]);
            p_back = 1e-20;
        }

        double prob = r_back*p_back + (1 - r_back) * (p_sym + AFB*p_asym);
        if(prob > 1) printf("Warning prob is too big %.2f %.2f %.2f %.2f %.2f  \n", r_back, AFB, p_back, p_sym, p_asym);
        if(print && p_sym < 1e-20){
            misses++;
            printf(" Warning p_sym is 0 or negative! for bin xf: %0.2f cost: %1.2f \n", v_elel_xF[i], v_elel_cost[i]);
            if(p_back < 1e-20) printf("p_back Is also 0 or negative! \n");
            if(prob < 1e-20) printf("Warning prob is also 0 or negative! \n");
            p_sym = 1e-20;
            prob = r_back*p_back + (1 - r_back) * (p_sym + AFB*p_asym);
        }
        //if(prob < 1e-20) printf(" Warning prob is 0 or negative! for bin xf: %0.2f cost: %1.2f \n", v_xF[i], v_cost[i]);
        prob = max(prob, 1e-20);
        if(prob >0.0) lnL += log(prob);
    }
    lnL_hists = get_hists_likelihood(h_elel_mc, h_elel_fitted_mc);
    //printf("lnL = %.3f  hists = %.3f \n", -2.0*lnL, lnL_hists);
    f += -2.0 * lnL + lnL_hists;
    if(print) {
        printf("%i misses out of %i events \n\n\n", misses, nElEl_DataEvents);
        print = false;
    }
    */

    update_templates(par, h_mumu_fitted_mc, h_mumu_sym, h_mumu_asym, h_mumu_dilu, 2);
    //print_hist(h_mumu_fitted_mc);
    for (int i=0; i<nMuMu_DataEvents; i++){
        Double_t p_sym = get_prob(v_mumu_xF[i], v_mumu_cost[i], h_mumu_sym);
        Double_t p_asym = get_prob(v_mumu_xF[i],  v_mumu_cost[i], h_mumu_asym);
        Double_t p_back = get_prob(v_mumu_xF[i], v_mumu_cost[i], h_mumu_back);



        if(p_back < 0){ 
            //Printf("P_back < 0! %.2e %.2f %.2f \n", p_back, v_xF[i], v_cost[i]);
            p_back = 1e-20;
        }

        double prob = r_back*p_back + (1 - r_back) * (p_sym + AFB*p_asym);
        if(prob > 1) printf("Warning prob is too big %.2f %.2f %.2f %.2f %.2f  \n", r_back, AFB, p_back, p_sym, p_asym);
        if(print && p_sym < 1e-20){
            misses++;
            printf(" Warning p_sym is 0 or negative! for bin xf: %0.2f cost: %1.2f \n", v_mumu_xF[i], v_mumu_cost[i]);
            if(p_back < 1e-20) printf("p_back Is also 0 or negative! \n");
            if(prob < 1e-20) printf("Warning prob is also 0 or negative! \n");
            p_sym = 1e-20;
            prob = r_back*p_back + (1 - r_back) * (p_sym + AFB*p_asym);
        }
        //if(prob < 1e-20) printf(" Warning prob is 0 or negative! for bin xf: %0.2f cost: %1.2f \n", v_xF[i], v_cost[i]);
        prob = max(prob, 1e-20);
        if(prob >0.0) lnL += log(prob);
    }
    lnL_hists = get_hists_likelihood(h_mumu_mc, h_mumu_fitted_mc);
    //printf("lnL = %.3f  hists = %.3f \n", -2.0*lnL, lnL_hists);
    f += -2.0 * lnL + lnL_hists;
    if(print) {
        printf("%i misses out of %i events \n\n\n", misses, nMuMu_DataEvents);
        print = false;
    }

}


void get_dilu_hist(TH2F *h_dilu, TH2F *h_mc_corr, TH2F* h_mc_inc){
    for(int i=1; i<=n_xf_bins; i++){    
        for(int j=1; j<=n_cost_bins; j++){
            Float_t bin_weight_corr = h_mc_corr->GetBinContent(i,j);
            Float_t bin_weight_inc = h_mc_inc->GetBinContent(i,j);
            Double_t xf_center = h_mc_inc->GetXaxis()->GetBinCenter(i);
            Double_t cost_center = h_mc_inc->GetYaxis()->GetBinCenter(j);
            h_dilu->Fill(xf_center, cost_center, bin_weight_corr/(bin_weight_corr + bin_weight_inc));
        }
    }
    return;
}





void setup(){
    //setup global variables
    //TH1::SetDefaultSumw2(kTRUE);
    /*
    h_elel_mc_count = new TH2F("h_elel_mc_count", "Events in bins for MC templates",
            n_xf_bins, xf_bins, n_cost_bins, cost_bins);
    h_elel_mc_count->SetDirectory(0);
    h_elel_mc = new TH2F("h_elel_mc", "Events in bins for MC templates",
            n_xf_bins, xf_bins, n_cost_bins, cost_bins);
    h_elel_mc->SetDirectory(0);
    h_elel_mc_errs = new TH2F("h_elel_mc_errs", "Events in bins for MC templates",
            n_xf_bins, xf_bins, n_cost_bins, cost_bins);
    h_elel_mc_errs->SetDirectory(0);
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

    TH2F *h_elel_mc_corr = (TH2F *)h_elel_mc->Clone("h_elel_mc_corr");
    TH2F *h_elel_mc_inc = (TH2F *) h_elel_mc->Clone("h_elel_mc_inc");
    h_elel_dilu = (TH2F *) h_elel_mc->Clone("h_elel_dilu");

    //gen_mc_template(t_elel_mc, alpha, h_elel_sym, h_elel_asym, h_elel_sym_count, m_low, m_high, FLAG_ELECTRONS);
    gen_finite_mc_template(t_elel_mc, h_elel_mc_corr, h_elel_mc_inc, h_elel_mc_errs, h_elel_mc_count, m_low, m_high, FLAG_ELECTRONS);
    h_elel_mc->Add(h_elel_mc_inc, h_elel_mc_corr);
    h_elel_fitted_mc = (TH2F *) h_elel_mc->Clone("h_elel_fitted_mc");
    get_dilu_hist(h_elel_dilu, h_elel_mc_corr, h_elel_mc_inc);
    gen_fit_templates_from_mc_and_dilu(h_elel_sym, h_elel_asym, h_elel_mc, h_elel_dilu);
    //gen_fit_templates_from_finite_mc(h_elel_sym, h_elel_asym, h_elel_mc_corr, h_elel_mc_inc);
    print_hist(h_elel_mc);
    print_hist(h_elel_mc_errs);
    print_hist(h_elel_mc_count);
    TTree *elel_ts[2] = {t_elel_back, t_elel_nosig};

    gen_fakes_template(t_elel_WJets, t_elel_QCD, t_elel_WJets_contam, t_elel_QCD_contam, h_elel_back, m_low, m_high, FLAG_ELECTRONS);
    gen_combined_background_template(2, elel_ts, h_elel_back, m_low, m_high, FLAG_ELECTRONS);

    nElEl_DataEvents = gen_data_template(t_elel_data, h_elel_data,  &v_elel_xF, &v_elel_cost, m_low, m_high, FLAG_ELECTRONS);
    */

    h_mumu_mc_count = new TH2F("h_mumu_mc_count", "Events in bins for MC templates",
            n_xf_bins, xf_bins, n_cost_bins, cost_bins);
    h_mumu_mc_count->SetDirectory(0);
    h_mumu_mc = new TH2F("h_mumu_mc", "Events in bins for MC templates",
            n_xf_bins, xf_bins, n_cost_bins, cost_bins);
    h_mumu_mc->SetDirectory(0);
    h_mumu_mc_errs = new TH2F("h_mumu_mc_errs", "Events in bins for MC templates",
            n_xf_bins, xf_bins, n_cost_bins, cost_bins);
    h_mumu_mc_errs->SetDirectory(0);
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

    TH2F *h_mumu_mc_corr = (TH2F *)h_mumu_mc->Clone("h_mumu_mc_corr");
    TH2F *h_mumu_mc_inc = (TH2F *) h_mumu_mc->Clone("h_mumu_mc_inc");
    h_mumu_dilu = (TH2F *) h_mumu_mc->Clone("h_mumu_dilu");

    gen_mc_template(t_mumu_mc, alpha, h_mumu_sym, h_mumu_asym, h_mumu_sym_count, m_low, m_high, FLAG_MUONS);
    //gen_finite_mc_template(t_mumu_mc, h_mumu_mc_corr, h_mumu_mc_inc, h_mumu_mc_errs, h_mumu_mc_count, m_low, m_high, FLAG_MUONS);
    //h_mumu_mc->Add(h_mumu_mc_inc, h_mumu_mc_corr);
    h_mumu_fitted_mc = (TH2F *) h_mumu_mc->Clone("h_mumu_fitted_mc");
    //get_dilu_hist(h_mumu_dilu, h_mumu_mc_corr, h_mumu_mc_inc);
    //gen_fit_templates_from_mc_and_dilu(h_mumu_sym, h_mumu_asym, h_mumu_mc, h_mumu_dilu);
    //gen_fit_templates_from_finite_mc(h_mumu_sym, h_mumu_asym, h_mumu_mc_corr, h_mumu_mc_inc);
    print_hist(h_mumu_sym);
    print_hist(h_mumu_asym);
    TTree *mumu_ts[2] = {t_mumu_back, t_mumu_nosig};

    gen_fakes_template(t_mumu_WJets, t_mumu_QCD, t_mumu_WJets_contam, t_mumu_QCD_contam, h_mumu_back, m_low, m_high, FLAG_MUONS);
    gen_combined_background_template(2, mumu_ts, h_mumu_back, m_low, m_high, FLAG_MUONS);

    nMuMu_DataEvents = gen_data_template(t_mumu_data, h_mumu_data,  &v_mumu_xF, &v_mumu_cost, m_low, m_high, FLAG_MUONS);
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
    v_elel_cost.clear();
    v_elel_xF.clear();
    v_mumu_cost.clear();
    v_mumu_xF.clear();
    printf("Finishing cleanup\n");
}

void finite_mc_fit(){
    Double_t AFB_fit[n_m_bins], AFB_fit_err[n_m_bins], r_elel_back_fit[n_m_bins], r_elel_back_fit_err[n_m_bins], 
             r_mumu_back_fit[n_m_bins], r_mumu_back_fit_err[n_m_bins];

    init();
    TTree *tout= new TTree("T_fit_res", "Tree with Fit Results");
    tout->SetDirectory(0);

    Double_t AFB, AFB_err, r_elel_back, r_elel_back_err, r_mumu_back, r_mumu_back_err;

    tout->Branch("var_low", &m_low);
    tout->Branch("var_high", &m_high);
    tout->Branch("nElElEvents", &nElEl_DataEvents);
    tout->Branch("nMuMuEvents", &nMuMu_DataEvents);
    tout->Branch("AFB", &AFB);
    tout->Branch("AFB_err", &AFB_err);
    tout->Branch("r_elel_back", &r_elel_back);
    tout->Branch("r_elel_back_err", &r_elel_back_err);
    tout->Branch("r_mumu_back", &r_mumu_back);
    tout->Branch("r_mumu_back_err", &r_mumu_back_err);

    unsigned int nElElEvents[n_m_bins];
    unsigned int nMuMuEvents[n_m_bins];

    Double_t r_mumu_back_starts[] = {0.124, 0.143, 0.16, 0.16, 0.17, 0.17, 0.17, 0.1, 0.1};
    Double_t r_elel_back_starts[] = {0.1, 0.16, 0.2, 0.2, 0.17, 0.16, 0.16, 0.13, 0.13};

    for(int i=0; i<1; i++){
        printf("Starting loop \n");
        m_low = m_bins[i];
        m_high = m_bins[i+1];
        alpha = alphas[i];

        setup();

        //float AFB_starts[6] = {0.64, 0.617, 0.589, 0.567, 0.535, 0.518};
        //float r_back_starts[6] = {0.087, 0.156, 0.211, 0.185, 0.171, 0.128};

        float AFB_starts[6] = {0.59, 0.608, 0.634, 0.642, 0.584, 0.605};
        float r_back_starts[6] = {0.110, 0.175, 0.183, 0.186, 0.166, 0.095};



        float AFB_start = AFB_starts[i];
        float AFB_start_error = 0.04;
        float AFB_max = 0.75;
        float r_back_start = r_back_starts[i];
        float r_back_start_error = 0.04;
        float r_back_max = 0.6;

        int n_params = 2+ (n_cost_bins * n_xf_bins);

        TVirtualFitter * minuit = TVirtualFitter::Fitter(0,n_params);
        minuit->SetFCN(fcn);
        minuit->SetParameter(0,"AFB", AFB_start, AFB_start_error, -AFB_max, AFB_max);
        minuit->SetParameter(1,"r_back", r_back_start, r_back_start_error, 0., 0.5);
        setHistParams(minuit, h_mumu_mc, 2);
        Double_t arglist[100];
        arglist[0] = 10000.;
        //minuit->ExecuteCommand("MIGRAD", arglist,0);
        Double_t up = 1.0;
        minuit->SetErrorDef(up);
        arglist[0] = 0.;
        //minuit->ExecuteCommand("MINOS", arglist, 0);



        AFB_fit[i] = minuit->GetParameter(0); 
        AFB_fit_err[i] = minuit->GetParError(0);
        r_elel_back_fit[i] = minuit->GetParameter(1); 
        r_elel_back_fit_err[i] = minuit->GetParError(1);

        nElElEvents[i] = nElEl_DataEvents;
        //nMuMuEvents[i] = nMuMu_DataEvents;

        AFB= AFB_fit[i];
        AFB_err = AFB_fit_err[i];

        //r_elel_back= r_elel_back_fit[i];
        //r_elel_back_err = r_elel_back_fit_err[i];

        r_mumu_back= r_mumu_back_fit[i];
        r_mumu_back_err = r_mumu_back_fit_err[i];
        tout->Fill();

        cleanup();
    }
    TFile *fout = TFile::Open(fout_name, "RECREATE");
    fout->cd();
    tout->Write();
    for(int i=0; i<1; i++){
        printf("\n Fit on M=[%.0f, %.0f], %i ElEl Events, %i MuMu Events: AFB = %0.3f +/- %0.3f r_elel_back = %0.3f +/- %0.3f r_mumu_back =%0.3f +/- %0.3f \n", 
                    m_bins[i], m_bins[i+1], nElElEvents[i], nMuMuEvents[i], AFB_fit[i], AFB_fit_err[i], r_elel_back_fit[i], r_elel_back_fit_err[i], r_mumu_back_fit[i], r_mumu_back_fit_err[i]);

    }
    //printf("fit results written to %s \n", fout_name.Data());

}



