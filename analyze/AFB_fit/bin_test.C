
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
//#include "../TemplateMaker_systematics.C"
#include "../TemplateMaker.C"
#include "root_files.h"



int n_xf_bins;
Float_t *xf_bins;
int n_xf_bins_v1 = 5;
Float_t xf_bins_v1[] = {0., 0.02, 0.04, 0.08, 0.12, 1.0};
int n_xf_bins_v2 = 6;
Float_t xf_bins_v2[] = {0., 0.02, 0.04, 0.08, 0.12, 0.25, 1.0};
//int n_cost_bins = 8;
//Float_t cost_bins[] = {-1.0, -.75, -.5, -.25, 0., 0.25, 0.5,  0.75, 1.0};
int n_cost_bins = 10;
Float_t cost_bins[] = {-1.0, -.8, -.6, -.4, -.2, 0., 0.2, 0.4, 0.6, 0.8, 1.0};
//int n_cost_bins = 12;
//Float_t cost_bins[] = {-1.0, -0.8333, -0.6667, -0.5, -0.3333, -0.1667, 0., 0.1667, 0.3333, 0.5, 0.6667, 0.8333, 1.0};
//int n_cost_bins = 14;
//Float_t cost_bins[] = {-1.0, -.857, -.714, -.571, -.429, -0.286, -.143,  0., 0.143, .286, 0.429, 0.571, 0.714, 0.857, 1.0};
int n_m_bins = 6;
Float_t m_bins[] = {150,200,   250,    350,    500,    700, 100000};
Double_t alphas[6] = {0.109, 0.078, 0.0762, 0.112, 0.065, 0.06};
Double_t alpha_unc[6] = {0.015, 0.015, 0.02, 0.03,   0.02, 0.02};
Double_t alpha;

const TString fout_name("AFB_fit/fit_results/m_bins/combined_fit_nominal_mar19.root");



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

        if(p_elel_back < 0){ 
            //Printf("P_back < 0! %.2e %.2f %.2f \n", p_back, v_xF[i], v_cost[i]);
            p_elel_back = 1e-20;
        }


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


        if(p_mumu_back < 0){ 
            //Printf("P_back < 0! %.2e %.2f %.2f \n", p_back, v_xF[i], v_cost[i]);
            p_mumu_back = 1e-20;
        }

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
    if(m_low <= m_bins[n_m_bins-2] ){
        n_xf_bins = n_xf_bins_v1;
        xf_bins = xf_bins_v1;
    }
    else{
        n_xf_bins = n_xf_bins_v2;
        xf_bins = xf_bins_v2;
    }

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

    nElEl_DataEvents = gen_data_template(t_elel_data, h_elel_data,  &v_elel_xF, &v_elel_cost, m_low, m_high, FLAG_ELECTRONS);

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

    nMuMu_DataEvents = gen_data_template(t_mumu_data, h_mumu_data, &v_mumu_xF, &v_mumu_cost, m_low, m_high, FLAG_MUONS);
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
void bin_test(){
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

    for(int i=0; i<n_m_bins; i++){
        printf("Starting loop \n");
        m_low = m_bins[i];
        m_high = m_bins[i+1];
        alpha = alphas[i];

        setup();

        for(int i=0; i<= n_xf_bins; i++){
            printf("\n");
            for(int j=0; j<= n_cost_bins; j++){
                printf("%.0f ",    h_elel_sym_count->GetBinContent(i,j));
            }
        }


        cleanup();
    }

}



