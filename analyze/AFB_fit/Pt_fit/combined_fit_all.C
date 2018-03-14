
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
//#include "../../TemplateMaker_systematics.C"
#include "../../TemplateMaker.C"
#include "../root_files.h"



int n_xf_bins = 4;
float xf_max = 1.0;
Float_t xf_bins[] = {0., 0.05, 0.1, 0.20, 1.0};
//int n_cost_bins = 8;
//Float_t cost_bins[] = {-1.0, -.75, -.5, -.25, 0., 0.25, 0.5,  0.75, 1.0};
int n_cost_bins = 10;
Float_t cost_bins[] = {-1.0, -.8, -.6, -.4, -.2, 0., 0.2, 0.4, 0.6, 0.8, 1.0};
//int n_cost_bins = 12;
//Float_t cost_bins[] = {-1.0, -.8333, -.6667, -.5, -.3333, -0.1667, 0., 0.1667, 0.3333, 0.5, 0.6667, 0.8333, 1.0};
//int n_cost_bins = 14;
//Float_t cost_bins[] = {-1.0, -.857, -.714, -.571, -.429, -0.286, -.143,  0., 0.143, .286, 0.429, 0.571, 0.714, 0.857, 1.0};
int n_m_bins = 6;
int n_pt_bins = 6;
Float_t m_bins[] = {150,200,   250,    350,    500,    700, 100000};
Double_t m_alphas[6] = {0.109, 0.078, 0.0762, 0.112, 0.065, 0.06};
Double_t m_alpha_unc[6] = {0.015, 0.015, 0.02, 0.03,   0.02, 0.02};
Float_t pt_bins[] =        {0.,25.,  50., 80.,   120.,   200., 10000.};
Double_t pt_alphas[6] =    {0.007, 0.136, 0.337, 0.546, 0.776, 0.945};
Double_t pt_alpha_unc[6] = {0.006, 0.015, 0.035, 0.05, 0.08, .15};
Double_t alpha;

const TString fout_name("AFB_fit/fit_results/pt_bins/combined_pt_fit_jan31_nominal.root");


int FLAG2 = FLAG_PT_BINS;
int n_bins = n_pt_bins;
//int FLAG2 = FLAG_M_BINS;
//int n_bins = n_m_bins;



float var_low;
float var_high;
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


    gen_mc_template(t_elel_mc, alpha, h_elel_sym, h_elel_asym, h_elel_sym_count, var_low, var_high, FLAG_ELECTRONS, FLAG2);
    TTree *elel_ts[2] = {t_elel_back, t_elel_nosig};

    gen_fakes_template(t_elel_WJets, t_elel_QCD, t_elel_WJets_contam, t_elel_QCD_contam, h_elel_back, var_low, var_high, FLAG_ELECTRONS,FLAG2);
    gen_combined_background_template(2, elel_ts, h_elel_back, var_low, var_high, FLAG_ELECTRONS, FLAG2);

    nElEl_DataEvents = gen_data_template(t_elel_data, h_elel_data,  &v_elel_xF, &v_elel_cost, var_low, var_high, FLAG_ELECTRONS,FLAG2);

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


    gen_mc_template(t_mumu_mc, alpha, h_mumu_sym, h_mumu_asym, h_mumu_sym_count, var_low, var_high, FLAG_MUONS,FLAG2);
    TTree *mumu_ts[2] = {t_mumu_back, t_mumu_nosig};
    gen_fakes_template(t_mumu_WJets, t_mumu_QCD, t_mumu_WJets_contam, t_mumu_QCD_contam, h_mumu_back, var_low, var_high, FLAG_MUONS,FLAG2);
    gen_combined_background_template(2, mumu_ts, h_mumu_back, var_low, var_high, FLAG_MUONS,FLAG2);

    nMuMu_DataEvents = gen_data_template(t_mumu_data, h_mumu_data, &v_mumu_xF, &v_mumu_cost, var_low, var_high, FLAG_MUONS,FLAG2);
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
void combined_fit_all(){
    Double_t AFB_fit[n_bins], AFB_fit_err[n_bins], r_elel_back_fit[n_bins], r_elel_back_fit_err[n_bins], 
             r_mumu_back_fit[n_bins], r_mumu_back_fit_err[n_bins];

    init();
    TTree *tout= new TTree("T_fit_res", "Tree with Fit Results");
    tout->SetDirectory(0);

    Double_t AFB, AFB_err, r_elel_back, r_elel_back_err, r_mumu_back, r_mumu_back_err;

    tout->Branch("var_low", &var_low);
    tout->Branch("var_high", &var_high);
    tout->Branch("nElElEvents", &nElEl_DataEvents);
    tout->Branch("nMuMuEvents", &nMuMu_DataEvents);
    tout->Branch("AFB", &AFB);
    tout->Branch("AFB_err", &AFB_err);
    tout->Branch("r_elel_back", &r_elel_back);
    tout->Branch("r_elel_back_err", &r_elel_back_err);
    tout->Branch("r_mumu_back", &r_mumu_back);
    tout->Branch("r_mumu_back_err", &r_mumu_back_err);

    unsigned int nElElEvents[n_bins];
    unsigned int nMuMuEvents[n_bins];

    Double_t r_mumu_back_starts[] = {0.124, 0.143, 0.16, 0.16, 0.17, 0.17, 0.17, 0.1, 0.1};
    Double_t r_elel_back_starts[] = {0.1, 0.16, 0.2, 0.2, 0.17, 0.16, 0.16, 0.13, 0.13};

    for(int i=0; i<n_bins; i++){
        printf("Starting loop \n");
        if (FLAG2 == FLAG_M_BINS){
            var_low = m_bins[i];
            var_high = m_bins[i+1];
            alpha = m_alphas[i];
        }
        else if (FLAG2 == FLAG_PT_BINS){
            var_low = pt_bins[i];
            var_high = pt_bins[i+1];
            alpha = pt_alphas[i];
        }

        setup();

        /*
           printf("ElEL: Integrals are %f %f %f %f  \n", h_elel_data->Integral(), h_elel_sym->Integral(), 
           h_elel_asym->Integral(), h_elel_back->Integral() );
           printf("MuMu: Integrals are %f %f %f %f  \n", h_mumu_data->Integral(), h_mumu_sym->Integral(), 
           h_mumu_asym->Integral(), h_mumu_back->Integral() );
           h_elel_sym->Print();
           h_mumu_sym->Print();
           */



        float AFB_start = 0.5;
        float AFB_start_error = 0.1;
        float AFB_max = 0.75;
        float r_back_start = 0.2;
        float r_back_start_error = 0.1;
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

        AFB= AFB_fit[i];
        AFB_err = AFB_fit_err[i];

        r_elel_back= r_elel_back_fit[i];
        r_elel_back_err = r_elel_back_fit_err[i];

        r_mumu_back= r_mumu_back_fit[i];
        r_mumu_back_err = r_mumu_back_fit_err[i];
        tout->Fill();

        cleanup();
    }
    TFile *fout = TFile::Open(fout_name, "RECREATE");
    fout->cd();
    tout->Write();
    for(int i=0; i<n_bins; i++){
        if(FLAG2 == FLAG_M_BINS){
            printf("\n Fit on M=[%.0f, %.0f], %i ElEl Events, %i MuMu Events: AFB = %0.3f +/- %0.3f r_elel_back = %0.3f +/- %0.3f r_mumu_back =%0.3f +/- %0.3f \n", 
                    m_bins[i], m_bins[i+1], nElElEvents[i], nMuMuEvents[i], AFB_fit[i], AFB_fit_err[i], r_elel_back_fit[i], r_elel_back_fit_err[i], r_mumu_back_fit[i], r_mumu_back_fit_err[i]);
        }
        else if(FLAG2 == FLAG_PT_BINS){
            printf("\n Fit on Pt=[%.0f, %.0f], %i ElEl Events, %i MuMu Events: AFB = %0.3f +/- %0.3f r_elel_back = %0.3f +/- %0.3f r_mumu_back =%0.3f +/- %0.3f \n", 
                    pt_bins[i], pt_bins[i+1], nElElEvents[i], nMuMuEvents[i], AFB_fit[i], AFB_fit_err[i], r_elel_back_fit[i], r_elel_back_fit_err[i], r_mumu_back_fit[i], r_mumu_back_fit_err[i]);
        }

    }
    printf("fit results written to %s \n", fout_name.Data());

}



