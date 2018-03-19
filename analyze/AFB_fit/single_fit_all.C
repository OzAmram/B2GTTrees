
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
#include "root_files.h"
//#include "../TemplateMaker.C"





//int FLAG = FLAG_ELECTRONS;
int FLAG = FLAG_MUONS;
bool do_both = true;
const TString mumu_fout_name("AFB_fit/fit_results/m_bins/MuMu_fit_mu_HLT_up_mar19.root");
const TString elel_fout_name("AFB_fit/fit_results/m_bins/ElEl_fit_el_HLT_up_mar19.root");


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






vector<double> v_xF;
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
        if(p_back < 0){ 
            //Printf("P_back < 0! %.2e %.2f %.2f \n", p_back, v_xF[i], v_cost[i]);
            p_back = 1e-20;
        }

        double prob = r_back*p_back + (1 - r_back) * (p_sym + AFB*p_asym);
        if(prob > 1) printf("Warning prob is too big %.2f %.2f %.2f %.2f %.2f  \n", r_back, AFB, p_back, p_sym, p_asym);
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

    if(FLAG == FLAG_MUONS){
        gen_mc_template(t_mumu_mc, alpha, h_sym, h_asym, h_sym_count, m_low, m_high, FLAG);
        TTree *ts[2] = {t_mumu_back, t_mumu_nosig};

        gen_fakes_template(t_mumu_WJets, t_mumu_QCD, t_mumu_WJets_contam, t_mumu_QCD_contam, h_back, m_low, m_high, FLAG);
        gen_combined_background_template(2, ts, h_back, m_low, m_high, FLAG);

        nDataEvents = gen_data_template(t_mumu_data, h_data, &v_xF, &v_cost, m_low, m_high, FLAG);
    }
    if(FLAG == FLAG_ELECTRONS){
        gen_mc_template(t_elel_mc, alpha, h_sym, h_asym, h_sym_count, m_low, m_high, FLAG);
        TTree *ts[2] = {t_elel_back, t_elel_nosig};

        gen_fakes_template(t_elel_WJets, t_elel_QCD, t_elel_WJets_contam, t_elel_QCD_contam, h_back, m_low, m_high, FLAG);
        gen_combined_background_template(2, ts, h_back, m_low, m_high, FLAG);

        nDataEvents = gen_data_template(t_elel_data, h_data, &v_xF, &v_cost, m_low, m_high, FLAG);
    }
    printf("\n\n\n Printing MC counts in each bin:\n");
    for(int i=1; i<=n_xf_bins; i++){
        for(int j=1; j<=n_cost_bins; j++){
            printf("%.0f ", h_sym_count->GetBinContent(i,j));
        }
        printf("\n");
    }

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
    v_cost.clear();
    v_xF.clear();
    printf("Finishing cleanup\n");
}
void single_fit_all(){
    init();
    Double_t AFB_fit[n_m_bins], AFB_fit_err[n_m_bins], r_back_fit[n_m_bins], r_back_fit_err[n_m_bins];
    unsigned int nEvents[n_m_bins];
    TTree *tout= new TTree("T_fit_res", "Tree with Fit Results");
    tout->SetDirectory(0);

    Double_t AFB, AFB_err, r_back, r_back_err;

    tout->Branch("var_low", &m_low);
    tout->Branch("var_high", &m_high);
    tout->Branch("nEvents", &nDataEvents);
    tout->Branch("AFB", &AFB);
    tout->Branch("AFB_err", &AFB_err);
    tout->Branch("r_back", &r_back);
    tout->Branch("r_back_err", &r_back_err);

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
        minuit->SetParameter(1,"r_back", r_back_start, r_back_start_error, 0., 0.5);
        Double_t arglist[100];
        arglist[0] = 10000.;
        minuit->ExecuteCommand("MIGRAD", arglist,0);
        Double_t up = 1.0;
        minuit->SetErrorDef(up);
        arglist[0] = 0.;
        minuit->ExecuteCommand("MINOS", arglist, 0);



        AFB_fit[i] = minuit->GetParameter(0); 
        AFB_fit_err[i] = minuit->GetParError(0);
        r_back_fit[i] = minuit->GetParameter(1); 
        r_back_fit_err[i] = minuit->GetParError(1);

        AFB = AFB_fit[i];
        AFB_err = AFB_fit_err[i];

        r_back = r_back_fit[i];
        r_back_err = r_back_fit_err[i];

        nEvents[i] = nDataEvents;
        tout->Fill();



        cleanup();
    }
    TFile *fout;
    if(FLAG == FLAG_MUONS) fout = TFile::Open(mumu_fout_name, "RECREATE");
    if(FLAG == FLAG_ELECTRONS) fout = TFile::Open(elel_fout_name, "RECREATE");
    fout->cd();
    tout->Write();
    if(do_both){
        
        do_both = false;
        FLAG = 1 - FLAG;
        printf("Starting 2nd run \n");
        single_fit_all();
        FLAG = 1 - FLAG;
    }

    if(FLAG == FLAG_MUONS){
        printf("MuMu fit results, written out to %s \n" , mumu_fout_name.Data());
    }
    if(FLAG == FLAG_ELECTRONS){
        printf("ElEl fit results, written out to %s \n" , elel_fout_name.Data());
    }
    for(int i=0; i<n_m_bins; i++){
        printf("\n Fit on M=[%.0f, %.0f], %i Events: AFB = %0.3f +/- %0.3f r_back = %0.3f +/- %0.3f \n", 
                    m_bins[i], m_bins[i+1], nEvents[i], AFB_fit[i], AFB_fit_err[i], r_back_fit[i], r_back_fit_err[i]);

    }

}



