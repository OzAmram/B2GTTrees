
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
//#include "../TemplateMaker_systematics.C"
#include "../../TemplateMaker.C"



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
Double_t m_alphas[6] = {0.095, 0.08695, 0.0762, 0.112, 0.065, 0.0605};
Double_t alpha_unc[6] = {0.015, 0.015, 0.015, 0.03,   0.02, 0.015};
Float_t pt_bins[] =        {0.,25.,  50., 80.,   120.,   200., 10000.};
Double_t pt_alphas[6] =    {0.007, 0.136, 0.337, 0.546, 0.776, 0.945};
Double_t pt_alpha_unc[6] = {0.006, 0.015, 0.035, 0.05, 0.08, .15};
Double_t alpha;


//int FLAG1 = FLAG_ELECTRONS;
int FLAG1 = FLAG_MUONS;
const TString fout_name("AFB_fit/fit_results/MuMu_pt_fit_jan24_nominal.root");

int FLAG2 = FLAG_PT_BINS;
int n_bins = n_pt_bins;
//int FLAG2 = FLAG_M_BINS;


float var_low;
float var_high;
//alpha = 0.0981;

bool print = true;


Double_t med_btag = 0.4432;

TH2F *h_asym, *h_sym, *h_back,  *h_data, *h_mc;
TH2F *h_mc_count, *h_sym_count;
//void *x_axp, *y_axp, *z_axp;

//define axis globally for convinience 
TAxis *x_ax, *y_ax, *z_ax;




TFile *f_mc, *f_back, *f_data, *f_QCD, *f_WJets, *f_WJets_contam, *f_QCD_contam;
TTree *t_mc, *t_back, *t_data, *t_QCD, *t_WJets, *t_WJets_contam, *t_QCD_contam, *t_nosig;


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

void init(){
    if(FLAG1 == FLAG_ELECTRONS){
         f_mc = (TFile*) TFile::Open("output_files/ElEl_DY_jan22.root");
         t_mc = (TTree *) f_mc ->Get("T_data");
         t_nosig = (TTree *) f_mc ->Get("T_back");
         f_back = (TFile*) TFile::Open("output_files/ElEl_combined_back_jan22.root");
         t_back = (TTree *) f_back ->Get("T_data");

         f_data = TFile::Open("output_files/SingleElectron_data_jan22.root");
         t_data = (TTree *)f_data->Get("T_data"); 

         f_QCD = TFile::Open("../analyze/output_files/ElEl_QCD_est_nov2.root");
         t_QCD = (TTree *)f_QCD->Get("T_data");

         f_WJets = TFile::Open("../analyze/output_files/ElEl_WJets_est_nov2.root");
         t_WJets = (TTree *)f_WJets->Get("T_data");

         f_WJets_contam = TFile::Open("../analyze/FakeRate/root_files/ElEl_fakerate_WJets_MC_dec4.root");
         t_WJets_contam = (TTree *)f_WJets_contam->Get("T_data");

         f_QCD_contam = TFile::Open("../analyze/FakeRate/root_files/ElEl_fakerate_QCD_MC_dec4.root");
         t_QCD_contam = (TTree *)f_QCD_contam->Get("T_data");
    }
    else{
         f_mc = (TFile*) TFile::Open("output_files/MuMu_DY_jan16.root");
         t_mc = (TTree *) f_mc ->Get("T_data");
         t_nosig = (TTree *) f_mc ->Get("T_back");
         f_back = (TFile*) TFile::Open("output_files/MuMu_combined_back_jan22.root");
         t_back = (TTree *) f_back ->Get("T_data");

         f_data = TFile::Open("output_files/SingleMuon_data_jan22.root");
         t_data = (TTree *)f_data->Get("T_data"); 

         f_QCD = TFile::Open("../analyze/output_files/MuMu_QCD_est_nov2.root");
         t_QCD = (TTree *)f_QCD->Get("T_data");

         f_WJets = TFile::Open("../analyze/output_files/MuMu_WJets_est_Nov2.root");
         t_WJets = (TTree *)f_WJets->Get("T_data");

         f_WJets_contam = TFile::Open("../analyze/FakeRate/root_files/MuMu_fakerate_Wjets_MC_dec4.root");
         t_WJets_contam = (TTree *)f_WJets_contam->Get("T_data");

         f_QCD_contam = TFile::Open("../analyze/FakeRate/root_files/MuMu_fakerate_QCD_MC_dec4.root");
         t_QCD_contam = (TTree *)f_QCD_contam->Get("T_data");
    }
    return;
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

    gen_mc_template(t_mc, alpha, h_sym, h_asym, h_sym_count, var_low, var_high, FLAG1, FLAG2);
    TTree *ts[2] = {t_back, t_nosig};

    gen_fakes_template(t_WJets, t_QCD, t_WJets_contam, t_QCD_contam, h_back, var_low, var_high, FLAG1, FLAG2);
    gen_combined_background_template(2, ts, h_back, var_low, var_high, FLAG1, FLAG2);

    nDataEvents = gen_data_template(t_data, h_data, &v_xF, &v_cost, var_low, var_high, FLAG1, FLAG2);
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
    Double_t AFB_fit[n_bins], AFB_fit_err[n_bins], r_back_fit[n_bins], r_back_fit_err[n_bins];
    unsigned int nEvents[n_bins];
    TTree *tout= new TTree("T_fit_res", "Tree with Fit Results");
    tout->SetDirectory(0);

    Double_t AFB, AFB_err, r_back, r_back_err;

    tout->Branch("var_low", &var_low);
    tout->Branch("var_high", &var_high);
    tout->Branch("nEvents", &nDataEvents);
    tout->Branch("AFB", &AFB);
    tout->Branch("AFB_err", &AFB_err);
    tout->Branch("r_back", &r_back);
    tout->Branch("r_back_err", &r_back_err);

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

        printf("Integrals are %f %f %f %f  \n", h_data->Integral(), h_sym->Integral(), 
                                               h_asym->Integral(), h_back->Integral() );
        h_sym->Print();



        float AFB_start = 0.5;
        float AFB_start_error = 0.1;
        float AFB_max = 0.75;
        float r_back_start = 0.2;
        float r_back_start_error = 0.1;
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
    TFile *fout = TFile::Open(fout_name, "RECREATE");
    fout->cd();
    tout->Write();
    for(int i=0; i<n_bins; i++){
        if (FLAG2 == FLAG_M_BINS){
            printf("\n Fit on M=[%.0f, %.0f], %i Events: AFB = %0.3f +/- %0.3f r_back = %0.3f +/- %0.3f \n", 
                        m_bins[i], m_bins[i+1], nEvents[i], AFB_fit[i], AFB_fit_err[i], r_back_fit[i], r_back_fit_err[i]);
        }
        if (FLAG2 == FLAG_PT_BINS){
            printf("\n Fit on pt=[%.0f, %.0f], %i Events: AFB = %0.3f +/- %0.3f r_back = %0.3f +/- %0.3f \n", 
                        pt_bins[i], pt_bins[i+1], nEvents[i], AFB_fit[i], AFB_fit_err[i], r_back_fit[i], r_back_fit_err[i]);
        }

    }
    printf("fit results written to %s \n", fout_name.Data());

}



