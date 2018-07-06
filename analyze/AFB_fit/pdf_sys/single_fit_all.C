
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
#include "../../TemplateMaker_systematics.C"
//#include "../../TemplateMaker.C"
#include "../root_files.h"




int FLAG1 = FLAG_ELECTRONS;
//int FLAG1 = FLAG_MUONS;
const TString mumu_fout_name("AFB_fit/fit_results/m_bins/MuMu_m_fit_july6_pdf_sys.root");
const TString elel_fout_name("AFB_fit/fit_results/m_bins/ElEl_m_fit_july6_pdf_sys.root");

bool do_both = true;
Int_t N_PDFS = 100;
//int FLAG2 = FLAG_PT_BINS;
//int n_bins = n_pt_bins;
int FLAG2 = FLAG_M_BINS;
int n_bins = n_m_bins;


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

    if(FLAG1 == FLAG_MUONS){
        gen_mc_template(t_mumu_mc, alpha, h_sym, h_asym, h_sym_count, var_low, var_high, FLAG1, FLAG2);
        TTree *ts[2] = {t_mumu_back, t_mumu_nosig};

        gen_fakes_template(t_mumu_WJets, t_mumu_QCD, t_mumu_WJets_contam, t_mumu_QCD_contam, h_back, var_low, var_high, FLAG1, FLAG2);
        gen_combined_background_template(2, ts, h_back, var_low, var_high, FLAG1, FLAG2);

        nDataEvents = gen_data_template(t_mumu_data, h_data, &v_xF, &v_cost, var_low, var_high, FLAG1, FLAG2);
    }
    if(FLAG1 == FLAG_ELECTRONS){
        gen_mc_template(t_elel_mc, alpha, h_sym, h_asym, h_sym_count, var_low, var_high, FLAG1, FLAG2);
        TTree *ts[2] = {t_elel_back, t_elel_nosig};

        gen_fakes_template(t_elel_WJets, t_elel_QCD, t_elel_WJets_contam, t_elel_QCD_contam, h_back, var_low, var_high, FLAG1, FLAG2);
        gen_combined_background_template(2, ts, h_back, var_low, var_high, FLAG1, FLAG2);

        nDataEvents = gen_data_template(t_elel_data, h_data, &v_xF, &v_cost, var_low, var_high, FLAG1, FLAG2);
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
    Float_t AFB_fit[N_PDFS], AFB_fit_err[N_PDFS], r_back_fit[N_PDFS], r_back_fit_err[N_PDFS];
    unsigned int nEvents[n_bins];
    TTree *tout= new TTree("T_fit_res", "Tree with Fit Results");
    tout->SetDirectory(0);

    Double_t AFB, AFB_err, r_back, r_back_err;

    tout->Branch("var_low", &var_low);
    tout->Branch("var_high", &var_high);
    tout->Branch("nEvents", &nDataEvents);
    tout->Branch("N_PDFS", &N_PDFS);
    tout->Branch("AFB", &AFB_fit, "AFB_fit[N_PDFS]/F");
    tout->Branch("AFB_err", &AFB_fit_err, "AFB_fit_err[N_PDFS]/F");
    tout->Branch("r_back", &r_back_fit, "r_back_fit[N_PDFS]/F");
    tout->Branch("r_back_err", &r_back_fit_err, "r_back_fit_err[N_PDFS]/F");

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

        for(int j = 0; j < N_PDFS; j++){

            printf("generating new MC template with new pdf \n");
            if(FLAG1 == FLAG_MUONS){
                gen_mc_template(t_mumu_mc, alpha, h_sym, h_asym, h_sym_count, var_low, var_high, FLAG1, FLAG2, j);
            }
            if(FLAG1 == FLAG_ELECTRONS){
                gen_mc_template(t_elel_mc, alpha, h_sym, h_asym, h_sym_count, var_low, var_high, FLAG1, FLAG2, j);
            }


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



            AFB_fit[j] = minuit->GetParameter(0); 
            AFB_fit_err[j] = minuit->GetParError(0);
            r_back_fit[j] = minuit->GetParameter(1); 
            r_back_fit_err[j] = minuit->GetParError(1);
        }

        tout->Fill();



        cleanup();
    }
    TFile *fout;
    if(FLAG1 == FLAG_MUONS) fout = TFile::Open(mumu_fout_name, "RECREATE");
    if(FLAG1 == FLAG_ELECTRONS) fout = TFile::Open(elel_fout_name, "RECREATE");
    fout->cd();
    tout->Write();
    if(do_both){
        
        do_both = false;
        FLAG1 = 1 - FLAG1;
        printf("Starting 2nd run \n");
        single_fit_all();
        FLAG1 = 1 - FLAG1;
    }
    if(FLAG1 == FLAG_MUONS){
        printf("MuMu fit results, written out to %s \n" , mumu_fout_name.Data());
    }
    if(FLAG1 == FLAG_ELECTRONS){
        printf("ElEl fit results, written out to %s \n" , elel_fout_name.Data());
    }

}



