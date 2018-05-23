
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
#include "root_files.h"



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

void print_hist(TH2F *h){
    for(int i=1; i<= n_xf_bins; i++){
        printf("\n");
        for(int j=1; j<= n_cost_bins; j++){
            printf("%.2e ",    h->GetBinContent(i,j));
        }
    }
    printf("\n\n");
    return;
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

    printf("Finishing setup \n");
    return;
}


void compute_afb_shift(){
    Double_t mumu_rbks[] = {0.109, 0.174, 0.185, 0.186, 0.166, 0.095};
    Double_t elel_rbks[] = {0.087, 0.157, 0.210, 0.185, 0.176, 0.128};
    Double_t AFBs[] = {0.611, 0.612, 0.614, 0.608, 0.559, 0.532};
    Double_t AFB_shifts[6][11];

    Double_t mumu_rbk, elel_rbk, AFB;
    TH1F *utype_dilu, *dtype_dilu, *mix_dilu;

    init();
    for(int k=0; k<6; k++){
        printf("Starting loop \n");
        m_low = m_bins[k];
        m_high = m_bins[k+1];
        alpha = alphas[k];
        mumu_rbk = mumu_rbks[k];
        elel_rbk = elel_rbks[k];
        AFB = AFBs[k];

        setup();

        char f_in[100];
        sprintf(f_in, "../setlimits/dilus/M%i_dilus.root", (int) m_low);
        printf("Getting dilus from file %s \n", f_in);
        auto f_dilu = (TFile*) TFile::Open(f_in);
        TH1F *utype_dilu = (TH1F *) f_dilu ->Get("utype_dilu");
        TH1F *dtype_dilu = (TH1F *) f_dilu ->Get("dtype_dilu");
        TH1F *mix_dilu = (TH1F *) f_dilu ->Get("mix_dilu");

        TH1F *h_dilu_ratio = new TH1F("dilu_ratio", "", n_xf_bins, xf_bins);
        for(int l=0; l<=10; l++){
            Float_t x = 0.1*((float) l);
            for(int i=1; i<=n_xf_bins; i++){
                Float_t u_dilu = utype_dilu->GetBinContent(i);
                Float_t d_dilu = dtype_dilu->GetBinContent(i);
                Float_t m_dilu = mix_dilu->GetBinContent(i);
                Float_t ratio = (u_dilu*x + d_dilu*(1-x))/m_dilu;
                h_dilu_ratio->SetBinContent(i,ratio);
            }

            //TH1F *h_dilu_ratio = (TH1F *)dtype_dilu->Clone("h_dilu_ratio");
            //h_dilu_ratio->Divide(mix_dilu);

            TVectorD a(3);
            TMatrixDSym m(3);
            for(int i=1; i<=n_xf_bins; i++){
                for(int j=1; j<=n_cost_bins; j++){
                    Double_t p_mumu_asym = h_mumu_asym->GetBinContent(i,j);
                    Double_t p_mumu_sym = h_mumu_sym->GetBinContent(i,j);
                    Double_t p_mumu_bk = h_mumu_back->GetBinContent(i,j);

                    Double_t p_elel_asym = h_elel_asym->GetBinContent(i,j);
                    Double_t p_elel_sym = h_elel_sym->GetBinContent(i,j);
                    Double_t p_elel_bk = h_elel_back->GetBinContent(i,j);

                    Double_t dilu_ratio = h_dilu_ratio->GetBinContent(i);

                    Double_t g = mumu_rbk*p_mumu_bk + (1.-mumu_rbk) *(p_mumu_sym + AFB*dilu_ratio*p_mumu_asym) +
                        elel_rbk*p_elel_bk + (1.-elel_rbk) *(p_elel_sym + AFB*dilu_ratio*p_elel_asym);
                    Double_t f = mumu_rbk*p_mumu_bk + (1.-mumu_rbk) *(p_mumu_sym + AFB*p_mumu_asym) +
                        elel_rbk*p_elel_bk + (1.-elel_rbk) *(p_elel_sym + AFB*p_elel_asym);


                    Double_t df_dAFB = ((1-mumu_rbk)*p_mumu_asym + (1-elel_rbk)*p_elel_asym);
                    Double_t res = (g/f)*df_dAFB;
                    a(0)+= res;
                    m(0,0) += (g/(f*f)) *pow(df_dAFB,2);
                    //printf("g=%.2e f=%.2e a=%.2e m=%.2e \n", g,f,a,m);
                    Double_t df_delel = (p_elel_bk -(p_elel_sym + AFB*p_elel_asym));
                    Double_t df_dmumu = (p_mumu_bk -(p_mumu_sym + AFB*p_mumu_asym));
                    a(1) += (g/f) * df_delel;
                    m(0,1) += (g/(f*f)) *df_dAFB*df_delel - (g/f) *p_elel_asym;
                    m(1,0) += (g/(f*f)) *df_dAFB*df_delel - (g/f) *p_elel_asym;
                    m(1,2) += (g/(f*f)) *df_dmumu*df_delel;
                    m(2,1) += (g/(f*f)) *df_dmumu*df_delel;
                    a(2) += (g/f) * df_dmumu;

                    m(0,2) += (g/(f*f)) *df_dAFB*df_dmumu - (g/f) *p_mumu_asym;
                    m(2,0) += (g/(f*f)) *df_dAFB*df_dmumu - (g/f) *p_mumu_asym;


                }
            }
            m.Invert();
            //Double_t delta_AFB = (a(0)/m(0,0));
            Double_t delta_AFB = m(0,0)*a(0) + m(0,1)*a(1) + m(0,2)*a(2);
            AFB_shifts[k][l] = delta_AFB;
        }
    }
    for(int k=0; k<6; k++){
        printf("\nFor M_low %.0f, shifts are: ", m_bins[k]);
        for(int l=0; l<=10; l++){
            printf("%.3f  ", AFB_shifts[k][l]);
        }


    }

}

