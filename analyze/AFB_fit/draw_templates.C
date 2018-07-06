



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
#include "TRatioPlot.h"
#include "TProfile.h"
#include "TStyle.h"
#include "TCanvas.h"
#include "TPostScript.h"
#include "TLorentzVector.h"
#include "TVector3.h"
#include "TFitter.h"
#include "TSystem.h"
#include "Math/Functor.h"
#include "../TemplateMaker_systematics.C"
#include "FitUtils.C"



float m_low;
float m_high;
int FLAG = FLAG_ELECTRONS;

bool print = true;


Double_t med_btag = 0.4432;

TH2F *h_asym, *h_sym, *h_back,  *h_data, *h_mc;
TH2F *h_mc_count, *h_sym_count, *h_qcd;
vector<double> v_xF;
vector<double> v_cost;
unsigned int nDataEvents;


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
    h_qcd = new TH2F("h_qcd", "Combined background template",
            n_xf_bins, xf_bins, n_cost_bins, cost_bins);

    h_qcd->SetDirectory(0);
    h_data = new TH2F("h_data", "Data template of (x_f, cost_r) xF > 0.15",
            n_xf_bins, xf_bins, n_cost_bins, cost_bins);

    h_data->SetDirectory(0);
    printf("Generating templates \n");

    if(FLAG == FLAG_MUONS){
        bool do_RC = false;
        int flag2 = FLAG_M_BINS;
        nDataEvents = gen_data_template(t_mumu_data, h_data, &v_xF, &v_cost, m_low, m_high, FLAG, flag2, do_RC);
        gen_mc_template(t_mumu_mc, alpha, h_sym, h_asym, h_sym_count, m_low, m_high, FLAG, flag2, do_RC);
        TTree *ts[2] = {t_mumu_back, t_mumu_nosig};

        gen_fakes_template(t_mumu_WJets, t_mumu_QCD, t_mumu_WJets_contam, t_mumu_QCD_contam, h_qcd, m_low, m_high, FLAG, flag2);
        if(h_qcd->Integral() > 0.) h_qcd->Scale(1./h_qcd->Integral());
        else h_qcd->Scale(0.);
        gen_combined_background_template(2, ts, h_back, m_low, m_high, FLAG, flag2, do_RC);

    }
    if(FLAG == FLAG_ELECTRONS){
        gen_mc_template(t_elel_mc, alpha, h_sym, h_asym, h_sym_count, m_low, m_high, FLAG);
        TTree *ts[2] = {t_elel_back, t_elel_nosig};

        gen_fakes_template(t_elel_WJets, t_elel_QCD, t_elel_WJets_contam, t_elel_QCD_contam, h_qcd, m_low, m_high, FLAG);
        if(h_qcd->Integral() > 0.) h_qcd->Scale(1./h_qcd->Integral());
        else h_qcd->Scale(0.);
        gen_combined_background_template(2, ts, h_back, m_low, m_high, FLAG);

        nDataEvents = gen_data_template(t_elel_data, h_data, &v_xF, &v_cost, m_low, m_high, FLAG);
    }
    /*
    printf("\n\n\n Printing MC counts in each bin:\n");
    for(int i=1; i<=n_xf_bins; i++){
        for(int j=1; j<=n_cost_bins; j++){
            printf("%.0f ", h_sym_count->GetBinContent(i,j));
        }
        printf("\n");
    }
    */

    //f_data->Close();
    //f_back->Close();
    //f_mc->Close();
    printf("Finishing setup \n");
    return;
}

void draw_template(char title[80], TH2F *h){
    TH1F *h_1d[n_xf_bins];
    for(int i=1; i<=n_xf_bins; i++){
        char title[10];
        sprintf(title, "h_%i", i);
        h_1d[i-1] = new TH1F(title, "", n_cost_bins, cost_bins);
        h_1d[i-1]->SetStats(0);
        for(int j=1; j<=n_cost_bins; j++){
            h_1d[i-1]->SetBinContent(j, h->GetBinContent(i,j));
        }
    }
    TCanvas *c1 = new TCanvas("c1", "", 200, 10, 900, 700);
    h_1d[0]->SetLineColor(kBlack);
    h_1d[0]->SetMaximum(0.15);
    h_1d[0]->SetMinimum(0.);
    h_1d[1]->SetLineColor(kBlue);
    h_1d[2]->SetLineColor(kRed);
    h_1d[3]->SetLineColor(kGreen);
    h_1d[4]->SetLineColor(kOrange);
    TLegend *leg = new TLegend(0.5, 0.7, 0.65, 0.85);
    for(int i=0; i<n_xf_bins; i++){

        h_1d[i]->SetLineWidth(3);
        h_1d[i]->Draw("same");
        char title[10];
        sprintf(title, "xf bin %i", i);
        leg->AddEntry(h_1d[i], title, "l");
    }
    leg->Draw();
    c1->Print(title);
}

void draw_fit(char title_base[80], TH2F *h_data, TH2F *h_sym, TH2F *h_asym, TH2F *h_back, TH2F *h_qcd, double AFB, double r_back, double r_qcd){
    TH1F *h_1d[n_xf_bins];
    TH1F *h_fits[n_xf_bins];
    for(int i=1; i<=n_xf_bins; i++){
        char title[10], title2[10];
        sprintf(title, "h_%i", i);
        sprintf(title2, "fit_%i", i);
        h_1d[i-1] = new TH1F(title, "", n_cost_bins, cost_bins);
        h_fits[i-1] = new TH1F(title2, "", n_cost_bins, cost_bins);
        for(int j=1; j<=n_cost_bins; j++){
            float sym = h_sym->GetBinContent(i,j);
            float asym = h_asym->GetBinContent(i,j);
            float back = h_back->GetBinContent(i,j);
            float qcd = h_qcd->GetBinContent(i,j);
    
            float fit = r_back*back + r_qcd*qcd + (1 - r_back - r_qcd) * (sym + AFB*asym);
            h_fits[i-1]->SetBinContent(j, fit);
            h_fits[i-1]->SetBinError(j, 0);
            h_1d[i-1]->SetBinContent(j, h_data->GetBinContent(i,j));
            h_1d[i-1]->SetBinError(j, h_data->GetBinError(i,j));
        }
        TCanvas *c1 = new TCanvas("c1", "", 200, 10, 900, 700);
        h_1d[i-1]->SetLineColor(kRed);
        h_1d[i-1]->SetStats(0);
        h_fits[i-1]->SetLineColor(kBlue);
        h_fits[i-1]->SetStats(0);
        auto rp = new TRatioPlot(h_1d[i-1], h_fits[i-1], "diffsig");
        rp->SetH1DrawOpt("P E");
        rp->Draw();
        c1->Update();
        char out_title[80];
        sprintf(out_title, "%s_%i.pdf", title_base, i);
        c1->Print(out_title);
    }
}


void draw_templates(){
    init();
    for(int i=1; i<2; i++){
        printf("Starting loop \n");
        m_low = m_bins[i];
        m_high = m_bins[i+1];
        alpha = alphas[i];

        setup();
        printf("\n\n\n Printing QCD template in each bin:\n");
        print_hist(h_qcd);
        printf("\n\n\n Print back template \n");
        print_hist(h_back);

        printf("Integrals are %f %f %f %f %f  \n", h_data->Integral(), h_sym->Integral(), 
                                               h_asym->Integral(), h_back->Integral(), h_qcd->Integral() );
        draw_template("ElEl_QCD_template.pdf", h_qcd);
        draw_template("ElEl_back_template.pdf", h_back);
        draw_fit("ElEl_fit_qcd_sep", h_data, h_sym, h_asym, h_back, h_qcd, 0.692, 0.29, 0.008);
    }
}
