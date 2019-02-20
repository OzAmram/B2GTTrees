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
#include "TemplateUtils.h"

TH1F *h1_elel_dy, *h1_mumu_dy;
TH2F *h_elel_asym, *h_elel_sym, *h_elel_back,  *h_elel_dy_gg, *h_elel_data, *h_elel_mc, *h_elel_qcd;
TH1F *h1_elel_mn, *h1_elel_pl, *h1_elel_back,  *h1_elel_dy_gg, *h1_elel_data, *h1_elel_mc, *h1_elel_qcd;
TH2F *h_elel_mc_count, *h_elel_sym_count;
TH2F *h_mumu_asym, *h_mumu_sym, *h_mumu_back,  *h_mumu_dy_gg, *h_mumu_data, *h_mumu_mc, *h_mumu_qcd;
TH1F *h1_mumu_mn, *h1_mumu_pl, *h1_mumu_back,  *h1_mumu_dy_gg, *h1_mumu_data, *h1_mumu_mc, *h1_mumu_qcd;


void make_mc_templates(){
    
        h_mumu_sym = new TH2F((string("mumu_sym")  ).c_str(), "Symmetric template of mc",
                n_xf_bins, xf_bins, n_cost_bins, cost_bins);
        h_mumu_sym->SetDirectory(0);
        h_mumu_asym = new TH2F((string("mumu_asym")  ).c_str(), "Asymmetric template of mc",
                n_xf_bins, xf_bins, n_cost_bins, cost_bins);
        h_mumu_asym->SetDirectory(0);
        h_mumu_back = new TH2F((string("mumu_bk")  ).c_str(), "Combined background template",
                n_xf_bins, xf_bins, n_cost_bins, cost_bins);
        h_mumu_back->SetDirectory(0);

        printf("making muon mc \n");
        gen_mc_template(t_mumu_mc, alpha, h_mumu_sym, h_mumu_asym, m_low, m_high, FLAG_MUONS, FLAG_M_BINS, do_RC);
        TTree *mumu_ts[1] = {t_mumu_mc};
        gen_combined_background_template(1, mumu_ts, h_mumu_back, m_low, m_high, FLAG_MUONS, FLAG_M_BINS, do_RC);

        h_elel_sym = new TH2F((string("ee_sym")  ).c_str(), "Symmetric template of mc",
                n_xf_bins, xf_bins, n_cost_bins, cost_bins);
        h_elel_sym->SetDirectory(0);
        h_elel_asym = new TH2F((string("ee_asym")  ).c_str(), "Asymmetric template of mc",
                n_xf_bins, xf_bins, n_cost_bins, cost_bins);
        h_elel_asym->SetDirectory(0);
        h_elel_back = new TH2F((string("ee_bk")  ).c_str(), "Combined background template",
                n_xf_bins, xf_bins, n_cost_bins, cost_bins);
        h_elel_back->SetDirectory(0);

        printf("making elel mc \n");
        gen_mc_template(t_elel_mc, alpha, h_elel_sym, h_elel_asym,  m_low, m_high, FLAG_ELECTRONS, FLAG_M_BINS, do_RC);
        TTree *elel_ts[1] = {t_elel_mc};
        gen_combined_background_template(1, elel_ts, h_elel_back, m_low, m_high, FLAG_ELECTRONS, FLAG_M_BINS,do_RC);

        h1_mumu_back = convert2d(h_mumu_back);
        h1_elel_back = convert2d(h_elel_back);

        TH2F *h_elel_dy, *h_mumu_dy;
        float afb = 0.61;

        printf("trying to add");
        h_elel_dy = (TH2F *) h_elel_sym->Clone("h_elel_dy");
        h_mumu_dy = (TH2F *) h_mumu_sym->Clone("h_mumu_dy");
        h_elel_dy->Add(h_elel_sym, h_elel_asym, 1., afb);
        h_mumu_dy->Add(h_mumu_sym, h_mumu_asym, 1., afb);

        printf("trying to convert \n");
        h1_elel_dy = convert2d(h_elel_dy);
        h1_mumu_dy = convert2d(h_mumu_dy);
        h1_elel_dy->Scale(0.5);
        h1_mumu_dy->Scale(0.5);
        printf("ElEl Back is %.2f, mc is %.2f sym is %.2f \n", h1_elel_back->Integral(), h1_elel_dy->Integral(), h_elel_sym->Integral());
        printf("MuMu Back is %.2f, mc is %.2f sym is %.2f \n", h1_mumu_back->Integral(), h1_mumu_dy->Integral(), h_mumu_sym->Integral());

        /*
        auto h_mumu_pl = *h_mumu_sym + *h_mumu_asym;
        auto h_mumu_mn = *h_mumu_sym - *h_mumu_asym;
        h_mumu_pl.Scale(0.5);
        h_mumu_mn.Scale(0.5);
        h1_mumu_pl = convert2d(&h_mumu_pl);
        h1_mumu_mn = convert2d(&h_mumu_mn);
        auto h_elel_pl = *h_elel_sym + *h_elel_asym;
        auto h_elel_mn = *h_elel_sym - *h_elel_asym;
        h_elel_pl.Scale(0.5);
        h_elel_mn.Scale(0.5);
        h1_elel_pl = convert2d(&h_elel_pl);
        h1_elel_mn = convert2d(&h_elel_mn);
        */

}




void template_check(){
    init();
    printf("init done \n");

    int m_bin = 1;
    m_low = m_bins[m_bin];
    m_high = m_bins[m_bin+1];
    alpha = alphas[m_bin];
    setup_all_SFs();
    make_mc_templates();


    h1_elel_back->SetLineColor(kGreen +3);
    h1_elel_dy->SetLineColor(kRed+1);

    h1_mumu_back->SetLineColor(kGreen +3);
    h1_mumu_dy->SetLineColor(kRed+1);

    TCanvas *c_elel = new TCanvas("c_elel", "Histograms", 200, 10, 900, 700);
    h1_elel_dy->Draw("hist");
    h1_elel_back->Draw("hist same");

    TCanvas *c_mumu = new TCanvas("c_mumu", "Histograms", 200, 10, 900, 700);
    h1_mumu_dy->Draw("hist");
    h1_mumu_back->Draw("hist same");

    c_elel->Print("elel_template_check.pdf");
    c_mumu->Print("mumu_template_check.pdf");

}

