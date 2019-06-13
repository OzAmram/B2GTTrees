
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
#include "TRatioPlot.h"
#include "TH1F.h"
#include "TProfile.h"
#include "TStyle.h"
#include "TCanvas.h"
#include "TPostScript.h"
#include "TLorentzVector.h"
#include "TVector3.h"
#include "TFitter.h"
#include "TSystem.h"
#include "fit_gen_cost.C"
#include "../../Utils/PlotUtils.C"




void draw_generator_cmp(){
    gStyle->SetOptStat(0);
    TFile *f_mad= TFile::Open("../generator_stuff/root_files/madgraph_m700_may23.root");
    TTree *t_mad = (TTree *)f_mad->Get("T_lhe");

    TFile *f_pwg = TFile::Open("../generator_stuff/root_files/powheg_m700_april30.root");
    TTree *t_pwg = (TTree *)f_pwg->Get("T_lhe");

    float m_low = 700;
    float m_high = 800.;
    bool phot_ind = false;

    char title[80];
    sprintf(title, "POWHEG vs. aMC@NLO (%.0f < M < %.0f)", m_low, m_high);
    TH1F *h_pwg_cost_st = new TH1F("h_pwg_cost_st", title, 20, -1., 1.);
    TH1F *h_pwg_cost_r = new TH1F("h_pwg_cost_r", title, 20, -1., 1.);
    TH1F *h_pwg_pt = new TH1F("h_pwg_pt", title, 20, 0., 300.);
    TH1F *h_pwg_xf = new TH1F("h_pwg_xf", title, 20, 0., 0.5);

    TH1F *h_mad_cost_st = new TH1F("h_mad_cost_st", title, 20, -1., 1.);
    TH1F *h_mad_cost_r = new TH1F("h_mad_cost_r", title, 20, -1., 1.);
    TH1F *h_mad_pt = new TH1F("h_mad_pt", title, 20, 0., 300.);
    TH1F *h_mad_xf = new TH1F("h_mad_xf", title, 20, 0.,.5);


    make_gen_cost(t_mad,  h_mad_cost_st, h_mad_cost_r, h_mad_pt, h_mad_xf, m_low, m_high, phot_ind);
    make_gen_cost(t_pwg,  h_pwg_cost_st, h_pwg_cost_r, h_pwg_pt, h_pwg_xf, m_low, m_high, phot_ind);

    h_mad_cost_st->Scale(1./h_mad_cost_st->Integral());
    h_pwg_cost_st->Scale(1./h_pwg_cost_st->Integral());

    h_mad_xf->Scale(1./h_mad_xf->Integral());
    h_pwg_xf->Scale(1./h_pwg_xf->Integral());

    h_mad_cost_r->Scale(1./h_mad_cost_r->Integral());
    h_pwg_cost_r->Scale(1./h_pwg_cost_r->Integral());


    h_mad_pt->Scale(1./h_mad_pt->Integral());
    h_pwg_pt->Scale(1./h_pwg_pt->Integral());


    h_mad_cost_r->Scale(1./h_mad_cost_r->Integral());
    h_pwg_cost_r->Scale(1./h_pwg_cost_r->Integral());


    
    make_ratio_plot("pwg_vs_mad_m700_cost_st_cmp.pdf", h_mad_cost_st, "amc@NLO",h_pwg_cost_st, "POWHEG", "aMC/POWHEG", "cos(#theta_{*})", false);
    make_ratio_plot("pwg_vs_mad_m700_pt_cmp.pdf", h_mad_pt, "amc@NLO",h_pwg_pt, "POWHEG", "aMC/POWHEG", "dilepton p_{T} (GeV)", false);
    make_ratio_plot("pwg_vs_mad_m700_xf_cmp.pdf", h_mad_xf, "amc@NLO",h_pwg_xf, "POWHEG", "aMC/POWHEG", "xF", false);
    make_ratio_plot("pwg_vs_mad_m700_cost_r_cmp.pdf", h_mad_cost_r, "amc@NLO",h_pwg_cost_r, "POWHEG", "aMC/POWHEG", "cos(#theta_{r})", false);

    return;
}

