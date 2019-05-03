
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


void unzero_bins(TH1 *h){
    int nBins = h->GetNbinsX();
    for(int i=1; i<= nBins; i++){
        float val = max(h->GetBinContent(i), 1e-8);
        h->SetBinContent(i,val);
    }
}


void make_ratio_plot(string title, TH1F* h1, char h1_label[80], TH1F* h2, char h2_label[80], char ratio_label[80], 
        char axis_label[80], bool logy=false, bool write_out = true){

    unzero_bins(h1);
    unzero_bins(h2);


    TCanvas *c = new TCanvas(title.c_str(), "Histograms", 200, 10, 900, 700);
    TPad *pad1 = new TPad((title+"p1").c_str(), "pad1", 0.,0.3,0.98,1.);
    pad1->SetBottomMargin(0);
    pad1->Draw();
    pad1->cd();
    if(logy) pad1->SetLogy();
    //h1->SetMinimum(0.1);
    h1->Draw("hist E");
    gStyle->SetEndErrorSize(4);
    h2->Draw("hist E same");


    gStyle->SetLegendBorderSize(0);
    TLegend *leg1 = new TLegend(0.15, 0.7, 0.35, 0.85);
    leg1->AddEntry(h1, h1_label, "l");
    leg1->AddEntry(h2, h2_label, "l");
    leg1->Draw();

    //gPad->BuildLegend();
    c->cd();
    TPad *pad2 = new TPad((title+"p2").c_str(), "pad2", 0.,0,.98,0.3);
    //pad2->SetTopMargin(0);
    pad2->SetBottomMargin(0.2);
    pad2->SetGridy();
    pad2->Draw();
    pad2->cd();

    auto ratio = (TH1F *) h1->Clone("h_ratio");
    ratio->SetMinimum(0.01);
    ratio->SetMaximum(2.0);
    ratio->Sumw2();
    ratio->SetStats(0);
    ratio->Divide(h2);
    ratio->SetMarkerStyle(21);
    ratio->SetLineColor(kBlack);
    ratio->Draw("ep");
    c->cd();

    ratio->SetTitle("");
    // Y axis ratio plot settings
    ratio->GetYaxis()->SetTitle(ratio_label);
    ratio->GetYaxis()->SetNdivisions(505);
    ratio->GetYaxis()->SetTitleSize(20);
    ratio->GetYaxis()->SetTitleFont(43);
    ratio->GetYaxis()->SetTitleOffset(1.2);
    ratio->GetYaxis()->SetLabelFont(43); // Absolute font size in pixel (precision 3)
    ratio->GetYaxis()->SetLabelSize(15);
    // X axis ratio plot settings
    ratio->GetXaxis()->SetTitle(axis_label);
    ratio->GetXaxis()->SetTitleSize(20);
    ratio->GetXaxis()->SetTitleFont(43);
    ratio->GetXaxis()->SetTitleOffset(3.);
    ratio->GetXaxis()->SetLabelFont(43); // Absolute font size in pixel (precision 3)
    ratio->GetXaxis()->SetLabelSize(20);

    //lumi_sqrtS = "";       // used with iPeriod = 0, e.g. for simulation-only plots (default is an empty string)
    //int iPeriod = 4; 
    //CMS_lumi(pad1, iPeriod, 33 );
    if(write_out) c->Print(title.c_str());
    return;
}

void draw_generator_cmp(){
    gStyle->SetOptStat(0);
    TFile *f_mad= TFile::Open("../generator_stuff/root_files/madgraph_m200_evts.root");
    TTree *t_mad = (TTree *)f_mad->Get("T_lhe");

    TFile *f_pwg = TFile::Open("../generator_stuff/root_files/powheg_m200_april30.root");
    TTree *t_pwg = (TTree *)f_pwg->Get("T_lhe");

    char title[80] = "POWHEG vs. aMC@NLO (200 < M < 400)";
    TH1F *h_pwg_cost_st = new TH1F("h_pwg_cost_st", title, 20, -1., 1.);
    TH1F *h_pwg_cost_r = new TH1F("h_pwg_cost_r", title, 20, -1., 1.);
    TH1F *h_pwg_pt = new TH1F("h_pwg_pt", title, 20, 0., 300.);
    TH1F *h_pwg_xf = new TH1F("h_pwg_xf", title, 20, 0., 1.);

    TH1F *h_mad_cost_st = new TH1F("h_mad_cost_st", title, 20, -1., 1.);
    TH1F *h_mad_cost_r = new TH1F("h_mad_cost_r", title, 20, -1., 1.);
    TH1F *h_mad_pt = new TH1F("h_mad_pt", title, 20, 0., 300.);
    TH1F *h_mad_xf = new TH1F("h_mad_xf", title, 20, 0., 1.);

    float m_low = 200;
    float m_high = 400.;
    bool phot_ind = false;

    make_gen_cost(t_mad,  h_mad_cost_st, h_mad_cost_r, h_mad_pt, h_mad_xf, m_low, m_high, phot_ind);
    make_gen_cost(t_pwg,  h_pwg_cost_st, h_pwg_cost_r, h_pwg_pt, h_pwg_xf, m_low, m_high, phot_ind);

    h_mad_cost_st->Scale(1./h_mad_cost_st->Integral());
    h_pwg_cost_st->Scale(1./h_pwg_cost_st->Integral());

    h_mad_cost_st->SetLineColor(kRed);
    h_pwg_cost_st->SetLineColor(kBlue);

    h_mad_cost_st->SetLineWidth(3);
    h_pwg_cost_st->SetLineWidth(3);

    h_mad_pt->Scale(1./h_mad_pt->Integral());
    h_pwg_pt->Scale(1./h_pwg_pt->Integral());

    h_mad_pt->SetLineColor(kRed);
    h_pwg_pt->SetLineColor(kBlue);

    h_mad_pt->SetLineWidth(3);
    h_pwg_pt->SetLineWidth(3);

    h_mad_cost_r->Scale(1./h_mad_cost_r->Integral());
    h_pwg_cost_r->Scale(1./h_pwg_cost_r->Integral());

    h_mad_cost_r->SetLineColor(kRed);
    h_pwg_cost_r->SetLineColor(kBlue);

    h_mad_cost_r->SetLineWidth(3);
    h_pwg_cost_r->SetLineWidth(3);

    make_ratio_plot("pwg_vs_mad_m200_cost_st_cmp.pdf", h_mad_cost_st, "amc@NLO",h_pwg_cost_st, "POWHEG", "aMC/POWHEG", "cos(#theta_{st})", false);
    make_ratio_plot("pwg_vs_mad_m200_pt_cmp.pdf", h_mad_pt, "amc@NLO",h_pwg_pt, "POWHEG", "aMC/POWHEG", "dilepton p_{T} (GeV)", false);

    return;
}

