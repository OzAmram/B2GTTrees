
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
#include "TCanvas.h"
#include "TPostScript.h"
#include "TLorentzVector.h"
#include "TVector3.h"
#include "TFitter.h"
#include "TSystem.h"
#include "Math/Functor.h"
#include "../../analyze/HistMaker.C"
#include "../tdrstyle.C"
#include "../CMS_lumi.C"
#include "root_files.h"

void make_ratio_plot(char title[80], TH1F* h1, char h1_label[80], TH1F* h2, char h2_label[80], char ratio_label[80], 
        char axis_label[80], bool logy=false){

    h1->SetLineColor(kRed);
    h2->SetLineColor(kBlue);
    h1->SetLineWidth(3);
    h2->SetLineWidth(3);

    TCanvas *c = new TCanvas("c_m", "Histograms", 200, 10, 900, 700);
    TPad *pad1 = new TPad("pad1", "pad1", 0.,0.3,0.98,1.);
    pad1->SetBottomMargin(0);
    pad1->Draw();
    pad1->cd();
    if(logy) pad1->SetLogy();
    h2->Draw("hist E");
    //m_stack->SetMaximum(65000);
    h1->SetMinimum(0.1);
    gStyle->SetEndErrorSize(4);
    h1->Draw("hist E same");


    gStyle->SetLegendBorderSize(0);
    TLegend *leg1 = new TLegend(0.3, 0.3);
    leg1->AddEntry(h1, h1_label, "l");
    leg1->AddEntry(h2, h2_label, "l");
    leg1->Draw();

    //gPad->BuildLegend();
    c->cd();
    TPad *pad2 = new TPad("pad2", "pad2", 0.,0,.98,0.3);
    //pad2->SetTopMargin(0);
    pad2->SetBottomMargin(0.2);
    pad2->SetGridy();
    pad2->Draw();
    pad2->cd();

    auto ratio = (TH1F *) h1->Clone("h_ratio");
    ratio->SetMinimum(0.8);
    ratio->SetMaximum(1.2);
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
    int iPeriod = 4; 
    CMS_lumi(pad1, iPeriod, 33 );
    c->Print(title);
    return;
}



void draw_roch_cmp(){

    TFile *f_data2 = TFile::Open("../analyze/output_files/SingleMuon_data_jan22.root");
    TFile *f_data1 = TFile::Open("../analyze/output_files/SingleMuon_data_slim_june25.root");
    //t_data = (TTree *)f_data->Get("T_data");
    //TFile * f_mc2 = TFile::Open("../analyze/output_files/MuMu_DY_mar19.root");
    //TFile *f_mc1 = TFile::Open("../analyze/output_files/MuMu_DY_slim_june25.root");
    //TFile *f_mc2 = TFile::Open("../analyze/output_files/MuMu_DY_rcoff_slim_june27.root");
    //f_mc = TFile::Open("../analyze/output_files/MuMu_DY_april9_unbinned.root");
    //t_mc_nosig = (TTree *)f_mc->Get("T_back");
    TTree *t1 = (TTree *)f_data1->Get("T_data");
    TTree *t2 = (TTree *)f_data2->Get("T_data");
    setTDRStyle();




    TH1F *h1_cost = new TH1F("data_cost", "Data", 40, -1.,1.);
    TH1F *h1_m = new TH1F("data_m", "Data Dimuon Mass Distribution", 30, 140, 400);
    TH1F *h1_pt = new TH1F("mc_pt", "MC signal", 40, 0, 1000);

    TH1F *h2_cost = new TH1F("h2_cost", "Data", 40, -1.,1.);
    TH1F *h2_m = new TH1F("h2_m", "Data Dimuon Mass Distribution", 30, 140, 400);
    TH1F *h2_pt = new TH1F("h2_pt", "MC signal", 40, 0, 1000);

    make_m_cost_pt_hist(t1, h1_m, h1_cost, h1_pt, true, FLAG_MUONS);
    make_m_cost_pt_hist(t2, h2_m, h2_cost, h2_pt, true, FLAG_MUONS);

    printf("ON integral is %.2f. OFF integral is %.2f \n", h1_m->Integral(), h2_m->Integral());

    make_ratio_plot("MuMu_data_rocc_cost_cmp.pdf", h2_cost, "RC OFF ",h1_cost, "RC ON", "OFF/On", "Cos(#theta)", false);
    make_ratio_plot("MuMu_data_rocc_m_cmp.pdf", h2_m, "RC OFF ",h1_m, "RC ON", "OFF/On", "M (GeV)", false);
}

