
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
#include "../analyze/TemplateMaker.C"





void draw_cmp(){
    TFile *f_data = TFile::Open("../analyze/output_files/DYToLL_data_2016_jun07.root");
    TTree *t_data = (TTree *)f_data->Get("T_data");
    TFile *f_mc = TFile::Open("../analyze/output_files/DYToLL_mc_2016_jun13.root");
    TTree *t_mc = (TTree *)f_mc->Get("T_data");
    TTree *t_mc_nosig = (TTree *)f_mc->Get("T_back");
    TFile *f_ttbar = TFile::Open("../analyze/output_files/ttbar_background_jun05.root");
    TTree *t_ttbar = (TTree *)f_ttbar->Get("T_data");

    TFile *f_ww = TFile::Open("../analyze/output_files/WW_background_jun06.root");
    TTree *t_ww = (TTree *)f_ww->Get("T_data");
    TFile *f_wz = TFile::Open("../analyze/output_files/WZ_background_jun06.root");
    TTree *t_wz = (TTree *)f_wz->Get("T_data");
    TFile *f_zz = TFile::Open("../analyze/output_files/ZZ_background_jun06.root");
    TTree *t_zz = (TTree *)f_zz->Get("T_data");
    TH1F *data_m = new TH1F("data_m", "Data Dimuon Mass Distribution", 30, 150, 1000);

    TH1F *mc_m = new TH1F("mc_m", "MC Signal (qqbar, qglu, qbarglu)", 30, 150, 1000);
    mc_m->SetFillColor(kRed);
    mc_m->SetMarkerStyle(21);
    mc_m->SetMarkerColor(kRed);
    TH1F *mc_nosig_m = new TH1F("mc_nosig_m", "MC no signal (qq, gluglu qbarqbar)", 30, 150, 1000);
    mc_nosig_m->SetFillColor(kYellow);
    mc_nosig_m->SetMarkerStyle(21);
    mc_nosig_m->SetMarkerColor(kYellow);
    TH1F *ttbar_m = new TH1F("ttbar_m", "TTBar Background", 30, 150, 1000);
    ttbar_m->SetFillColor(kBlue);
    ttbar_m->SetMarkerStyle(21);
    ttbar_m->SetMarkerColor(kBlue);


    TH1F *data_cost = new TH1F("data_cost", "Data", 40, -1.,1.);
    TH1F *mc_cost = new TH1F("mc_cost", "MC Signal (qqbar, qglu, qbarglu)", 40, -1,1);
    mc_cost->SetFillColor(kRed);
    mc_cost->SetMarkerStyle(21);
    mc_cost->SetMarkerColor(kRed);
    TH1F *mc_nosig_cost = new TH1F("mc_nosig_cost", "MC no signal (qq, gluglu qbarqbar)", 40, -1.,1.);
    mc_nosig_cost->SetFillColor(kYellow);
    mc_nosig_cost->SetMarkerStyle(21);
    mc_nosig_cost->SetMarkerColor(kYellow);
    TH1F *ttbar_cost = new TH1F("ttbar_cost", "TTbar Background", 40, -1.,1.);
    ttbar_cost->SetFillColor(kBlue);
    ttbar_cost->SetMarkerStyle(21);
    ttbar_cost->SetMarkerColor(kBlue);



    TH1F *ww_m = new TH1F("ww_m", "MC Signal (qqbar, qglu, qbarglu)", 30, 150, 1000);
    TH1F *ww_cost = new TH1F("ww_cost", "MC Signal (qqbar, qglu, qbarglu)", 40, -1,1);
    TH1F *wz_m = new TH1F("wz_m", "MC Signal (qqbar, qglu, qbarglu)", 30, 150, 1000);
    TH1F *wz_cost = new TH1F("wz_cost", "MC Signal (qqbar, qglu, qbarglu)", 40, -1,1);
    TH1F *zz_m = new TH1F("zz_m", "MC Signal (qqbar, qglu, qbarglu)", 30, 150, 1000);
    TH1F *zz_cost = new TH1F("zz_cost", "MC Signal (qqbar, qglu, qbarglu)", 40, -1,1);

    TH1F *diboson_m = new TH1F("diboson_m", "DiBoson (WW, WZ, ZZ)", 30, 150, 1000);
    TH1F *diboson_cost = new TH1F("diboson_cost", "DiBoson (WW, WZ,ZZ)", 40, -1,1);

    make_m_cost_hist(t_data, data_m, data_cost, true);
    make_m_cost_hist(t_mc, mc_m, mc_cost, false);
    make_m_cost_hist(t_mc_nosig, mc_nosig_m, mc_nosig_cost, false);
    make_m_cost_hist(t_ttbar, ttbar_m, ttbar_cost, false);
    ttbar_m->Scale(1.24);
    ttbar_cost->Scale(1.24);

    make_m_cost_hist(t_ww, ww_m, ww_cost, false);
    make_m_cost_hist(t_wz, wz_m, wz_cost, false);
    make_m_cost_hist(t_zz, zz_m, zz_cost, false);

    diboson_m->Add(ww_m);
    diboson_m->Add(wz_m);
    diboson_m->Add(zz_m);

    diboson_cost->Add(ww_cost);
    diboson_cost->Add(wz_cost);
    diboson_cost->Add(zz_cost);

    diboson_m->SetFillColor(kGreen);
    diboson_cost->SetFillColor(kGreen);
    
    TCanvas *c1 = new TCanvas("c1", "Histograms", 200, 10, 900, 700);
    //mc_m->Draw();
    

    //#lumi is fb^-1, convert to pb^-1
    /*
    float lumi = 35.867;

    mc_m->Scale(lumi*1000);
    mc_nosig_m->Scale(lumi*1000);
    mc_cost->Scale(lumi*1000);
   
    ttbar_m->Scale(lumi*1000);
    ttbar_cost->Scale(lumi*1000);
    */
    
    /*
    float scale1 = mc_cost->Integral();
    printf("DY has %f integral \n", scale1);
    float scale2 = mc_nosig_cost->Integral();
    float scale3 = ttbar_cost->Integral();
    printf("TTbar has %f integral \n", scale3);

    float tot_scale = scale1 + scale2 + scale3;

    //mc_m->Scale(1./tot_scale);
    //mc_nosig_m->Scale(1./tot_scale);
    mc_cost->Scale(1./tot_scale);
    mc_nosig_cost->Scale(1./tot_scale);
    ttbar_cost->Scale(1./tot_scale);
    */



    

    THStack *m_stack = new THStack("m_stack", "MuMu Mass Distribution: Data vs MC (All SF's Applied, EMu weighted ttbar); (GeV)");
    m_stack->Add(ttbar_m);
    m_stack->Add(diboson_m);
    m_stack->Add(mc_nosig_m);
    m_stack->Add(mc_m);


    THStack *cost_stack = new THStack("cost_stack", "Cos(#theta) Distribution: Data vs MC (All SF's Applied, EMu weighted ttbar);Cos(#theta)");
    cost_stack->Add(ttbar_cost);
    cost_stack->Add(diboson_cost);
    cost_stack->Add(mc_nosig_cost);
    cost_stack->Add(mc_cost);

    TCanvas *c_m = new TCanvas("c_m", "Histograms", 200, 10, 900, 700);
    //c_m->SetLogy();
    c_m->cd();
    m_stack->Draw("hist");
    gStyle->SetEndErrorSize(4);
    data_m->SetMarkerStyle(kFullCircle);
    data_m->SetMarkerColor(1);
    data_m->Draw("E1 same");
    c_m->Update();
    gPad->BuildLegend();


    TCanvas *c_cost = new TCanvas("c_cost", "Histograms", 200, 10, 900, 700);
    //c_cost->SetLogy();
    c_cost->cd();
    cost_stack->Draw("hist");
    data_cost->SetMarkerStyle(kFullCircle);
    data_cost->SetMarkerColor(1);
    //cost_stack->SetMaximum(3000);
    data_cost->Draw("P E same");
    c_cost->Update();
    gPad->BuildLegend();

    /*
    leg = new TLegend(0.1,0.6,0.5,0.9);
    leg->AddEntry(h_cost1, "No extra cuts", "f");
    leg->AddEntry(h_cost_cut, "met_pt < 80, mu1_pt > 40, mu2_pt > 20", "f");
    leg->Draw();
    */

    //TCanvas *c_cost_cut = new TCanvas("c_cost_cut", "Histograms", 200, 10, 900, 700);
    //c_cost_cut->cd();

/*
    TCanvas *c_met_pt = new TCanvas("c_met_pt", "Histograms", 200, 10, 900, 700);
    c_met_pt->cd();
    h_met_pt->SetFillColor(35);
    h_met_pt->Draw();
    c_met_pt->Update();
    TH1F *met_pt_int = h_integral(h_met_pt, "met_pt_int");


*/


 
}

    
    
