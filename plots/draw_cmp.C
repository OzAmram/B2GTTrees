
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
    TFile *f_data = TFile::Open("../analyze/output_files/DYToLL_data_2016_may9.root");
    TTree *t_data = (TTree *)f_data->Get("T_data");
    TFile *f_mc = TFile::Open("../analyze/output_files/DYToLL_mc_2016_may26.root");
    TTree *t_mc = (TTree *)f_mc->Get("T_data");
    TTree *t_mc_nosig = (TTree *)f_mc->Get("T_back");
    TFile *f_back = TFile::Open("../analyze/output_files/ttbar_background_may26.root");
    TTree *t_back = (TTree *)f_back->Get("T_data");

    TH1F *data_m = new TH1F("data_m", "Data Dimuon Mass Distribution", 30, 150, 1000);

    TH1F *mc_m = new TH1F("mc_m", "MC Signal (qqbar, qglu, qbarglu)", 30, 150, 1000);
    mc_m->SetFillColor(kRed);
    mc_m->SetMarkerStyle(21);
    mc_m->SetMarkerColor(kRed);
    TH1F *mc_nosig_m = new TH1F("mc_nosig_m", "MC no signal (qq, gluglu qbarqbar)", 30, 150, 1000);
    mc_nosig_m->SetFillColor(kYellow);
    mc_nosig_m->SetMarkerStyle(21);
    mc_nosig_m->SetMarkerColor(kYellow);
    TH1F *back_m = new TH1F("h_m", "TTBar Background", 30, 150, 1000);
    back_m->SetFillColor(kBlue);
    back_m->SetMarkerStyle(21);
    back_m->SetMarkerColor(kBlue);


    TH1F *data_cost = new TH1F("data_cost", "Data", 40, -1.,1.);
    TH1F *mc_cost = new TH1F("mc_cost", "MC Signal (qqbar, qglu, qbarglu)", 40, -1,1);
    mc_cost->SetFillColor(kRed);
    mc_cost->SetMarkerStyle(21);
    mc_cost->SetMarkerColor(kRed);
    TH1F *mc_nosig_cost = new TH1F("mc_nosig_cost", "MC no signal (qq, gluglu qbarqbar)", 40, -1.,1.);
    mc_nosig_cost->SetFillColor(kYellow);
    mc_nosig_cost->SetMarkerStyle(21);
    mc_nosig_cost->SetMarkerColor(kYellow);
    TH1F *back_cost = new TH1F("back_cost", "TTbar Background", 40, -1.,1.);
    back_cost->SetFillColor(kBlue);
    back_cost->SetMarkerStyle(21);
    back_cost->SetMarkerColor(kBlue);

    make_m_cost_hist(t_data, data_m, data_cost, true);
    make_m_cost_hist(t_mc, mc_m, mc_cost, false);
    make_m_cost_hist(t_mc_nosig, mc_nosig_m, mc_nosig_cost, false);
    make_m_cost_hist(t_back, back_m, back_cost, false);

    
    TCanvas *c1 = new TCanvas("c1", "Histograms", 200, 10, 900, 700);
    mc_m->Draw();
    

    //#lumi is fb^-1, convert to pb^-1
    /*
    float lumi = 35.867;

    mc_m->Scale(lumi*1000);
    mc_nosig_m->Scale(lumi*1000);
    mc_cost->Scale(lumi*1000);
   
    back_m->Scale(lumi*1000);
    back_cost->Scale(lumi*1000);
    */
    
    /*
    float scale1 = mc_cost->Integral();
    printf("DY has %f integral \n", scale1);
    float scale2 = mc_nosig_cost->Integral();
    float scale3 = back_cost->Integral();
    printf("TTbar has %f integral \n", scale3);

    float tot_scale = scale1 + scale2 + scale3;

    //mc_m->Scale(1./tot_scale);
    //mc_nosig_m->Scale(1./tot_scale);
    mc_cost->Scale(1./tot_scale);
    mc_nosig_cost->Scale(1./tot_scale);
    back_cost->Scale(1./tot_scale);
    */



    

    THStack *m_stack = new THStack("m_stack", "MuMu Mass Distribution: Data vs MC (All SF's Applied); (GeV)");
    m_stack->Add(back_m);
    m_stack->Add(mc_nosig_m);
    m_stack->Add(mc_m);


    THStack *cost_stack = new THStack("cost_stack", "Cos(#theta) Distribution: Data vs MC (All SF's Applied);Cos(#theta)");
    cost_stack->Add(back_cost);
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
    data_cost->Draw("P same");
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

    
    
