
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
#include "Math/Functor.h"
#include "../../analyze/TemplateMaker.C"
#include "../tdrstyle.C"
#include "../CMS_lumi.C"

const int type = FLAG_MUONS;



void draw_cmp(){
    setTDRStyle();
    TFile *f_data = TFile::Open("../analyze/output_files/SingleMuon_data_aug28.root");
    TTree *t_data = (TTree *)f_data->Get("T_data");
    TFile *f_mc = TFile::Open("../analyze/output_files/MuMu_DY_aug30.root");
    TTree *t_mc = (TTree *)f_mc->Get("T_data");
    TTree *t_mc_nosig = (TTree *)f_mc->Get("T_back");
    TFile *f_ttbar = TFile::Open("../analyze/output_files/MuMu_TTbar_aug30.root");
    TTree *t_ttbar = (TTree *)f_ttbar->Get("T_data");

    TFile *f_QCD = TFile::Open("../analyze/output_files/MuMu_QCD_est_sep20.root");
    TTree *t_QCD = (TTree *)f_QCD->Get("T_data");

    TFile *f_WJets = TFile::Open("../analyze/output_files/MuMu_Wjets_est_sep20.root");
    TTree *t_WJets = (TTree *)f_WJets->Get("T_data");

    TFile *f_diboson = TFile::Open("../analyze/output_files/MuMu_diboson_aug30.root");
    TTree *t_diboson = (TTree *)f_diboson->Get("T_data");

    TFile *f_wt = TFile::Open("../analyze/output_files/MuMu_WT_aug30.root");
    TTree *t_wt = (TTree *)f_wt->Get("T_data");
    TH1F *data_m = new TH1F("data_m", "Data Dimuon Mass Distribution", 30, 150, 2000);

    TH1F *mc_m = new TH1F("mc_m", "MC Signal (qqbar, qglu, qbarglu)", 30, 150, 2000);
    mc_m->SetFillColor(kRed+1);
    mc_m->SetMarkerStyle(21);
    mc_m->SetMarkerColor(kRed+1);
    TH1F *mc_nosig_m = new TH1F("mc_nosig_m", "MC no signal (qq, gluglu qbarqbar)", 30, 150, 2000);
    mc_nosig_m->SetFillColor(kMagenta);
    mc_nosig_m->SetMarkerStyle(21);
    mc_nosig_m->SetMarkerColor(kMagenta);
    TH1F *ttbar_m = new TH1F("ttbar_m", "TTBar Background", 30, 150, 2000);
    ttbar_m->SetFillColor(kBlue);
    ttbar_m->SetMarkerStyle(21);
    ttbar_m->SetMarkerColor(kBlue);


    TH1F *data_cost = new TH1F("data_cost", "Data", 40, -1.,1.);
    TH1F *mc_cost = new TH1F("mc_cost", "MC Signal (qqbar, qglu, qbarglu)", 40, -1,1);
    mc_cost->SetFillColor(kRed+1);
    mc_cost->SetMarkerStyle(21);
    mc_cost->SetMarkerColor(kRed+1);
    TH1F *mc_nosig_cost = new TH1F("mc_nosig_cost", "MC no signal (qq, gluglu qbarqbar)", 40, -1.,1.);
    mc_nosig_cost->SetFillColor(kMagenta);
    mc_nosig_cost->SetMarkerStyle(21);
    mc_nosig_cost->SetMarkerColor(kMagenta);
    TH1F *ttbar_cost = new TH1F("ttbar_cost", "TTbar Background", 40, -1.,1.);
    ttbar_cost->SetFillColor(kBlue);
    ttbar_cost->SetMarkerStyle(21);
    ttbar_cost->SetMarkerColor(kBlue);



    /*
    TH1F *ww_m = new TH1F("ww_m", "MC Signal (qqbar, qglu, qbarglu)", 30, 150, 2000);
    TH1F *ww_cost = new TH1F("ww_cost", "MC Signal (qqbar, qglu, qbarglu)", 40, -1,1);
    TH1F *wz_m = new TH1F("wz_m", "MC Signal (qqbar, qglu, qbarglu)", 30, 150, 2000);
    TH1F *wz_cost = new TH1F("wz_cost", "MC Signal (qqbar, qglu, qbarglu)", 40, -1,1);
    TH1F *zz_m = new TH1F("zz_m", "MC Signal (qqbar, qglu, qbarglu)", 30, 150, 2000);
    TH1F *zz_cost = new TH1F("zz_cost", "MC Signal (qqbar, qglu, qbarglu)", 40, -1,1);
    TH1F *diboson_m = new TH1F("diboson_m", "DiBoson (WW, WZ, ZZ)", 30, 150, 2000);
    TH1F *diboson_cost = new TH1F("diboson_cost", "DiBoson (WW, WZ,ZZ)", 40, -1,1);
    make_m_cost_hist(t_ww, ww_m, ww_cost, false);
    make_m_cost_hist(t_wz, wz_m, wz_cost, false);
    make_m_cost_hist(t_zz, zz_m, zz_cost, false);
    diboson_m->Add(ww_m);
    diboson_m->Add(wz_m);
    diboson_m->Add(zz_m);

    diboson_cost->Add(ww_cost);
    diboson_cost->Add(wz_cost);
    diboson_cost->Add(zz_cost);
    */

    TH1F *diboson_m = new TH1F("diboson_m", "DiBoson (WW, WZ, ZZ)", 30, 150, 2000);
    TH1F *diboson_cost = new TH1F("diboson_cost", "DiBoson (WW, WZ,ZZ)", 40, -1,1);

    TH1F *QCD_m = new TH1F("QCD_m", "QCD", 30, 150, 2000);
    TH1F *QCD_cost = new TH1F("QCD_cost", "QCD", 40, -1,1);

    TH1F *WJets_m = new TH1F("WJets_m", "WJets", 30, 150, 2000);
    TH1F *WJets_cost = new TH1F("WJets_cost", "WJets", 40, -1,1);

    TH1F *wt_m = new TH1F("wt_m", "tw + #bar{t}w", 30, 150, 2000);
    TH1F *wt_cost = new TH1F("wt_cost", "tw + #bar{t}w", 40, -1,1);


    wt_m->SetFillColor(kOrange+7); 
    wt_cost->SetFillColor(kOrange+7); 

    make_m_cost_hist(t_data, data_m, data_cost, true);
    make_m_cost_hist(t_mc, mc_m, mc_cost, false);
    make_m_cost_hist(t_mc_nosig, mc_nosig_m, mc_nosig_cost, false);
    make_m_cost_hist(t_ttbar, ttbar_m, ttbar_cost, false);

    make_m_cost_hist(t_QCD, QCD_m, QCD_cost, true, FLAG_QCD);
    
    make_m_cost_hist(t_WJets, WJets_m, WJets_cost, true, FLAG_WJETS);


    make_m_cost_hist(t_diboson, diboson_m, diboson_cost, false, FLAG_MUONS);

    make_m_cost_hist(t_wt, wt_m, wt_cost, false);


    ttbar_m->Scale(1.05);
    ttbar_cost->Scale(1.05);
    diboson_m->Scale(1.05);
    diboson_cost->Scale(1.05);
    wt_m->Scale(1.05);
    wt_cost->Scale(1.05);

    diboson_m->SetFillColor(kGreen+3);
    diboson_cost->SetFillColor(kGreen + 3);
    
    QCD_m->SetFillColor(kRed -7);
    QCD_cost->SetFillColor(kRed -7);

    QCD_m->Add(WJets_m);
    QCD_cost->Add(WJets_cost);

    //mc_m->Draw();
    
    //TCanvas *c1 = new TCanvas("c1", "Histograms", 200, 10, 900, 700);

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



    

    THStack *m_stack = new THStack("m_stack", "MuMu Mass Distribution: Data vs MC ; m_{#mu^{+}#mu^{-}} (GeV)");
    m_stack->Add(ttbar_m);
    m_stack->Add(QCD_m);
    m_stack->Add(wt_m);
    m_stack->Add(diboson_m);
    m_stack->Add(mc_nosig_m);
    m_stack->Add(mc_m);


    THStack *cost_stack = new THStack("cost_stack", "Cos(#theta) Distribution: Data vs MC;Cos(#theta)_{r}");
    cost_stack->Add(ttbar_cost);
    cost_stack->Add(QCD_cost);
    cost_stack->Add(wt_cost);
    cost_stack->Add(diboson_cost);
    cost_stack->Add(mc_nosig_cost);
    cost_stack->Add(mc_cost);

    TCanvas *c_m = new TCanvas("c_m", "Histograms", 200, 10, 900, 700);
    TPad *pad1 = new TPad("pad1", "pad1", 0.,0.3,0.98,1.);
    pad1->SetBottomMargin(0);
    pad1->Draw();
    pad1->cd();
    pad1->SetLogy();
    m_stack->Draw("hist");
    m_stack->SetMaximum(65000);
    gStyle->SetEndErrorSize(4);
    data_m->SetMarkerStyle(kFullCircle);
    data_m->SetMarkerColor(1);
    data_m->DrawCopy("P E same");


    gStyle->SetLegendBorderSize(0);
    TLegend *leg1 = new TLegend(0.5, 0.65, 0.75, 0.8);
    leg1->AddEntry(data_m, "data", "p");
    leg1->AddEntry(mc_m, "DY (q#bar{q}, qg #bar{q}g)", "f");
    leg1->AddEntry(mc_nosig_m, "DY no asymmety(gg, qq, #bar{q}#bar{q})", "f");
    leg1->AddEntry(diboson_m, "WW + WZ + ZZ", "f");
    leg1->AddEntry(wt_m, "tW + #bar{t}W", "f");
    leg1->AddEntry(QCD_m, "QCD + WJets", "f");
    leg1->AddEntry(ttbar_m, "t#bar{t}", "f");
    leg1->Draw();

    //gPad->BuildLegend();
    c_m->cd();
    TPad *pad2 = new TPad("pad2", "pad2", 0.,0,.98,0.3);
    //pad2->SetTopMargin(0);
    pad2->SetBottomMargin(0.2);
    pad2->Draw();
    pad2->cd();



    TList *stackHists = m_stack->GetHists();
    TH1* m_mc_sum = (TH1*)stackHists->At(0)->Clone();
    m_mc_sum->Reset();

    for (int i=0;i<stackHists->GetSize();++i) {
      m_mc_sum->Add((TH1*)stackHists->At(i));
    }
    auto ratio = (TH1F *) data_m->Clone("h_ratio");
    ratio->SetMinimum(0.75);
    ratio->SetMaximum(1.25);
    ratio->Sumw2();
    ratio->SetStats(0);
    ratio->Divide(m_mc_sum);
    ratio->SetMarkerStyle(21);
    ratio->Draw("ep");
    TLine *l1 = new TLine(150,1,2000,1);
    l1->SetLineStyle(2);
    l1->Draw();
    c_m->cd();

    ratio->SetTitle("");
    // Y axis ratio plot settings
   ratio->GetYaxis()->SetTitle("Data/MC");
   ratio->GetYaxis()->SetNdivisions(505);
   ratio->GetYaxis()->SetTitleSize(20);
   ratio->GetYaxis()->SetTitleFont(43);
   ratio->GetYaxis()->SetTitleOffset(1.2);
   ratio->GetYaxis()->SetLabelFont(43); // Absolute font size in pixel (precision 3)
   ratio->GetYaxis()->SetLabelSize(15);
   // X axis ratio plot settings
   ratio->GetXaxis()->SetTitle("M_{#mu#mu} (GeV)");
   ratio->GetXaxis()->SetTitleSize(20);
   ratio->GetXaxis()->SetTitleFont(43);
   ratio->GetXaxis()->SetTitleOffset(3.);
   ratio->GetXaxis()->SetLabelFont(43); // Absolute font size in pixel (precision 3)
   ratio->GetXaxis()->SetLabelSize(20);
 
    writeExtraText = true;
    extraText = "Preliminary";
    //lumi_sqrtS = "";       // used with iPeriod = 0, e.g. for simulation-only plots (default is an empty string)
    int iPeriod = 4; 
    CMS_lumi(pad1, iPeriod, 33 );

    /*
    TCanvas *c_m2 = new TCanvas("c_m2", "Ratio test", 900, 700);
    auto rp = new TRatioPlot(data_m, m_mc_sum);
    rp->fHistDrawProxy = m_stack;
    c_m2->SetTicks(0,1);
    rp->Draw();
    c1->Update();
    */


    TCanvas *c_cost = new TCanvas("c_cost", "Histograms", 200, 10, 900, 700);
    //c_cost->SetLogy();
    c_cost->cd();
    cost_stack->Draw("hist");
    data_cost->SetMarkerStyle(kFullCircle);
    data_cost->SetMarkerColor(1);
    //cost_stack->SetMaximum(3000);
    data_cost->Draw("P E same");
    c_cost->Update();
    TLegend *leg2 = new TLegend(0.5, 0.65, 0.75, 0.8);
    leg2->AddEntry(data_m, "data", "p");
    leg2->AddEntry(mc_m, "DY (q#bar{q}, qg #bar{q}g)", "f");
    leg2->AddEntry(mc_nosig_m, "DY no asymmety(gg, qq, #bar{q}#bar{q})", "f");
    leg2->AddEntry(diboson_m, "WW + WZ + ZZ", "f");
    leg2->AddEntry(wt_m, "tW + #bar{t}W", "f");
    leg2->AddEntry(QCD_m, "QCD + WJets", "f");
    leg2->AddEntry(ttbar_m, "t#bar{t}", "f");
    leg2->Draw();

    CMS_lumi(c_cost, iPeriod, 11 );
    c_cost->Update();
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

    
    
