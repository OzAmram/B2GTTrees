

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
#include "../../analyze/HistMaker.C"
#include "../tdrstyle.C"
#include "../CMS_lumi.C"





void draw_emu_samesign(){
    setTDRStyle();
    TFile *f_data = TFile::Open("../analyze/output_files/EMu_data_samesign_jan17.root");
    TTree *t_data = (TTree *)f_data->Get("T_data");

                                
    TFile *f_ttbar = TFile::Open("../analyze/output_files/EMu_ttbar_wt_samesign_jan17.root");
    TTree *t_ttbar = (TTree *)f_ttbar->Get("T_data");

    TFile *f_DYToLL = TFile::Open("../analyze/output_files/EMu_DY_samesign_jan17.root");
    TTree *t_dy = (TTree *)f_DYToLL->Get("T_data");

    TFile *f_diboson = TFile::Open("../analyze/output_files/EMu_diboson_samesign_jan18.root");
    TTree *t_diboson = (TTree *)f_diboson->Get("T_data");



    TH1F *data_m = new TH1F("data_m", "MC Signal (qqbar, qglu, qbarglu)", 30, 150, 1000);
    TH1F *ttbar_m = new TH1F("ttbar_m", "MC Signal (qqbar, qglu, qbarglu)", 30, 150, 1000);
    TH1F *diboson_m = new TH1F("diboson_m", "MC Signal (qqbar, qglu, qbarglu)", 30, 150, 1000);
    TH1F *dy_m = new TH1F("dy_m", "MC Signal (qqbar, qglu, qbarglu)", 30, 150, 1000);

    int cost_nbins = 20;
    TH1F *dy_cost = new TH1F("mc_cost", "MC signal", cost_nbins, -1, 1);
    TH1F *data_cost = new TH1F("data_cost", "MC signal", cost_nbins, -1., 1.);
    TH1F *ttbar_cost = new TH1F("ttbar_cost", "MC signal", cost_nbins, -1., 1.);
    TH1F *diboson_cost = new TH1F("diboson_cost", "MC signal", cost_nbins, -1., 1.);
    dy_cost->SetFillColor(kRed+1);
    dy_cost->SetMarkerColor(kRed+1);
    ttbar_cost->SetFillColor(kBlue);
    ttbar_cost->SetMarkerStyle(21);
    ttbar_cost->SetMarkerColor(kBlue);
    diboson_cost->SetFillColor(kGreen+3);

    int xf_nbins = 16;
    TH1F *dy_xf = new TH1F("mc_xf", "MC signal", xf_nbins, 0, 0.8);
    TH1F *data_xf = new TH1F("data_xf", "MC signal", xf_nbins, 0, 0.8);
    TH1F *ttbar_xf = new TH1F("ttbar_xf", "MC signal", xf_nbins, 0, 0.8);
    TH1F *diboson_xf = new TH1F("diboson_xf", "MC signal", xf_nbins, 0, 0.8);
    dy_xf->SetFillColor(kRed+1);
    dy_xf->SetMarkerColor(kRed+1);
    ttbar_xf->SetFillColor(kBlue);
    ttbar_xf->SetMarkerStyle(21);
    ttbar_xf->SetMarkerColor(kBlue);
    diboson_xf->SetFillColor(kGreen+3);

    Double_t m_low = 150;
    Double_t m_high = 10000;

    int type  = FLAG_MUONS;

    make_emu_m_cost_xf_hist(t_data, data_m, data_cost, data_xf, true, type);
    make_emu_m_cost_xf_hist(t_ttbar, ttbar_m, ttbar_cost, ttbar_xf, false, type);
    make_emu_m_cost_xf_hist(t_diboson, diboson_m, diboson_cost, diboson_xf, false, type);
    make_emu_m_cost_xf_hist(t_dy, dy_m, dy_cost, dy_xf, false, type);

    //correct for wrong ttbar xsec
    //ttbar_m->Scale(831.76/730.6);

    Double_t data_count = data_m->Integral();
    Double_t mc_count = ttbar_m->Integral() + diboson_m->Integral() + dy_m->Integral();

    Double_t mc_unc = sqrt(mc_count);

    printf("Data count %.0f \n", data_count);
    printf("MC count %.0f \n", mc_count);
    Double_t ratio = mc_count / data_count ;
    printf("MC is %1.3f of observed events \n", ratio);



    dy_m->SetFillColor(kRed+1);
    ttbar_m->SetFillColor(kBlue);
    diboson_m->SetFillColor(kGreen+3);

    THStack *m_stack = new THStack("m_stack", "Samesign EMu Mass Distribution: Data vs MC ; m_{e#mu} (GeV)");
    m_stack->Add(dy_m);
    m_stack->Add(diboson_m);
    m_stack->Add(ttbar_m);

    THStack *cost_stack = new THStack("cost_stack", "Samesign EMu cos Distribution (symmeterized): Data vs MC ; cos(#theta_r)");
    cost_stack->Add(dy_cost);
    cost_stack->Add(diboson_cost);
    cost_stack->Add(ttbar_cost);

    THStack *xf_stack = new THStack("xf_stack", "Samesign EMu x_F Distribution (symmeterized): Data vs MC ; x_F");
    xf_stack->Add(dy_xf);
    xf_stack->Add(diboson_xf);
    xf_stack->Add(ttbar_xf);



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
    leg1->AddEntry(ttbar_m, "t#bar{t} + Wt", "f");
    leg1->AddEntry(diboson_m, "WW + WZ + ZZ", "f");
    leg1->AddEntry(dy_m, "DY #rightarrow #tau#tau", "f");
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
    auto h_ratio = (TH1F *) data_m->Clone("h_ratio");
    h_ratio->SetMinimum(0.);
    h_ratio->SetMaximum(10.);
    h_ratio->Sumw2();
    h_ratio->SetStats(0);
    h_ratio->Divide(m_mc_sum);
    h_ratio->SetMarkerStyle(21);
    h_ratio->Draw("ep");
    TLine *l1 = new TLine(150,1,1000,1);
    l1->SetLineStyle(2);
    l1->Draw();
    c_m->cd();

    h_ratio->SetTitle("");
    // Y axis h_ratio plot settings
   h_ratio->GetYaxis()->SetTitle("Data/MC");
   h_ratio->GetYaxis()->SetNdivisions(505);
   h_ratio->GetYaxis()->SetTitleSize(20);
   h_ratio->GetYaxis()->SetTitleFont(43);
   h_ratio->GetYaxis()->SetTitleOffset(1.2);
   h_ratio->GetYaxis()->SetLabelFont(43); // Absolute font size in pixel (precision 3)
   h_ratio->GetYaxis()->SetLabelSize(15);
   // X axis h_ratio plot settings
   h_ratio->GetXaxis()->SetTitle("M_{e#mu} (GeV)");
   h_ratio->GetXaxis()->SetTitleSize(20);
   h_ratio->GetXaxis()->SetTitleFont(43);
   h_ratio->GetXaxis()->SetTitleOffset(3.);
   h_ratio->GetXaxis()->SetLabelFont(43); // Absolute font size in pixel (precision 3)
   h_ratio->GetXaxis()->SetLabelSize(20);
 
    writeExtraText = true;
    extraText = "Preliminary";
    //lumi_sqrtS = "";       // used with iPeriod = 0, e.g. for simulation-only plots (default is an empty string)
    int iPeriod = 4; 
    CMS_lumi(pad1, iPeriod, 11 );




    TCanvas *c_cost = new TCanvas("c_cost", "Histograms", 200, 10, 900, 700);
    TPad *cost_pad1 = new TPad("pad1c", "pad1", 0.,0.3,0.98,1.);
    cost_pad1->SetBottomMargin(0);
    cost_pad1->SetLogy();
    cost_pad1->Draw();
    cost_pad1->cd();
    cost_stack->Draw("hist");
    data_cost->SetMarkerStyle(kFullCircle);
    data_cost->SetMarkerColor(1);
    cost_stack->SetMinimum(1);
    cost_stack->SetMaximum(1100);
    data_cost->Draw("P E same");
    cost_pad1->Update();
    leg1->Draw();
    c_cost->Update();

    c_cost->cd();
    TPad *cost_pad2 = new TPad("cost_pad2", "pad2", 0.,0,.98,0.3);
    //pad2->SetTopMargin(0);
    cost_pad2->SetBottomMargin(0.2);
    cost_pad2->SetGridy();
    cost_pad2->Draw();
    cost_pad2->cd();
    TList *cost_stackHists = cost_stack->GetHists();
    TH1* cost_mc_sum = (TH1*)cost_stackHists->At(0)->Clone();
    cost_mc_sum->Reset();

    for (int i=0;i<cost_stackHists->GetSize();++i) {
      cost_mc_sum->Add((TH1*)cost_stackHists->At(i));
    }
    auto cost_ratio = (TH1F *) data_cost->Clone("h_cost_ratio");
    cost_ratio->SetMinimum(0.);
    cost_ratio->SetMaximum(10.);
    cost_ratio->Sumw2();
    cost_ratio->SetStats(0);
    cost_ratio->Divide(cost_mc_sum);
    cost_ratio->SetMarkerStyle(21);
    cost_ratio->Draw("ep");
    TLine *l2 = new TLine(0,1,2000,1);
    l2->SetLineStyle(2);
    l2->Draw();
    c_cost->cd();

    cost_ratio->SetTitle("");
    // Y axis cost_ratio plot settings
   cost_ratio->GetYaxis()->SetTitle("Data/MC");
   cost_ratio->GetYaxis()->SetNdivisions(505);
   cost_ratio->GetYaxis()->SetTitleSize(20);
   cost_ratio->GetYaxis()->SetTitleFont(43);
   cost_ratio->GetYaxis()->SetTitleOffset(1.2);
   cost_ratio->GetYaxis()->SetLabelFont(43); // Absolute font size in pixel (precision 3)
   cost_ratio->GetYaxis()->SetLabelSize(15);
   // X axis cost_ratio plot settings
   cost_ratio->GetXaxis()->SetTitle("dimuon Cos(#theta)_{r} (GeV)");
   cost_ratio->GetXaxis()->SetTitleSize(20);
   cost_ratio->GetXaxis()->SetTitleFont(43);
   cost_ratio->GetXaxis()->SetTitleOffset(3.);
   cost_ratio->GetXaxis()->SetLabelFont(43); // Absolute font size in pixel (precision 3)
   cost_ratio->GetXaxis()->SetLabelSize(20);

    CMS_lumi(cost_pad1, iPeriod, 11);



    TCanvas *c_xf = new TCanvas("c_xf", "Histograms", 200, 10, 900, 700);
    TPad *xf_pad1 = new TPad("pad1", "pad1", 0.,0.3,0.98,1.);
    xf_pad1->SetBottomMargin(0);
    xf_pad1->Draw();
    xf_pad1->cd();
    xf_stack->Draw("hist");
    data_xf->SetMarkerStyle(kFullCircle);
    data_xf->SetMarkerColor(1);
    xf_stack->SetMinimum(1);
    xf_stack->SetMaximum(10000);
    data_xf->SetMinimum(1);
    data_xf->SetMaximum(10000);
    data_xf->Draw("P E same");
    xf_pad1->SetLogy();
    c_xf->Update();
    leg1->Draw();

    c_xf->cd();
    TPad *xf_pad2 = new TPad("xf_pad2", "pad2", 0.,0,.98,0.3);
    //pad2->SetTopMargin(0);
    xf_pad2->SetBottomMargin(0.2);
    xf_pad2->SetGridy();
    xf_pad2->Draw();
    xf_pad2->cd();
    TList *xf_stackHists = xf_stack->GetHists();
    TH1* xf_mc_sum = (TH1*)xf_stackHists->At(0)->Clone();
    xf_mc_sum->Reset();

    for (int i=0;i<xf_stackHists->GetSize();++i) {
      xf_mc_sum->Add((TH1*)xf_stackHists->At(i));
    }
    auto xf_ratio = (TH1F *) data_xf->Clone("h_xf_ratio");
    xf_ratio->SetMinimum(0.);
    xf_ratio->SetMaximum(10.);
    xf_ratio->Sumw2();
    xf_ratio->SetStats(0);
    xf_ratio->Divide(xf_mc_sum);
    xf_ratio->SetMarkerStyle(21);
    xf_ratio->Draw("ep");
    c_xf->cd();

    xf_ratio->SetTitle("");
    // Y axis xf_ratio plot settings
   xf_ratio->GetYaxis()->SetTitle("Data/MC");
   xf_ratio->GetYaxis()->SetNdivisions(505);
   xf_ratio->GetYaxis()->SetTitleSize(20);
   xf_ratio->GetYaxis()->SetTitleFont(43);
   xf_ratio->GetYaxis()->SetTitleOffset(1.2);
   xf_ratio->GetYaxis()->SetLabelFont(43); // Absolute font size in pixel (precision 3)
   xf_ratio->GetYaxis()->SetLabelSize(15);
   // X axis xf_ratio plot settings
   xf_ratio->GetXaxis()->SetTitle("dimuon xf (GeV)");
   xf_ratio->GetXaxis()->SetTitleSize(20);
   xf_ratio->GetXaxis()->SetTitleFont(43);
   xf_ratio->GetXaxis()->SetTitleOffset(3.);
   xf_ratio->GetXaxis()->SetLabelFont(43); // Absolute font size in pixel (precision 3)
   xf_ratio->GetXaxis()->SetLabelSize(20);
    CMS_lumi(xf_pad1, iPeriod, 11 );
    c_xf->Update();
}








