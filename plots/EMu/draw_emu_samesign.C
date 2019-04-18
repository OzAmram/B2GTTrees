

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
#include "../../Utils/HistMaker.C"
#include "../../Utils/root_files.h"
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
    TH1F *qcd_m = new TH1F("qcd_m", "MC Signal (qqbar, qglu, qbarglu)", 30, 150, 1000);
    dy_m->SetFillColor(kRed+1);
    ttbar_m->SetFillColor(kBlue);
    diboson_m->SetFillColor(kGreen+3);
    qcd_m->SetFillColor(kRed-7);

    int cost_nbins = 10;
    TH1F *dy_cost = new TH1F("mc_cost", "MC signal", cost_nbins, 0, 1);
    TH1F *data_cost = new TH1F("data_cost", "MC signal", cost_nbins, 0., 1.);
    TH1F *ttbar_cost = new TH1F("ttbar_cost", "MC signal", cost_nbins, 0., 1.);
    TH1F *diboson_cost = new TH1F("diboson_cost", "MC signal", cost_nbins, 0., 1.);
    TH1F *qcd_cost = new TH1F("qcd_cost", "MC signal", cost_nbins, 0., 1.);
    dy_cost->SetFillColor(kRed+1);
    dy_cost->SetMarkerColor(kRed+1);
    ttbar_cost->SetFillColor(kBlue);
    ttbar_cost->SetMarkerStyle(21);
    ttbar_cost->SetMarkerColor(kBlue);
    diboson_cost->SetFillColor(kGreen+3);
    qcd_cost->SetFillColor(kRed-7);

    int xf_nbins = 16;
    TH1F *dy_xf = new TH1F("mc_xf", "MC signal", xf_nbins, 0, 0.8);
    TH1F *data_xf = new TH1F("data_xf", "MC signal", xf_nbins, 0, 0.8);
    TH1F *ttbar_xf = new TH1F("ttbar_xf", "MC signal", xf_nbins, 0, 0.8);
    TH1F *diboson_xf = new TH1F("diboson_xf", "MC signal", xf_nbins, 0, 0.8);
    TH1F *qcd_xf = new TH1F("qcd_xf", "MC signal", xf_nbins, 0, 0.8);
    dy_xf->SetFillColor(kRed+1);
    dy_xf->SetMarkerColor(kRed+1);
    ttbar_xf->SetFillColor(kBlue);
    ttbar_xf->SetMarkerStyle(21);
    ttbar_xf->SetMarkerColor(kBlue);
    diboson_xf->SetFillColor(kGreen+3);
    qcd_xf->SetFillColor(kRed-7);

    int pt_nbins = 10;
    TH1F *dy_pt = new TH1F("mc_pt", "MC signal", xf_nbins, 0, 600);
    TH1F *data_pt = new TH1F("data_pt", "MC signal", xf_nbins, 0, 600);
    TH1F *ttbar_pt = new TH1F("ttbar_pt", "MC signal", xf_nbins, 0, 600);
    TH1F *diboson_pt = new TH1F("diboson_pt", "MC signal", xf_nbins, 0, 600);
    TH1F *qcd_pt = new TH1F("qcd_pt", "MC signal", xf_nbins, 0, 600);
    dy_pt->SetFillColor(kRed+1);
    dy_pt->SetMarkerColor(kRed+1);
    ttbar_pt->SetFillColor(kBlue);
    ttbar_pt->SetMarkerStyle(21);
    ttbar_pt->SetMarkerColor(kBlue);
    diboson_pt->SetFillColor(kGreen+3);
    qcd_pt->SetFillColor(kRed-7);

    Double_t m_low = 150;
    Double_t m_high = 10000;
    bool ss= true;

    int type  = FLAG_MUONS;

    make_emu_m_cost_pt_xf_hist(t_data, data_m, data_cost, data_pt, data_xf, true,  m_low, m_high, ss);
    make_emu_m_cost_pt_xf_hist(t_ttbar, ttbar_m, ttbar_cost, ttbar_pt, ttbar_xf, false,  m_low, m_high, ss);
    make_emu_m_cost_pt_xf_hist(t_diboson, diboson_m, diboson_cost, diboson_pt, diboson_xf, false,  m_low, m_high, ss);
    make_emu_m_cost_pt_xf_hist(t_dy, dy_m, dy_cost, dy_xf, dy_pt, false,  m_low, m_high, ss);

    make_qcd_from_emu_m_cost_pt_xf_hist(t_data, t_ttbar, t_diboson, t_dy, qcd_m, qcd_cost, qcd_pt, qcd_xf, m_low, m_high);

    //correct for wrong ttbar xsec
    //ttbar_m->Scale(831.76/730.6);

    Double_t data_count = data_m->Integral();
    Double_t mc_count = ttbar_m->Integral() + diboson_m->Integral() + dy_m->Integral();

    Double_t mc_unc = sqrt(mc_count);

    printf("Data count %.0f \n", data_count);
    printf("MC count %.0f \n", mc_count);
    Double_t ratio = mc_count / data_count ;
    printf("MC is %1.3f of observed events \n", ratio);




    THStack *m_stack = new THStack("m_stack", "Samesign EMu Mass Distribution: Data vs MC ; m_{e#mu} (GeV)");
    m_stack->Add(dy_m);
    m_stack->Add(diboson_m);
    m_stack->Add(ttbar_m);
    m_stack->Add(qcd_m);

    THStack *cost_stack = new THStack("cost_stack", "Samesign EMu cos Distribution (symmeterized): Data vs MC ; cos(#theta_r)");
    cost_stack->Add(dy_cost);
    cost_stack->Add(diboson_cost);
    cost_stack->Add(ttbar_cost);
    cost_stack->Add(qcd_cost);

    THStack *xf_stack = new THStack("xf_stack", "Samesign EMu x_F Distribution (symmeterized): Data vs MC ; x_F");
    xf_stack->Add(dy_xf);
    xf_stack->Add(diboson_xf);
    xf_stack->Add(ttbar_xf);
    xf_stack->Add(qcd_xf);


    THStack *pt_stack = new THStack("pt_stack", "Samesign EMu Dilepton pt; p_T (GeV)");
    pt_stack->Add(dy_pt);
    pt_stack->Add(diboson_pt);
    pt_stack->Add(ttbar_pt);
    pt_stack->Add(qcd_pt);

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
    leg1->AddEntry(qcd_m, "WJets + QCD (inferred)", "f");
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
   h_ratio->GetYaxis()->SetTitle("Obs./Exp.");
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
    //cost_pad1->SetLogy();
    cost_pad1->Draw();
    cost_pad1->cd();
    cost_stack->Draw("hist");
    data_cost->SetMarkerStyle(kFullCircle);
    data_cost->SetMarkerColor(1);
    cost_stack->SetMinimum(1);
    cost_stack->SetMaximum(3000);
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
   cost_ratio->GetYaxis()->SetTitle("Obs./Exp.");
   cost_ratio->GetYaxis()->SetNdivisions(505);
   cost_ratio->GetYaxis()->SetTitleSize(20);
   cost_ratio->GetYaxis()->SetTitleFont(43);
   cost_ratio->GetYaxis()->SetTitleOffset(1.2);
   cost_ratio->GetYaxis()->SetLabelFont(43); // Absolute font size in pixel (precision 3)
   cost_ratio->GetYaxis()->SetLabelSize(15);
   // X axis cost_ratio plot settings
   cost_ratio->GetXaxis()->SetTitle("e#mu |Cos(#theta)_{r}|");
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
   xf_ratio->GetYaxis()->SetTitle("Obs./Exp.");
   xf_ratio->GetYaxis()->SetNdivisions(505);
   xf_ratio->GetYaxis()->SetTitleSize(20);
   xf_ratio->GetYaxis()->SetTitleFont(43);
   xf_ratio->GetYaxis()->SetTitleOffset(1.2);
   xf_ratio->GetYaxis()->SetLabelFont(43); // Absolute font size in pixel (precision 3)
   xf_ratio->GetYaxis()->SetLabelSize(15);
   // X axis xf_ratio plot settings
   xf_ratio->GetXaxis()->SetTitle("e#mu xf (GeV)");
   xf_ratio->GetXaxis()->SetTitleSize(20);
   xf_ratio->GetXaxis()->SetTitleFont(43);
   xf_ratio->GetXaxis()->SetTitleOffset(3.);
   xf_ratio->GetXaxis()->SetLabelFont(43); // Absolute font size in pixel (precision 3)
   xf_ratio->GetXaxis()->SetLabelSize(20);
    CMS_lumi(xf_pad1, iPeriod, 11 );
    c_xf->Update();




    TCanvas *c_pt = new TCanvas("c_pt", "Histograms", 200, 10, 900, 700);
    TPad *pt_pad1 = new TPad("pad1", "pad1", 0.,0.3,0.98,1.);
    pt_pad1->SetBottomMargin(0);
    pt_pad1->Draw();
    pt_pad1->cd();
    pt_stack->Draw("hist");
    data_pt->SetMarkerStyle(kFullCircle);
    data_pt->SetMarkerColor(1);
    pt_stack->SetMinimum(1);
    pt_stack->SetMaximum(10000);
    data_pt->SetMinimum(1);
    data_pt->SetMaximum(10000);
    data_pt->Draw("P E same");
    pt_pad1->SetLogy();
    c_pt->Update();
    leg1->Draw();

    c_pt->cd();
    TPad *pt_pad2 = new TPad("pt_pad2", "pad2", 0.,0,.98,0.3);
    //pad2->SetTopMargin(0);
    pt_pad2->SetBottomMargin(0.2);
    pt_pad2->SetGridy();
    pt_pad2->Draw();
    pt_pad2->cd();
    TList *pt_stackHists = pt_stack->GetHists();
    TH1* pt_mc_sum = (TH1*)pt_stackHists->At(0)->Clone();
    pt_mc_sum->Reset();

    for (int i=0;i<pt_stackHists->GetSize();++i) {
      pt_mc_sum->Add((TH1*)pt_stackHists->At(i));
    }
    auto pt_ratio = (TH1F *) data_pt->Clone("h_pt_ratio");
    pt_ratio->SetMinimum(0.);
    pt_ratio->SetMaximum(10.);
    pt_ratio->Sumw2();
    pt_ratio->SetStats(0);
    pt_ratio->Divide(pt_mc_sum);
    pt_ratio->SetMarkerStyle(21);
    pt_ratio->Draw("ep");
    c_pt->cd();

    pt_ratio->SetTitle("");
    // Y axis pt_ratio plot settings
   pt_ratio->GetYaxis()->SetTitle("Obs./Exp.");
   pt_ratio->GetYaxis()->SetNdivisions(505);
   pt_ratio->GetYaxis()->SetTitleSize(20);
   pt_ratio->GetYaxis()->SetTitleFont(43);
   pt_ratio->GetYaxis()->SetTitleOffset(1.2);
   pt_ratio->GetYaxis()->SetLabelFont(43); // Absolute font size in pixel (precision 3)
   pt_ratio->GetYaxis()->SetLabelSize(15);
   // X axis pt_ratio plot settings
   pt_ratio->GetXaxis()->SetTitle("e#mu pt (GeV)");
   pt_ratio->GetXaxis()->SetTitleSize(20);
   pt_ratio->GetXaxis()->SetTitleFont(43);
   pt_ratio->GetXaxis()->SetTitleOffset(3.);
   pt_ratio->GetXaxis()->SetLabelFont(43); // Absolute font size in pixel (precision 3)
   pt_ratio->GetXaxis()->SetLabelSize(20);
   CMS_lumi(pt_pad1, iPeriod, 11 );
   c_pt->Update();
}









