
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
#include "../tdrstyle.C"
#include "../CMS_lumi.C"
#include "../../Utils/HistMaker.C"
#include "../../Utils/root_files.h"





void draw_emu_new(){
    setTDRStyle();
    init_emu();

                                
    //TFile *f_ttbar = TFile::Open("../analyze/output_files/EMu_TTbar_Mu_june29.root");
    TFile *f_ttbar = TFile::Open("../analyze/output_files/EMu_ttbar_mar7.root");
    TTree *t_ttbar = (TTree *)f_ttbar->Get("T_data");

    TFile *f_diboson = TFile::Open("../analyze/output_files/EMu_diboson_mar19.root");
    TTree *t_diboson = (TTree *)f_diboson->Get("T_data");

    TFile *f_wt = TFile::Open("../analyze/output_files/EMu_wt_back_mar18.root");
    TTree *t_wt = (TTree *)f_wt->Get("T_data");

    t_ttbar->Print(); 
    t_diboson->Print();
    t_wt->Print();

    TH1F *data_m = new TH1F("data_m", "MC Signal (qqbar, qglu, qbarglu)", 30, 150, 1000);
    TH1F *ttbar_m = new TH1F("ttbar_m", "MC Signal (qqbar, qglu, qbarglu)", 30, 150, 1000);
    TH1F *diboson_m = new TH1F("diboson_m", "MC Signal (qqbar, qglu, qbarglu)", 30, 150, 1000);
    TH1F *wt_m = new TH1F("wt_m", "MC Signal (qqbar, qglu, qbarglu)", 30, 150, 1000);
    TH1F *dy_m = new TH1F("dy_m", "MC Signal (qqbar, qglu, qbarglu)", 30, 150, 1000);
    TH1F *qcd_m = new TH1F("qcd_m", "MC Signal (qqbar, qglu, qbarglu)", 30, 150, 1000);

    TH1F *h_dummy = new TH1F("h_dummy", "", 100, 0, 100.);

    Double_t m_low = 150;
    Double_t m_high = 10000;

    int type  = FLAG_MUONS;

    make_emu_m_cost_pt_xf_hist(t_emu_data, data_m, h_dummy, h_dummy, h_dummy, true, type);
    make_emu_m_cost_pt_xf_hist(t_ttbar, ttbar_m, h_dummy, h_dummy, h_dummy, false, type);
    make_emu_m_cost_pt_xf_hist(t_diboson, diboson_m, h_dummy, h_dummy, h_dummy, false, type);
    make_emu_m_cost_pt_xf_hist(t_wt, wt_m, h_dummy, h_dummy, h_dummy, false, type);
    make_emu_m_cost_pt_xf_hist(t_emu_dy, dy_m, h_dummy, h_dummy, h_dummy, false, type);
    Fakerate_est_emu(t_emu_WJets, t_emu_QCD,t_emu_WJets_contam, qcd_m, type);

    //correct for wrong ttbar xsec
    //ttbar_m->Scale(831.76/730.6);

    Double_t data_count = data_m->Integral();
    Double_t fake_count = qcd_m->Integral();
    Double_t mc_count = ttbar_m->Integral() + diboson_m->Integral() + wt_m->Integral() + dy_m->Integral();

    Double_t fake_unc = 0.35 * fake_count;
    Double_t mc_unc = sqrt(mc_count);

    printf("Data count %.0f \n", data_count);
    printf("MC count %.0f \n", mc_count);
    printf("Fake count %.0f \n", fake_count);
    Double_t ratio = data_count / (mc_count + fake_count);
    Double_t unc = sqrt( (data_count/(mc_count + fake_count)/(mc_count + fake_count)) +  
                          pow(data_count/(mc_count + fake_count)/(mc_count + fake_count), 2) * (fake_unc*fake_unc + mc_unc*mc_unc));
    printf("Ratio is %1.3f +/- %1.3f \n", ratio, unc);



    dy_m->SetFillColor(kRed+1);
    ttbar_m->SetFillColor(kBlue);
    wt_m->SetFillColor(kOrange+7); 
    diboson_m->SetFillColor(kGreen+3);
    qcd_m->SetFillColor(kRed -7);

    THStack *m_stack = new THStack("m_stack", "EMu Mass Distribution: Data vs MC ; m_{e#mu} (GeV)");
    m_stack->Add(dy_m);
    m_stack->Add(diboson_m);
    m_stack->Add(wt_m);
    m_stack->Add(ttbar_m);
    m_stack->Add(qcd_m);



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
    leg1->AddEntry(ttbar_m, "t#bar{t}", "f");
    leg1->AddEntry(qcd_m, "QCD and W+Jets", "f");
    leg1->AddEntry(wt_m, "tW + #bar{t}W", "f");
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
    h_ratio->SetMinimum(0.75);
    h_ratio->SetMaximum(1.25);
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
}









