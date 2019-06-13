
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
#include "root_files.h"

const int type = FLAG_ELECTRONS;



void draw_cmp(){
    setTDRStyle();
    init();
    ee_init();
    init_gamgam();

    int n_pt_bins = 40;
    TH1F *mc_pt = new TH1F("mc_pt", "MC signal", n_pt_bins, 0, 1000);
    TH1F *mc_nosig_pt = new TH1F("mc_nosig_pt", "MC signal", n_pt_bins, 0, 1000);
    TH1F *data_pt = new TH1F("data_pt", "MC signal", n_pt_bins, 0, 1000);
    TH1F *ttbar_pt = new TH1F("ttbar_pt", "MC signal", n_pt_bins, 0, 1000);
    TH1F *diboson_pt = new TH1F("diboson_pt", "MC signal", n_pt_bins, 0, 1000);
    TH1F *wt_pt = new TH1F("wt_pt", "MC signal", n_pt_bins, 0, 1000);
    TH1F *QCD_pt = new TH1F("QCD_pt", "MC signal", n_pt_bins, 0, 1000);
    TH1F *gg_pt = new TH1F("gg_pt", "MC signal", n_pt_bins, 0, 1000);
    mc_pt->SetFillColor(kRed+1);
    mc_pt->SetMarkerColor(kRed+1);
    mc_nosig_pt->SetFillColor(kMagenta);
    ttbar_pt->SetFillColor(kBlue);
    ttbar_pt->SetMarkerStyle(21);
    ttbar_pt->SetMarkerColor(kBlue);
    diboson_pt->SetFillColor(kGreen+3);
    QCD_pt->SetFillColor(kRed -7);
    wt_pt->SetFillColor(kOrange+7); 

    int xf_nbins = 10;
    TH1F *mc_xf = new TH1F("mc_xf", "MC signal", xf_nbins, 0, 0.5);
    TH1F *mc_nosig_xf = new TH1F("mc_nosig_xf", "MC signal", xf_nbins, 0, 0.5);
    TH1F *data_xf = new TH1F("data_xf", "MC signal", xf_nbins, 0, 0.5);
    TH1F *ttbar_xf = new TH1F("ttbar_xf", "MC signal", xf_nbins, 0, 0.5);
    TH1F *diboson_xf = new TH1F("diboson_xf", "MC signal", xf_nbins, 0, 0.5);
    TH1F *wt_xf = new TH1F("wt_xf", "MC signal", xf_nbins, 0, 0.5);
    TH1F *QCD_xf = new TH1F("QCD_xf", "MC signal", xf_nbins, 0, 0.5);
    TH1F *gg_xf = new TH1F("gg_xf", "MC signal", xf_nbins, 0, 0.5);
    mc_xf->SetFillColor(kRed+1);
    mc_xf->SetMarkerColor(kRed+1);
    mc_nosig_xf->SetFillColor(kMagenta);
    ttbar_xf->SetFillColor(kBlue);
    ttbar_xf->SetMarkerStyle(21);
    ttbar_xf->SetMarkerColor(kBlue);
    diboson_xf->SetFillColor(kGreen+3);
    QCD_xf->SetFillColor(kRed -7);
    wt_xf->SetFillColor(kOrange+7); 

    int n_m_bins = 30;
    TH1F *data_m = new TH1F("data_m", "Data Dimuon Mass Distribution", n_m_bins, 150, 2000);
    TH1F *mc_m = new TH1F("mc_m", "MC Signal (qqbar, qglu, qbarglu)", n_m_bins, 150, 2000);
    TH1F *mc_nosig_m = new TH1F("mc_nosig_m", "MC no signal (qq, gluglu qbarqbar)", n_m_bins, 150, 2000);
    TH1F *ttbar_m = new TH1F("ttbar_m", "TTBar Background", n_m_bins, 150, 2000);
    TH1F *diboson_m = new TH1F("diboson_m", "DiBoson (WW, WZ, ZZ)", n_m_bins, 150, 2000);
    TH1F *QCD_m = new TH1F("QCD_m", "QCD", n_m_bins, 150, 2000);
    TH1F *gg_m = new TH1F("gg_m", "QCD", n_m_bins, 150, 2000);
    TH1F *WJets_m = new TH1F("WJets_m", "WJets", n_m_bins, 150, 2000);
    TH1F *wt_m = new TH1F("wt_m", "tw + #bar{t}w", n_m_bins, 150, 2000);
    mc_m->SetFillColor(kRed+1);
    mc_m->SetMarkerStyle(21);
    mc_m->SetMarkerColor(kRed+1);
    mc_nosig_m->SetFillColor(kMagenta);
    mc_nosig_m->SetMarkerStyle(21);
    mc_nosig_m->SetMarkerColor(kMagenta);
    ttbar_m->SetFillColor(kBlue);
    ttbar_m->SetMarkerStyle(21);
    ttbar_m->SetMarkerColor(kBlue);


    int num_cost_bins = 10;
    TH1F *data_cost = new TH1F("data_cost", "Data", n_cost_bins, -1.,1.);
    TH1F *mc_cost = new TH1F("mc_cost", "MC Signal (qqbar, qglu, qbarglu)", n_cost_bins, -1,1);
    TH1F *mc_nosig_cost = new TH1F("mc_nosig_cost", "MC no signal (qq, gluglu qbarqbar)", n_cost_bins, -1.,1.);
    TH1F *ttbar_cost = new TH1F("ttbar_cost", "TTbar Background", n_cost_bins, -1.,1.);
    TH1F *diboson_cost = new TH1F("diboson_cost", "DiBoson (WW, WZ,ZZ)", n_cost_bins, -1,1);
    TH1F *QCD_cost = new TH1F("QCD_cost", "QCD", n_cost_bins, -1,1);
    TH1F *gg_cost = new TH1F("gg_cost", "QCD", n_cost_bins, -1,1);
    TH1F *WJets_cost = new TH1F("WJets_cost", "WJets", n_cost_bins, -1,1);
    TH1F *wt_cost = new TH1F("wt_cost", "tw + #bar{t}w", n_cost_bins, -1,1);

    mc_cost->SetFillColor(kRed+1);
    mc_cost->SetMarkerStyle(21);
    mc_cost->SetMarkerColor(kRed+1);
    mc_nosig_cost->SetFillColor(kMagenta);
    mc_nosig_cost->SetMarkerStyle(21);
    mc_nosig_cost->SetMarkerColor(kMagenta);
    ttbar_cost->SetFillColor(kBlue);
    ttbar_cost->SetMarkerStyle(21);
    ttbar_cost->SetMarkerColor(kBlue);


    wt_m->SetFillColor(kOrange+7); 
    wt_cost->SetFillColor(kOrange+7); 
    diboson_m->SetFillColor(kGreen+3);
    diboson_cost->SetFillColor(kGreen + 3);
    QCD_m->SetFillColor(kRed -7);
    QCD_cost->SetFillColor(kRed -7);

    gg_cost->SetFillColor(kOrange);
    gg_m->SetFillColor(kOrange);
    gg_pt->SetFillColor(kOrange);
    gg_xf->SetFillColor(kOrange);


    float m_low = 150.;
    float m_high = 100000.;

    bool do_RC = false;
    make_m_cost_pt_xf_hist(t_elel_data, data_m, data_cost, data_pt, data_xf, true, type, do_RC, m_low, m_high);
    make_m_cost_pt_xf_hist(t_elel_mc, mc_m, mc_cost, mc_pt, mc_xf, false, type,  do_RC, m_low, m_high);
    make_m_cost_pt_xf_hist(t_elel_nosig, mc_nosig_m, mc_nosig_cost, mc_nosig_pt, mc_nosig_xf, false, type, do_RC, m_low, m_high);
    make_m_cost_pt_xf_hist(t_ttbar, ttbar_m, ttbar_cost, ttbar_pt, ttbar_xf, false, type, do_RC, m_low, m_high);
    make_m_cost_pt_xf_hist(t_wt, wt_m, wt_cost, wt_pt, wt_xf, false, type, do_RC, m_low, m_high);
    make_m_cost_pt_xf_hist(t_elel_gamgam, gg_m, gg_cost, gg_pt, gg_xf, false, type, do_RC, m_low, m_high);
    make_m_cost_pt_xf_hist(t_diboson, diboson_m, diboson_cost, diboson_pt, diboson_xf, false, type,  do_RC, m_low, m_high);


    symmetrize1d(gg_cost);

    //gg_cost->Scale(0.);
    //gg_m->Scale(0.);

    bool ss_qcd = true;
    bool in_os_region = true;
    Fakerate_est_el(t_elel_WJets, t_elel_QCD, t_elel_WJets_contam, t_elel_QCD_contam, QCD_m, QCD_cost, QCD_pt, QCD_xf, m_low, m_high, ss_qcd, in_os_region);


    printf("Data integral is %.2f \n", data_cost->Integral());
    printf("DY integral is %.2f \n", mc_cost->Integral());
    printf("ttbar integral is %.2f \n", ttbar_cost->Integral());






    bool scale_qcd_error=true;
    float qcd_err = 0.4;
    setHistError(QCD_m, qcd_err);
    setHistError(QCD_cost, qcd_err);
    setHistError(QCD_pt, qcd_err);
    setHistError(QCD_xf, qcd_err);

    float gam_err = 0.4;
    setHistError(gg_m, gam_err);
    setHistError(gg_cost, gam_err);
    setHistError(gg_pt, gam_err);
    setHistError(gg_xf, gam_err);



    

    THStack *m_stack = new THStack("m_stack", "MuMu Mass Distribution: Data vs MC ; m_{#mu^{+}#mu^{-}} (GeV)");
    m_stack->Add(diboson_m);
    m_stack->Add(QCD_m);
    m_stack->Add(wt_m);
    m_stack->Add(ttbar_m);
    m_stack->Add(gg_m);
    m_stack->Add(mc_nosig_m);
    m_stack->Add(mc_m);


    THStack *cost_stack = new THStack("cost_stack", "Cos(#theta) Distribution: Data vs MC; MuMu Cos(#theta)_{r}");
    cost_stack->Add(diboson_cost);
    cost_stack->Add(QCD_cost);
    cost_stack->Add(wt_cost);
    cost_stack->Add(ttbar_cost);
    cost_stack->Add(gg_cost);
    cost_stack->Add(mc_nosig_cost);
    cost_stack->Add(mc_cost);

    THStack *pt_stack = new THStack("pt_stack", "Dimuon Pt Distribution: Data vs MC; Dimuon Pt (GeV)");
    pt_stack->Add(diboson_pt);
    pt_stack->Add(QCD_pt);
    pt_stack->Add(wt_pt);
    pt_stack->Add(ttbar_pt);
    pt_stack->Add(gg_pt);
    pt_stack->Add(mc_nosig_pt);
    pt_stack->Add(mc_pt);

    THStack *xf_stack = new THStack("xf_stack", "Di-electron x_F Distribution: Data vs MC; x_F");
    xf_stack->Add(diboson_xf);
    xf_stack->Add(QCD_xf);
    xf_stack->Add(wt_xf);
    xf_stack->Add(ttbar_xf);
    xf_stack->Add(gg_xf);
    xf_stack->Add(mc_nosig_xf);
    xf_stack->Add(mc_xf);

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
    TLegend *leg1 = new TLegend(0.3, 0.3);
    leg1->AddEntry(data_m, "data", "p");
    leg1->AddEntry(mc_m, "DY (q#bar{q}, qg #bar{q}g)", "f");
    leg1->AddEntry(mc_nosig_m, "DY no asymmety(gg, qq, #bar{q}#bar{q})", "f");
    leg1->AddEntry(ttbar_m, "t#bar{t}", "f");
    leg1->AddEntry(gg_m, "#gamma#gamma to ee", "f");
    leg1->AddEntry(wt_m, "tW + #bar{t}W", "f");
    leg1->AddEntry(QCD_m, "QCD + WJets", "f");
    leg1->AddEntry(diboson_m, "WW + WZ + ZZ", "f");
    leg1->Draw();

    //gPad->BuildLegend();
    c_m->cd();
    TPad *pad2 = new TPad("pad2", "pad2", 0.,0,.98,0.3);
    //pad2->SetTopMargin(0);
    pad2->SetBottomMargin(0.2);
    pad2->SetGridy();
    pad2->Draw();
    pad2->cd();



    TList *stackHists = m_stack->GetHists();
    TH1* m_mc_sum = (TH1*)stackHists->At(0)->Clone();
    m_mc_sum->Reset();

    for (int i=0;i<stackHists->GetSize();++i) {
      m_mc_sum->Add((TH1*)stackHists->At(i));
    }
    auto m_ratio = (TH1F *) data_m->Clone("h_m_ratio");
    m_ratio->SetMinimum(0.75);
    m_ratio->SetMaximum(1.25);
    m_ratio->Sumw2();
    m_ratio->SetStats(0);
    m_ratio->Divide(m_mc_sum);
    m_ratio->SetMarkerStyle(21);
    m_ratio->Draw("ep");
    TLine *l1 = new TLine(150,1,2000,1);
    l1->SetLineStyle(2);
    l1->Draw();
    c_m->cd();

    m_ratio->SetTitle("");
    // Y axis m_ratio plot settings
   m_ratio->GetYaxis()->SetTitle("Obs/Exp");
   m_ratio->GetYaxis()->SetNdivisions(505);
   m_ratio->GetYaxis()->SetTitleSize(20);
   m_ratio->GetYaxis()->SetTitleFont(43);
   m_ratio->GetYaxis()->SetTitleOffset(1.2);
   m_ratio->GetYaxis()->SetLabelFont(43); // Absolute font size in pixel (precision 3)
   m_ratio->GetYaxis()->SetLabelSize(15);
   // X axis m_ratio plot settings
   m_ratio->GetXaxis()->SetTitle("M_{ee} (GeV)");
   m_ratio->GetXaxis()->SetTitleSize(20);
   m_ratio->GetXaxis()->SetTitleFont(43);
   m_ratio->GetXaxis()->SetTitleOffset(3.);
   m_ratio->GetXaxis()->SetLabelFont(43); // Absolute font size in pixel (precision 3)
   m_ratio->GetXaxis()->SetLabelSize(20);
 
    //lumi_sqrtS = "";       // used with iPeriod = 0, e.g. for simulation-only plots (default is an empty string)
    int iPeriod = 4; 
    //CMS_lumi(pad1, iPeriod, 33 );



    TCanvas *c_cost = new TCanvas("c_cost", "Histograms", 200, 10, 900, 700);
    TPad *cost_pad1 = new TPad("pad1c", "pad1", 0.,0.3,0.98,1.);
    cost_pad1->SetBottomMargin(0);
    //c_cost->SetLogy();
    cost_pad1->Draw();
    cost_pad1->cd();
    cost_stack->Draw("hist");
    data_cost->SetMarkerStyle(kFullCircle);
    data_cost->SetMarkerColor(1);
    cost_stack->SetMinimum(1);
    data_cost->Draw("P E same");
    cost_pad1->Update();
    TLegend *leg2 = (TLegend *) leg1->Clone("leg2");
    leg2->Draw();

    //CMS_lumi(cost_pad1, iPeriod, 11 );
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
    cost_ratio->SetMinimum(0.7);
    cost_ratio->SetMaximum(1.3);
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
   cost_ratio->GetYaxis()->SetTitle("Obs/Exp");
   cost_ratio->GetYaxis()->SetNdivisions(505);
   cost_ratio->GetYaxis()->SetTitleSize(20);
   cost_ratio->GetYaxis()->SetTitleFont(43);
   cost_ratio->GetYaxis()->SetTitleOffset(1.2);
   cost_ratio->GetYaxis()->SetLabelFont(43); // Absolute font size in pixel (precision 3)
   cost_ratio->GetYaxis()->SetLabelSize(15);
   // X axis cost_ratio plot settings
   cost_ratio->GetXaxis()->SetTitle("dielectron Cos(#theta)_{r} (GeV)");
   cost_ratio->GetXaxis()->SetTitleSize(20);
   cost_ratio->GetXaxis()->SetTitleFont(43);
   cost_ratio->GetXaxis()->SetTitleOffset(3.);
   cost_ratio->GetXaxis()->SetLabelFont(43); // Absolute font size in pixel (precision 3)
   cost_ratio->GetXaxis()->SetLabelSize(20);

    //TCanvas *c_cost_cut = new TCanvas("c_cost_cut", "Histograms", 200, 10, 900, 700);
    //c_cost_cut->cd();

    TCanvas *c_pt = new TCanvas("c_pt", "Histograms", 200, 10, 900, 700);
    TPad *pt_pad1 = new TPad("pad1", "pad1", 0.,0.3,0.98,1.);
    pt_pad1->SetBottomMargin(0);
    pt_pad1->Draw();
    pt_pad1->cd();
    pt_stack->Draw("hist");
    data_pt->SetMarkerStyle(kFullCircle);
    data_pt->SetMarkerColor(1);
    pt_stack->SetMinimum(1);
    pt_stack->SetMaximum(100000);
    data_pt->SetMinimum(1);
    data_pt->SetMaximum(100000);
    data_pt->Draw("P E same");
    pt_pad1->SetLogy();
    c_pt->Update();
    TLegend *leg3 = (TLegend *) leg1->Clone("leg3");
    leg3->Draw();

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
    pt_ratio->SetMinimum(0.7);
    pt_ratio->SetMaximum(1.3);
    pt_ratio->Sumw2();
    pt_ratio->SetStats(0);
    pt_ratio->Divide(pt_mc_sum);
    pt_ratio->SetMarkerStyle(21);
    pt_ratio->Draw("ep");
    c_pt->cd();

    pt_ratio->SetTitle("");
    // Y axis pt_ratio plot settings
   pt_ratio->GetYaxis()->SetTitle("Obs/Exp");
   pt_ratio->GetYaxis()->SetNdivisions(505);
   pt_ratio->GetYaxis()->SetTitleSize(20);
   pt_ratio->GetYaxis()->SetTitleFont(43);
   pt_ratio->GetYaxis()->SetTitleOffset(1.2);
   pt_ratio->GetYaxis()->SetLabelFont(43); // Absolute font size in pixel (precision 3)
   pt_ratio->GetYaxis()->SetLabelSize(15);
   // X axis pt_ratio plot settings
   pt_ratio->GetXaxis()->SetTitle("dielectron pt (GeV)");
   pt_ratio->GetXaxis()->SetTitleSize(20);
   pt_ratio->GetXaxis()->SetTitleFont(43);
   pt_ratio->GetXaxis()->SetTitleOffset(3.);
   pt_ratio->GetXaxis()->SetLabelFont(43); // Absolute font size in pixel (precision 3)
   pt_ratio->GetXaxis()->SetLabelSize(20);
    //CMS_lumi(pt_pad1, iPeriod, 11 );
    c_pt->Update();




    TCanvas *c_xf = new TCanvas("c_xf", "Histograms", 200, 10, 900, 700);
    TPad *xf_pad1 = new TPad("pad1", "pad1", 0.,0.3,0.98,1.);
    xf_pad1->SetBottomMargin(0);
    xf_pad1->Draw();
    xf_pad1->cd();
    xf_stack->Draw("hist");
    data_xf->SetMarkerStyle(kFullCircle);
    data_xf->SetMarkerColor(1);
    xf_stack->SetMinimum(1);
    xf_stack->SetMaximum(100000);
    data_xf->SetMinimum(1);
    data_xf->SetMaximum(100000);
    data_xf->Draw("P E same");
    xf_pad1->SetLogy();
    c_xf->Update();
    TLegend *leg4 = (TLegend *) leg1->Clone("leg4");
    leg4->Draw();

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
    xf_ratio->SetMinimum(0.7);
    xf_ratio->SetMaximum(1.3);
    xf_ratio->Sumw2();
    xf_ratio->SetStats(0);
    xf_ratio->Divide(xf_mc_sum);
    xf_ratio->SetMarkerStyle(21);
    xf_ratio->Draw("ep");
    c_xf->cd();

    xf_ratio->SetTitle("");
    // Y axis xf_ratio plot settings
   xf_ratio->GetYaxis()->SetTitle("Obs/Exp");
   xf_ratio->GetYaxis()->SetNdivisions(505);
   xf_ratio->GetYaxis()->SetTitleSize(20);
   xf_ratio->GetYaxis()->SetTitleFont(43);
   xf_ratio->GetYaxis()->SetTitleOffset(1.2);
   xf_ratio->GetYaxis()->SetLabelFont(43); // Absolute font size in pixel (precision 3)
   xf_ratio->GetYaxis()->SetLabelSize(15);
   // X axis xf_ratio plot settings
   xf_ratio->GetXaxis()->SetTitle(" x_F");
   xf_ratio->GetXaxis()->SetTitleSize(20);
   xf_ratio->GetXaxis()->SetTitleFont(43);
   xf_ratio->GetXaxis()->SetTitleOffset(3.);
   xf_ratio->GetXaxis()->SetLabelFont(43); // Absolute font size in pixel (precision 3)
   xf_ratio->GetXaxis()->SetLabelSize(20);
    //CMS_lumi(xf_pad1, iPeriod, 11 );
    c_xf->Update();

    float m_chi2 = computeChi2(m_ratio);
    float cost_chi2 = computeChi2(cost_ratio);

    printf("M-bins chi2/dof = %.1f/%i \n", m_chi2, n_m_bins);
    printf("cost-bins chi2/dof = %.1f/%i \n", cost_chi2, n_cost_bins);
 
}

    
    
