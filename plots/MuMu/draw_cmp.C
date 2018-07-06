
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
#include "root_files.h"

const int type = FLAG_MUONS;



void draw_cmp(){
    setTDRStyle();
    init();

    TH1F *data_m = new TH1F("data_m", "Data Dimuon Mass Distribution", 30, 150, 2000);
    TH1F *mc_pt = new TH1F("mc_pt", "MC signal", 40, 0, 1000);
    TH1F *mc_nosig_pt = new TH1F("mc_nosig_pt", "MC signal", 40, 0, 1000);
    TH1F *data_pt = new TH1F("data_pt", "MC signal", 40, 0, 1000);
    TH1F *ttbar_pt = new TH1F("ttbar_pt", "MC signal", 40, 0, 1000);
    TH1F *diboson_pt = new TH1F("diboson_pt", "MC signal", 40, 0, 1000);
    TH1F *wt_pt = new TH1F("wt_pt", "MC signal", 40, 0, 1000);
    TH1F *QCD_pt = new TH1F("QCD_pt", "MC signal", 40, 0, 1000);
    mc_pt->SetFillColor(kRed+1);
    mc_pt->SetMarkerColor(kRed+1);
    mc_nosig_pt->SetFillColor(kMagenta);
    ttbar_pt->SetFillColor(kBlue);
    ttbar_pt->SetMarkerStyle(21);
    ttbar_pt->SetMarkerColor(kBlue);
    diboson_pt->SetFillColor(kGreen+3);
    QCD_pt->SetFillColor(kRed -7);
    wt_pt->SetFillColor(kOrange+7); 

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


    TH1F *data_cost = new TH1F("data_cost", "Data", 20, -1.,1.);
    TH1F *mc_cost = new TH1F("mc_cost", "MC Signal (qqbar, qglu, qbarglu)", 20, -1,1);
    mc_cost->SetFillColor(kRed+1);
    mc_cost->SetMarkerStyle(21);
    mc_cost->SetMarkerColor(kRed+1);
    TH1F *mc_nosig_cost = new TH1F("mc_nosig_cost", "MC no signal (qq, gluglu qbarqbar)", 20, -1.,1.);
    mc_nosig_cost->SetFillColor(kMagenta);
    mc_nosig_cost->SetMarkerStyle(21);
    mc_nosig_cost->SetMarkerColor(kMagenta);
    TH1F *ttbar_cost = new TH1F("ttbar_cost", "TTbar Background", 20, -1.,1.);
    ttbar_cost->SetFillColor(kBlue);
    ttbar_cost->SetMarkerStyle(21);
    ttbar_cost->SetMarkerColor(kBlue);




    TH1F *diboson_m = new TH1F("diboson_m", "DiBoson (WW, WZ, ZZ)", 30, 150, 2000);
    TH1F *diboson_cost = new TH1F("diboson_cost", "DiBoson (WW, WZ,ZZ)", 20, -1,1);

    TH1F *QCD_m = new TH1F("QCD_m", "QCD", 30, 150, 2000);
    TH1F *QCD_cost = new TH1F("QCD_cost", "QCD", 20, -1,1);

    TH1F *WJets_m = new TH1F("WJets_m", "WJets", 30, 150, 2000);
    TH1F *WJets_cost = new TH1F("WJets_cost", "WJets", 20, -1,1);

    TH1F *wt_m = new TH1F("wt_m", "tw + #bar{t}w", 30, 150, 2000);
    TH1F *wt_cost = new TH1F("wt_cost", "tw + #bar{t}w", 20, -1,1);


    wt_m->SetFillColor(kOrange+7); 
    wt_cost->SetFillColor(kOrange+7); 
    diboson_m->SetFillColor(kGreen+3);
    diboson_cost->SetFillColor(kGreen + 3);
    QCD_m->SetFillColor(kRed -7);
    QCD_cost->SetFillColor(kRed -7);

    bool do_RC = false;

    make_m_cost_pt_hist(t_data, data_m, data_cost, data_pt, true, type, do_RC);
    make_m_cost_pt_hist(t_mc, mc_m, mc_cost, mc_pt, false, type, do_RC);
    make_m_cost_pt_hist(t_mc_nosig, mc_nosig_m, mc_nosig_cost, mc_nosig_pt, false, type, do_RC);
    make_m_cost_pt_hist(t_ttbar, ttbar_m, ttbar_cost, ttbar_pt, false, type, do_RC);
    make_m_cost_pt_hist(t_wt, wt_m, wt_cost, wt_pt, false, type, do_RC);
    make_m_cost_pt_hist(t_diboson, diboson_m, diboson_cost, diboson_pt, false, type, do_RC);

    Fakerate_est_mu(t_WJets, t_QCD, t_WJets_mc, t_QCD_mc, QCD_m, QCD_cost, QCD_pt);




    Double_t EMu_ratio= 1.05;
    ttbar_m->Scale(EMu_ratio);
    ttbar_cost->Scale(EMu_ratio);
    diboson_m->Scale(EMu_ratio);
    diboson_cost->Scale(EMu_ratio);
    wt_m->Scale(EMu_ratio);
    wt_cost->Scale(EMu_ratio);


    int nBins_x = QCD_m->GetXaxis()->GetNbins();
    int nBins_y = QCD_cost->GetYaxis()->GetNbins();
    //printf("Get size %i \n", nBins);
    for (int i=1; i <= nBins_x; i++){
        for (int j=1; j <= nBins_y; j++){

            Double_t m_val = QCD_m->GetBinContent(i,j);
            Double_t cost_val = QCD_cost->GetBinContent(i,j);

            QCD_m->SetBinError(i,j, 0.2*m_val);
            QCD_cost->SetBinError(i,j, 0.2*cost_val);
        }
    }


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


    THStack *cost_stack = new THStack("cost_stack", "Cos(#theta) Distribution: Data vs MC; MuMu Cos(#theta)_{r}");
    cost_stack->Add(ttbar_cost);
    cost_stack->Add(QCD_cost);
    cost_stack->Add(wt_cost);
    cost_stack->Add(diboson_cost);
    cost_stack->Add(mc_nosig_cost);
    cost_stack->Add(mc_cost);

    THStack *pt_stack = new THStack("pt_stack", "Dimuon Pt Distribution: Data vs MC; Dimuon Pt (GeV)");
    pt_stack->Add(ttbar_pt);
    pt_stack->Add(QCD_pt);
    pt_stack->Add(wt_pt);
    pt_stack->Add(diboson_pt);
    pt_stack->Add(mc_nosig_pt);
    pt_stack->Add(mc_pt);

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
    pad2->SetGridy();
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
 
    //lumi_sqrtS = "";       // used with iPeriod = 0, e.g. for simulation-only plots (default is an empty string)
    int iPeriod = 4; 
    CMS_lumi(pad1, iPeriod, 33 );



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
    TLegend *leg2 = new TLegend(0.5, 0.65, 0.75, 0.8);
    leg2->AddEntry(data_m, "data", "p");
    leg2->AddEntry(mc_m, "DY (q#bar{q}, qg #bar{q}g)", "f");
    leg2->AddEntry(mc_nosig_m, "DY no asymmety(gg, qq, #bar{q}#bar{q})", "f");
    leg2->AddEntry(diboson_m, "WW + WZ + ZZ", "f");
    leg2->AddEntry(wt_m, "tW + #bar{t}W", "f");
    leg2->AddEntry(QCD_m, "QCD + WJets", "f");
    leg2->AddEntry(ttbar_m, "t#bar{t}", "f");
    leg2->Draw();

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
    /*
    leg = new TLegend(0.1,0.6,0.5,0.9);
    leg->AddEntry(h_cost1, "No extra cuts", "f");
    leg->AddEntry(h_cost_cut, "met_pt < 80, mu1_pt > 40, mu2_pt > 20", "f");
    leg->Draw();
    */

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
    TLegend *leg3 = new TLegend(0.5, 0.65, 0.75, 0.8);
    leg3->AddEntry(data_pt, "data", "p");
    leg3->AddEntry(mc_pt, "DY (q#bar{q}, qg #bar{q}g)", "f");
    leg3->AddEntry(mc_nosig_pt, "DY no asymmety(gg, qq, #bar{q}#bar{q})", "f");
    leg3->AddEntry(diboson_pt, "WW + WZ + ZZ", "f");
    leg3->AddEntry(wt_pt, "tW + #bar{t}W", "f");
    leg3->AddEntry(QCD_pt, "QCD + WJets", "f");
    leg3->AddEntry(ttbar_pt, "t#bar{t}", "f");
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
   pt_ratio->GetYaxis()->SetTitle("Data/MC");
   pt_ratio->GetYaxis()->SetNdivisions(505);
   pt_ratio->GetYaxis()->SetTitleSize(20);
   pt_ratio->GetYaxis()->SetTitleFont(43);
   pt_ratio->GetYaxis()->SetTitleOffset(1.2);
   pt_ratio->GetYaxis()->SetLabelFont(43); // Absolute font size in pixel (precision 3)
   pt_ratio->GetYaxis()->SetLabelSize(15);
   // X axis pt_ratio plot settings
   pt_ratio->GetXaxis()->SetTitle("dimuon pt (GeV)");
   pt_ratio->GetXaxis()->SetTitleSize(20);
   pt_ratio->GetXaxis()->SetTitleFont(43);
   pt_ratio->GetXaxis()->SetTitleOffset(3.);
   pt_ratio->GetXaxis()->SetLabelFont(43); // Absolute font size in pixel (precision 3)
   pt_ratio->GetXaxis()->SetLabelSize(20);
    CMS_lumi(pt_pad1, iPeriod, 11 );
    c_pt->Update();




 
}

    
    
