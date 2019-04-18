


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

const int type = FLAG_ELECTRONS;


void draw_samesign_cmp(){
    init();
    init_ss();
    init_emu_ss();




    setTDRStyle();

    TH1F *data_m = new TH1F("data_m", "Data Dimuon Mass Distribution", 30, 150, 2000);


    int pt_bins = 20.;
    TH1F *data_pt = new TH1F("data_pt", "MC signal", pt_bins, 0, 800);
    TH1F *back_pt = new TH1F("back_pt", "MC signal", pt_bins, 0, 800);
    TH1F *QCD_pt = new TH1F("QCD_pt", "MC signal", pt_bins, 0, 800);
    TH1F *DY_pt = new TH1F("DY_pt", "MC signal", pt_bins, 0, 800);

    int xf_nbins = 16;
    TH1F *data_xf = new TH1F("data_xf", "MC signal", xf_nbins, 0, 0.8);
    TH1F *QCD_xf = new TH1F("QCD_xf", "MC signal", xf_nbins, 0, 0.8);







    TH1F *back_m = new TH1F("back_m", "back (WW, WZ, ZZ)", 30, 150, 2000);
    TH1F *back_xf = new TH1F("back_xf", "MC signal", xf_nbins, 0, 0.8);

    TH1F *DY_m = new TH1F("DY_m", "DY (WW, WZ, ZZ)", 30, 150, 2000);
    TH1F *DY_xf = new TH1F("DY_xf", "MC signal", xf_nbins, 0, 0.8);

    TH1F *QCD_m = new TH1F("QCD_m", "QCD", 30, 150, 2000);

    TH1F *WJets_m = new TH1F("WJets_m", "WJets", 30, 150, 2000);

    int n_cost_bins = 10;
    TH1F *data_cost = new TH1F("data_cost", "Data", n_cost_bins, -1.,1.);
    TH1F *back_cost = new TH1F("back_cost", "back (WW, WZ,ZZ)", n_cost_bins, -1.,1);
    TH1F *DY_cost = new TH1F("DY_cost", "DY (WW, WZ,ZZ)", n_cost_bins, -1.,1);
    TH1F *QCD_cost = new TH1F("QCD_cost", "QCD", n_cost_bins, -1.,1);
    TH1F *WJets_cost = new TH1F("WJets_cost", "WJets", n_cost_bins, -1.,1);

    DY_pt->SetFillColor(kRed+1);
    back_pt->SetFillColor(kGreen+3);
    QCD_pt->SetFillColor(kRed -7);

    DY_xf->SetFillColor(kRed+1);
    back_xf->SetFillColor(kGreen+3);
    QCD_xf->SetFillColor(kRed -7);


    DY_m->SetFillColor(kRed+1);
    DY_cost->SetFillColor(kRed+1);
    back_m->SetFillColor(kGreen+3);
    back_cost->SetFillColor(kGreen + 3);
    QCD_m->SetFillColor(kRed -7);
    QCD_cost->SetFillColor(kRed -7);

    bool do_RC = false;
    bool qcd_from_emu = false;
    float m_low = 150.;
    float m_high = 10000.;
    bool ss = true;
    bool in_os_region = false;

    make_m_cost_pt_xf_hist(t_elel_ss_data, data_m, data_cost, data_pt, data_xf, true, type,  do_RC, m_low, m_high, ss);
    make_m_cost_pt_xf_hist(t_elel_ss_back, back_m, back_cost, back_pt, back_xf, false, type,  do_RC, m_low, m_high, ss);
    make_m_cost_pt_xf_hist(t_elel_ss_dy, DY_m, DY_cost, DY_pt, DY_xf, false, type,  do_RC, m_low, m_high, ss);

    
    if(qcd_from_emu) make_qcd_from_emu_m_cost_pt_xf_hist(t_emu_ss_data, t_emu_ss_ttbar, t_emu_ss_diboson, t_emu_ss_dy, QCD_m, QCD_cost, QCD_pt, QCD_xf, m_low, m_high);
    

    else
        Fakerate_est_el(t_elel_WJets, t_elel_QCD, t_elel_WJets_contam, t_elel_QCD_contam, QCD_m, QCD_cost, QCD_pt, QCD_xf, m_low, m_high, ss, in_os_region);
   

    printf("Integrals of data, QCD, back, DY are %.2f %.2f %.2f %.2f \n", data_m->Integral(), QCD_m->Integral(), back_m->Integral(), DY_m->Integral());



    bool normalize = false;
    bool from_fit = false;
    
    if(normalize){
        Double_t n_data = data_m->Integral();
        Double_t n_mc = back_m->Integral() +  DY_m->Integral();
        Double_t n_QCD = QCD_m->Integral();
        Double_t qcd_ratio = (n_data - n_mc) / n_QCD;
        printf("Ratio of obs to expected QCD is %.2f \n", qcd_ratio);


        QCD_m->Scale(qcd_ratio);
        QCD_pt->Scale(qcd_ratio);
        QCD_cost->Scale(qcd_ratio);
        QCD_xf->Scale(qcd_ratio);
    }
    if(from_fit){
        float qcd_ratio, dy_ratio;
        if(qcd_from_emu){
            //qcd_ratio = 0.51;
            qcd_ratio = 0.8;
            dy_ratio = 1.4;
        }
        else{
            qcd_ratio = 0.85;
            dy_ratio = 1.15;
        }

        QCD_m->Scale(qcd_ratio);
        QCD_pt->Scale(qcd_ratio);
        QCD_cost->Scale(qcd_ratio);
        QCD_xf->Scale(qcd_ratio);

        DY_m->Scale(qcd_ratio);
        DY_pt->Scale(qcd_ratio);
        DY_cost->Scale(qcd_ratio);
        DY_xf->Scale(qcd_ratio);
    }


    bool scale_qcd_error=true;
    float qcd_err = 0.4;
    if(scale_qcd_error){
        int nBins_x = QCD_m->GetXaxis()->GetNbins();
        int nBins_y = QCD_cost->GetYaxis()->GetNbins();
        //printf("Get size %i \n", nBins);
        for (int i=1; i <= nBins_x; i++){
            for (int j=1; j <= nBins_y; j++){

                Double_t m_val = QCD_m->GetBinContent(i,j);
                Double_t cost_val = QCD_cost->GetBinContent(i,j);
                Double_t xf_val = QCD_xf->GetBinContent(i,j);

                QCD_m->SetBinError(i,j, qcd_err*m_val);
                QCD_cost->SetBinError(i,j, qcd_err*cost_val);
                QCD_xf->SetBinError(i,j, qcd_err*cost_val);
            }
        }
    }







    int nBins_x = QCD_m->GetXaxis()->GetNbins();
    int nBins_y = QCD_cost->GetYaxis()->GetNbins();
    //printf("Get size %i \n", nBins);
    /*
    for (int i=1; i <= nBins_x; i++){
        for (int j=1; j <= nBins_y; j++){

            Double_t m_val = QCD_m->GetBinContent(i,j);
            Double_t cost_val = QCD_cost->GetBinContent(i,j);

            QCD_m->SetBinError(i,j, 0.2*m_val);
            QCD_cost->SetBinError(i,j, 0.2*cost_val);
        }
    }
    */





    

    THStack *m_stack = new THStack("m_stack", "MuMu Mass Distribution: Data vs MC ; m_{e^{+}e^{-}} (GeV)");
    m_stack->Add(back_m);
    m_stack->Add(QCD_m);
    m_stack->Add(DY_m);


    THStack *cost_stack = new THStack("cost_stack", "Cos(#theta) Distribution: Data vs MC; ee Cos(#theta)_{r}");
    cost_stack->Add(back_cost);
    cost_stack->Add(QCD_cost);
    cost_stack->Add(DY_cost);

    THStack *pt_stack = new THStack("pt_stack", "Dimuon Pt Distribution: Data vs MC; Dielectron Pt (GeV)");
    pt_stack->Add(back_pt);
    pt_stack->Add(QCD_pt);
    pt_stack->Add(DY_pt);

    THStack *xf_stack = new THStack("xf_stack", "Dimuon x_F Distribution: Data vs MC; x_F");
    xf_stack->Add(back_xf);
    xf_stack->Add(QCD_xf);
    xf_stack->Add(DY_xf);

    TCanvas *c_m = new TCanvas("c_m", "Histograms", 200, 10, 900, 700);
    TPad *pad1 = new TPad("pad1", "pad1", 0.,0.3,0.98,1.);
    pad1->SetBottomMargin(0);
    pad1->Draw();
    pad1->cd();
    pad1->SetLogy();
    m_stack->Draw("hist");
    //m_stack->SetMaximum(65000);
    gStyle->SetEndErrorSize(4);
    data_m->SetMarkerStyle(kFullCircle);
    data_m->SetMarkerColor(1);
    data_m->DrawCopy("P E same");


    gStyle->SetLegendBorderSize(0);
    TLegend *leg1 = new TLegend(0.25, 0.25);
    leg1->AddEntry(data_m, "data", "p");
    leg1->AddEntry(DY_m, "DY (miss-sign)", "f");
    leg1->AddEntry(QCD_m, "QCD + WJets", "f");
    leg1->AddEntry(back_m, "t#bar{t} + wt + WW + WZ + ZZ", "f");
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
    ratio->SetMinimum(0.5);
    ratio->SetMaximum(1.5);
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
   ratio->GetYaxis()->SetTitle("Obs/Exp");
   ratio->GetYaxis()->SetNdivisions(505);
   ratio->GetYaxis()->SetTitleSize(20);
   ratio->GetYaxis()->SetTitleFont(43);
   ratio->GetYaxis()->SetTitleOffset(1.2);
   ratio->GetYaxis()->SetLabelFont(43); // Absolute font size in pixel (precision 3)
   ratio->GetYaxis()->SetLabelSize(15);
   // X axis ratio plot settings
   ratio->GetXaxis()->SetTitle("M_{ee} (GeV)");
   ratio->GetXaxis()->SetTitleSize(20);
   ratio->GetXaxis()->SetTitleFont(43);
   ratio->GetXaxis()->SetTitleOffset(3.);
   ratio->GetXaxis()->SetLabelFont(43); // Absolute font size in pixel (precision 3)
   ratio->GetXaxis()->SetLabelSize(20);
 
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
    cost_stack->SetMaximum(2000);
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
    cost_ratio->SetMinimum(0.5);
    cost_ratio->SetMaximum(1.5);
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
   cost_ratio->GetXaxis()->SetTitle("dielectron Cos(#theta)_{r}");
   cost_ratio->GetXaxis()->SetTitleSize(20);
   cost_ratio->GetXaxis()->SetTitleFont(43);
   cost_ratio->GetXaxis()->SetTitleOffset(3.);
   cost_ratio->GetXaxis()->SetLabelFont(43); // Absolute font size in pixel (precision 3)
   cost_ratio->GetXaxis()->SetLabelSize(20);

    //CMS_lumi(cost_pad1, iPeriod, 11);

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
    //pt_stack->SetMinimum(1);
    //pt_stack->SetMaximum(100000);
    //data_pt->SetMinimum(1);
    //data_pt->SetMaximum(100000);
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
    pt_ratio->SetMinimum(0.5);
    pt_ratio->SetMaximum(1.5);
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
    //xf_stack->SetMinimum(1);
    //xf_stack->SetMaximum(100000);
    //data_xf->SetMinimum(1);
    //data_xf->SetMaximum(100000);
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
    xf_ratio->SetMinimum(0.5);
    xf_ratio->SetMaximum(1.5);
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
   xf_ratio->GetXaxis()->SetTitle("dielectron xf (GeV)");
   xf_ratio->GetXaxis()->SetTitleSize(20);
   xf_ratio->GetXaxis()->SetTitleFont(43);
   xf_ratio->GetXaxis()->SetTitleOffset(3.);
   xf_ratio->GetXaxis()->SetLabelFont(43); // Absolute font size in pixel (precision 3)
   xf_ratio->GetXaxis()->SetLabelSize(20);
    //CMS_lumi(xf_pad1, iPeriod, 11 );
    c_xf->Update();
 
}

    
    
