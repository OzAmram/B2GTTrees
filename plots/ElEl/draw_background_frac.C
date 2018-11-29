

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



void draw_background_frac(){
    setTDRStyle();
    init();
    int nBins = 6;

    Double_t m_bins[] = {150,200,250,350,500,700, 10000};
    TH1F *mc_m = new TH1F("mc_m", "MC Signal (qqbar, qglu, qbarglu)", nBins, m_bins);
    TH1F *mc_nosig_m = new TH1F("mc_nosig_m", "MC no Asymmetry production (qq, qbarqbar, gluglu)", nBins, m_bins);
    TH1F *mc_cost = new TH1F("mc_cost", "MC Signal (qqbar, qglu, qbarglu)", 40, -1,1);
    TH1F *mc_nosig_cost = new TH1F("mc_nosig_cost", "MC no Asymmetry production (qq, qbarqbar, gluglu)", 40, -1,1);

    TH1F *ttbar_m = new TH1F("h_m", "TTBar Background", nBins, m_bins);
    TH1F *ttbar_cost = new TH1F("back_cost", "TTbar", 40, -1,1);

    TH1F *diboson_m = new TH1F("diboson_m", "MC Signal (qqbar, qglu, qbarglu)", nBins, m_bins);
    TH1F *QCD_m = new TH1F("QCD_m", "MC Signal (qqbar, qglu, qbarglu)", nBins, m_bins);
    TH1F *diboson_cost = new TH1F("diboson_cost", "MC Signal (qqbar, qglu, qbarglu)", 40, -1,1);
    TH1F *QCD_cost = new TH1F("QCD_cost", "MC Signal (qqbar, qglu, qbarglu)", 40, -1,1);


    TH1F *mc_pt = new TH1F("mc_pt", "MC signal", 40, 0, 1000);
    TH1F *mc_nosig_pt = new TH1F("mc_nosig_pt", "MC signal", 40, 0, 1000);
    TH1F *data_pt = new TH1F("data_pt", "MC signal", 40, 0, 1000);
    TH1F *ttbar_pt = new TH1F("ttbar_pt", "MC signal", 40, 0, 1000);
    TH1F *diboson_pt = new TH1F("diboson_pt", "MC signal", 40, 0, 1000);
    TH1F *wt_pt = new TH1F("wt_pt", "MC signal", 40, 0, 1000);
    TH1F *QCD_pt = new TH1F("QCD_pt", "MC signal", 40, 0, 1000);

    TH1F *mc_xf = new TH1F("mc_xf", "MC signal", 40, 0, 1000);
    TH1F *mc_nosig_xf = new TH1F("mc_nosig_xf", "MC signal", 40, 0, 1000);
    TH1F *data_xf = new TH1F("data_xf", "MC signal", 40, 0, 1000);
    TH1F *ttbar_xf = new TH1F("ttbar_xf", "MC signal", 40, 0, 1000);
    TH1F *diboson_xf = new TH1F("diboson_xf", "MC signal", 40, 0, 1000);
    TH1F *wt_xf = new TH1F("wt_xf", "MC signal", 40, 0, 1000);
    TH1F *QCD_xf = new TH1F("QCD_xf", "MC signal", 40, 0, 1000);
    
    /*
    TH1F *ww_m = new TH1F("ww_m", "MC Signal (qqbar, qglu, qbarglu)", nBins, m_bins);
    TH1F *ww_cost = new TH1F("ww_cost", "MC Signal (qqbar, qglu, qbarglu)", 40, -1,1);
    TH1F *wz_m = new TH1F("wz_m", "MC Signal (qqbar, qglu, qbarglu)", nBins, m_bins);
    TH1F *wz_cost = new TH1F("wz_cost", "MC Signal (qqbar, qglu, qbarglu)", 40, -1,1);
    TH1F *zz_m = new TH1F("zz_m", "MC Signal (qqbar, qglu, qbarglu)", nBins, m_bins);
    TH1F *zz_cost = new TH1F("zz_cost", "MC Signal (qqbar, qglu, qbarglu)", 40, -1,1);
    */

    TH1F *wt_m = new TH1F("wt_m", "W top", nBins, m_bins);
    TH1F *wt_cost = new TH1F("wt_cost", "W top", 40, -1,1);

    int type = FLAG_ELECTRONS;
    make_m_cost_pt_xf_hist(t_mc, mc_m, mc_cost, mc_pt, mc_xf, false,type, false);
    make_m_cost_pt_xf_hist(t_mc_nosig, mc_nosig_m, mc_nosig_pt, mc_nosig_cost, mc_nosig_xf, false,type, false);
    make_m_cost_pt_xf_hist(t_ttbar, ttbar_m, ttbar_cost, ttbar_pt, ttbar_xf, false,type, true);
    make_m_cost_pt_xf_hist(t_diboson, diboson_m, diboson_cost, diboson_pt, diboson_xf, false,type, true);
    make_m_cost_pt_xf_hist(t_wt, wt_m, wt_cost, wt_pt, wt_xf, false,type, true);
    

    Fakerate_est_el(t_WJets, t_QCD, t_WJets_mc, t_QCD_mc, QCD_m, QCD_cost, QCD_pt, QCD_xf);


    /*
     
    make_m_cost_hist(t_ww, ww_m, ww_cost, false);
    make_m_cost_hist(t_wz, wz_m, wz_cost, false);
    make_m_cost_hist(t_zz, zz_m, zz_cost, false);
    make_m_cost_hist(t_wt, wt_m, wt_cost, false);
    */
        
    

    Double_t ttbar_frac[nBins], ttbar_frac_unc[nBins], diboson_frac[nBins], diboson_frac_unc[nBins], bin_center[nBins];
    Double_t back_frac[nBins], back_frac_unc[nBins], nosig_frac[nBins], nosig_frac_unc[nBins];
    Double_t wt_frac[nBins], wt_frac_unc[nBins], QCD_frac[nBins], QCD_frac_unc[nBins];
    for (int i=1; i <= nBins; i++){
        Double_t N_mc = mc_m->GetBinContent(i);
        Double_t N_mc_nosig = mc_nosig_m->GetBinContent(i);
        Double_t N_ttbar =  ttbar_m->GetBinContent(i);
        Double_t N_diboson = diboson_m->GetBinContent(i);
        Double_t N_QCD = QCD_m->GetBinContent(i);
        Double_t N_wt = wt_m->GetBinContent(i);
        Double_t denom = N_ttbar + N_diboson + N_mc + N_wt + N_mc_nosig + N_QCD;
        bin_center[i-1] = mc_m->GetBinCenter(i);
        printf("bin center %f \n", bin_center[i-1]);
        nosig_frac[i-1] = N_mc_nosig/denom;
        //nosig_frac_unc[i-1] = std::sqrt(  N_mc_nosig*pow((1/denom - N_mc_nosig/pow(denom,2)), 2) +
                                          //(denom - N_mc_nosig) * pow(1/denom, 4));
        diboson_frac[i-1] = N_diboson/(denom);
        //diboson_frac_unc[i-1] = std::sqrt(  N_diboson*pow((1/denom - N_diboson/pow(denom,2)), 2) +
                                          //(denom - N_diboson) * pow(1/denom, 4));
        ttbar_frac[i-1] =  (N_ttbar)/(denom);
        //ttbar_frac_unc[i-1] = std::sqrt(  N_ttbar*pow((1/denom - N_ttbar/pow(denom,2)), 2) +
                                          //(denom - N_ttbar) * pow(1/denom, 4));
        wt_frac[i-1] =  (N_wt)/(denom);
        //wt_frac_unc[i-1] = std::sqrt(  N_wt*pow((1/denom - N_wt/pow(denom,2)), 2) +
                                          //(denom - N_wt) * pow(1/denom, 4));
        QCD_frac[i-1] =  (N_QCD)/(denom);
        QCD_frac_unc[i-1]  = 0.3 * (N_QCD) / (denom);
                                          
        back_frac[i-1]  = diboson_frac[i-1] + ttbar_frac[i-1] + nosig_frac[i-1] + wt_frac[i-1] + QCD_frac[i-1];
        back_frac_unc[i-1] = sqrt(pow(0.02*(back_frac[i-1]-QCD_frac[i-1]),2) + pow(QCD_frac_unc[i-1],2));
        //back_frac_unc[i-1] = std::sqrt(pow(ttbar_frac_unc[i-1],2) + pow(diboson_frac_unc[i-1], 2) + 
                                       //pow(nosig_frac_unc[i-1],2) + pow(wt_frac_unc[i-1], 2) +
                                       //pow(QCD_frac_unc[i-1],2) );
    }
    bin_center[nBins-1] = 800.;
    Double_t fit_res[] = {0.079, 0.147, 0.212, 0.175, 0.170, 0.137};
    Double_t fit_errs[] = {0.005, 0.008, 0.010, 0.013, 0.022, 0.035};

    Double_t comb_fit_res[] = {0.078, 0.147, 0.212, 0.176, 0.170, 0.139};
    Double_t comb_fit_errs[] = {0.005, 0.008, 0.010, 0.013, 0.022, 0.035};

    TGraphErrors *mc_nosig_frac = new TGraphErrors(nBins, bin_center, nosig_frac, 0, 0);
    mc_nosig_frac->SetTitle("MC no asym events (qq, gluglu, qbarqbar)");

    TGraphErrors *ttbar_mc_frac = new TGraphErrors(nBins, bin_center, ttbar_frac, 0, 0);
    ttbar_mc_frac->SetTitle("ttbar fraction from MC (scaled with e#mu) ");

    TGraphErrors *wt_mc_frac = new TGraphErrors(nBins, bin_center, wt_frac, 0, 0);
    wt_mc_frac->SetTitle("W top fraction from MC");

    TGraphErrors *diboson_mc_frac = new TGraphErrors(nBins, bin_center, diboson_frac, 0, 0);
    diboson_mc_frac->SetTitle("diboson fraction from MC ");

    TGraphErrors *QCD_est_frac = new TGraphErrors(nBins, bin_center, QCD_frac, 0, QCD_frac_unc);
    diboson_mc_frac->SetTitle("QCD fraction from Fakerate Est.");

    TGraphErrors *back_mc_frac = new TGraphErrors(nBins, bin_center, back_frac, 0, back_frac_unc);
    back_mc_frac->SetTitle("Total background fraction from MC ");

    TGraphErrors *fit_frac = new TGraphErrors(nBins, bin_center, fit_res, 0, fit_errs);
    fit_frac->SetTitle("Fraction of background events from ee fit");

    TGraphErrors *comb_fit_frac = new TGraphErrors(nBins, bin_center, comb_fit_res, 0, comb_fit_errs);
    comb_fit_frac->SetTitle("Fraction of background events from combined fit");

    fit_frac->SetMaximum(0.4);
    fit_frac->SetMinimum(0.0);
    comb_fit_frac->SetMaximum(0.4);
    comb_fit_frac->SetMinimum(0.0);
    ttbar_mc_frac->SetMaximum(0.4);
    ttbar_mc_frac->SetMinimum(0.0);
    wt_mc_frac->SetMaximum(0.4);
    wt_mc_frac->SetMinimum(0.0);
    diboson_mc_frac->SetMaximum(0.4);
    diboson_mc_frac->SetMinimum(0.0);
    mc_nosig_frac->SetMaximum(0.4);
    mc_nosig_frac->SetMinimum(0.0);
    back_mc_frac->SetMaximum(0.4);
    back_mc_frac->SetMinimum(0.0);

    QCD_est_frac->SetMaximum(0.4);
    QCD_est_frac->SetMinimum(0.0);

    TCanvas *c3 = new TCanvas("c3", "canvas", 200,10, 900,700);
    fit_frac->SetMarkerStyle(kFullSquare);
    //fit_frac->SetLineWidth(2);

    comb_fit_frac->SetMarkerStyle(kFullTriangleUp);
    comb_fit_frac->SetMarkerColor(kRed);
    //comb_fit_frac->SetLineWidth(2);
    //comb_fit_frac->SetLineColor(2);

    c3->Update();
    ttbar_mc_frac->SetMarkerStyle(kFullSquare);
    ttbar_mc_frac->SetMarkerColor(kBlue);
    ttbar_mc_frac->SetLineWidth(3);
    ttbar_mc_frac->SetLineColor(kBlue);


    wt_mc_frac->SetMarkerStyle(kFullSquare);
    wt_mc_frac->SetMarkerColor(kOrange+7);
    wt_mc_frac->SetLineWidth(3);
    wt_mc_frac->SetLineColor(kOrange+7);

    diboson_mc_frac->SetMarkerStyle(kFullSquare);
    diboson_mc_frac->SetMarkerColor(kGreen +3);
    diboson_mc_frac->SetLineWidth(3);
    diboson_mc_frac->SetLineColor(kGreen +3);

    back_mc_frac->SetMarkerStyle(kFullSquare);
    back_mc_frac->SetMarkerColor(kMagenta +3);
    back_mc_frac->SetLineWidth(3);
    back_mc_frac->SetLineColor(kMagenta +3);

    mc_nosig_frac->SetMarkerStyle(kFullSquare);
    mc_nosig_frac->SetMarkerColor(kMagenta);
    mc_nosig_frac->SetLineWidth(3);
    mc_nosig_frac->SetLineColor(kMagenta);

    QCD_est_frac->SetMarkerStyle(kFullSquare);
    QCD_est_frac->SetMarkerColor(kRed - 7);
    QCD_est_frac->SetLineWidth(3);
    QCD_est_frac->SetLineColor(kRed -7);

    TMultiGraph *mg = new TMultiGraph();
    mg->Add(ttbar_mc_frac);
    mg->Add(wt_mc_frac);
    mg->Add(diboson_mc_frac);
    mg->Add(QCD_est_frac);
    mg->Add(back_mc_frac);
    mg->Add(mc_nosig_frac);
    mg->Add(fit_frac);
    mg->Add(comb_fit_frac);

    mg->SetTitle("Fraction of background events (Nominal e#mu scaling)");
    

    mg->Draw("A C P");

    mg->GetXaxis()->SetTitle("M_{ee} (GeV)");
    mg->GetYaxis()->SetTitle("Fraction of selected events");
    mg->GetYaxis()->SetRangeUser(0, 0.4);

    gStyle->SetLegendBorderSize(0);
    TLegend *leg2 = new TLegend(0.5, 0.65, 0.75, 0.8);
    leg2->AddEntry(fit_frac, "R_{bk} value from ee fit", "p");
    leg2->AddEntry(comb_fit_frac, "R_{bk} value from combined fit", "p");
    leg2->AddEntry(back_mc_frac, "Total background fraction from MC + QCD + WJets", "l");
    leg2->AddEntry(QCD_est_frac, "QCD + WJets", "l");
    leg2->AddEntry(ttbar_mc_frac, "t#bar{t}", "l");
    leg2->AddEntry(mc_nosig_frac, "DY no asymmety(gg, qq, #bar{q}#bar{q})", "l");
    leg2->AddEntry(diboson_mc_frac, "WW + WZ + ZZ", "l");
    leg2->AddEntry(wt_mc_frac, "tW + #bar{t}W", "l");
    leg2->Draw();

    writeExtraText = true;
    //lumi_sqrtS = "";       // used with iPeriod = 0, e.g. for simulation-only plots (default is an empty string)
    int iPeriod = 4; 
    CMS_lumi(c3, iPeriod, 33 );
    c3->Update();

    //gPad->BuildLegend();
    

    return;
}


