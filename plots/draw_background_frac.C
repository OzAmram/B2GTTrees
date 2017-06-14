

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



void draw_background_frac(){
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

    int nBins = 6;

    Double_t m_bins[] = {150,200,250,350,500,700,10000};
    TH1F *mc_m = new TH1F("mc_m", "MC Signal (qqbar, qglu, qbarglu)", nBins, m_bins);
    TH1F *mc_nosig_m = new TH1F("mc_nosig_m", "MC no Asymmetry production (qq, qbarqbar, gluglu)", nBins, m_bins);
    TH1F *mc_cost = new TH1F("mc_cost", "MC Signal (qqbar, qglu, qbarglu)", 40, -1,1);
    TH1F *mc_nosig_cost = new TH1F("mc_nosig_cost", "MC no Asymmetry production (qq, qbarqbar, gluglu)", 40, -1,1);

    TH1F *ttbar_m = new TH1F("h_m", "TTBar Background", nBins, m_bins);
    TH1F *ttbar_cost = new TH1F("back_cost", "TTbar", 40, -1,1);

    TH1F *ww_m = new TH1F("ww_m", "MC Signal (qqbar, qglu, qbarglu)", nBins, m_bins);
    TH1F *ww_cost = new TH1F("ww_cost", "MC Signal (qqbar, qglu, qbarglu)", 40, -1,1);
    TH1F *wz_m = new TH1F("wz_m", "MC Signal (qqbar, qglu, qbarglu)", nBins, m_bins);
    TH1F *wz_cost = new TH1F("wz_cost", "MC Signal (qqbar, qglu, qbarglu)", 40, -1,1);
    TH1F *zz_m = new TH1F("zz_m", "MC Signal (qqbar, qglu, qbarglu)", nBins, m_bins);
    TH1F *zz_cost = new TH1F("zz_cost", "MC Signal (qqbar, qglu, qbarglu)", 40, -1,1);


    make_m_cost_hist(t_mc, mc_m, mc_cost, false);
    make_m_cost_hist(t_mc_nosig, mc_nosig_m, mc_nosig_cost, false);
    make_m_cost_hist(t_ttbar, ttbar_m, ttbar_cost, false);
    make_m_cost_hist(t_ww, ww_m, ww_cost, false);
    make_m_cost_hist(t_wz, wz_m, wz_cost, false);
    make_m_cost_hist(t_zz, zz_m, zz_cost, false);

    Double_t ttbar_frac[nBins], ttbar_frac_unc[nBins], diboson_frac[nBins], diboson_frac_unc[nBins], bin_center[nBins];
    Double_t back_frac[nBins], back_frac_unc[nBins], nosig_frac[nBins], nosig_frac_unc[nBins];
    for (int i=1; i <= nBins; i++){
        Double_t N_mc = mc_m->GetBinContent(i);
        Double_t N_mc_nosig = mc_nosig_m->GetBinContent(i);
        Double_t N_ttbar = ttbar_m->GetBinContent(i);
        Double_t N_ww = ww_m->GetBinContent(i);
        Double_t N_wz = wz_m->GetBinContent(i);
        Double_t N_zz = wz_m->GetBinContent(i);
        Double_t N_diboson = N_ww+ + N_wz+ N_zz;
        Double_t N_back = N_ttbar + N_diboson;
        Double_t denom = N_ttbar + N_diboson + N_mc + N_mc_nosig;
        bin_center[i-1] = mc_m->GetBinCenter(i);
        printf("bin center %f \n", bin_center[i-1]);
        Double_t scaling = 1.24;
        nosig_frac[i-1] = N_mc_nosig/denom;
        nosig_frac_unc[i-1] = std::sqrt(  N_mc_nosig*pow((1/denom - N_mc_nosig/pow(denom,2)), 2) +
                                          (N_mc + N_ttbar + N_diboson) * pow(1/denom, 4));
        diboson_frac[i-1] = N_diboson/(denom);
        diboson_frac_unc[i-1] = std::sqrt(  N_diboson*pow((1/denom - N_diboson/pow(denom,2)), 2) +
                                          (N_mc + N_ttbar + N_mc_nosig) * pow(1/denom, 4));
        ttbar_frac[i-1] =  scaling*(N_ttbar)/(denom);
        ttbar_frac_unc[i-1] = scaling*std::sqrt(  N_ttbar*pow((1/denom - N_ttbar/pow(denom,2)), 2) +
                                          (N_mc + N_diboson + N_mc_nosig) * pow(1/denom, 4));
        back_frac[i-1]  = diboson_frac[i-1] + ttbar_frac[i-1];
        back_frac_unc[i-1] = std::sqrt(pow(ttbar_frac_unc[i-1],2) + pow(diboson_frac_unc[i-1], 2));
    }
    bin_center[nBins-1] = 850;
    Double_t fit_res[] = {0.192, 0.158, 0.195, 0.162, 0.152, 0.095};
    Double_t fit_errs[] = {0.008, 0.011, 0.009, 0.013, 0.024, 0.029};
    TGraphErrors *mc_nosig_frac = new TGraphErrors(nBins, bin_center, nosig_frac, 0, nosig_frac_unc);
    mc_nosig_frac->SetTitle("MC no asym events (qq, gluglu, qbarqbar)");
    TGraphErrors *ttbar_mc_frac = new TGraphErrors(nBins, bin_center, ttbar_frac, 0, ttbar_frac_unc);
    ttbar_mc_frac->SetTitle("ttbar fraction from MC (scaled with EMu) ");
    TGraphErrors *diboson_mc_frac = new TGraphErrors(nBins, bin_center, diboson_frac, 0, diboson_frac_unc);
    diboson_mc_frac->SetTitle("diboson fraction from MC ");
    TGraphErrors *back_mc_frac = new TGraphErrors(nBins, bin_center, back_frac, 0, back_frac_unc);
    back_mc_frac->SetTitle("Total background fraction from MC ");
    TGraphErrors *fit_frac = new TGraphErrors(nBins, bin_center, fit_res, 0, fit_errs);
    fit_frac->SetTitle("Fraction of background events from fit results");
    fit_frac->SetMaximum(0.2);
    fit_frac->SetMinimum(0.0);
    ttbar_mc_frac->SetMaximum(0.2);
    ttbar_mc_frac->SetMinimum(0.0);
    diboson_mc_frac->SetMaximum(0.2);
    diboson_mc_frac->SetMinimum(0.0);
    mc_nosig_frac->SetMaximum(0.2);
    mc_nosig_frac->SetMinimum(0.0);
    back_mc_frac->SetMaximum(0.2);
    back_mc_frac->SetMinimum(0.0);
    TCanvas *c3 = new TCanvas("c3", "canvas", 200,10, 900,700);
    fit_frac->SetMarkerStyle(kFullSquare);
    fit_frac->SetLineWidth(2);
    c3->Update();
    ttbar_mc_frac->SetMarkerStyle(kFullSquare);
    ttbar_mc_frac->SetMarkerColor(kBlue);
    ttbar_mc_frac->SetLineWidth(3);
    ttbar_mc_frac->SetLineColor(kBlue);


    diboson_mc_frac->SetMarkerStyle(kFullSquare);
    diboson_mc_frac->SetMarkerColor(kGreen);
    diboson_mc_frac->SetLineWidth(3);
    diboson_mc_frac->SetLineColor(kGreen);

    back_mc_frac->SetMarkerStyle(kFullSquare);
    back_mc_frac->SetMarkerColor(kRed);
    back_mc_frac->SetLineWidth(3);
    back_mc_frac->SetLineColor(kRed);

    mc_nosig_frac->SetMarkerStyle(kFullSquare);
    mc_nosig_frac->SetMarkerColor(kYellow);
    mc_nosig_frac->SetLineWidth(3);
    mc_nosig_frac->SetLineColor(kYellow);
    TMultiGraph *mg = new TMultiGraph();
    mg->Add(ttbar_mc_frac);
    mg->Add(diboson_mc_frac);
    mg->Add(back_mc_frac);
    mg->Add(mc_nosig_frac);
    mg->Add(fit_frac);

    mg->SetTitle("Fraction of background events");
    

    mg->Draw("A C P");

    mg->GetXaxis()->SetTitle("M (GeV)");
    mg->GetXaxis()->CenterTitle();
    mg->GetYaxis()->SetTitle("R_{ttbar}");
    mg->GetYaxis()->CenterTitle();


    c3->Update();

    gPad->BuildLegend();
    

    return;
}


