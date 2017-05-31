

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



void draw_ttbar_frac(){
    TFile *f_data = TFile::Open("../analyze/output_files/DYToLL_data_2016_may9.root");
    TTree *t_data = (TTree *)f_data->Get("T_data");
    TFile *f_mc = TFile::Open("../analyze/output_files/DYToLL_mc_2016_may26.root");
    TTree *t_mc = (TTree *)f_mc->Get("T_data");
    TTree *t_mc_nosig = (TTree *)f_mc->Get("T_back");
    TFile *f_back = TFile::Open("../analyze/output_files/ttbar_background_may26.root");
    TTree *t_back = (TTree *)f_back->Get("T_data");



    int nBins = 6;

    Double_t m_bins[] = {150,200,250,350,500,700,10000};
    TH1F *mc_m = new TH1F("mc_m", "MC Signal (qqbar, qglu, qbarglu)", nBins, m_bins);
    TH1F *mc_cost = new TH1F("mc_cost", "MC Signal (qqbar, qglu, qbarglu)", 40, -1,1);

    TH1F *back_m = new TH1F("h_m", "TTBar Background", nBins, m_bins);
    TH1F *back_cost = new TH1F("back_cost", "TTbar", 40, -1,1);


    make_m_cost_hist(t_mc, mc_m, mc_cost);
    make_m_cost_hist(t_back, back_m, back_cost);

    Double_t ttbar_frac[nBins], bin_center[nBins];
    for (int i=1; i <= nBins; i++){
        Double_t N_m = mc_m->GetBinContent(i);
        Double_t N_ttbar = back_m->GetBinContent(i);
        bin_center[i-1] = mc_m->GetBinCenter(i);
        printf("bin center %f \n", bin_center[i-1]);
        Double_t scaling = 1.0;
        ttbar_frac[i-1] = scaling*N_ttbar/(N_m + N_ttbar);
    }
    bin_center[nBins-1] = 850;
    Double_t fit_res[] = {0.170, 0.148, 0.173, 0.157, 0.137, 0.096};
    Double_t fit_errs[] = {0.008, 0.012, 0.011, 0.015, 0.024, 0.032};
    TGraph *mc_frac = new TGraph(nBins, bin_center, ttbar_frac);
    mc_frac->SetTitle("Fraction from MC");
    TGraphErrors *fit_frac = new TGraphErrors(nBins, bin_center, fit_res, 0, fit_errs);
    fit_frac->SetTitle("Fraction from fit results");
    fit_frac->SetMaximum(0.2);
    fit_frac->SetMinimum(0.0);
    mc_frac->SetMaximum(0.2);
    mc_frac->SetMinimum(0.0);
    TCanvas *c3 = new TCanvas("c3", "canvas", 200,10, 900,700);
    fit_frac->SetMarkerStyle(kFullSquare);
    fit_frac->SetLineWidth(2);
    c3->Update();
    mc_frac->SetMarkerStyle(kFullSquare);
    mc_frac->SetMarkerColor(kBlue);
    mc_frac->SetLineWidth(3);
    mc_frac->SetLineColor(kBlue);


    TMultiGraph *mg = new TMultiGraph();
    mg->Add(mc_frac);
    mg->Add(fit_frac);

    mg->SetTitle("Fraction of ttbar events");
    

    mg->Draw("A C P");

    mg->GetXaxis()->SetTitle("M (GeV)");
    mg->GetXaxis()->CenterTitle();
    mg->GetYaxis()->SetTitle("R_{ttbar}");
    mg->GetYaxis()->CenterTitle();


    c3->Update();

    gPad->BuildLegend();
    

    return;
}


