

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

bool has_no_bjets(Int_t nJets, Double_t jet1_pt, Double_t jet2_pt, 
        Double_t jet1_cmva, Double_t jet2_cmva){
    Double_t med_btag = 0.4432;
    if(nJets ==0) return true;
    else if(nJets == 1){
        if(jet1_pt < 20.) return true;
        else return jet1_cmva < med_btag;
    }
    else{
        return (jet1_pt < 20. || jet1_cmva < med_btag) && (jet2_pt < 20. || jet2_cmva < med_btag);
    }
}

void read_tree(TTree *t1, TH1F *h_m,  bool is_data=false){
    //read event data
    Long64_t size  =  t1->GetEntries();
    Double_t m, xF, cost, mu1_pt, mu2_pt, jet1_cmva, jet2_cmva, gen_weight;
    Double_t bcdef_HLT_SF, bcdef_iso_SF, bcdef_id_SF;
    Double_t jet1_pt, jet2_pt, jet1_b_weight, jet2_b_weight;
    Float_t cost_pt, met_pt;
    Int_t nJets;
    nJets = 2;
    t1->SetBranchAddress("m", &m);
    t1->SetBranchAddress("xF", &xF);
    t1->SetBranchAddress("cost", &cost);
    t1->SetBranchAddress("met_pt", &met_pt);
    t1->SetBranchAddress("jet2_CMVA", &jet2_cmva);
    t1->SetBranchAddress("jet1_CMVA", &jet1_cmva);
    t1->SetBranchAddress("jet1_pt", &jet1_pt);
    t1->SetBranchAddress("jet2_pt", &jet2_pt);
    if(!is_data){
        t1->SetBranchAddress("nJets", &nJets);
        t1->SetBranchAddress("gen_weight", &gen_weight);
        t1->SetBranchAddress("bcdef_HLT_SF", &bcdef_HLT_SF);
        t1->SetBranchAddress("bcdef_iso_SF", &bcdef_iso_SF);
        t1->SetBranchAddress("bcdef_id_SF", &bcdef_id_SF);
        t1->SetBranchAddress("jet1_b_weight", &jet1_b_weight);
        t1->SetBranchAddress("jet2_b_weight", &jet2_b_weight);
    }

    for (int i=0; i<size; i++) {
        t1->GetEntry(i);
        bool no_bjets = has_no_bjets(nJets, jet1_pt, jet2_pt, jet1_cmva, jet2_cmva);

        if(m >= 150. && met_pt < 50. && no_bjets){
            if(is_data){
                h_m->Fill(m);
            }
            else{
                Double_t weight = gen_weight * bcdef_HLT_SF * bcdef_iso_SF * bcdef_id_SF;
                if (nJets >= 1) weight *= jet1_b_weight;
                if (nJets >= 2) weight *= jet2_b_weight;
                //Double_t weight = gen_weight;
                h_m->Fill(m,weight);
            }

            
        }
    }
}

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

    TH1F *back_m = new TH1F("h_m", "TTBar Background", nBins, m_bins);


    read_tree(t_mc, mc_m);
    read_tree(t_back, back_m);

    Double_t ttbar_frac[nBins], bin_center[nBins];
    for (int i=1; i <= nBins; i++){
        Double_t N_m = mc_m->GetBinContent(i);
        Double_t N_ttbar = back_m->GetBinContent(i);
        bin_center[i-1] = mc_m->GetBinCenter(i);
        printf("bin center %f \n", bin_center[i-1]);
        ttbar_frac[i-1] = N_ttbar/(N_m + N_ttbar);
    }
    bin_center[nBins-1] = 850;
    Double_t fit_res[] = {0.122, 0.125, 0.158, 0.153, 0.144, 0.096};
    TGraph *mc_frac = new TGraph(nBins, bin_center, ttbar_frac);
    mc_frac->SetTitle("Fraction from MC");
    TGraph *fit_frac = new TGraph(nBins, bin_center, fit_res);
    fit_frac->SetTitle("Fraction from fit results");
    fit_frac->SetMaximum(0.2);
    fit_frac->SetMinimum(0.0);
    mc_frac->SetMaximum(0.2);
    mc_frac->SetMinimum(0.0);
    TCanvas *c3 = new TCanvas("c3", "canvas", 200,10, 900,700);
    fit_frac->SetMarkerStyle(kFullSquare);
    fit_frac->SetLineWidth(3);
    c3->Update();
    mc_frac->SetMarkerStyle(kFullSquare);
    mc_frac->SetMarkerColor(kBlue);
    mc_frac->SetLineWidth(3);


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


