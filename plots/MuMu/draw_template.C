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
#include "../analyze/HistMaker.C"


int n_xf_bins = 5;
float xf_max = 1.0;
Float_t xf_bins[] = {0., 0.05, 0.1, 0.15, 0.25, 1.0};
int n_m_bins = 1;
float m_max = 1000.;
int n_cost_bins = 10;
Float_t cost_bins[] = {-1.0, -.8, -.6, -.4, -.2, 0., 0.2, 0.4, 0.6, 0.8, 1.0};

int gen_mc_template_p(TTree *t1, TH2F* h_sym, TH2F *h_asym, TH2F *h_count, Double_t m_low, Double_t m_high){
    Long64_t nEntries  =  t1->GetEntries();
    //printf("size is %i \n", nEntries);
    Double_t m, xF, cost, gen_weight, reweight, jet1_cmva, jet2_cmva, cost_st;
    Double_t bcdef_HLT_SF, bcdef_iso_SF, bcdef_id_SF;
    Double_t gh_HLT_SF, gh_iso_SF, gh_id_SF;
    Double_t jet1_pt, jet2_pt, jet1_b_weight, jet2_b_weight;
    Float_t cost_pt, met_pt;
    Int_t nJets;
    t1->SetBranchAddress("m", &m);
    t1->SetBranchAddress("xF", &xF);
    t1->SetBranchAddress("cost", &cost);
    t1->SetBranchAddress("cost_st", &cost_st);
    t1->SetBranchAddress("jet1_CMVA", &jet1_cmva);
    t1->SetBranchAddress("jet2_CMVA", &jet2_cmva);
    t1->SetBranchAddress("met_pt", &met_pt);
    t1->SetBranchAddress("jet1_pt", &jet1_pt);
    t1->SetBranchAddress("jet2_pt", &jet2_pt);
    t1->SetBranchAddress("nJets", &nJets);
    t1->SetBranchAddress("gen_weight", &gen_weight);
    t1->SetBranchAddress("bcdef_HLT_SF", &bcdef_HLT_SF);
    t1->SetBranchAddress("bcdef_iso_SF", &bcdef_iso_SF);
    t1->SetBranchAddress("bcdef_id_SF", &bcdef_id_SF);
    t1->SetBranchAddress("gh_HLT_SF", &gh_HLT_SF);
    t1->SetBranchAddress("gh_iso_SF", &gh_iso_SF);
    t1->SetBranchAddress("gh_id_SF", &gh_id_SF);
    t1->SetBranchAddress("jet1_b_weight", &jet1_b_weight);
    t1->SetBranchAddress("jet2_b_weight", &jet2_b_weight);

    h_sym->Sumw2();
    h_asym->Sumw2();
    TH2D *h_sym_bcdef = (TH2D *) h_sym->Clone("h_sym_bcdef");
    TH2D *h_sym_gh = (TH2D *)h_sym->Clone("h_sym_gh");
    TH2D *h_asym_bcdef = (TH2D *)h_asym->Clone("h_asym_bcdef");
    TH2D *h_asym_gh = (TH2D *)h_asym->Clone("h_asym_gh");
    int n = 0;
    for (int i=0; i<nEntries; i++) {
        t1->GetEntry(i);
        bool no_bjets = has_no_bjets(nJets, jet1_pt, jet2_pt, jet1_cmva, jet2_cmva);
        if(m >= m_low && m <= m_high && met_pt < 50.  && no_bjets && xF > 0.25){
            reweight = (4./3.)*cost_st*(2. + alpha)/
                (1. + cost_st*cost_st + alpha*(1.- cost_st*cost_st));
            n++;


            Double_t bcdef_weight = gen_weight * bcdef_HLT_SF * bcdef_iso_SF * bcdef_id_SF;
            Double_t gh_weight = gen_weight * gh_HLT_SF * gh_iso_SF * gh_id_SF;
            if (nJets >= 1){
                bcdef_weight *= jet1_b_weight;
                gh_weight *= jet1_b_weight;
            }
            if (nJets >= 2){
                bcdef_weight *= jet2_b_weight;
                gh_weight *= jet2_b_weight;
            }


            h_sym_bcdef->Fill(xF, cost, bcdef_weight); 
            h_sym_bcdef->Fill(xF, -cost, bcdef_weight); 
            h_sym_gh->Fill(xF, cost, gh_weight); 
            h_sym_gh->Fill(xF, -cost, gh_weight); 

            h_asym_bcdef->Fill(xF, cost, reweight * bcdef_weight);
            h_asym_bcdef->Fill(xF, -cost, -reweight * bcdef_weight);
            h_asym_gh->Fill(xF, cost, reweight * gh_weight);
            h_asym_gh->Fill(xF, -cost, -reweight * gh_weight);

            h_count->Fill(xF, cost, 1);
        }
    }

    h_sym_bcdef->Scale(1000*bcdef_lumi);
    h_sym_gh->Scale(1000*gh_lumi);
    h_asym_bcdef->Scale(1000*bcdef_lumi);
    h_asym_gh->Scale(1000*gh_lumi);

    h_sym->Add(h_sym_bcdef, h_sym_gh);
    h_asym->Add(h_asym_bcdef, h_asym_gh);
    printf("N sym is %i \n", n);
    float norm = h_sym -> Integral();
    h_sym->Scale(1./norm);
    h_asym->Scale(1./norm);
    t1->ResetBranchAddresses();
    return 0;
}

void draw_template(){
    TFile* f_mc = (TFile*) TFile::Open("../analyze/output_files/DYToLL_mc_2016_jun13.root");
    TTree *t_mc = (TTree *) f_mc ->Get("T_data");

    TH2F * h_sym = new TH2F("h_sym", "Symmetric template of mc (xF, cost_r) xF>0.15",
            n_xf_bins, xf_bins, n_cost_bins, cost_bins);
    TH2F *h_asym = new TH2F("h_asym", "Asymmetric template of mc (xF cost_r) xF > 0.15",
            n_xf_bins, xf_bins, n_cost_bins, cost_bins);

    TH2F *h_sym_count = new TH2F("h_sym_count", "Events in bins for MC templates",
            n_xf_bins, xf_bins, n_cost_bins, cost_bins);

  
    Double_t m_low = 200;
    Double_t m_high = 250;
   
    gen_mc_template_p(t_mc, h_sym, h_asym, h_sym_count, m_low, m_high);
    h_sym->Print();

    TH1D *h_sym_cost = h_sym->ProjectionY("h_sym_cost");
    TH1D *h_asym_cost = h_asym->ProjectionY("h_asym_cost");
    TH1D *h_count_cost = h_sym_count->ProjectionY("h_count_cost");

    TCanvas *c1 = new TCanvas("c1", "ZZ back M", 100,200, 900, 700);
    c1->cd();
    h_sym_cost->Print();
    h_sym_cost->Draw("hist");
    h_sym_cost->SetFillColor(kBlue);
    h_sym_cost->SetMarkerStyle(21);
    h_sym_cost->SetMarkerColor(kBlue);
    c1->Update();

    TCanvas *c2 = new TCanvas("c2", "ZZ back cost", 900, 700);
    c2->cd();
    h_asym_cost->Print();
    h_asym_cost->Draw("hist");
    h_asym_cost->SetFillColor(kBlue);
    h_asym_cost->SetMarkerStyle(21);
    h_asym_cost->SetMarkerColor(kBlue);
    c2->Update();

    TCanvas *c3 = new TCanvas("c3", "ZZ back cost", 900, 700);
    c3->cd();
    h_count_cost->Print();
    h_count_cost->Draw("hist");
    h_count_cost->SetFillColor(kBlue);
    h_count_cost->SetMarkerStyle(21);
    h_count_cost->SetMarkerColor(kBlue);
    c3->Update();
}
