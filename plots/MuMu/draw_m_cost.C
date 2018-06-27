
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



int n_cost_bins = 10;
Float_t cost_bins[] = {-1.0, -.8, -.6, -.4, -.2, 0., 0.2, 0.4, 0.6, 0.8, 1.0};

void make_tau_m_cost_hist(TTree *t1, TH1F *h_m, TH1F *h_cost, bool is_data=false){
    //read event data
    Long64_t size  =  t1->GetEntries();
    Double_t m, xF, cost, mu1_pt, mu2_pt, jet1_cmva, jet2_cmva, gen_weight;
    Double_t bcdef_HLT_SF, bcdef_iso_SF, bcdef_id_SF;
    Double_t gh_HLT_SF, gh_iso_SF, gh_id_SF, pu_SF;
    Double_t jet1_pt, jet2_pt, jet1_b_weight, jet2_b_weight;
    Float_t met_pt;
    Bool_t double_muon_trig;
    Int_t nJets;
    nJets = 2;
    TH2D *h_m_bcdef = (TH2D *)h_m->Clone("h_m_bcdef");
    TH2D *h_m_gh = (TH2D *)h_m->Clone("h_m_gh");
    TH2D *h_cost_bcdef = (TH2D *)h_cost->Clone("h_cost_bcdef");
    TH2D *h_cost_gh = (TH2D *)h_cost->Clone("h_cost_gh");
    is_tau_event = false;
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
        t1->SetBranchAddress("gh_HLT_SF", &gh_HLT_SF);
        t1->SetBranchAddress("gh_iso_SF", &gh_iso_SF);
        t1->SetBranchAddress("gh_id_SF", &gh_id_SF);
        t1->SetBranchAddress("pu_SF", &pu_SF);
        t1->SetBranchAddress("jet1_b_weight", &jet1_b_weight);
        t1->SetBranchAddress("jet2_b_weight", &jet2_b_weight);
        t1->SetBranchAddress("double_muon_trig", &double_muon_trig);
    }

    for (int i=0; i<size; i++) {
        t1->GetEntry(i);
        bool no_bjets = has_no_bjets(nJets, jet1_pt, jet2_pt, jet1_cmva, jet2_cmva);

        if(m >= 150. && met_pt < 50. && no_bjets && !double_muon_trig){
            if(is_data){
                h_m->Fill(m);
                h_cost->Fill(cost);
            }
            else{
                Double_t bcdef_weight = gen_weight * bcdef_HLT_SF * bcdef_iso_SF * bcdef_id_SF * pu_SF;
                Double_t gh_weight = gen_weight * gh_HLT_SF * gh_iso_SF * gh_id_SF* pu_SF;
                if (nJets >= 1){
                    bcdef_weight *= jet1_b_weight;
                    gh_weight *= jet1_b_weight;
                }
                if (nJets >= 2){
                    bcdef_weight *= jet2_b_weight;
                    gh_weight *= jet2_b_weight;
                }
                //Double_t weight = gen_weight;
                h_m_bcdef->Fill(m,bcdef_weight);
                h_m_gh->Fill(m,gh_weight);
                h_cost_bcdef->Fill(cost,bcdef_weight);
                h_cost_gh->Fill(cost,gh_weight);
            }


        }
    }
    if(!is_data){
        /*
        h_m_bcdef ->Scale(bcdef_lumi * 1000);
        h_cost_bcdef ->Scale(bcdef_lumi * 1000);
        h_m_gh ->Scale(gh_lumi * 1000);
        h_cost_gh ->Scale(gh_lumi * 1000);
        */
        h_m->Add(h_m_bcdef, h_m_gh);
        h_cost->Add(h_cost_bcdef, h_cost_gh);
    }
    t1->ResetBranchAddresses();
    printf("Tot xsection %.2e \n", h_m->Integral());
}

void draw_m_cost(){
    TFile *f = TFile::Open("../analyze/output_files/MuMu_fakerate_non_Wjets_MC_sep5.root");
    TTree *t = (TTree *)f->Get("T_data");
    int n_m_bins = 6; 
    Double_t m_bins[] = {150,200,250,350,500,700,1000}; 
    TH1F *h_m_mc = new TH1F("h_m_mc", "DYtoTauTau, dimuon Mass distribution; M_{#mu#mu} (GeV)", n_m_bins, m_bins);

    TH1F *h_cost_mc = new TH1F("h_cost_mc", "Colins-Soper Angular Distribution; c_{*}", n_cost_bins, cost_bins);


    TH1F *h_m_mc_nosig = new TH1F("h_m_mcnosig", "DYtoTauTau, dimuon Mass distribution; M_{#mu#mu} (GeV)", n_m_bins, m_bins);

    TH1F *h_cost_mc_nosig = new TH1F("h_cost_mcnosig", "DYtoTauTau, ditau angular distribution; c_{*}", n_cost_bins, cost_bins);
    //TH1F *h_cost = new TH1F("h_cost", "DYtoTauTau, dimuon angular distribution; Cos(#theta)", 20, -1, 1);


    make_tau_m_cost_hist(t, h_m_mc, h_cost_mc, false);

    h_cost_mc->Scale(1./h_cost_mc->Integral());
    //h_cost_mc_nosig->Scale(1./h_cost_mc_nosig->Integral());
    //
    TCanvas *c1 = new TCanvas("c1", "X", 900, 700);
    h_m_mc->Draw("hist");


    TCanvas *c2 = new TCanvas("c2", "ZZ back cost", 900, 700);
    c2->cd();
    //h_cost_mc->Print();
    h_cost_mc->Draw("hist");
    h_cost_mc->SetMaximum(0.22);
    h_cost_mc->SetStats(kFALSE);
    h_cost_mc_nosig->SetMaximum(0.22);
    h_cost_mc->SetLineColor(kBlue);
    //h_cost_mc->SetMarkerStyle(21);
    h_cost_mc->SetMarkerColor(kBlue);

    c2->Update();

}
