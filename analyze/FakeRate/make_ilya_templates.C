
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <iostream>
#include <cstring>
#include <algorithm>
#include <vector>

#include "TROOT.h"
#include "TFile.h"
#include "TTree.h"
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
#include "TH2D.h"
#include "Math/Functor.h"
#include "../TemplateMaker_systematics.C"




void gen_templates(TTree *t_WJets, TTree *t_QCD, TTree *t_WJets_contam, 
        TTree* t_QCD_contam, TH2F *h, int flag1 = FLAG_MUONS){
    bool ss = true;
    h->Sumw2();
    TH2D *h_err;
    TLorentzVector *lep_p=0;
    TLorentzVector *lep_m=0;
    Double_t pt;
    FakeRate FR;
    if(flag1 == FLAG_MUONS){
        //TH2D *FR;
        setup_new_mu_fakerate(&FR);
        h_err = (TH2D *) FR.h->Clone("h_err");
        h_err->Reset();
        for (int l=0; l<=3; l++){
            TTree *t;
            if (l==0) t = t_WJets;
            if (l==1) t = t_QCD;
            if (l==2) t = t_WJets_contam;
            if (l==3) t = t_QCD_contam;
            Double_t m, xF, cost, jet1_cmva, jet2_cmva, gen_weight;
            Double_t bcdef_HLT_SF, bcdef_iso_SF, bcdef_id_SF;
            Double_t gh_HLT_SF, gh_iso_SF, gh_id_SF;
            Double_t jet1_pt, jet2_pt, jet1_b_weight, jet2_b_weight, pu_SF;
            Double_t evt_fakerate, mu1_fakerate, mu2_fakerate, mu1_eta, mu1_pt, mu2_eta, mu2_pt;
            Int_t iso_mu;
            Bool_t double_muon_trig;
            Float_t met_pt;
            Int_t nJets;
            nJets = 2;
            pu_SF=1;
            t->SetBranchAddress("m", &m);
            t->SetBranchAddress("xF", &xF);
            t->SetBranchAddress("cost", &cost);
            t->SetBranchAddress("met_pt", &met_pt);
            t->SetBranchAddress("jet2_CMVA", &jet2_cmva);
            t->SetBranchAddress("jet1_CMVA", &jet1_cmva);
            t->SetBranchAddress("jet1_pt", &jet1_pt);
            t->SetBranchAddress("jet2_pt", &jet2_pt);
            //t1->SetBranchAddress("evt_fakerate", &evt_fakerate);
            //t1->SetBranchAddress("mu_fakerate", &mu1_fakerate);
            t->SetBranchAddress("mu1_pt", &mu1_pt);
            t->SetBranchAddress("mu2_pt", &mu2_pt);
            t->SetBranchAddress("mu1_eta", &mu1_eta);
            t->SetBranchAddress("mu2_eta", &mu2_eta);
            t->SetBranchAddress("nJets", &nJets);
            t->SetBranchAddress("mu_p", &lep_p);
            t->SetBranchAddress("mu_m", &lep_m);
            if(l==0 || l==2 ){
                t->SetBranchAddress("iso_muon", &iso_mu);
            }
            if(l==2 || l==3){
                t->SetBranchAddress("gen_weight", &gen_weight);
                t->SetBranchAddress("jet1_b_weight", &jet1_b_weight);
                t->SetBranchAddress("jet2_b_weight", &jet2_b_weight);
                t->SetBranchAddress("pu_SF", &pu_SF);
                t->SetBranchAddress("bcdef_HLT_SF", &bcdef_HLT_SF);
                t->SetBranchAddress("bcdef_iso_SF", &bcdef_iso_SF);
                t->SetBranchAddress("bcdef_id_SF", &bcdef_id_SF);
                t->SetBranchAddress("gh_HLT_SF", &gh_HLT_SF);
                t->SetBranchAddress("gh_iso_SF", &gh_iso_SF);
                t->SetBranchAddress("gh_id_SF", &gh_id_SF);
            }

            Long64_t size  =  t->GetEntries();
            for (int i=0; i<size; i++) {
                t->GetEntry(i);
                bool no_bjets = has_no_bjets(nJets, jet1_pt, jet2_pt, jet1_cmva, jet2_cmva);
                bool not_cosmic = notCosmic(*lep_p, *lep_m);
                if(l==0){
                    if(iso_mu ==1) mu1_fakerate = get_new_fakerate_prob(mu1_pt, mu1_eta, FR.h);
                    if(iso_mu ==0) mu1_fakerate = get_new_fakerate_prob(mu2_pt, mu2_eta, FR.h);
                    evt_fakerate = mu1_fakerate/(1-mu1_fakerate);
                }
                if(l==1){
                    mu1_fakerate = get_new_fakerate_prob(mu1_pt, mu1_eta, FR.h);
                    mu2_fakerate = get_new_fakerate_prob(mu2_pt, mu2_eta, FR.h);
                    evt_fakerate = -(mu1_fakerate/(1-mu1_fakerate)) * (mu2_fakerate/(1-mu2_fakerate));
                }
                if(l==2){
                    Double_t bcdef_weight = gen_weight * bcdef_HLT_SF *  bcdef_id_SF * bcdef_iso_SF;
                    Double_t gh_weight = gen_weight * gh_HLT_SF * gh_id_SF * gh_iso_SF;
                    Double_t mc_weight = (bcdef_weight *bcdef_lumi + gh_weight * gh_lumi) * 1000;
                    if(iso_mu ==1) mu1_fakerate = get_new_fakerate_prob(mu1_pt, mu1_eta, FR.h);
                    if(iso_mu ==0) mu1_fakerate = get_new_fakerate_prob(mu2_pt, mu2_eta, FR.h);
                    evt_fakerate = -(mu1_fakerate * mc_weight)/(1-mu1_fakerate);
                }
                if(l==3){
                    Double_t bcdef_weight = gen_weight * bcdef_HLT_SF *  bcdef_id_SF;
                    Double_t gh_weight = gen_weight * gh_HLT_SF * gh_id_SF;
                    Double_t mc_weight = (bcdef_weight *bcdef_lumi + gh_weight * gh_lumi) * 1000;
                    mu1_fakerate = get_new_fakerate_prob(mu1_pt, mu1_eta, FR.h);
                    mu2_fakerate = get_new_fakerate_prob(mu2_pt, mu2_eta, FR.h);
                    evt_fakerate = mc_weight * (mu1_fakerate/(1-mu1_fakerate)) * (mu2_fakerate/(1-mu2_fakerate));
                }
                cost = get_cost(*lep_p, *lep_m);
                TLorentzVector cm = *lep_p + *lep_m;
                pt = cm.Pt();
                xF = compute_xF(cm); 
                float y = abs(cm.Rapidity());
                bool pass =  met_pt < 50.  && no_bjets && not_cosmic;
                m = cm.M();

                if(pass){
                    if(l==0 && iso_mu ==1) h_err->Fill(min(abs(mu1_eta), 2.3), min(mu1_pt, 150.));
                    if(l==0 && iso_mu ==0) h_err->Fill(min(abs(mu2_eta), 2.3), min(mu2_pt, 150.));
                    h->Fill(y, m, evt_fakerate);
                }


            }
            printf("After iter %i current fakerate est is %.0f \n", l, h->Integral());
        }
    }
    else{

        //TH2D *FR;
        setup_new_el_fakerate(&FR);
        h_err = (TH2D *) FR.h->Clone("h_err");
        h_err->Reset();
        for (int l=0; l<=3; l++){
            TTree *t;
            if (l==0) t = t_WJets;
            if (l==1) t = t_QCD;
            if (l==2) t = t_WJets_contam;
            if (l==3) t = t_QCD_contam;
            Double_t m, xF, cost, jet1_cmva, jet2_cmva, gen_weight;
            Double_t jet1_pt, jet2_pt, jet1_b_weight, jet2_b_weight, pu_SF;
            Double_t el_id_SF, el_reco_SF;
            Double_t evt_fakerate, el1_fakerate, el2_fakerate, el1_eta, el1_pt, el2_eta, el2_pt;
            Int_t iso_el;
            Float_t met_pt;
            Int_t nJets;
            nJets = 2;
            pu_SF=1;
            t->SetBranchAddress("m", &m);
            t->SetBranchAddress("xF", &xF);
            t->SetBranchAddress("cost", &cost);
            t->SetBranchAddress("met_pt", &met_pt);
            t->SetBranchAddress("jet2_CMVA", &jet2_cmva);
            t->SetBranchAddress("jet1_CMVA", &jet1_cmva);
            t->SetBranchAddress("jet1_pt", &jet1_pt);
            t->SetBranchAddress("jet2_pt", &jet2_pt);
            //t1->SetBranchAddress("evt_fakerate", &evt_fakerate);
            //t1->SetBranchAddress("el_fakerate", &el1_fakerate);
            t->SetBranchAddress("el1_pt", &el1_pt);
            t->SetBranchAddress("el2_pt", &el2_pt);
            t->SetBranchAddress("el1_eta", &el1_eta);
            t->SetBranchAddress("el2_eta", &el2_eta);
            t->SetBranchAddress("nJets", &nJets);
            t->SetBranchAddress("el_p", &lep_p);
            t->SetBranchAddress("el_m", &lep_m);
            if(l==0 || l==2 ){
                t->SetBranchAddress("iso_el", &iso_el);
            }
            if(l==2 || l==3){
                t->SetBranchAddress("el_id_SF", &el_id_SF);
                t->SetBranchAddress("el_reco_SF", &el_reco_SF);
                t->SetBranchAddress("gen_weight", &gen_weight);
                t->SetBranchAddress("jet1_b_weight", &jet1_b_weight);
                t->SetBranchAddress("jet2_b_weight", &jet2_b_weight);
                t->SetBranchAddress("pu_SF", &pu_SF);
            }

            Long64_t size  =  t->GetEntries();
            for (int i=0; i<size; i++) {
                t->GetEntry(i);
                bool no_bjets = has_no_bjets(nJets, jet1_pt, jet2_pt, jet1_cmva, jet2_cmva);
                bool not_cosmic = notCosmic(*lep_p, *lep_m);
                if(l==0){
                    if(iso_el ==1) el1_fakerate = get_new_fakerate_prob(el1_pt, el1_eta, FR.h);
                    if(iso_el ==0) el1_fakerate = get_new_fakerate_prob(el2_pt, el2_eta, FR.h);
                    evt_fakerate = el1_fakerate/(1-el1_fakerate);
                }
                if(l==1){
                    el1_fakerate = get_new_fakerate_prob(el1_pt, el1_eta, FR.h);
                    el2_fakerate = get_new_fakerate_prob(el2_pt, el2_eta, FR.h);
                    evt_fakerate = -(el1_fakerate/(1-el1_fakerate)) * (el2_fakerate/(1-el2_fakerate));
                }
                if(l==2){

                    Double_t mc_weight = gen_weight * el_id_SF * el_reco_SF  * 1000. * el_lumi;
                    if(iso_el ==1) el1_fakerate = get_new_fakerate_prob(el1_pt, el1_eta, FR.h);
                    if(iso_el ==0) el1_fakerate = get_new_fakerate_prob(el2_pt, el2_eta, FR.h);
                    evt_fakerate = -(el1_fakerate * mc_weight)/(1-el1_fakerate);
                }
                if(l==3){
                    Double_t mc_weight = gen_weight * el_id_SF * el_reco_SF * 1000. * el_lumi;

                    el1_fakerate = get_new_fakerate_prob(el1_pt, el1_eta, FR.h);
                    el2_fakerate = get_new_fakerate_prob(el2_pt, el2_eta, FR.h);
                    evt_fakerate = mc_weight * (el1_fakerate/(1-el1_fakerate)) * (el2_fakerate/(1-el2_fakerate));
                }


                cost = get_cost(*lep_p, *lep_m);
                TLorentzVector cm = *lep_p + *lep_m;
                pt = cm.Pt();
                xF = compute_xF(cm); 
                float y = abs(cm.Rapidity());
                m = cm.M();
                bool pass =  met_pt < 50.  && no_bjets && not_cosmic;
                if(pass){
                    if(l==0 && iso_el ==1) h_err->Fill(min(abs(el1_eta), 2.3), min(el1_pt, 150.));
                    if(l==0 && iso_el ==0) h_err->Fill(min(abs(el2_eta), 2.3), min(el2_pt, 150.));
                    //if(l==3) printf("Evt fr %.2e \n", evt_fakerate);
                    //if(l==3) printf("cost, fr %.2f %.2e \n", cost, evt_fakerate);
                    h->Fill(y, m, evt_fakerate);
                }
            }

            printf("After iter %i current fakerate est is %.0f \n", l, h->Integral());
        }
    }
    printf("Performing fakes cleanup (removing neg. bins) \n");
    cleanup_template(h);
    set_fakerate_errors(h_err, FR.h, h);
    printf("Total Fakerate weight Weight is %.2f \n", h->Integral());
    if(h->Integral() < 0.){
        h->Scale(0.);
        printf("zeroing Fakes template \n");
    }
    return;
}

void make_ilya_templates(){

    auto f_mumu_QCD = TFile::Open("../analyze/output_files/MuMu_QCD_est_mlow_feb21.root");
    auto t_mumu_QCD = (TTree *)f_mumu_QCD->Get("T_data");

    auto f_mumu_WJets = TFile::Open("../analyze/output_files/MuMu_WJets_est_mlow_feb21.root");
    auto t_mumu_WJets = (TTree *)f_mumu_WJets->Get("T_data");

    auto f_mumu_WJets_contam = TFile::Open("../analyze/output_files/MuMu_WJets_mc_mlow_comb_mar6.root");
    auto t_mumu_WJets_contam = (TTree *)f_mumu_WJets_contam->Get("T_data");

    auto f_mumu_QCD_contam = TFile::Open("../analyze/output_files/MuMu_QCD_mc_mlow_comb_mar6.root");
    auto t_mumu_QCD_contam = (TTree *)f_mumu_QCD_contam->Get("T_data");



    auto f_elel_QCD = TFile::Open("../analyze/output_files/ElEl_QCD_est_mlow_feb21.root");
    auto t_elel_QCD = (TTree *)f_elel_QCD->Get("T_data");

    auto f_elel_WJets = TFile::Open("../analyze/output_files/ElEl_WJets_est_mlow_feb21.root");
    auto t_elel_WJets = (TTree *)f_elel_WJets->Get("T_data");

    auto f_elel_WJets_contam = TFile::Open("../analyze/output_files/ElEl_WJets_mc_mlow_comb_mar6.root");
    auto t_elel_WJets_contam = (TTree *)f_elel_WJets_contam->Get("T_data");

    auto f_elel_QCD_contam = TFile::Open("../analyze/output_files/ElEl_QCD_mc_mlow_comb_mar6.root");
    auto t_elel_QCD_contam = (TTree *)f_elel_QCD_contam->Get("T_data");

    float eta_bins[5] = {0., 1.0, 1.25, 1.5, 2.4};
    float mass_bins[46] = {15, 20, 25, 30, 35, 40, 45, 50, 55, 60, 64, 68, 72, 76, 81, 86, 91, 96, 101, 106, 
        110, 115, 120, 126, 133, 141, 150, 160, 171, 185, 200, 220, 243, 273, 320, 380, 440,
        510, 600, 700, 830, 1000, 1200, 1500, 2000, 3000}; 
    int n_mbins = 45;
    int n_etabins = 4;

    TH2F *h_elel = new TH2F("h_fakes_elel", "Electron pairs fakes est", n_etabins, eta_bins, n_mbins, mass_bins);
    TH2F *h_mumu = new TH2F("h_fakes_mumu", "Muon pairs fakes est", n_etabins, eta_bins, n_mbins, mass_bins);

    gen_templates(t_mumu_WJets, t_mumu_QCD, t_mumu_WJets_contam, t_mumu_QCD_contam, h_mumu, FLAG_MUONS);
    gen_templates(t_elel_WJets, t_elel_QCD, t_elel_WJets_contam, t_elel_QCD_contam, h_elel, FLAG_ELECTRONS);

    h_elel->Scale(0.5);
    float mu_scaling = 1./(1. + R_mu_ss_os);
    h_mumu->Scale(mu_scaling);

    TFile *fout = new TFile("lepton_pair_fakes_est.root", "RECREATE");
    fout->cd();
    h_elel->Write();
    h_mumu->Write();
    fout->Close();
}

    
