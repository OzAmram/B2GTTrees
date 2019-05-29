


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
#include "../../Utils/HistMaker.C"
#include "../../Utils/root_files.h"
#include "../../Utils/PlotUtils.C"
#include "../tdrstyle.C"
#include "../CMS_lumi.C"


TFile *f_data, *f_mc, *f_mc_nosig, *f_ttbar, *f_QCD, *f_WJets, *f_WJets_mc, *f_QCD_mc, *f_diboson, *f_wt;
TTree *t_data, *t_mc, *t_mc_nosig, *t_ttbar, *t_QCD, *t_WJets, *t_WJets_mc, *t_QCD_mc, *t_diboson, *t_wt;

void make_cost_hists(TTree *t1, TH1F *h_cost1, TH1F* h_cost2, TH1F* h_delta){
    //read event data
    Long64_t size  =  t1->GetEntries();
    Double_t m, xF, cost, mu1_pt, mu2_pt, jet1_cmva, jet2_cmva, gen_weight, cost_st;
    Double_t bcdef_HLT_SF, bcdef_iso_SF, bcdef_id_SF;
    Double_t gh_HLT_SF, gh_iso_SF, gh_id_SF, el_id_SF, el_reco_SF, el_HLT_SF;
    Double_t bcdef_trk_SF, gh_trk_SF;
    Double_t jet1_pt, jet2_pt, jet1_b_weight, jet2_b_weight, pu_SF;
    Double_t lep1_pt_corr, lep2_pt_corr, lep1_pt, lep2_pt;
    TLorentzVector *gen_mu_p = 0;
    TLorentzVector *gen_mu_m = 0;
    TLorentzVector *gen_el_p = 0;
    TLorentzVector *gen_el_m = 0;
    TLorentzVector *lep_p = 0;
    TLorentzVector *lep_m = 0;
    TLorentzVector cm;
    Float_t met_pt;
    Int_t nJets;
    nJets = 2;
    pu_SF=1;
    t1->SetBranchAddress("m", &m);
    t1->SetBranchAddress("xF", &xF);
    t1->SetBranchAddress("cost", &cost);
    t1->SetBranchAddress("cost_st", &cost_st);
    t1->SetBranchAddress("met_pt", &met_pt);
    t1->SetBranchAddress("jet2_CMVA", &jet2_cmva);
    t1->SetBranchAddress("jet1_CMVA", &jet1_cmva);
    t1->SetBranchAddress("jet1_pt", &jet1_pt);
    t1->SetBranchAddress("jet2_pt", &jet2_pt);
    t1->SetBranchAddress("nJets", &nJets);
    t1->SetBranchAddress("gen_weight", &gen_weight);
    //t1->SetBranchAddress("jet1_b_weight", &jet1_b_weight);
    //t1->SetBranchAddress("jet2_b_weight", &jet2_b_weight);
    t1->SetBranchAddress("pu_SF", &pu_SF);
    t1->SetBranchAddress("mu_p", &lep_p);
    t1->SetBranchAddress("mu_m", &lep_m);
    t1->SetBranchAddress("bcdef_HLT_SF", &bcdef_HLT_SF);
    t1->SetBranchAddress("bcdef_iso_SF", &bcdef_iso_SF);
    t1->SetBranchAddress("bcdef_id_SF", &bcdef_id_SF);
    t1->SetBranchAddress("gh_HLT_SF", &gh_HLT_SF);
    t1->SetBranchAddress("gh_iso_SF", &gh_iso_SF);
    t1->SetBranchAddress("gh_id_SF", &gh_id_SF);
    t1->SetBranchAddress("gh_trk_SF", &gh_trk_SF);
    t1->SetBranchAddress("bcdef_trk_SF", &bcdef_trk_SF);
    t1->SetBranchAddress("mu1_pt", &lep1_pt);
    t1->SetBranchAddress("mu2_pt", &lep2_pt);
    const double root2 = sqrt(2);
    for (int i=0; i<size; i++) {
        t1->GetEntry(i);
        bool no_bjets = true;

        if(m >= 150. && met_pt < 50. && no_bjets){
            TLorentzVector mu_p, mu_m;
            double cost1 = get_cost(*lep_p, *lep_m); 
            double cost2 = get_cost_v2(*lep_p, *lep_m); 
            double evt_weight = gen_weight * 1000. * mu_lumi; 
            h_cost1->Fill(cost1, evt_weight);
            h_cost2->Fill(cost2, evt_weight);


            h_delta->Fill(cost1-cost2, evt_weight);
        }



    }
    //h_cost1->Scale(1./h_cost1->Integral());
    //h_cost2->Scale(1./h_cost2->Integral());
    //h_delta->Scale(1./h_delta->Integral());

    t1->ResetBranchAddresses();
}


void draw_cost_calcs_cmp(){

    init();
    setTDRStyle();
    TH1F *cost1 = new TH1F("bin_cost", "Binned MC", 10, -1.,1.);
    TH1F *cost2 = new TH1F("cost2", "Binned MC", 10, -1.,1.);
    TH1F *delta = new TH1F("delta", "Binned MC", 500, -0.2,0.2);
    
    cost1->SetLineColor(kBlue);
    cost1->SetLineWidth(3);
    cost2->SetLineColor(kRed);
    cost2->SetLineWidth(3);

    make_cost_hists(t_mumu_mc, cost1, cost2, delta);
    printf("integrals are %.1f %.1f \n", cost1->Integral(), cost2->Integral());

    
    make_ratio_plot("cost_calcs_cmp.pdf", cost1, "v1 ",cost2, "v2", "v1/v2", "cost", false, false);

    TCanvas *c2 = new TCanvas("c2", "", 0,0, 800, 800);
    delta->Draw("hist");

}


