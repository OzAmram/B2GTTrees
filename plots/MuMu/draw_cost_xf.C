
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



int n_cost_bins = 10;
Float_t cost_bins[] = {-1.0, -.8, -.6, -.4, -.2, 0., 0.2, 0.4, 0.6, 0.8, 1.0};

void make_m_cost_xf_hist(TTree *t1, TH1F *h_cost_st, TH1F *h_cost, TH1F *h_xF, bool is_data=false){
    //read event data
    Long64_t size  =  t1->GetEntries();
    Double_t m, xF, cost, cost_st, mu1_pt, mu2_pt, jet1_cmva, jet2_cmva, gen_weight;
    Double_t bcdef_HLT_SF, bcdef_iso_SF, bcdef_id_SF;
    Double_t gh_HLT_SF, gh_iso_SF, gh_id_SF;
    Double_t jet1_pt, jet2_pt, jet1_b_weight=1., jet2_b_weight=1.;
    Float_t met_pt;
    Bool_t is_tau_event;
    Int_t nJets;
    nJets = 2;
    TH2D *h_cost_st_bcdef = (TH2D *)h_cost_st->Clone("h_cost_st_bcdef");
    TH2D *h_cost_st_gh = (TH2D *)h_cost_st->Clone("h_cost_st_gh");
    TH2D *h_cost_bcdef = (TH2D *)h_cost->Clone("h_cost_bcdef");
    TH2D *h_cost_gh = (TH2D *)h_cost->Clone("h_cost_gh");
    TH2D *h_xF_bcdef = (TH2D *)h_xF->Clone("h_xF_bcdef");
    TH2D *h_xF_gh = (TH2D *)h_xF->Clone("h_xF_gh");
    is_tau_event = false;
    t1->SetBranchAddress("m", &m);
    t1->SetBranchAddress("xF", &xF);
    t1->SetBranchAddress("cost", &cost);
    t1->SetBranchAddress("cost_st", &cost_st);
    t1->SetBranchAddress("met_pt", &met_pt);
    t1->SetBranchAddress("jet2_CMVA", &jet2_cmva);
    t1->SetBranchAddress("jet1_CMVA", &jet1_cmva);
    t1->SetBranchAddress("jet1_pt", &jet1_pt);
    t1->SetBranchAddress("jet2_pt", &jet2_pt);
    t1->SetBranchAddress("is_tau_event", &is_tau_event);
    if(!is_data){
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
    }

    for (int i=0; i<size; i++) {
        t1->GetEntry(i);
        bool no_bjets = has_no_bjets(nJets, jet1_pt, jet2_pt, jet1_cmva, jet2_cmva);

        if(!is_tau_event && m >= 150 && met_pt < 50. && no_bjets){
            if(is_data){
                h_cost_st->Fill(cost_st);
                h_cost->Fill(cost);
                h_xF->Fill(xF);
            }
            else{
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
                //Double_t weight = gen_weight;
                h_cost_st_bcdef->Fill(cost_st,bcdef_weight);
                h_cost_st_gh->Fill(cost_st,gh_weight);
                h_cost_bcdef->Fill(cost,bcdef_weight);
                h_cost_gh->Fill(cost,gh_weight);
                h_xF_bcdef->Fill(xF,bcdef_weight);
                h_xF_gh->Fill(xF,gh_weight);
            }


        }
    }
    if(!is_data){
        h_cost_st_bcdef ->Scale(bcdef_lumi * 1000);
        h_cost_bcdef ->Scale(bcdef_lumi * 1000);
        h_xF_bcdef ->Scale(bcdef_lumi * 1000);
        h_cost_st_gh ->Scale(gh_lumi * 1000);
        h_cost_gh ->Scale(gh_lumi * 1000);
        h_xF_gh ->Scale(gh_lumi * 1000);
        h_cost_st->Add(h_cost_st_bcdef, h_cost_st_gh);
        h_cost->Add(h_cost_bcdef, h_cost_gh);
        h_xF->Add(h_xF_bcdef, h_xF_gh);
    }
    t1->ResetBranchAddresses();
}

void draw_cost_xf(){
    mumu_init();
    setTDRStyle();
    int n_m_bins = 6; 
    Double_t m_bins[] = {150,200,250,350,500,700,1000000}; 

    TH1F *h_xF_mc = new TH1F("h_xF_mc", "#mu^{+}#mu^{-} Feynman-x; |xF|", 20, 0, 0.5);

    TH1F *h_cost_mc = new TH1F("h_cost_mc", "Colins-Soper Angular Distribution", n_cost_bins, cost_bins);
    TH1F *h_cost_st_mc = new TH1F("h_cost_st_mc", "Colins-Soper True Angular Distribution; c_{*}", n_cost_bins, cost_bins);


    TH1F *h_cost_st_mc_nosig = new TH1F("h_cost_st_mc_nosig", "Colins-Soper Angular Distribution; c_{*}", n_cost_bins, cost_bins);

    TH1F *h_xF_mc_nosig = new TH1F("h_xF_mc_nosig", "#mu^{+}#mu^{-} Feynman-x; |xF|", 20, 0, 0.5);

    TH1F *h_cost_mc_nosig = new TH1F("h_cost_mcnosig", "DYtoTauTau, ditau angular distribution; c_{r}", n_cost_bins, cost_bins);
    //TH1F *h_cost = new TH1F("h_cost", "DYtoTauTau, dimuon angular distribution; Cos(#theta)", 20, -1, 1);


    make_m_cost_xf_hist(t_mc, h_cost_st_mc, h_cost_mc, h_xF_mc, false);
    make_m_cost_xf_hist(t_mc_nosig, h_cost_st_mc_nosig, h_cost_mc_nosig, h_xF_mc_nosig, false);

    h_cost_mc->Scale(1./h_cost_mc->Integral());
    h_cost_st_mc->Scale(1./h_cost_st_mc->Integral());
    h_cost_mc_nosig->Scale(1./h_cost_mc_nosig->Integral());

    h_xF_mc->Scale(1./h_xF_mc->Integral());
    h_xF_mc_nosig->Scale(1./h_xF_mc_nosig->Integral());

    TCanvas *c2 = new TCanvas("c2", "ZZ back cost", 900, 700);
    c2->cd();
    //h_cost_mc->Print();
    h_cost_mc->Draw("hist");
    h_cost_mc->GetXaxis()->SetTitle("DY to #mu#mu Colins-Soper Cos(#theta)");
    h_cost_mc->SetLineWidth(3);
    h_cost_mc->SetMaximum(0.22);
    h_cost_mc->SetStats(kFALSE);
    h_cost_mc_nosig->SetMaximum(0.22);
    h_cost_mc->SetLineColor(kBlue);
    //h_cost_mc->SetMarkerStyle(21);
    h_cost_mc->SetMarkerColor(kBlue);
    h_cost_mc->SetMinimum(0);

    h_cost_st_mc->Draw("hist same");
    h_cost_st_mc->SetLineWidth(3);
    h_cost_st_mc->SetLineColor(kRed);
    //h_cost_mc->SetMarkerStyle(21);
    h_cost_st_mc->SetMarkerColor(kRed);
    h_cost_st_mc->SetMinimum(0);

    gStyle->SetLegendBorderSize(0);
    TLegend *leg1 = new TLegend(0.5, 0.65, 0.75, 0.8);
    leg1->AddEntry(h_cost_mc, 
            //"#mu^{+}#mu^{-} Reconstructed Angle #scale[1.5]{#font[22]{(c_{r})}}   ", "f");
            "Reconstructed Angle #scale[1.5]{#font[22]{(c_{r})}}   ", "f");
    leg1->AddEntry(h_cost_st_mc, 
            //"#tau^{+}#tau^{-} Correct Angle #scale[1.5]{#font[22]{(c_{*})}} ", "f");
            "Correct Angle #scale[1.5]{#font[22]{(c_{*})}} ", "f");
    leg1->Draw();
    c2->Update();



    TCanvas *c3 = new TCanvas("c3", "ZZ back xF", 900, 700);
    c3->cd();
    //h_xF_mc->Print();
    h_xF_mc->Draw("hist");
    h_xF_mc->SetMaximum(0.3);
    h_xF_mc->SetStats(kFALSE);
    h_xF_mc_nosig->SetMaximum(0.3);
    h_xF_mc->SetLineColor(kBlue);
    h_xF_mc->SetLineWidth(3);
    h_xF_mc->SetMarkerStyle(21);
    h_xF_mc->SetMarkerColor(kBlue);

    h_xF_mc_nosig->Draw("hist same");
    h_xF_mc_nosig->SetLineColor(kRed);
    h_xF_mc_nosig->SetLineWidth(3);
    //h_xF_mc->SetMarkerStyle(21);
    h_xF_mc_nosig->SetMarkerColor(kRed);
    TLegend *leg2 = new TLegend(0.5, 0.65, 0.75, 0.8);
    leg2->AddEntry(h_xF_mc, "Events with an asymmetry (q #bar{q} + q g + #bar{q} g)  ", "f");
    leg2->AddEntry(h_xF_mc_nosig, "Events without an asymmetry (q q + #bar{q} #bar{q} + g g)  ", "f");
    leg2->Draw();
    c3->Update();


    /*
    lumiTextSize     = 0.2;
    lumiTextOffset   = 0.2;
    cmsTextSize      = 0.35;
    cmsTextOffset    = 0.1;  // only used in outOfFrame version
    */

    writeExtraText = true;
    extraText = "Preliminary Simulation";
    lumi_sqrtS = "";       // used with iPeriod = 0, e.g. for simulation-only plots (default is an empty string)
    int iPeriod = 0; 
    CMS_lumi( c3, iPeriod, 33 );
    CMS_lumi( c2, iPeriod, 11 );

}
