


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
//#include "../../analyze/HistMaker.C"
#include "../tdrstyle.C"
#include "../CMS_lumi.C"


Double_t bcdef_lumi = 5.746 + 2.572 + 4.242 + 4.024 + 3.104;
// adding new Hv2 data set to get full 2016 luminosity
Double_t gh_lumi =  7.573 + 0.215 + 8.434;
Double_t tot_lumi = 35.9;
Double_t emu_scaling = 1.05;
Double_t emu_unc = 0.05;
TFile *f_data, *f_mc, *f_mc_nosig, *f_ttbar, *f_QCD, *f_WJets, *f_WJets_mc, *f_QCD_mc, *f_diboson, *f_wt;
TTree *t_data, *t_mc, *t_mc_nosig, *t_ttbar, *t_QCD, *t_WJets, *t_WJets_mc, *t_QCD_mc, *t_diboson, *t_wt;

double get_cost(TLorentzVector lep_p, TLorentzVector lep_m){

    TLorentzVector cm = lep_p + lep_m;
    double root2 = sqrt(2);
    double lep_p_pls = (lep_p.E()+lep_p.Pz())/root2;
    double lep_p_min = (lep_p.E()-lep_p.Pz())/root2;
    double lep_m_pls = (lep_m.E()+lep_m.Pz())/root2;
    double lep_m_min = (lep_m.E()-lep_m.Pz())/root2;
    double qt2 = cm.Pt() * cm.Pt();
    double cm_m2 = cm.M2();
    //cost_p = cos(theta)_r (reconstructed collins soper angle, sign
    //may be 'wrong' if lepton pair direction is not the same as inital
    //quark direction)
    double cost = 2*(lep_m_pls*lep_p_min - lep_m_min*lep_p_pls)/sqrt(cm_m2*(cm_m2 + qt2));
    if(cm.Pz() <0.) cost = -cost;
    return cost;
}

double get_cost_v2(TLorentzVector lep_p, TLorentzVector lep_m){

    double Ebeam = 6500.;
    double Pbeam = sqrt(Ebeam*Ebeam - 0.938*0.938);
    TLorentzVector cm = lep_p + lep_m;
    TLorentzVector p1(0., 0., Pbeam, Ebeam);
    TLorentzVector p2(0., 0., -Pbeam, Ebeam);

    if(cm.Pz() < 0. ){
        TLorentzVector p = p1;
        p1 = p2;
        p2 = p;
    }

    TVector3 beta = -cm.BoostVector();
    lep_m.Boost(beta);
    lep_p.Boost(beta);
    p1.Boost(beta);
    p2.Boost(beta);

    // Now calculate the direction of the new z azis

    TVector3 p1u = p1.Vect();
    p1u.SetMag(1.0);
    TVector3 p2u = p2.Vect();
    p2u.SetMag(1.0);
    TVector3 pzu = p1u - p2u;
    pzu.SetMag(1.0);
    lep_m.RotateUz(pzu); 
    double cost = lep_m.CosTheta();
    return cost;
}

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
    t1->SetBranchAddress("mu1_pt_corr", &lep1_pt_corr);
    t1->SetBranchAddress("mu2_pt_corr", &lep2_pt_corr);
    t1->SetBranchAddress("mu1_pt", &lep1_pt);
    t1->SetBranchAddress("mu2_pt", &lep2_pt);
    const double root2 = sqrt(2);
    for (int i=0; i<size; i++) {
        t1->GetEntry(i);
        bool no_bjets = true;

        if(m >= 150. && met_pt < 50. && no_bjets){
            TLorentzVector mu_p, mu_m;
            if((lep_p->Pt() > lep_m->Pt() && lep1_pt_corr > lep2_pt_corr) ||
                    (lep_p->Pt() < lep_m->Pt() && lep1_pt_corr < lep2_pt_corr)){
                mu_p.SetPtEtaPhiE(lep1_pt, lep_p->Eta(), lep_p->Phi(), lep_p->E());
                mu_m.SetPtEtaPhiE(lep2_pt, lep_m->Eta(), lep_m->Phi(), lep_m->E());
            }
            else{
                mu_p.SetPtEtaPhiE(lep2_pt, lep_p->Eta(), lep_p->Phi(), lep_p->E());
                mu_m.SetPtEtaPhiE(lep1_pt, lep_m->Eta(), lep_m->Phi(), lep_m->E());
            }
            double cost1 = get_cost_v2(*lep_p, *lep_m); 
            double cost2 = get_cost_v2(mu_p, mu_m); 
            h_cost1->Fill(cost1, gen_weight);
            h_cost2->Fill(cost2, gen_weight);
            h_delta->Fill(cost1-cost2, gen_weight);
        }



    }
    h_cost1->Scale(1./h_cost1->Integral());
    h_cost2->Scale(1./h_cost2->Integral());
    h_delta->Scale(1./h_delta->Integral());

    t1->ResetBranchAddresses();
}

void make_ratio_plot(char title[80], TH1F* h1, char h1_label[80], TH1F* h2, char h2_label[80], char ratio_label[80], 
        char axis_label[80], bool logy=false){

    TCanvas *c = new TCanvas("c_m", "Histograms", 200, 10, 900, 700);
    TPad *pad1 = new TPad("pad1", "pad1", 0.,0.3,0.98,1.);
    pad1->SetBottomMargin(0);
    pad1->Draw();
    pad1->cd();
    if(logy) pad1->SetLogy();
    h1->Draw("hist E");
    //m_stack->SetMaximum(65000);
    h1->SetMinimum(0.1);
    gStyle->SetEndErrorSize(4);
    h2->Draw("hist E same");


    gStyle->SetLegendBorderSize(0);
    TLegend *leg1 = new TLegend(0.3, 0.3);
    leg1->AddEntry(h1, h1_label, "l");
    leg1->AddEntry(h2, h2_label, "l");
    leg1->Draw();

    //gPad->BuildLegend();
    c->cd();
    TPad *pad2 = new TPad("pad2", "pad2", 0.,0,.98,0.3);
    //pad2->SetTopMargin(0);
    pad2->SetBottomMargin(0.2);
    pad2->SetGridy();
    pad2->Draw();
    pad2->cd();

    auto ratio = (TH1F *) h1->Clone("h_ratio");
    ratio->SetMinimum(0.01);
    ratio->SetMaximum(2.0);
    ratio->Sumw2();
    ratio->SetStats(0);
    ratio->Divide(h2);
    ratio->SetMarkerStyle(21);
    ratio->SetLineColor(kBlack);
    ratio->Draw("ep");
    c->cd();

    ratio->SetTitle("");
    // Y axis ratio plot settings
    ratio->GetYaxis()->SetTitle(ratio_label);
    ratio->GetYaxis()->SetNdivisions(505);
    ratio->GetYaxis()->SetTitleSize(20);
    ratio->GetYaxis()->SetTitleFont(43);
    ratio->GetYaxis()->SetTitleOffset(1.2);
    ratio->GetYaxis()->SetLabelFont(43); // Absolute font size in pixel (precision 3)
    ratio->GetYaxis()->SetLabelSize(15);
    // X axis ratio plot settings
    ratio->GetXaxis()->SetTitle(axis_label);
    ratio->GetXaxis()->SetTitleSize(20);
    ratio->GetXaxis()->SetTitleFont(43);
    ratio->GetXaxis()->SetTitleOffset(3.);
    ratio->GetXaxis()->SetLabelFont(43); // Absolute font size in pixel (precision 3)
    ratio->GetXaxis()->SetLabelSize(20);

    //lumi_sqrtS = "";       // used with iPeriod = 0, e.g. for simulation-only plots (default is an empty string)
    int iPeriod = 4; 
    CMS_lumi(pad1, iPeriod, 33 );
    c->Print(title);
    return;
}


void draw_cost_calcs_cmp(){

    TFile *f_mc= TFile::Open("../analyze/output_files/MuMu_DY_slim_june25.root");
    TTree *t_mc = (TTree *)f_mc->Get("T_data");

    setTDRStyle();
    TH1F *cost1 = new TH1F("bin_cost", "Binned MC", 40, -1.,1.);
    TH1F *cost2 = new TH1F("cost2", "Binned MC", 40, -1.,1.);
    TH1F *delta = new TH1F("delta", "Binned MC", 40, -0.2,0.2);
    
    cost1->SetLineColor(kBlue);
    cost1->SetLineWidth(3);
    cost2->SetLineColor(kRed);
    cost2->SetLineWidth(3);

    make_cost_hists(t_mc, cost1, cost2, delta);

    
    make_ratio_plot("cost_calcs_cmp.pdf", cost1, "v1 ",cost2, "v2", "Binned/Unbinned", "cost", false);

    TCanvas *c2 = new TCanvas("c2", "", 0,0, 800, 800);
    delta->Draw("hist");

}


