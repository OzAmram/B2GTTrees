



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
#include "TH2D.h"
#include "TLegend.h"
#include "Math/Functor.h"
#include "../../analyze/HistMaker.C"
#include "../tdrstyle.C"
#include "../CMS_lumi.C"
#include "../../generator_stuff/gen_level/GenLoader.cc"

const int type = FLAG_MUONS;

TFile *f_data, *f_mc, *f_mc_nosig, *f_ttbar, *f_QCD, *f_WJets, *f_WJets_mc, *f_QCD_mc, *f_diboson, *f_wt;
TTree *t_data, *t_mc, *t_mc_nosig, *t_ttbar, *t_QCD, *t_WJets, *t_WJets_mc, *t_QCD_mc, *t_diboson, *t_wt;

void make_m_cost_pt_gen_hist(TTree *t1, TH1F *h_m, TH1F *h_cost, TH1F *h_pt){
    printf("making gen_Hist \n");
    //read event data
    Long64_t size  =  t1->GetEntries();
    printf("Size is %i \n", (int) size);
    Double_t m, xF, cost, mu1_pt, mu2_pt, jet1_cmva, jet2_cmva, gen_weight, cost_st;
    Double_t bcdef_HLT_SF, bcdef_iso_SF, bcdef_id_SF;
    Double_t gh_HLT_SF, gh_iso_SF, gh_id_SF, el_id_SF, el_reco_SF, el_HLT_SF;
    Double_t bcdef_trk_SF, gh_trk_SF;
    Double_t jet1_pt, jet2_pt, jet1_b_weight, jet2_b_weight, pu_SF;
    TLorentzVector *gen_mu_p = 0;
    TLorentzVector *gen_mu_m = 0;
    TLorentzVector *gen_el_p = 0;
    TLorentzVector *gen_el_m = 0;
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
    t1->SetBranchAddress("gen_mu_p", &gen_mu_p);
    t1->SetBranchAddress("gen_mu_m", &gen_mu_m);
    const double root2 = sqrt(2);
    for (int i=0; i<size; i++) {
        //printf("%i \n", i);
        t1->GetEntry(i);
        bool no_bjets = has_no_bjets(nJets, jet1_pt, jet2_pt, jet1_cmva, jet2_cmva);

        if(m >= 150. && met_pt < 50. && no_bjets){
            cm = *gen_mu_p + *gen_mu_m;
            Double_t pt = cm.Pt();
            double mu_p_pls = (gen_mu_p->E()+gen_mu_p->Pz())/root2;
            double mu_p_min = (gen_mu_p->E()-gen_mu_p->Pz())/root2;
            double mu_m_pls = (gen_mu_m->E()+gen_mu_m->Pz())/root2;
            double mu_m_min = (gen_mu_m->E()-gen_mu_m->Pz())/root2;
            double qt2 = cm.Px()*cm.Px()+cm.Py()*cm.Py();
            double cm_m2 = cm.M2();
            //cost_p = cos(theta)_r (reconstructed collins soper angle, sign
            //may be 'wrong' if lepton pair direction is not the same as inital
            //quark direction)
            double gen_cost = 2*(mu_m_pls*mu_p_min - mu_m_min*mu_p_pls)/sqrt(cm_m2*(cm_m2 + qt2));
            if(cost > 0) gen_cost = fabs(gen_cost);
            else gen_cost = -fabs(gen_cost);

            double final_weight = 1000.*tot_lumi*gen_weight;
            //Double_t weight = gen_weight;
            h_m->Fill(m,final_weight);
            h_cost->Fill(gen_cost, final_weight);
            h_pt->Fill(pt, final_weight);


        }
    }
    h_cost->Scale(1./h_cost->Integral());

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
    TLegend *leg1 = new TLegend(0.25, 0.05, 0.4, 0.2);
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

void fill_gen_hist(TTree *t, TH1F *h_m, TH1F *h_pt){
    Long64_t size = t->GetEntries();
    GenLoader *g = new GenLoader(t);
    bool got_lep1, got_lep2;
    TLorentzVector lep1, lep2;
    for(int i =0; i<size; i++){
        got_lep1 = got_lep2 = false;
        g->load(i);
        for (int j=0; j<g->fGens->GetEntriesFast(); ++j) {
            TGenParticle *p1 = (TGenParticle*)((*(g->fGens))[j]);
            if((abs(p1->pdgId) == 11 || abs(p1->pdgId) == 13) && p1->status == 23) {
                if(!got_lep1){
                    got_lep1 = true;
                    lep1.SetPtEtaPhiM(p1->pt, p1->eta, p1->phi, p1->mass);
                }
                else if(!got_lep2){
                    got_lep2 = true;
                    lep2.SetPtEtaPhiM(p1->pt, p1->eta, p1->phi, p1->mass);
                }
                else{printf("Extra lepton! \n");}
            }
        }
        if(got_lep1 && got_lep2){
            TLorentzVector cm = lep1 + lep2;
            if(cm.M() > 150.){
                h_m->Fill(cm.M(), g->fGenInfo->weight);
                h_pt->Fill(cm.Pt(), g->fGenInfo->weight);
            }
        }
    }
    h_m->Scale(1./h_m->Integral());
    h_pt->Scale(1./h_pt->Integral());
    return;
}

void draw_4way_cmp(){
    TFile *f_mc_unbinned = TFile::Open("../analyze/output_files/MuMu_DY_april9_unbinned.root");
    TTree *t_mc_unbinned = (TTree *)f_mc_unbinned->Get("T_data");

    TFile *f_mc_binned = TFile::Open("../analyze/output_files/MuMu_DY_zpeak_june14.root");
    TTree *t_mc_binned = (TTree *)f_mc_binned->Get("T_data");
    setTDRStyle();


    

    TH1F *binned_m = new TH1F("bin_m", "Binned MC", 40, 150.,1000.);
    TH1F *binned_cost = new TH1F("bin_cost", "Binned MC", 40, -1.,1.);
    TH1F *binned_pt = new TH1F("bin_pt", "Binned MC", 20, 0.,400.);


    TH1F *unbinned_m = new TH1F("unbin_m", "unbinned MC", 40, 150.,1000.);
    TH1F *unbinned_cost = new TH1F("unbin_cost", "unbinned MC", 40, -1.,1.);
    TH1F *unbinned_pt = new TH1F("unbin_pt", "unbinned MC", 20, 0.,400.);

    printf("making full hists \n");
    make_m_cost_pt_gen_hist(t_mc_binned, binned_m, binned_cost, binned_pt);
    make_m_cost_pt_gen_hist(t_mc_unbinned, unbinned_m, unbinned_cost, unbinned_pt);
    printf("bakc from hists \n");
    binned_pt->Scale(1./binned_pt->Integral());
    unbinned_pt->Scale(1./unbinned_pt->Integral());
    TFile *f_unbinned = TFile::Open("../generator_stuff/gen_level/DY_M50_gen.root");
    TTree *t_unbinned = (TTree *)f_unbinned->Get("Events");

    TFile *f_binned = TFile::Open("../generator_stuff/gen_level/DY_M100_gen.root");
    TTree *t_binned = (TTree *)f_binned->Get("Events");

    TH1F *h_unbin_m = new TH1F("h_unbin_m", "", 20, 100, 200);
    TH1F *h_unbin_pt = new TH1F("h_unbin_pt", "", 20, 0, 400);

    printf("making gen hists \n");
    TH1F *h_bin_m = new TH1F("h_bin_m", "", 20, 100, 200);
    TH1F *h_bin_pt = new TH1F("h_bin_pt", "", 20, 0, 400);
    fill_gen_hist(t_unbinned, h_unbin_m, h_unbin_pt);
    fill_gen_hist(t_binned, h_bin_m, h_bin_pt);


    binned_pt->SetLineColor(kRed);
    binned_pt->SetLineWidth(3);
    unbinned_pt->SetLineColor(kBlue);
    unbinned_pt->SetLineWidth(3);

    h_bin_pt->SetLineColor(kGreen);
    h_bin_pt->SetLineWidth(3);
    h_unbin_pt->SetLineColor(kBlack);
    h_unbin_pt->SetLineWidth(3);
    


    printf("drawing\n");
    TCanvas *c = new TCanvas("c_m", "Histograms", 200, 10, 900, 700);
    c->SetLogy();
    binned_pt->Draw("hist E");
    //m_stack->SetMaximum(65000);
    gStyle->SetEndErrorSize(4);
    unbinned_pt->Draw("hist E same");
    h_unbin_pt->Draw("hist E same");
    h_bin_pt->Draw("hist E same");
    h_unbin_pt->Draw("hist E same");


    gStyle->SetLegendBorderSize(0);
    TLegend *leg1 = new TLegend(0.25, 0.2, 0.4, 0.4);
    leg1->AddEntry(unbinned_pt, "DY-M50 FULL", "l");
    leg1->AddEntry(binned_pt, "DY-M100 FULL", "l");
    leg1->AddEntry(h_bin_pt, "DY-M50 GEN", "l");
    leg1->AddEntry(h_unbin_pt, "DY-M100 GEN", "l");
    leg1->Draw();

    c->Print("4way_cmp.pdf");


    TCanvas *c2 = new TCanvas("c_m", "Histograms", 200, 10, 900, 700);
    c2->SetLogy();
    binned_m->Draw("hist E");
    //m_stack->SetMaximum(65000);
    gStyle->SetEndErrorSize(4);
    unbinned_m->Draw("hist E same");
    h_unbin_m->Draw("hist E same");
    h_bin_m->Draw("hist E same");
    h_unbin_m->Draw("hist E same");


    gStyle->SetLegendBorderSize(0);
    TLegend *leg2 = new TLegend(0.25, 0.2, 0.4, 0.4);
    leg2->AddEntry(unbinned_pt, "DY-M50 FULL", "l");
    leg2->AddEntry(binned_pt, "DY-M100 FULL", "l");
    leg2->AddEntry(h_bin_pt, "DY-M50 GEN", "l");
    leg2->AddEntry(h_unbin_pt, "DY-M100 GEN", "l");
    leg2->Draw();

    c2->Print("4way_m_cmp.pdf");



}



