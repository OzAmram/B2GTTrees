


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

const int type = FLAG_MUONS;

TFile *f_data, *f_mc, *f_mc_nosig, *f_ttbar, *f_QCD, *f_WJets, *f_WJets_mc, *f_QCD_mc, *f_diboson, *f_wt;
TTree *t_data, *t_mc, *t_mc_nosig, *t_ttbar, *t_QCD, *t_WJets, *t_WJets_mc, *t_QCD_mc, *t_diboson, *t_wt;

void make_m_cost_pt_gen_hist(TTree *t1, TH1F *h_m, TH1F *h_cost, TH1F *h_pt){
    //read event data
    Long64_t size  =  t1->GetEntries();
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
    t1->SetBranchAddress("jet1_b_weight", &jet1_b_weight);
    t1->SetBranchAddress("jet2_b_weight", &jet2_b_weight);
    t1->SetBranchAddress("pu_SF", &pu_SF);
    t1->SetBranchAddress("gen_el_p", &gen_el_p);
    t1->SetBranchAddress("gen_el_m", &gen_el_m);
    t1->SetBranchAddress("el_id_SF", &el_id_SF);
    t1->SetBranchAddress("el_reco_SF", &el_reco_SF);
    t1->SetBranchAddress("el_HLT_SF", &el_HLT_SF);
    const double root2 = sqrt(2);
    for (int i=0; i<size; i++) {
        t1->GetEntry(i);
        bool no_bjets = has_no_bjets(nJets, jet1_pt, jet2_pt, jet1_cmva, jet2_cmva);
        if(m >= 150 && met_pt < 50.  && no_bjets){
            cm = *gen_el_p + *gen_el_m;
            Double_t pt = cm.Pt();

            double el_p_pls = (gen_el_p->E()+gen_el_p->Pz())/root2;
            double el_p_min = (gen_el_p->E()-gen_el_p->Pz())/root2;
            double el_m_pls = (gen_el_m->E()+gen_el_m->Pz())/root2;
            double el_m_min = (gen_el_m->E()-gen_el_m->Pz())/root2;
            double qt2 = cm.Px()*cm.Px()+cm.Py()*cm.Py();
            double cm_m2 = cm.M2();
            //cost_p = cos(theta)_r (reconstructed collins soper angle, sign
            //may be 'wrong' if lepton pair direction is not the same as inital
            //quark direction)
            double gen_cost = 2*(el_m_pls*el_p_min - el_m_min*el_p_pls)/sqrt(cm_m2*(cm_m2 + qt2));
            if(cost > 0) gen_cost = fabs(gen_cost);
            else gen_cost = -fabs(gen_cost);

            Double_t evt_weight = gen_weight *pu_SF * el_id_SF * el_reco_SF * el_HLT_SF * el_lumi *1000.;
            h_m->Fill(m, evt_weight);
            h_cost->Fill(cost, evt_weight);
            h_pt->Fill(pt, evt_weight);
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


void draw_gen_cmp(){
    TFile *f_mc_unbinned = TFile::Open("../analyze/output_files/ElEl_DY_slim_nov2.root");
    TTree *t_mc_unbinned = (TTree *)f_mc_unbinned->Get("T_data");


    TFile *f_mc_binned = TFile::Open("../analyze/output_files/ElEl_DY_unbinned_slim_sep11.root");
    TTree *t_mc_binned = (TTree *)f_mc_binned->Get("T_data");
    setTDRStyle();




    TH1F *binned_m = new TH1F("bin_m", "Binned MC", 40, 150.,1000.);
    TH1F *binned_cost = new TH1F("bin_cost", "Binned MC", 40, -1.,1.);
    TH1F *binned_pt = new TH1F("bin_pt", "Binned MC", 25, 0.,400.);

    TH1F *unbinned_m = new TH1F("unbin_m", "unbinned MC", 40, 150.,1000.);
    TH1F *unbinned_cost = new TH1F("unbin_cost", "unbinned MC", 40, -1.,1.);
    TH1F *unbinned_pt = new TH1F("unbin_pt", "unbinned MC", 25, 0.,400.);

    binned_m->SetLineColor(kBlue);
    binned_m->SetLineWidth(3);
    unbinned_m->SetLineColor(kRed);
    unbinned_m->SetLineWidth(3);

    binned_cost->SetLineColor(kBlue);
    binned_cost->SetLineWidth(3);
    unbinned_cost->SetLineColor(kRed);
    unbinned_cost->SetLineWidth(3);

    binned_pt->SetLineColor(kBlue);
    binned_pt->SetLineWidth(3);
    unbinned_pt->SetLineColor(kRed);
    unbinned_pt->SetLineWidth(3);
    
    make_m_cost_pt_gen_hist(t_mc_binned, binned_m, binned_cost, binned_pt);
    make_m_cost_pt_gen_hist(t_mc_unbinned, unbinned_m, unbinned_cost, unbinned_pt);



    make_ratio_plot("ElEl_gen_m_cmp.pdf", binned_m, "Binned MC",unbinned_m, "Unbinned MC", "Binned/Unbinned", "M (GeV)", true);
    make_ratio_plot("ElEl_gen_cost_cmp.pdf", binned_cost, "Binned MC",unbinned_cost, "Unbinned MC", "Binned/Unbinned", "Cos(#theta_{r})", false);
    make_ratio_plot("ElEl_gen_pt_cmp.pdf", binned_pt, "Binned MC",unbinned_pt, "Unbinned MC", "Binned/Unbinned", "Dilepton Pt (GeV)", true);




}



