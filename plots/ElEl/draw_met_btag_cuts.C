


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
#include "../tdrstyle.C"
#include "../CMS_lumi.C"
#include "../../Utils/HistMaker.C"
#include "../../Utils/root_files.h"

const int type = FLAG_ELECTRONS;


void make_cut_hists(TTree *t1, TH1F *h_both, TH1F* h_metonly, TH1F *h_btagonly ){
    //read event data
    Long64_t size  =  t1->GetEntries();
    Double_t m, xF, cost, mu1_pt, mu2_pt, jet1_cmva, jet2_cmva, gen_weight, cost_st;
    Double_t bcdef_HLT_SF, bcdef_iso_SF, bcdef_id_SF;
    Double_t gh_HLT_SF, gh_iso_SF, gh_id_SF, el_id_SF, el_reco_SF, el_HLT_SF;
    Double_t bcdef_trk_SF, gh_trk_SF;
    Double_t jet1_pt, jet2_pt, jet1_b_weight, jet2_b_weight, pu_SF;
    jet1_b_weight = jet2_b_weight =1.;
    TLorentzVector *el_p = 0;
    TLorentzVector *el_m = 0;
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
    //t1->SetBranchAddress("jet1_b_weight", &jet1_b_weight);
    //t1->SetBranchAddress("jet2_b_weight", &jet2_b_weight);
    t1->SetBranchAddress("pu_SF", &pu_SF);
    t1->SetBranchAddress("el_p", &el_p);
    t1->SetBranchAddress("el_m", &el_m);
    t1->SetBranchAddress("el_id_SF", &el_id_SF);
    t1->SetBranchAddress("el_reco_SF", &el_reco_SF);
    t1->SetBranchAddress("el_HLT_SF", &el_HLT_SF);
    const double root2 = sqrt(2);
    int nEvent=0;
    for (int i=0; i<size; i++) {
        t1->GetEntry(i);
        bool no_bjets = has_no_bjets(nJets, jet1_pt, jet2_pt, jet1_cmva, jet2_cmva);

        if(m >= 150.){
            nEvent++;
            cost = get_cost_v2(*el_p, *el_m);

            Double_t evt_weight = 1000. * el_lumi* gen_weight *pu_SF * el_id_SF * el_reco_SF * el_HLT_SF;
            if (nJets >= 1){
                evt_weight *= jet1_b_weight;
            }
            if (nJets >= 2){
                evt_weight *= jet2_b_weight;
            }

            if(no_bjets && met_pt < 50. ){
                h_both->Fill(cost,evt_weight);
            }

            if(met_pt < 50.){
                h_metonly->Fill(cost, evt_weight);
            }
            if(no_bjets){
                h_btagonly->Fill(cost, evt_weight);
            }


        }
    }
    printf("%i events \n", nEvent);
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
    ratio->SetMinimum(0.6);
    ratio->SetMaximum(1.4);
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


void draw_met_btag_cuts(){
    init();
    setTDRStyle();


    

    TH1F *h_both = new TH1F("h1", "Binned MC", 20, -1.,1.);
    TH1F *h_metonly = new TH1F("h2", "Binned MC", 20, -1.,1.);
    TH1F *h_btagonly = new TH1F("h3", "Binned MC", 20, -1.,1.);


    h_both->SetLineColor(kBlue);
    h_both->SetLineWidth(3);

    h_metonly->SetLineColor(kRed);
    h_metonly->SetLineWidth(3);

    h_btagonly->SetLineColor(kRed);
    h_btagonly->SetLineWidth(3);

    
    make_cut_hists(t_elel_mc, h_both, h_metonly, h_btagonly);




    make_ratio_plot("ElEl_met_cut.pdf", h_btagonly, "Before",h_both, "After MET Cut", "Before/After", "cos(#theta_{r})", false);
    make_ratio_plot("ElEl_btag_cut.pdf", h_metonly, "Before",h_both, "After anti-b-tag Cut", "Before/After", "cos(#theta_{r})", false);




}



