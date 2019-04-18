
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

void make_pt_cmp(TTree *t1, TH1F *h_def_ratio, TH1F *h_alt_ratio, bool is_data=false){
    //read event data
    Long64_t size  =  t1->GetEntries();
    Double_t m, xF, cost, mu1_pt, mu2_pt, jet1_cmva, jet2_cmva, gen_weight;
    Double_t bcdef_HLT_SF, bcdef_iso_SF, bcdef_id_SF;
    Double_t gh_HLT_SF, gh_iso_SF, gh_id_SF, el_id_SF, el_reco_SF, el_HLT_SF;
    Double_t bcdef_trk_SF, gh_trk_SF;
    Double_t jet1_pt, jet2_pt, jet1_b_weight, jet2_b_weight, pu_SF;
    Double_t mu1_pt_corr, mu2_pt_corr, mu1_pt_alt, mu2_pt_alt;
    jet1_b_weight = jet2_b_weight =1.;
    TLorentzVector *mu_p = 0;
    TLorentzVector *mu_m = 0;
    TLorentzVector *el_p = 0;
    TLorentzVector *el_m = 0;
    TLorentzVector cm;
    Float_t met_pt;
    Int_t nJets;
    nJets = 2;
    pu_SF=1;
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
        //t1->SetBranchAddress("jet1_b_weight", &jet1_b_weight);
        //t1->SetBranchAddress("jet2_b_weight", &jet2_b_weight);
        t1->SetBranchAddress("pu_SF", &pu_SF);
    }
    t1->SetBranchAddress("mu_p", &mu_p);
    t1->SetBranchAddress("mu_m", &mu_m);
    t1->SetBranchAddress("mu1_pt", &mu1_pt);
    t1->SetBranchAddress("mu2_pt", &mu2_pt);
    if(!is_data){
        t1->SetBranchAddress("bcdef_HLT_SF", &bcdef_HLT_SF);
        t1->SetBranchAddress("bcdef_iso_SF", &bcdef_iso_SF);
        t1->SetBranchAddress("bcdef_id_SF", &bcdef_id_SF);
        t1->SetBranchAddress("gh_HLT_SF", &gh_HLT_SF);
        t1->SetBranchAddress("gh_iso_SF", &gh_iso_SF);
        t1->SetBranchAddress("gh_id_SF", &gh_id_SF);
        t1->SetBranchAddress("gh_trk_SF", &gh_trk_SF);
        t1->SetBranchAddress("bcdef_trk_SF", &bcdef_trk_SF);
        t1->SetBranchAddress("mu1_pt_corr", &mu1_pt_corr);
        t1->SetBranchAddress("mu2_pt_corr", &mu2_pt_corr);
        t1->SetBranchAddress("mu1_pt_alt", &mu1_pt_alt);
        t1->SetBranchAddress("mu2_pt_alt", &mu2_pt_alt);
    }
    for (int i=0; i<size; i++) {
        t1->GetEntry(i);
        bool no_bjets = has_no_bjets(nJets, jet1_pt, jet2_pt, jet1_cmva, jet2_cmva);

        cm = *mu_p + *mu_m;
        Double_t pt = cm.Pt();
        cost = get_cost_v2(*mu_p, *mu_m);
        if(m >= 400. && met_pt < 50. && no_bjets){

            TLorentzVector mu_p_new, mu_m_new;
            if((mu_p->Pt() > mu_m->Pt() && mu1_pt_corr > mu2_pt_corr) ||
                    (mu_p->Pt() < mu_m->Pt() && mu1_pt_corr < mu2_pt_corr)){
                mu_p_new.SetPtEtaPhiE(mu1_pt, mu_p->Eta(), mu_p->Phi(), mu_p->E());
                mu_m_new.SetPtEtaPhiE(mu2_pt, mu_m->Eta(), mu_m->Phi(), mu_m->E());
            }
            else{
                mu_p_new.SetPtEtaPhiE(mu2_pt, mu_p->Eta(), mu_p->Phi(), mu_p->E());
                mu_m_new.SetPtEtaPhiE(mu1_pt, mu_m->Eta(), mu_m->Phi(), mu_m->E());
            }
            double cost_orig = get_cost_v2(mu_p_new, mu_m_new);
            cm = mu_p_new + mu_m_new;
            //if(cm.M() < 150.) continue;
            double m_orig = cm.M();
            double pt_orig = cm.Pt();

            TLorentzVector mu_p_alt, mu_m_alt;
            if((mu_p->Pt() > mu_m->Pt() && mu1_pt_corr > mu2_pt_corr) ||
                    (mu_p->Pt() < mu_m->Pt() && mu1_pt_corr < mu2_pt_corr)){
                mu_p_alt.SetPtEtaPhiE(mu1_pt_alt, mu_p->Eta(), mu_p->Phi(), mu_p->E());
                mu_m_alt.SetPtEtaPhiE(mu2_pt_alt, mu_m->Eta(), mu_m->Phi(), mu_m->E());
            }
            else{
                mu_p_alt.SetPtEtaPhiE(mu2_pt_alt, mu_p->Eta(), mu_p->Phi(), mu_p->E());
                mu_m_alt.SetPtEtaPhiE(mu1_pt_alt, mu_m->Eta(), mu_m->Phi(), mu_m->E());
            }
            double cost_alt = get_cost_v2(mu_p_alt, mu_m_alt);
            cm = mu_p_alt + mu_m_alt;
            //if(cm.M() < 150.) continue;
            double m_alt = cm.M();
            double pt_alt = cm.Pt();
                   
            double def_ratio1 = mu1_pt_corr/mu1_pt;
            double def_ratio2 = mu2_pt_corr/mu2_pt;

            double alt_ratio1 = mu1_pt_alt/mu1_pt;
            double alt_ratio2 = mu2_pt_alt/mu2_pt;

            double m_alt_ratio = m/m_orig;
            double cost_alt_ratio = cost/cost_orig;

            if(is_data){
                h_def_ratio->Fill(def_ratio1);
                h_def_ratio->Fill(def_ratio2);

                h_alt_ratio->Fill(alt_ratio1);
                h_alt_ratio->Fill(alt_ratio2);
            }
            else{
                Double_t bcdef_weight = gen_weight *pu_SF * bcdef_HLT_SF * bcdef_iso_SF * bcdef_id_SF * bcdef_trk_SF;
                Double_t gh_weight = gen_weight *pu_SF * gh_HLT_SF * gh_iso_SF * gh_id_SF * gh_trk_SF;
                if (nJets >= 1){
                    bcdef_weight *= jet1_b_weight;
                    gh_weight *= jet1_b_weight;
                }
                if (nJets >= 2){
                    bcdef_weight *= jet2_b_weight;
                    gh_weight *= jet2_b_weight;
                }
                Double_t evt_weight = 1000*(bcdef_lumi * bcdef_weight + gh_lumi *gh_weight);
                //h_def_ratio->Fill(def_ratio1, evt_weight);
                //h_def_ratio->Fill(def_ratio2, evt_weight);

                //h_alt_ratio->Fill(alt_ratio1, evt_weight);
                //h_alt_ratio->Fill(alt_ratio2, evt_weight);
                h_alt_ratio->Fill(cost_alt, evt_weight);
                h_def_ratio->Fill(cost, evt_weight);
            }


        }
    }


    t1->ResetBranchAddresses();
}
void make_ratio_plot(char title[80], TH1F* h1, char h1_label[80], TH1F* h2, char h2_label[80], char ratio_label[80], 
        char axis_label[80], bool logy=false){

    h1->SetLineColor(kRed);
    h2->SetLineColor(kBlue);
    h1->SetLineWidth(3);
    h2->SetLineWidth(3);

    TCanvas *c = new TCanvas("c_m", "Histograms", 200, 10, 900, 700);
    TPad *pad1 = new TPad("pad1", "pad1", 0.,0.3,0.98,1.);
    pad1->SetBottomMargin(0);
    pad1->Draw();
    pad1->cd();
    if(logy) pad1->SetLogy();
    h2->Draw("hist E");
    //m_stack->SetMaximum(65000);
    h1->SetMinimum(0.1);
    gStyle->SetEndErrorSize(4);
    h1->Draw("hist E same");


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
    ratio->SetMinimum(0.8);
    ratio->SetMaximum(1.2);
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



void draw_roch_cmp(){

    init();
    setTDRStyle();




    TH1F *h1_cost = new TH1F("data_cost", "Data", 40, -1.,1.);
    TH1F *h1_m = new TH1F("data_m", "Data Dimuon Mass Distribution", 60, 140, 2000);
    TH1F *h1_pt = new TH1F("mc_pt", "MC signal", 40, 0, 1000);
    TH1F *h1_xf = new TH1F("msdjf_pt", "MC signal", 40, 0, 1000);

    TH1F *h2_cost = new TH1F("h2_cost", "Data", 40, -1.,1.);
    TH1F *h2_m = new TH1F("h2_m", "Data Dimuon Mass Distribution", 60, 140, 2000);
    TH1F *h2_pt = new TH1F("h2_pt", "MC signal", 40, 0, 1000);
    TH1F *h2_xf = new TH1F("hx_pt", "MC signal", 40, 0, 1000);

    TH1F *h_def_ratio = new TH1F("d", "Roch default ratio", 40, 0.5, 1.5);
    TH1F *h_alt_ratio = new TH1F("a", "Roch alt ratio", 40, 0.5, 1.5);

    bool is_data = false;
    int type = FLAG_MUONS;
    make_m_cost_pt_xf_hist(t_mumu_mc, h1_m, h1_cost, h1_pt, h1_xf, is_data, type, true);
    make_m_cost_pt_xf_hist(t_mumu_mc, h2_m, h2_cost, h2_pt, h2_xf, is_data, type, false);

    //make_pt_cmp(t1, h1_cost, h2_cost);
    //TCanvas *c1 = new TCanvas("c1", "", 200, 10, 900, 700);
    //h_def_ratio->Draw("hist");
    //TCanvas *c2 = new TCanvas("c2", "", 200, 10, 900, 700);
    //h_alt_ratio->Draw("hist");
    //make_ratio_plot("MuMu_mc_rocc_alt_cost_cmp.pdf", h1_cost, "RC def ",h2_cost, "RC Alt", "OFF/On", "Cos(#theta)", false);

    printf("ON integral is %.2f. OFF integral is %.2f \n", h1_m->Integral(), h2_m->Integral());

    make_ratio_plot("MuMu_mc_roch_cost_cmp.pdf", h2_cost, "RC OFF ",h1_cost, "RC ON", "OFF/On", "Cos(#theta)", false);
    make_ratio_plot("MuMu_mc_roch_m_cmp.pdf", h2_m, "RC OFF ",h1_m, "RC ON", "OFF/On", "M (GeV)", true);
}

