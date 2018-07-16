

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

const int type = FLAG_MUONS;



void Fakerate_est_zpeak_mu(TTree *t_WJets, TTree *t_QCD, TTree *t_WJets_contam, TTree *t_QCD_contam, TH1F *h_m){
    FakeRate FR;
    //TH2D *FR;
    setup_new_mu_fakerate(&FR);
    //FR.h->Print();
    for (int l=0; l<=3; l++){
        printf("l=%i\n", l);
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
        TLorentzVector *mu_p = 0;
        TLorentzVector *mu_m = 0;
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
        t->SetBranchAddress("mu_p", &mu_p);
        t->SetBranchAddress("mu_m", &mu_m);
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
            if(l==0){
                if(iso_mu ==0) mu1_fakerate = get_new_fakerate_prob(mu1_pt, mu1_eta, FR.h);
                if(iso_mu ==1) mu1_fakerate = get_new_fakerate_prob(mu2_pt, mu2_eta, FR.h);
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
                if(iso_mu ==0) mu1_fakerate = get_new_fakerate_prob(mu1_pt, mu1_eta, FR.h);
                if(iso_mu ==1) mu1_fakerate = get_new_fakerate_prob(mu2_pt, mu2_eta, FR.h);
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



            if(m>= 105. && met_pt < 50.  && no_bjets){
                //if(l==3) printf("Evt rate %.2e \n", evt_fakerate);
                TLorentzVector cm = *mu_p + *mu_m;

                h_m->Fill(m, evt_fakerate);
            }
        }
        printf("After iter %i current fakerate est is %.0f \n", l, h_m->Integral());
    }
    printf("Total fakerate est is %.0f \n", h_m->Integral());
}

void draw_zpeak(){
    init();
    f_data = TFile::Open("../analyze/output_files/SingleMuon_data_july10.root");
    t_data = (TTree *)f_data->Get("T_data");
    f_mc = TFile::Open("../analyze/output_files/MuMu_DY_unbinned_july10.root");
    //f_mc = TFile::Open("../analyze/output_files/MuMu_DY_unbinned_backbatch_june28.root");
    t_mc = (TTree *)f_mc->Get("T_data");
    t_mc_nosig = (TTree *)f_mc->Get("T_back");
    f_ttbar = TFile::Open("../analyze/output_files/MuMu_TTbar_july10.root");
    t_ttbar = (TTree *)f_ttbar->Get("T_data");
    setTDRStyle();



    
    int n_bins = 23;
    float bins[] = {50., 55., 60., 64., 68., 72., 76., 81., 86., 91., 96., 101., 106., 110., 115., 120., 126., 133., 141., 150., 160., 171., 185., 200};

    //int n_bins = 50;
    //double start = 50;
    //double end = 2000.;
    TH1F *data_m = new TH1F("data_m", "Data Dimuon Mass Distribution", n_bins, bins);
    TH1F *mc_m = new TH1F("mc_m", "MC Signal (qqbar, qglu, qbarglu)", n_bins, bins);
    TH1F *diboson_m = new TH1F("diboson_m", "DiBoson (WW, WZ, ZZ)", n_bins, bins);
    TH1F *QCD_m = new TH1F("QCD_m", "QCD", n_bins, bins);
    TH1F *WJets_m = new TH1F("WJets_m", "WJets", n_bins, bins);
    TH1F *wt_m = new TH1F("wt_m", "tw + #bar{t}w", n_bins, bins);
    TH1F *mc_nosig_m = new TH1F("mc_nosig_m", "MC no signal (qq, gluglu qbarqbar)", n_bins, bins);
    TH1F *ttbar_m = new TH1F("ttbar_m", "TTBar Background", n_bins, bins);
    mc_m->SetFillColor(kRed+1);
    mc_m->SetMarkerStyle(21);
    mc_m->SetMarkerColor(kRed+1);
    mc_nosig_m->SetFillColor(kMagenta);
    mc_nosig_m->SetMarkerStyle(21);
    mc_nosig_m->SetMarkerColor(kMagenta);
    ttbar_m->SetFillColor(kBlue);
    ttbar_m->SetMarkerStyle(21);
    ttbar_m->SetMarkerColor(kBlue);


    TH1F *data_cost = new TH1F("data_cost", "Data Dimuon Mass Distribution", n_bins, bins);
    TH1F *mc_cost = new TH1F("mc_cost", "MC Signal (qqbar, qglu, qbarglu)", n_bins, bins);
    TH1F *diboson_cost = new TH1F("diboson_cost", "DiBoson (WW, WZ, ZZ)", n_bins, bins);
    TH1F *QCD_cost = new TH1F("QCD_cost", "QCD", n_bins, bins);
    TH1F *WJets_cost = new TH1F("WJets_cost", "WJets", n_bins, bins);
    TH1F *wt_cost = new TH1F("wt_cost", "tw + #bar{t}w", n_bins, bins);
    TH1F *mc_nosig_cost = new TH1F("mc_nosig_cost", "MC no signal (qq, gluglu qbarqbar)", n_bins, bins);
    TH1F *ttbar_cost = new TH1F("ttbar_cost", "TTBar Background", n_bins, bins);


    TH1F *data_pt = new TH1F("data_pt", "Data Dimuon Mass Distribution", n_bins, bins);
    TH1F *mc_pt = new TH1F("mc_pt", "MC Signal (qqbar, qglu, qbarglu)", n_bins, bins);
    TH1F *diboson_pt = new TH1F("diboson_pt", "DiBoson (WW, WZ, ZZ)", n_bins, bins);
    TH1F *QCD_pt = new TH1F("QCD_pt", "QCD", n_bins, bins);
    TH1F *WJets_pt = new TH1F("WJets_pt", "WJets", n_bins, bins);
    TH1F *wt_pt = new TH1F("wt_pt", "tw + #bar{t}w", n_bins, bins);
    TH1F *mc_nosig_pt = new TH1F("mc_nosig_pt", "MC no signal (qq, gluglu qbarqbar)", n_bins, bins);
    TH1F *ttbar_pt = new TH1F("ttbar_pt", "TTBar Background", n_bins, bins);


    wt_m->SetFillColor(kOrange+7); 
    diboson_m->SetFillColor(kGreen+3);
    QCD_m->SetFillColor(kRed -7);

    printf("data \n");
    make_m_cost_pt_hist(t_data, data_m, data_cost, data_pt, true, type, true, 50.);
    printf("mc \n");
    make_m_cost_pt_hist(t_mc, mc_m, mc_cost, mc_pt, false, type, true, 50.);
    make_m_cost_pt_hist(t_mc_nosig, mc_nosig_m, mc_nosig_cost, mc_nosig_pt, false, type, true, 50.);
    printf("ttbar \n");
    make_m_cost_pt_hist(t_ttbar, ttbar_m, ttbar_cost, ttbar_pt, false, type, true, 50.);
    printf("wt + diboson \n");
    make_m_cost_pt_hist(t_wt, wt_m, wt_cost, wt_pt, false, type, false, 50.);
    make_m_cost_pt_hist(t_diboson, diboson_m, diboson_cost, diboson_pt, false, type, false, 50.);

    //Fakerate_est_zpeak_mu(t_WJets, t_QCD, t_WJets_mc, t_QCD_mc, QCD_m);




    Double_t EMu_ratio= 0.97;
    ttbar_m->Scale(EMu_ratio);
    diboson_m->Scale(EMu_ratio);
    wt_m->Scale(EMu_ratio);


    int nBins_x = QCD_m->GetXaxis()->GetNbins();
    int nBins_y = QCD_m->GetYaxis()->GetNbins();
    //printf("Get size %i \n", nBins);
    for (int i=1; i <= nBins_x; i++){
        for (int j=1; j <= nBins_y; j++){

            Double_t m_val = QCD_m->GetBinContent(i,j);

            QCD_m->SetBinError(i,j, 0.2*m_val);
        }
    }

    

    THStack *m_stack = new THStack("m_stack", "MuMu Mass Distribution: Data vs MC ; m_{#mu^{+}#mu^{-}} (GeV)");
    m_stack->Add(diboson_m);
    m_stack->Add(wt_m);
    //m_stack->Add(QCD_m);
    m_stack->Add(ttbar_m);
    m_stack->Add(mc_nosig_m);
    m_stack->Add(mc_m);



    TCanvas *c_m = new TCanvas("c_m", "Histograms", 200, 10, 900, 700);
    TPad *pad1 = new TPad("pad1", "pad1", 0.,0.3,0.98,1.);
    pad1->SetBottomMargin(0);
    pad1->Draw();
    pad1->cd();
    pad1->SetLogy();
    m_stack->Draw("hist");
    //m_stack->SetMaximum(65000);
    m_stack->SetMinimum(0.1);
    gStyle->SetEndErrorSize(4);
    data_m->SetMarkerStyle(kFullCircle);
    data_m->SetMarkerColor(1);
    data_m->DrawCopy("P E same");


    gStyle->SetLegendBorderSize(0);
    TLegend *leg1 = new TLegend(0.5, 0.25, 0.75, 0.4);
    leg1->AddEntry(data_m, "data", "p");
    leg1->AddEntry(mc_m, "DY (q#bar{q}, qg #bar{q}g)", "f");
    leg1->AddEntry(mc_nosig_m, "DY no asymmety(gg, qq, #bar{q}#bar{q})", "f");
    leg1->AddEntry(diboson_m, "WW + WZ + ZZ", "f");
    leg1->AddEntry(wt_m, "tW + #bar{t}W", "f");
    //leg1->AddEntry(QCD_m, "QCD + WJets", "f");
    leg1->AddEntry(ttbar_m, "t#bar{t}", "f");
    leg1->Draw();

    //gPad->BuildLegend();
    c_m->cd();
    TPad *pad2 = new TPad("pad2", "pad2", 0.,0,.98,0.3);
    //pad2->SetTopMargin(0);
    pad2->SetBottomMargin(0.2);
    pad2->SetGridy();
    pad2->Draw();
    pad2->cd();



    TList *stackHists = m_stack->GetHists();
    TH1* m_mc_sum = (TH1*)stackHists->At(0)->Clone();
    m_mc_sum->Reset();

    for (int i=0;i<stackHists->GetSize();++i) {
      m_mc_sum->Add((TH1*)stackHists->At(i));
    }
    auto ratio = (TH1F *) data_m->Clone("h_ratio");
    ratio->SetMinimum(0.75);
    ratio->SetMaximum(1.25);
    ratio->Sumw2();
    ratio->SetStats(0);
    ratio->Divide(m_mc_sum);
    ratio->SetMarkerStyle(21);
    ratio->Draw("ep");
    TLine *l1 = new TLine(50,1,2000,1);
    l1->SetLineStyle(2);
    l1->Draw();
    c_m->cd();

    ratio->SetTitle("");
    // Y axis ratio plot settings
   ratio->GetYaxis()->SetTitle("Data/MC");
   ratio->GetYaxis()->SetNdivisions(505);
   ratio->GetYaxis()->SetTitleSize(20);
   ratio->GetYaxis()->SetTitleFont(43);
   ratio->GetYaxis()->SetTitleOffset(1.2);
   ratio->GetYaxis()->SetLabelFont(43); // Absolute font size in pixel (precision 3)
   ratio->GetYaxis()->SetLabelSize(15);
   // X axis ratio plot settings
   ratio->GetXaxis()->SetTitle("M_{#mu#mu} (GeV)");
   ratio->GetXaxis()->SetTitleSize(20);
   ratio->GetXaxis()->SetTitleFont(43);
   ratio->GetXaxis()->SetTitleOffset(3.);
   ratio->GetXaxis()->SetLabelFont(43); // Absolute font size in pixel (precision 3)
   ratio->GetXaxis()->SetLabelSize(20);
 
    //lumi_sqrtS = "";       // used with iPeriod = 0, e.g. for simulation-only plots (default is an empty string)
    int iPeriod = 4; 
    CMS_lumi(pad1, iPeriod, 33 );
    c_m->Print("MuMu_RC_on_zpeak.pdf");


 
}

    
    
