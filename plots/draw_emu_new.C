
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
#include "../analyze/TemplateMaker.C"
#include "tdrstyle.C"
#include "CMS_lumi.C"



Double_t count_tree(TTree *t1,  TH1F* h, Double_t m_low, Double_t m_high, bool is_data=false){
    //count events in the tree
    Long64_t size  =  t1->GetEntries();
    Int_t nJets;
    Double_t m, xF, cost, mu1_pt, mu2_pt, jet1_cmva, jet2_cmva, gen_weight;
    Double_t jet1_pt, jet2_pt;
    Double_t bcdef_HLT_SF, bcdef_iso_SF, bcdef_id_SF;
    Double_t gh_HLT_SF, gh_iso_SF, gh_id_SF, el_id_SF;
    Double_t jet1_b_weight, jet2_b_weight;
    Float_t cost_pt, met_pt;
    TLorentzVector *el=0; 
    TLorentzVector *mu=0;
    TLorentzVector cm;
    t1->SetBranchAddress("met_pt", &met_pt);
    t1->SetBranchAddress("nJets", &nJets);
    t1->SetBranchAddress("jet2_CMVA", &jet2_cmva);
    t1->SetBranchAddress("jet1_CMVA", &jet1_cmva);
    t1->SetBranchAddress("jet1_pt", &jet1_pt);
    t1->SetBranchAddress("jet2_pt", &jet2_pt);
    t1->SetBranchAddress("el", &el);
    t1->SetBranchAddress("mu", &mu);
    if(!is_data){
        t1->SetBranchAddress("gen_weight", &gen_weight);
        t1->SetBranchAddress("bcdef_HLT_SF", &bcdef_HLT_SF);
        t1->SetBranchAddress("bcdef_iso_SF", &bcdef_iso_SF);
        t1->SetBranchAddress("bcdef_id_SF", &bcdef_id_SF);
        t1->SetBranchAddress("gh_HLT_SF", &gh_HLT_SF);
        t1->SetBranchAddress("gh_iso_SF", &gh_iso_SF);
        t1->SetBranchAddress("gh_id_SF", &gh_id_SF);
        t1->SetBranchAddress("el_id_SF", &el_id_SF);
        t1->SetBranchAddress("jet1_b_weight", &jet1_b_weight);
        t1->SetBranchAddress("jet2_b_weight", &jet2_b_weight);
    }
    jet1_cmva = -1.;
    jet2_cmva = -1.;
    Double_t count = 0;
    Double_t bcdef_count = 0;
    Double_t gh_count = 0;

    Double_t med_btag = 0.4432;
    for (int i=0; i<size; i++) {
        t1->GetEntry(i);
        cm = *el + *mu;
        m = cm.M();
        bool no_bjets = has_no_bjets(nJets, jet1_pt, jet2_pt, jet1_cmva, jet2_cmva);
        if(m > m_low && m < m_high && no_bjets && met_pt < 50.){
            if(is_data){
                count += 1;
                h->Fill(m);
            }
            else{
                Double_t bcdef_weight = gen_weight * bcdef_HLT_SF * bcdef_iso_SF * bcdef_id_SF
                                        *el_id_SF;
                Double_t gh_weight = gen_weight * gh_HLT_SF * gh_iso_SF * gh_id_SF
                                    *el_id_SF;
                if (nJets >= 1){
                    bcdef_weight *= jet1_b_weight;
                    gh_weight *= jet1_b_weight;
                }
                if (nJets >= 2){
                    bcdef_weight *= jet2_b_weight;
                    gh_weight *= jet2_b_weight;
                }
                //printf("%.2e %.2e \n", bcdef_weight, gh_weight);
                bcdef_count += bcdef_weight;
                gh_count += gh_weight;
                h->Fill(m, bcdef_weight*bcdef_lumi*1000 + gh_weight*gh_lumi*1000);
            }

            
        }
    }
    if(!is_data){
        bcdef_count *= bcdef_lumi*1000;
        gh_count *= gh_lumi*1000;
        count = bcdef_count + gh_count;
        //count = count * (bcdef_lumi + gh_lumi)*1000;
    }

    //printf("%f \n", count);
    return count;
}


void draw_emu_new(){
    setTDRStyle();
    TFile *f_data = TFile::Open("../analyze/output_files/EMu_data_jun21.root");
    TTree *t_data = (TTree *)f_data->Get("T_data");

                                
    TFile *f_ttbar = TFile::Open("../analyze/output_files/EMu_ttbar_jun26.root");
    TTree *t_ttbar = (TTree *)f_ttbar->Get("T_data");

    TFile *f_DYToLL = TFile::Open("../analyze/output_files/EMu_DY_jun26.root");
    TTree *t_DYToLL = (TTree *)f_DYToLL->Get("T_data");

    TFile *f_diboson = TFile::Open("../analyze/output_files/EMu_diboson_jun26.root");
    TTree *t_diboson = (TTree *)f_diboson->Get("T_data");

    TFile *f_wt = TFile::Open("../analyze/output_files/EMu_WT_jun26.root");
    TTree *t_wt = (TTree *)f_wt->Get("T_data");

    TH1F *data_m = new TH1F("data_m", "MC Signal (qqbar, qglu, qbarglu)", 30, 150, 1000);
    TH1F *ttbar_m = new TH1F("ttbar_m", "MC Signal (qqbar, qglu, qbarglu)", 30, 150, 1000);
    TH1F *diboson_m = new TH1F("diboson_m", "MC Signal (qqbar, qglu, qbarglu)", 30, 150, 1000);
    TH1F *wt_m = new TH1F("wt_m", "MC Signal (qqbar, qglu, qbarglu)", 30, 150, 1000);
    TH1F *dy_m = new TH1F("dy_m", "MC Signal (qqbar, qglu, qbarglu)", 30, 150, 1000);

    Double_t m_low = 150;
    Double_t m_high = 10000;
    Double_t data_count = count_tree(t_data, data_m, m_low, m_high, true);
    Double_t ttbar_count = count_tree(t_ttbar, ttbar_m, m_low, m_high);
    Double_t diboson_count = count_tree(t_diboson, diboson_m, m_low, m_high) ;
    Double_t wt_count = count_tree(t_wt, wt_m, m_low, m_high);
    Double_t DYToLL_count = count_tree(t_DYToLL, dy_m, m_low, m_high);

    printf("M low = %.0f, M high = %.0f \n", m_low, m_high);
    printf("Data count %.0f \n", data_count);
    printf("TTbar count %.0f \n", ttbar_count);
    printf("Diboson count %.0f \n", diboson_count);
    printf("wt count %.0f \n", wt_count);
    printf("DYToLL count %.0f \n", DYToLL_count);
    //Double_t ratio = (data_count - diboson_count - DYToLL_count - wt_count )/ttbar_count;
    //Double_t unc = sqrt((data_count + DYToLL_count + wt_count) *(1/ttbar_count/ttbar_count) 
            //+ ratio*ratio/ttbar_count);
    Double_t denom = ttbar_count + diboson_count + wt_count;
    Double_t num = data_count - DYToLL_count;
    Double_t ratio = (data_count - DYToLL_count)/(ttbar_count + diboson_count + wt_count);
    Double_t unc = sqrt((data_count + DYToLL_count) *pow(1/denom,2) +
            (ttbar_count + diboson_count + wt_count)*pow(num/denom/denom,2));
    printf("Ratio=(data - DYToLL)/(TTbar + tW + diboson) is %1.3f +/- %1.3f \n", ratio, unc);



    dy_m->SetFillColor(kRed+1);
    ttbar_m->SetFillColor(kBlue);
    wt_m->SetFillColor(kOrange+7); 
    diboson_m->SetFillColor(kGreen+3);

    THStack *m_stack = new THStack("m_stack", "EMu Mass Distribution: Data vs MC ; m_{e#mu} (GeV)");
    m_stack->Add(ttbar_m);
    m_stack->Add(wt_m);
    m_stack->Add(diboson_m);
    m_stack->Add(dy_m);



    TCanvas *c_m = new TCanvas("c_m", "Histograms", 200, 10, 900, 700);
    TPad *pad1 = new TPad("pad1", "pad1", 0.,0.3,0.98,1.);
    pad1->SetBottomMargin(0);
    pad1->Draw();
    pad1->cd();
    pad1->SetLogy();
    m_stack->Draw("hist");
    m_stack->SetMaximum(65000);
    gStyle->SetEndErrorSize(4);
    data_m->SetMarkerStyle(kFullCircle);
    data_m->SetMarkerColor(1);
    data_m->DrawCopy("P E same");


    gStyle->SetLegendBorderSize(0);
    TLegend *leg1 = new TLegend(0.5, 0.65, 0.75, 0.8);
    leg1->AddEntry(data_m, "data", "p");
    leg1->AddEntry(dy_m, "DY #rightarrow #tau#tau", "f");
    leg1->AddEntry(diboson_m, "WW + WZ + ZZ", "f");
    leg1->AddEntry(wt_m, "tW + #bar{t}W", "f");
    leg1->AddEntry(ttbar_m, "t#bar{t}", "f");
    leg1->Draw();

    //gPad->BuildLegend();
    c_m->cd();
    TPad *pad2 = new TPad("pad2", "pad2", 0.,0,.98,0.3);
    //pad2->SetTopMargin(0);
    pad2->SetBottomMargin(0.2);
    pad2->Draw();
    pad2->cd();



    TList *stackHists = m_stack->GetHists();
    TH1* m_mc_sum = (TH1*)stackHists->At(0)->Clone();
    m_mc_sum->Reset();

    for (int i=0;i<stackHists->GetSize();++i) {
      m_mc_sum->Add((TH1*)stackHists->At(i));
    }
    auto h_ratio = (TH1F *) data_m->Clone("h_ratio");
    h_ratio->SetMinimum(0.75);
    h_ratio->SetMaximum(1.25);
    h_ratio->Sumw2();
    h_ratio->SetStats(0);
    h_ratio->Divide(m_mc_sum);
    h_ratio->SetMarkerStyle(21);
    h_ratio->Draw("ep");
    TLine *l1 = new TLine(150,1,1000,1);
    l1->SetLineStyle(2);
    l1->Draw();
    c_m->cd();

    h_ratio->SetTitle("");
    // Y axis h_ratio plot settings
   h_ratio->GetYaxis()->SetTitle("Data/MC");
   h_ratio->GetYaxis()->SetNdivisions(505);
   h_ratio->GetYaxis()->SetTitleSize(20);
   h_ratio->GetYaxis()->SetTitleFont(43);
   h_ratio->GetYaxis()->SetTitleOffset(1.2);
   h_ratio->GetYaxis()->SetLabelFont(43); // Absolute font size in pixel (precision 3)
   h_ratio->GetYaxis()->SetLabelSize(15);
   // X axis h_ratio plot settings
   h_ratio->GetXaxis()->SetTitle("M_{e#mu} (GeV)");
   h_ratio->GetXaxis()->SetTitleSize(20);
   h_ratio->GetXaxis()->SetTitleFont(43);
   h_ratio->GetXaxis()->SetTitleOffset(3.);
   h_ratio->GetXaxis()->SetLabelFont(43); // Absolute font size in pixel (precision 3)
   h_ratio->GetXaxis()->SetLabelSize(20);
 
    writeExtraText = true;
    extraText = "Preliminary";
    //lumi_sqrtS = "";       // used with iPeriod = 0, e.g. for simulation-only plots (default is an empty string)
    int iPeriod = 4; 
    CMS_lumi(pad1, iPeriod, 11 );
}









