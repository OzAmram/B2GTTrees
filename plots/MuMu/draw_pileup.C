
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
#include "TRatioPlot.h"
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
#include "../tdrstyle.C"
#include "../CMS_lumi.C"
#include "root_files.h"

const int type = FLAG_MUONS;

void make_pileup_hist(TTree *t1, TH1F *h_before, TH1F *h_after, bool is_data=false, int flag1 = FLAG_MUONS){
    //read event data
    Long64_t size  =  t1->GetEntries();
    Double_t m, xF, cost, mu1_pt, mu2_pt, jet1_cmva, jet2_cmva, gen_weight;
    Double_t bcdef_HLT_SF, bcdef_iso_SF, bcdef_id_SF;
    Double_t gh_HLT_SF, gh_iso_SF, gh_id_SF, el_id_SF, el_reco_SF;
    Double_t jet1_pt, jet2_pt, jet1_b_weight, jet2_b_weight, pu_SF;
    jet1_b_weight = jet2_b_weight = 1.0;
    TLorentzVector *mu_p = 0;
    TLorentzVector *mu_m = 0;
    TLorentzVector *el_p = 0;
    TLorentzVector *el_m = 0;
    TLorentzVector cm;
    Float_t met_pt;
    Int_t nJets;
    Int_t pu_NtrueInt;
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
    t1->SetBranchAddress("pu_NtrueInt", &pu_NtrueInt);
    if(!is_data){
        t1->SetBranchAddress("nJets", &nJets);
        t1->SetBranchAddress("gen_weight", &gen_weight);
        t1->SetBranchAddress("jet1_b_weight", &jet1_b_weight);
        t1->SetBranchAddress("jet2_b_weight", &jet2_b_weight);
        t1->SetBranchAddress("pu_SF", &pu_SF);
    }
    if(flag1 == FLAG_MUONS){
        t1->SetBranchAddress("mu_p", &mu_p);
        t1->SetBranchAddress("mu_m", &mu_m);
        if(!is_data){
            t1->SetBranchAddress("bcdef_HLT_SF", &bcdef_HLT_SF);
            t1->SetBranchAddress("bcdef_iso_SF", &bcdef_iso_SF);
            t1->SetBranchAddress("bcdef_id_SF", &bcdef_id_SF);
            t1->SetBranchAddress("gh_HLT_SF", &gh_HLT_SF);
            t1->SetBranchAddress("gh_iso_SF", &gh_iso_SF);
            t1->SetBranchAddress("gh_id_SF", &gh_id_SF);
        }
        for (int i=0; i<size; i++) {
            t1->GetEntry(i);
            bool no_bjets = has_no_bjets(nJets, jet1_pt, jet2_pt, jet1_cmva, jet2_cmva);

            if(m >= 150. && met_pt < 50. && no_bjets){
                cm = *mu_p + *mu_m;
                Double_t pt = cm.Pt();
                if(is_data){
                    h_before->Fill(pu_NtrueInt);
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
                    Double_t evt_weight = (bcdef_weight*bcdef_lumi + gh_weight*gh_lumi) *1000;
                    h_before->Fill(pu_NtrueInt, evt_weight);
                    h_after->Fill(pu_NtrueInt, evt_weight *pu_SF);

                }


            }
        }
    }

    if(is_data) printf("N data events is %.0f \n", h_before->Integral());
    t1->ResetBranchAddresses();
}


void draw_pileup(){
    setTDRStyle();

    TFile *f7 = TFile::Open("../analyze/SFs/DataPileupHistogram_69200.root");
    TH1D *data_pu = (TH1D *) f7->Get("pileup")->Clone();
    data_pu->Scale(1./data_pu->Integral());
    data_pu->SetDirectory(0);

    init();
    mumu_init();
    TH1F *mc_pu_before = new TH1F("mc_pu_before", "MC signal", 100, 0, 100);
    TH1F *mc_nosig_pu_before = new TH1F("mc_nosig_pu_before", "MC signal", 100, 0, 100);
    TH1F *ttbar_pu_before = new TH1F("ttbar_pu_before", "MC signal", 100, 0, 100);
    TH1F *diboson_pu_before = new TH1F("diboson_pu_before", "MC signal", 100, 0, 100);
    TH1F *wt_pu_before = new TH1F("wt_pu_before", "MC signal", 100, 0, 100);

    TH1F *mc_pu_after = new TH1F("mc_pu_after", "MC signal", 100, 0, 100);
    TH1F *mc_nosig_pu_after = new TH1F("mc_nosig_pu_after", "MC signal", 100, 0, 100);
    TH1F *ttbar_pu_after = new TH1F("ttbar_pu_after", "MC signal", 100, 0, 100);
    TH1F *diboson_pu_after = new TH1F("diboson_pu_after", "MC signal", 100, 0, 100);
    TH1F *wt_pu_after = new TH1F("wt_pu_after", "MC signal", 100, 0, 100);

    mc_pu_before->SetFillColor(kRed+1);
    mc_pu_before->SetMarkerColor(kRed+1);
    mc_nosig_pu_before->SetFillColor(kMagenta);
    ttbar_pu_before->SetFillColor(kBlue);
    ttbar_pu_before->SetMarkerStyle(21);
    ttbar_pu_before->SetMarkerColor(kBlue);
    diboson_pu_before->SetFillColor(kGreen+3);
    wt_pu_before->SetFillColor(kOrange+7); 

    mc_pu_after->SetFillColor(kRed+1);
    mc_pu_after->SetMarkerColor(kRed+1);
    mc_nosig_pu_after->SetFillColor(kMagenta);
    ttbar_pu_after->SetFillColor(kBlue);
    ttbar_pu_after->SetMarkerStyle(21);
    ttbar_pu_after->SetMarkerColor(kBlue);
    diboson_pu_after->SetFillColor(kGreen+3);
    wt_pu_after->SetFillColor(kOrange+7); 







    make_pileup_hist(t_mumu_mc, mc_pu_before, mc_pu_after, false, type);
    make_pileup_hist(t_mumu_nosig, mc_nosig_pu_before, mc_nosig_pu_after, false, type);
    make_pileup_hist(t_ttbar, ttbar_pu_before, ttbar_pu_after, false, type);
    make_pileup_hist(t_wt, wt_pu_before, wt_pu_after, false);
    make_pileup_hist(t_diboson, diboson_pu_before, diboson_pu_after, false, type);









    Double_t pu_before_tot = ttbar_pu_before->Integral() + wt_pu_before->Integral() + diboson_pu_before->Integral() + 
                            mc_nosig_pu_before->Integral() + mc_pu_before->Integral();

    ttbar_pu_before->Scale(1./pu_before_tot);
    wt_pu_before->Scale(1./pu_before_tot);
    diboson_pu_before->Scale(1./pu_before_tot);
    mc_nosig_pu_before->Scale(1./pu_before_tot);
    mc_pu_before->Scale(1./pu_before_tot);

    Double_t pu_after_tot = ttbar_pu_after->Integral() + wt_pu_after->Integral() + diboson_pu_after->Integral() + 
                            mc_nosig_pu_after->Integral() + mc_pu_after->Integral();

    ttbar_pu_after->Scale(1./pu_after_tot);
    wt_pu_after->Scale(1./pu_after_tot);
    diboson_pu_after->Scale(1./pu_after_tot);
    mc_nosig_pu_after->Scale(1./pu_after_tot);
    mc_pu_after->Scale(1./pu_after_tot);

    printf("%.3e %.3e %.3e \n", mc_nosig_pu_after, mc_pu_after, ttbar_pu_after);

    THStack *pu_before_stack = new THStack("pu_before_stack", "DiMuon Pileup; # Primary Vertices");
    pu_before_stack->Add(ttbar_pu_before);
    pu_before_stack->Add(wt_pu_before);
    pu_before_stack->Add(diboson_pu_before);
    pu_before_stack->Add(mc_nosig_pu_before);
    pu_before_stack->Add(mc_pu_before);

    THStack *pu_after_stack = new THStack("pu_after_stack", "Dimuon Pileup; # Primary Vertices");
    pu_after_stack->Add(ttbar_pu_after);
    pu_after_stack->Add(wt_pu_after);
    pu_after_stack->Add(diboson_pu_after);
    pu_after_stack->Add(mc_nosig_pu_after);
    pu_after_stack->Add(mc_pu_after);


    TCanvas *c_pu_before = new TCanvas("c_pu_before", "Histograms", 200, 10, 900, 700);
    c_pu_before->cd();
    TPad *before_pad1 = new TPad("pad1", "pad1", 0.,0.3,0.98,1.);
    before_pad1->SetBottomMargin(0);
    before_pad1->Draw();
    before_pad1->cd();
    pu_before_stack->Draw("hist");
    pu_before_stack->GetXaxis()->SetRangeUser(0.,50.);
    data_pu->SetMarkerStyle(kFullCircle);
    data_pu->SetMarkerColor(1);
    //pu_before_stack->SetMinimum(1);
    //pu_before_stack->SetMaximum(100000);
    //data_pu->SetMinimum(1);
    //data_pu->SetMaximum(100000);
    data_pu->Draw("P same");
    before_pad1->SetLogy();
    before_pad1->Update();
    TLegend *leg3 = new TLegend(0.5, 0.65, 0.75, 0.8);
    leg3->AddEntry(data_pu, "data", "p");
    leg3->AddEntry(mc_pu_before, "DY (q#bar{q}, qg #bar{q}g)", "f");
    leg3->AddEntry(mc_nosig_pu_before, "DY no asymmety(gg, qq, #bar{q}#bar{q})", "f");
    leg3->AddEntry(diboson_pu_before, "WW + WZ + ZZ", "f");
    leg3->AddEntry(wt_pu_before, "tW + #bar{t}W", "f");
    leg3->AddEntry(ttbar_pu_before, "t#bar{t}", "f");
    leg3->Draw();
    c_pu_before->cd();

    TPad *before_pad2 = new TPad("before_pad2", "pad2", 0.,0,.98,0.3);
    //pad2->SetTopMargin(0);
    before_pad2->SetBottomMargin(0.2);
    before_pad2->SetGridy();
    before_pad2->Draw();
    before_pad2->cd();
    TList *before_stackHists = pu_before_stack->GetHists();
    TH1* before_mc_sum = (TH1*) before_stackHists->At(0)->Clone();
    before_mc_sum->Reset();

    for (int i=0;i<before_stackHists->GetSize();++i) {
      before_mc_sum->Add((TH1*)before_stackHists->At(i));
    }
    auto before_ratio = (TH1F *) data_pu->Clone("h_before_ratio");
    before_ratio->SetMinimum(0.7);
    before_ratio->SetMaximum(1.3);
    before_ratio->Sumw2();
    before_ratio->SetStats(0);
    before_ratio->Divide(before_mc_sum);
    before_ratio->SetMarkerStyle(21);
    before_ratio->Draw("ep");
    c_pu_before->cd();

    before_ratio->SetTitle("");
    // Y axis before_ratio plot settings
   before_ratio->GetYaxis()->SetTitle("Data/MC");
   before_ratio->GetYaxis()->SetNdivisions(505);
   before_ratio->GetYaxis()->SetTitleSize(20);
   before_ratio->GetYaxis()->SetTitleFont(43);
   before_ratio->GetYaxis()->SetTitleOffset(1.2);
   before_ratio->GetYaxis()->SetLabelFont(43); // Absolute font size in pixel (precision 3)
   before_ratio->GetYaxis()->SetLabelSize(15);
   // X axis before_ratio plot settings
   before_ratio->GetXaxis()->SetTitle("Number of vertices before reweighting");
   before_ratio->GetXaxis()->SetRangeUser(0.,50.);
   before_ratio->GetXaxis()->SetTitleSize(20);
   before_ratio->GetXaxis()->SetTitleFont(43);
   before_ratio->GetXaxis()->SetTitleOffset(3.);
   before_ratio->GetXaxis()->SetLabelFont(43); // Absolute font size in pixel (precision 3)
   before_ratio->GetXaxis()->SetLabelSize(20);

    int iPeriod = 4; 
    CMS_lumi(c_pu_before, iPeriod, 11 );
    c_pu_before->Update();

    TCanvas *c_pu_after = new TCanvas("c_pu_after", "Histograms", 200, 10, 900, 700);
    c_pu_after->cd();
    TPad *after_pad1 = new TPad("pad1", "pad1", 0.,0.3,0.98,1.);
    after_pad1->SetBottomMargin(0);
    after_pad1->Draw();
    after_pad1->cd();
    pu_after_stack->Draw("hist");
    pu_after_stack->GetXaxis()->SetRangeUser(0.,50.);
    data_pu->SetMarkerStyle(kFullCircle);
    data_pu->SetMarkerColor(1);
    //pu_after_stack->SetMinimum(1);
    //pu_after_stack->SetMaximum(100000);
    //data_pu->SetMinimum(1);
    //data_pu->SetMaximum(100000);
    data_pu->Draw("P  same");
    after_pad1->SetLogy();
    after_pad1->Update();
    TLegend *leg4 = new TLegend(0.5, 0.65, 0.75, 0.8);
    leg4->AddEntry(data_pu, "data", "p");
    leg4->AddEntry(mc_pu_after, "DY (q#bar{q}, qg #bar{q}g)", "f");
    leg4->AddEntry(mc_nosig_pu_after, "DY no asymmety(gg, qq, #bar{q}#bar{q})", "f");
    leg4->AddEntry(diboson_pu_after, "WW + WZ + ZZ", "f");
    leg4->AddEntry(wt_pu_after, "tW + #bar{t}W", "f");
    leg4->AddEntry(ttbar_pu_after, "t#bar{t}", "f");
    leg4->Draw();


    CMS_lumi(c_pu_after, iPeriod, 11 );
    c_pu_after->cd();


    TPad *after_pad2 = new TPad("after_pad2", "pad2", 0.,0,.98,0.3);
    //pad2->SetTopMargin(0);
    after_pad2->SetBottomMargin(0.2);
    after_pad2->SetGridy();
    after_pad2->Draw();
    after_pad2->cd();
    TList *after_stackHists = pu_after_stack->GetHists();
    TH1* after_mc_sum = (TH1*)after_stackHists->At(0)->Clone();
    after_mc_sum->Reset();

    for (int i=0;i<after_stackHists->GetSize();++i) {
      after_mc_sum->Add((TH1*)after_stackHists->At(i));
    }
    auto after_ratio = (TH1F *) data_pu->Clone("h_after_ratio");
    after_ratio->SetMinimum(0.7);
    after_ratio->SetMaximum(1.3);
    after_ratio->Sumw2();
    after_ratio->SetStats(0);
    after_ratio->Divide(after_mc_sum);
    after_ratio->SetMarkerStyle(21);
    after_ratio->Draw("ep");
    c_pu_after->cd();

    after_ratio->SetTitle("Number of Vertices after reweighting");
    // Y axis after_ratio plot settings
   after_ratio->GetYaxis()->SetTitle("Data/MC");
   after_ratio->GetYaxis()->SetNdivisions(505);
   after_ratio->GetYaxis()->SetTitleSize(20);
   after_ratio->GetYaxis()->SetTitleFont(43);
   after_ratio->GetYaxis()->SetTitleOffset(1.2);
   after_ratio->GetYaxis()->SetLabelFont(43); // Absolute font size in pixel (precision 3)
   after_ratio->GetYaxis()->SetLabelSize(15);
   // X axis after_ratio plot settings
   after_ratio->GetXaxis()->SetTitle("Pileup after reweighting");
   after_ratio->GetXaxis()->SetRangeUser(0.,50.);
   after_ratio->GetXaxis()->SetTitleSize(20);
   after_ratio->GetXaxis()->SetTitleFont(43);
   after_ratio->GetXaxis()->SetTitleOffset(3.);
   after_ratio->GetXaxis()->SetLabelFont(43); // Absolute font size in pixel (precision 3)
   after_ratio->GetXaxis()->SetLabelSize(20);



}



