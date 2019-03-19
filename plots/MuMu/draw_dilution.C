
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





void draw_dilution(){
    setTDRStyle();
    const int nBins = 10;
    TH1F *h_xF = new TH1F("h_xf", "xF distribution, M>150", nBins, 0, 0.5);
    TH1F *h_Nc = new TH1F("h_Nc", "Number Correct; |xF|", nBins, 0, 0.5);
    TH1F *h_Ni = new TH1F("h_Ni", "Number Incorrect", nBins, 0, 0.5);


    //read event data
    mumu_init();

    Long64_t size  =  t_mc->GetEntries();

    Double_t m, xF, cost, mu1_pt, mu2_pt, jet1_cmva, jet2_cmva, gen_weight, cost_st;
    Double_t bcdef_HLT_SF, bcdef_iso_SF, bcdef_id_SF;
    Double_t gh_HLT_SF, gh_iso_SF, gh_id_SF;
    Double_t jet1_pt, jet2_pt, jet1_b_weight=1., jet2_b_weight=1.;
    Float_t met_pt;
    Int_t nJets;
    nJets = 2;

    t_mc->SetBranchAddress("m", &m);
    t_mc->SetBranchAddress("xF", &xF);
    t_mc->SetBranchAddress("cost", &cost);
    t_mc->SetBranchAddress("cost_st", &cost_st);
    t_mc->SetBranchAddress("met_pt", &met_pt);
    t_mc->SetBranchAddress("jet2_CMVA", &jet2_cmva);
    t_mc->SetBranchAddress("jet1_CMVA", &jet1_cmva);
    t_mc->SetBranchAddress("jet1_pt", &jet1_pt);
    t_mc->SetBranchAddress("jet2_pt", &jet2_pt);
    t_mc->SetBranchAddress("nJets", &nJets);
    t_mc->SetBranchAddress("gen_weight", &gen_weight);
    t_mc->SetBranchAddress("bcdef_HLT_SF", &bcdef_HLT_SF);
    t_mc->SetBranchAddress("bcdef_iso_SF", &bcdef_iso_SF);
    t_mc->SetBranchAddress("bcdef_id_SF", &bcdef_id_SF);
    t_mc->SetBranchAddress("gh_HLT_SF", &gh_HLT_SF);
    t_mc->SetBranchAddress("gh_iso_SF", &gh_iso_SF);
    t_mc->SetBranchAddress("gh_id_SF", &gh_id_SF);
    t_mc->SetBranchAddress("jet1_b_weight", &jet1_b_weight);
    t_mc->SetBranchAddress("jet2_b_weight", &jet2_b_weight);

    for (int i=0; i<size; i++) {
        t_mc->GetEntry(i);
        bool no_bjets = has_no_bjets(nJets, jet1_pt, jet2_pt, jet1_cmva, jet2_cmva);

        if(m >= 150. && met_pt < 50. && no_bjets){
            //printf("%0.2f \n", xF);
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
            Double_t final_weight = (bcdef_lumi*bcdef_weight + gh_lumi * gh_weight)/(bcdef_weight + gh_weight);

            h_xF->Fill(xF, final_weight);
            Double_t ratio = cost_st/cost;
            if(ratio > 0) h_Nc->Fill(xF, final_weight);
            if(ratio < 0) h_Ni->Fill(xF, final_weight);


        }
    }

    TCanvas *c1 = new TCanvas("c1", "canvas", 200,10, 900,700);
    h_xF->SetFillColor(kBlue);
    h_xF->Draw();
    c1->Update();

    TCanvas *c2 = new TCanvas("c2", "canva", 100,100, 700,700);
    h_Nc ->SetLineColor(kBlue);
    h_Nc -> SetTitle("Guessing Lepton Pair Direction as Incident Quark Direction");
    h_Nc->SetStats(kFALSE);
    h_Nc ->SetLineWidth(2);
    h_Ni ->SetLineColor(kRed);
    h_Ni ->SetLineWidth(2);
    h_Nc ->Draw("hist");
    h_Ni ->Draw("hist same");

    TLegend *leg1 = new TLegend(0.5, 0.65, 0.75, 0.8);
    leg1->AddEntry(h_Nc, "c_{r} Correct Sign", "f");
    leg1->AddEntry(h_Ni, "c_{r} Incorrect Sign", "f");
    leg1->Draw();
    c2->Update();

    Double_t Nc, Ni, dilu[nBins], dilu_error[nBins], xF_center[nBins]; 

    for (int i=1; i < nBins; i++){
        Nc = h_Nc->GetBinContent(i);
        Ni = h_Ni->GetBinContent(i);
        
        xF_center[i] = h_Nc->GetBinCenter(i);
        dilu[i] = (Nc - Ni)/(Nc + Ni);
        dilu_error[i] = sqrt(pow(1./(Nc+Ni) - Nc/(Nc+Ni)/(Nc+Ni), 2)*Nc + 
                pow(1./(Nc+Ni) - Ni/(Nc+Ni)/(Nc+Ni), 2)*Ni);
        printf("Num events %.1f, xf %.2f Dilu error %1.3e \n", Nc+Ni, xF_center[i],
                dilu_error[i]);
    }

    TGraphErrors *g_dillu = new TGraphErrors(nBins, xF_center, dilu, 0, dilu_error);

    TCanvas *c3 = new TCanvas("c3", "canvas", 200,10, 900,700);
    g_dillu->Draw("A C P");
    g_dillu->SetMarkerStyle(20);
    g_dillu->SetTitle("Dilution Effect");
    g_dillu->GetXaxis()->SetTitle("|xF|");
    g_dillu->GetYaxis()->SetTitle("Dilution Factor");

    c3->Update();


    /*
    lumiTextSize     = 0.2;
    lumiTextOffset   = 0.2;
    cmsTextSize      = 0.35;
    cmsTextOffset    = 0.2;  // only used in outOfFrame version
    */
    writeExtraText = true;
    extraText = "Simulation";
    lumi_sqrtS = "";       // used with iPeriod = 0, e.g. for simulation-only plots (default is an empty string)
    int iPeriod = 0; 
    CMS_lumi( c3, iPeriod, 11 );
    CMS_lumi( c2, iPeriod, 33 );




}
