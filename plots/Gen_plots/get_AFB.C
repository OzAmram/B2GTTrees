
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



void make_hists_from_tree(TTree *t1, TH1F *h_m, TH1F *h_cost, float m_low=150., float m_high = 200.){
    //read event data
    h_m->Sumw2();
    h_cost->Sumw2();
    Long64_t size  =  t1->GetEntries();
    int nSelected=0;
    double root2 = sqrt(2);
    
    Int_t lep1_id, lep2_id;
    TLorentzVector *lep_pls = 0;
    TLorentzVector *lep_mns = 0;
    TLorentzVector cm;
    t1->SetBranchAddress("lep_pls", &lep_pls);
    t1->SetBranchAddress("lep_mns", &lep_mns);
    t1->SetBranchAddress("lep1_id", &lep1_id);
    t1->SetBranchAddress("lep2_id", &lep2_id);
    for (int i=0; i<size; i++) {
        t1->GetEntry(i);
        cm = *lep_pls + *lep_mns;
        //printf("M= %.2f \n", cm.M());
        if(cm.M() >= m_low && cm.M() <= m_high && (lep1_id < 14) && (lep2_id < 14)){
            nSelected++;
            h_m->Fill(cm.M());

            double mu_p_pls = (lep_pls->E()+lep_pls->Pz())/root2;
            double mu_p_min = (lep_pls->E()-lep_pls->Pz())/root2;
            double mu_m_pls = (lep_mns->E()+lep_mns->Pz())/root2;
            double mu_m_min = (lep_mns->E()-lep_mns->Pz())/root2;
            double qt2 = cm.Px()*cm.Px()+cm.Py()*cm.Py();
            double cm_m2 = cm.M2();
            double gen_cost = 2*(mu_m_pls*mu_p_min - mu_m_min*mu_p_pls)/sqrt(cm_m2*(cm_m2 + qt2));
            //pls and mns leptons are backwards in ttree, so flip sign
            gen_cost = -gen_cost;

            //flip sign for reconstruction
            if(cm.Pz() < 0.) gen_cost = -gen_cost;
            h_cost->Fill(gen_cost);
        }
        



    }
    printf("Selected %i events \n", nSelected);
    //h_cost->Scale(1./h_cost->Integral());
    //h_m->Scale(1./h_m->Integral());


    t1->ResetBranchAddresses();
    return;
}


void get_AFB(){
    gStyle->SetOptStat(0);
    //TFile *f_mad = TFile::Open("../generator_stuff/mass_binned_200k.root");
    TFile *f_mad = TFile::Open("../generator_stuff/madgraph_m500_evts.root");
    TTree *t_mad = (TTree *)f_mad->Get("T_lhe");


    TH1F *h_mad_m = new TH1F("h_mad_m", "", 10, 150, 1000);
    TH1F *h_mad_cost = new TH1F("h_mad_cost", "", 2, -1., 1.);

    float m_low = 500.;
    float m_high = 700.;

    make_hists_from_tree(t_mad, h_mad_m, h_mad_cost, m_low, m_high);

    float nF = h_mad_cost->GetBinContent(2);
    float nB = h_mad_cost->GetBinContent(1);
    float AFB = ((nF - nB))/((nF+nB));
    float dAFB = (1.-AFB*AFB)/sqrt((nF+nB));

    printf("nF, nB %.2f %.2f \n", nF, nB);
    printf("AFB is %.3f +/- %.3f \n", AFB, dAFB);
}
    

