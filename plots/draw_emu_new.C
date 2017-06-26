
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



Double_t count_tree(TTree *t1,  Double_t m_low, Double_t m_high, bool is_data=false){
    //count events in the tree
    Long64_t size  =  t1->GetEntries();
    Int_t nJets;
    Double_t m, xF, cost, mu1_pt, mu2_pt, jet1_cmva, jet2_cmva, gen_weight;
    Double_t jet1_pt, jet2_pt;
    Double_t bcdef_HLT_SF, bcdef_iso_SF, bcdef_id_SF;
    Double_t gh_HLT_SF, gh_iso_SF, gh_id_SF;
    Double_t el_id_SF, btag_weight;
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
        t1->SetBranchAddress("btag_weight", &btag_weight);
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
            }
            else{
                Double_t bcdef_weight = gen_weight * bcdef_HLT_SF * 
                    bcdef_iso_SF * bcdef_id_SF *el_id_SF * btag_weight;
                Double_t gh_weight = gen_weight * gh_HLT_SF * 
                    gh_iso_SF * gh_id_SF* el_id_SF *btag_weight;
                //printf("%.2e %.2e \n", bcdef_weight, gh_weight);
                bcdef_count += bcdef_weight;
                gh_count += gh_weight;
                //count += gen_weight;
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
    TFile *f_data = TFile::Open("../analyze/output_files/EMu_data_jun21.root");
    TTree *t_data = (TTree *)f_data->Get("T_data");

                                
    TFile *f_ttbar = TFile::Open("../analyze/output_files/EMu_ttbar_jun22.root");
    TTree *t_ttbar = (TTree *)f_ttbar->Get("T_data");

    TFile *f_DYToLL = TFile::Open("../analyze/output_files/EMu_DY_jun22.root");
    TTree *t_DYToLL = (TTree *)f_DYToLL->Get("T_data");

    TFile *f_diboson = TFile::Open("../analyze/output_files/EMu_diboson_jun22.root");
    TTree *t_diboson = (TTree *)f_diboson->Get("T_data");

    TFile *f_wt = TFile::Open("../analyze/output_files/EMu_WT_jun21.root");
    TTree *t_wt = (TTree *)f_wt->Get("T_data");


    Double_t m_low = 150;
    Double_t m_high = 200;
    Double_t data_count = count_tree(t_data, m_low, m_high, true);
    Double_t ttbar_count = count_tree(t_ttbar, m_low, m_high);
    Double_t diboson_count = count_tree(t_diboson, m_low, m_high) ;
    Double_t wt_count = count_tree(t_wt, m_low, m_high);
    Double_t DYToLL_count = count_tree(t_DYToLL, m_low, m_high);

    printf("Data count %.0f \n", data_count);
    printf("TTbar count %.0f \n", ttbar_count);
    printf("Diboson count %.0f \n", diboson_count);
    printf("wt count %.0f \n", wt_count);
    printf("DYToLL count %.0f \n", DYToLL_count);
    Double_t ratio = (data_count - diboson_count - DYToLL_count - wt_count )/ttbar_count;
    Double_t unc = ratio * sqrt((1/data_count) + (1/ttbar_count) + (1/wt_count));
    printf("Ratio (data - diboson - DYToLL)/TTbar is %1.2f +/- %1.2f \n", ratio, unc);
    return;
}









