//construct templates from reconstructed events
//
//
//

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


#define FLAG_MUONS  0
#define FLAG_ELECTRONS  1
#define FLAG_QCD  2
#define FLAG_WJETS  3
#define FLAG_NORMAL 4

Double_t bcdef_lumi = 5.746 + 2.572 + 4.242 + 4.024 + 3.104;
// adding new Hv2 data set to get full 2016 luminosity
Double_t gh_lumi =  7.573 + 0.215 + 8.434;
Double_t tot_lumi = 35.9;

bool has_no_bjets(Int_t nJets, Double_t jet1_pt, Double_t jet2_pt, 
        Double_t jet1_cmva, Double_t jet2_cmva){
    Double_t med_btag = 0.4432;
    if(nJets ==0) return true;
    else if(nJets == 1){
        if(jet1_pt < 20.) return true;
        else return jet1_cmva < med_btag;
    }
    else{
        return (jet1_pt < 20. || jet1_cmva < med_btag) && (jet2_pt < 20. || jet2_cmva < med_btag);
    }
}

int gen_data_template(TTree *t1, TH2F* h, vector<double> *v_m, vector<double> *v_xF, vector<double> *v_cost, Double_t m_low, Double_t m_high){
    Long64_t nEntries  =  t1->GetEntries();
    //printf("size is %i \n", nEntries);
    Double_t m, xF, cost, gen_weight, jet1_cmva, jet2_cmva,
             jet1_pt, jet2_pt;
    Float_t met_pt;
    Int_t nJets;
    t1->SetBranchAddress("m", &m);
    t1->SetBranchAddress("xF", &xF);
    t1->SetBranchAddress("cost", &cost);
    t1->SetBranchAddress("jet1_CMVA", &jet1_cmva);
    t1->SetBranchAddress("jet2_CMVA", &jet2_cmva);
    t1->SetBranchAddress("met_pt", &met_pt);
    t1->SetBranchAddress("jet1_pt", &jet1_pt);
    t1->SetBranchAddress("jet2_pt", &jet2_pt);
    //t1->SetBranchAddress("nJets", &nJets);
    nJets =2;
    int nEvents = 0;
    int n=0;
    for (int i=0; i<nEntries; i++) {
        t1->GetEntry(i);
        bool no_bjets = has_no_bjets(nJets, jet1_pt, jet2_pt, jet1_cmva, jet2_cmva);
        if(m >= m_low && m <= m_high && met_pt < 50. && no_bjets){
            n++;
            h->Fill(xF, cost, 1); 
            v_m->push_back(m);
            v_xF->push_back(xF);
            v_cost->push_back(cost);
            nEvents++;
        }
    }

    h->Scale(1./h->Integral());
    t1->ResetBranchAddresses();
    return nEvents;
}


int gen_mc_template(TTree *t1, Double_t alpha, TH2F* h_sym, TH2F *h_asym, TH2F *h_count, 
        Double_t m_low, Double_t m_high, int flag = FLAG_MUONS){
    Long64_t nEntries  =  t1->GetEntries();
    //printf("size is %i \n", nEntries);

    h_sym->Sumw2();
    h_asym->Sumw2();

    Double_t m, xF, cost, gen_weight, reweight, jet1_cmva, jet2_cmva, cost_st;
    Double_t bcdef_HLT_SF, bcdef_iso_SF, bcdef_id_SF;
    Double_t gh_HLT_SF, gh_iso_SF, gh_id_SF;
    Double_t el_id_SF, el_reco_SF, pu_SF;
    Double_t jet1_pt, jet2_pt, jet1_b_weight, jet2_b_weight;
    Float_t cost_pt, met_pt;
    Int_t nJets;
    t1->SetBranchAddress("m", &m);
    t1->SetBranchAddress("xF", &xF);
    t1->SetBranchAddress("cost", &cost);
    t1->SetBranchAddress("cost_st", &cost_st);
    t1->SetBranchAddress("jet1_CMVA", &jet1_cmva);
    t1->SetBranchAddress("jet2_CMVA", &jet2_cmva);
    t1->SetBranchAddress("met_pt", &met_pt);
    t1->SetBranchAddress("jet1_pt", &jet1_pt);
    t1->SetBranchAddress("jet2_pt", &jet2_pt);
    t1->SetBranchAddress("nJets", &nJets);
    t1->SetBranchAddress("gen_weight", &gen_weight);
    t1->SetBranchAddress("jet1_b_weight", &jet1_b_weight);
    t1->SetBranchAddress("jet2_b_weight", &jet2_b_weight);
    t1->SetBranchAddress("pu_SF", &pu_SF);
    int n = 0;
    if(flag == FLAG_MUONS){
        t1->SetBranchAddress("bcdef_HLT_SF", &bcdef_HLT_SF);
        t1->SetBranchAddress("bcdef_iso_SF", &bcdef_iso_SF);
        t1->SetBranchAddress("bcdef_id_SF", &bcdef_id_SF);
        t1->SetBranchAddress("gh_HLT_SF", &gh_HLT_SF);
        t1->SetBranchAddress("gh_iso_SF", &gh_iso_SF);
        t1->SetBranchAddress("gh_id_SF", &gh_id_SF);
        TH2D *h_sym_bcdef = (TH2D *) h_sym->Clone("h_sym_bcdef");
        TH2D *h_sym_gh = (TH2D *)h_sym->Clone("h_sym_gh");
        TH2D *h_asym_bcdef = (TH2D *)h_asym->Clone("h_asym_bcdef");
        TH2D *h_asym_gh = (TH2D *)h_asym->Clone("h_asym_gh");
        for (int i=0; i<nEntries; i++) {
            t1->GetEntry(i);
            bool no_bjets = has_no_bjets(nJets, jet1_pt, jet2_pt, jet1_cmva, jet2_cmva);
            if(m >= m_low && m <= m_high && met_pt < 50.
                    && no_bjets){
                reweight = (4./3.)*cost_st*(2. + alpha)/
                    (1. + cost_st*cost_st + alpha*(1.- cost_st*cost_st));
                n++;


                Double_t bcdef_weight = gen_weight * pu_SF * bcdef_HLT_SF * bcdef_iso_SF * bcdef_id_SF;
                Double_t gh_weight = gen_weight * pu_SF * gh_HLT_SF * gh_iso_SF * gh_id_SF;
                if (nJets >= 1){
                    bcdef_weight *= jet1_b_weight;
                    gh_weight *= jet1_b_weight;
                }
                if (nJets >= 2){
                    bcdef_weight *= jet2_b_weight;
                    gh_weight *= jet2_b_weight;
                }


                h_sym_bcdef->Fill(xF, cost, bcdef_weight); 
                h_sym_bcdef->Fill(xF, -cost, bcdef_weight); 
                h_sym_gh->Fill(xF, cost, gh_weight); 
                h_sym_gh->Fill(xF, -cost, gh_weight); 

                h_asym_bcdef->Fill(xF, cost, reweight * bcdef_weight);
                h_asym_bcdef->Fill(xF, -cost, -reweight * bcdef_weight);
                h_asym_gh->Fill(xF, cost, reweight * gh_weight);
                h_asym_gh->Fill(xF, -cost, -reweight * gh_weight);

                h_count->Fill(xF, cost, 1);
                h_count->Fill(xF, -cost, 1);
            }
        }

        h_sym_bcdef->Scale(1000*bcdef_lumi);
        h_sym_gh->Scale(1000*gh_lumi);
        h_asym_bcdef->Scale(1000*bcdef_lumi);
        h_asym_gh->Scale(1000*gh_lumi);

        h_sym->Add(h_sym_bcdef, h_sym_gh);
        h_asym->Add(h_asym_bcdef, h_asym_gh);
    }
    else if (flag == FLAG_ELECTRONS) {
        t1->SetBranchAddress("el_id_SF", &el_id_SF);
        t1->SetBranchAddress("el_reco_SF", &el_reco_SF);
        for (int i=0; i<nEntries; i++) {
            t1->GetEntry(i);
            bool no_bjets = has_no_bjets(nJets, jet1_pt, jet2_pt, jet1_cmva, jet2_cmva);
            if(m >= m_low && m <= m_high && met_pt < 50.
                    && no_bjets){
                reweight = (4./3.)*cost_st*(2. + alpha)/
                    (1. + cost_st*cost_st + alpha*(1.- cost_st*cost_st));
                n++;


                Double_t evt_weight = gen_weight * el_id_SF *el_reco_SF * pu_SF;
                if (nJets >= 1){
                    evt_weight *= jet1_b_weight;
                }
                if (nJets >= 2){
                    evt_weight *= jet2_b_weight;
                }


                h_sym->Fill(xF, cost, evt_weight); 
                h_sym->Fill(xF, -cost, evt_weight); 

                h_asym->Fill(xF, cost, reweight * evt_weight);
                h_asym->Fill(xF, -cost, -reweight * evt_weight);

                h_count->Fill(xF, cost, 1);
                h_count->Fill(xF, -cost, 1);
            }
        }
        h_sym->Scale(1000*(bcdef_lumi + gh_lumi));
        h_asym->Scale(1000*(bcdef_lumi + gh_lumi));

    }

    printf("N sym is %i \n", n);
    float norm = h_sym -> Integral();
    h_sym->Scale(1./norm);
    h_asym->Scale(1./norm);
    t1->ResetBranchAddresses();
    printf("MC templates generated from %i events \n \n", n);
    return 0;
}



int gen_background_template(TTree *t1, TH2F* h, TH2F* h_count, 
        Double_t m_low, Double_t m_high, int flag = FLAG_MUONS){
    Long64_t nEntries  =  t1->GetEntries();
    h->Sumw2();

    Double_t m, xF, cost, mu1_pt, mu2_pt, jet1_cmva, jet2_cmva, gen_weight;
    Double_t bcdef_HLT_SF, bcdef_iso_SF, bcdef_id_SF;
    Double_t gh_HLT_SF, gh_iso_SF, gh_id_SF, el_id_SF, el_reco_SF, pu_SF;
    Double_t jet1_pt, jet2_pt, jet1_b_weight, jet2_b_weight;
    Float_t cost_pt, met_pt;
    Int_t nJets;

    int nEvents = 0;

    t1->SetBranchAddress("m", &m);
    t1->SetBranchAddress("xF", &xF);
    t1->SetBranchAddress("cost", &cost);
    t1->SetBranchAddress("jet1_CMVA", &jet1_cmva);
    t1->SetBranchAddress("jet2_CMVA", &jet2_cmva);
    t1->SetBranchAddress("jet1_pt", &jet1_pt);
    t1->SetBranchAddress("jet2_pt", &jet2_pt);
    t1->SetBranchAddress("met_pt", &met_pt);
    t1->SetBranchAddress("nJets", &nJets);
    t1->SetBranchAddress("gen_weight", &gen_weight);
    t1->SetBranchAddress("jet1_b_weight", &jet1_b_weight);
    t1->SetBranchAddress("jet2_b_weight", &jet2_b_weight);
    t1->SetBranchAddress("pu_SF", &pu_SF);
    if(flag == FLAG_MUONS){

        t1->SetBranchAddress("bcdef_HLT_SF", &bcdef_HLT_SF);
        t1->SetBranchAddress("bcdef_iso_SF", &bcdef_iso_SF);
        t1->SetBranchAddress("bcdef_id_SF", &bcdef_id_SF);
        t1->SetBranchAddress("gh_HLT_SF", &gh_HLT_SF);
        t1->SetBranchAddress("gh_iso_SF", &gh_iso_SF);
        t1->SetBranchAddress("gh_id_SF", &gh_id_SF);

        TH2D *h_bcdef = (TH2D *)h->Clone("h_back_bcdef");
        TH2D *h_gh = (TH2D *)h->Clone("h_back_gh");

        for (int i=0; i<nEntries; i++) {
            t1->GetEntry(i);
            bool no_bjets = has_no_bjets(nJets, jet1_pt, jet2_pt, jet1_cmva, jet2_cmva);
            if(m >= m_low && m <= m_high && met_pt < 50. 
                    && no_bjets){

                Double_t bcdef_weight = gen_weight *pu_SF * bcdef_HLT_SF * bcdef_iso_SF * bcdef_id_SF;
                Double_t gh_weight = gen_weight *pu_SF * gh_HLT_SF * gh_iso_SF * gh_id_SF;
                if (nJets >= 1){
                    bcdef_weight *= jet1_b_weight;
                    gh_weight *= jet1_b_weight;
                }
                if (nJets >= 2){
                    bcdef_weight *= jet2_b_weight;
                    gh_weight *= jet2_b_weight;
                }
                h_bcdef->Fill(xF, cost, bcdef_weight);
                h_gh->Fill(xF, cost, gh_weight);
                h_count ->Fill(xF, cost, 1);
                nEvents++;
            }
        }
        printf("N ttbar events %i \n", nEvents);
        h_bcdef->Scale(bcdef_lumi*1000);
        h_gh->Scale(gh_lumi*1000);
        h->Add(h_bcdef, h_gh);
    }
    else if (flag == FLAG_ELECTRONS) {
        t1->SetBranchAddress("el_id_SF", &el_id_SF);
        t1->SetBranchAddress("el_reco_SF", &el_reco_SF);
        for (int i=0; i<nEntries; i++) {
            t1->GetEntry(i);
            bool no_bjets = has_no_bjets(nJets, jet1_pt, jet2_pt, jet1_cmva, jet2_cmva);
            if(m >= m_low && m <= m_high && met_pt < 50.  && no_bjets){


                Double_t evt_weight = gen_weight * pu_SF * el_id_SF * el_reco_SF;
                if (nJets >= 1){
                    evt_weight *= jet1_b_weight;
                }
                if (nJets >= 2){
                    evt_weight *= jet2_b_weight;
                }
                h->Fill(xF, cost, evt_weight);
                h_count ->Fill(xF, cost, 1);
                nEvents++;
            }
        }
        //h->Scale(1000*(bcdef_lumi + gh_lumi));
    }

    printf("back norm %f \n", h->Integral());
    h->Scale(1./h->Integral());
    t1->ResetBranchAddresses();
    return 0;
}

int gen_combined_background_template(int nTrees, TTree **ts, TH2F* h,  
        Double_t m_low, Double_t m_high, int flag = FLAG_MUONS){
    h->Sumw2();
    for(int i=0; i<nTrees; i++){
        TTree *t1 = ts[i];
        Long64_t nEntries  =  t1->GetEntries();

        Double_t m, xF, cost, mu1_pt, mu2_pt, jet1_cmva, jet2_cmva, gen_weight;
        Double_t bcdef_HLT_SF, bcdef_iso_SF, bcdef_id_SF;
        Double_t gh_HLT_SF, gh_iso_SF, gh_id_SF, el_id_SF, el_reco_SF;
        Double_t jet1_pt, jet2_pt, jet1_b_weight, jet2_b_weight, pu_SF;
        Float_t cost_pt, met_pt;
        Int_t nJets;

        int nEvents = 0;

        t1->SetBranchAddress("m", &m);
        t1->SetBranchAddress("xF", &xF);
        t1->SetBranchAddress("cost", &cost);
        t1->SetBranchAddress("jet1_CMVA", &jet1_cmva);
        t1->SetBranchAddress("jet2_CMVA", &jet2_cmva);
        t1->SetBranchAddress("jet1_pt", &jet1_pt);
        t1->SetBranchAddress("jet2_pt", &jet2_pt);
        t1->SetBranchAddress("met_pt", &met_pt);
        t1->SetBranchAddress("nJets", &nJets);
        t1->SetBranchAddress("gen_weight", &gen_weight);
        t1->SetBranchAddress("jet1_b_weight", &jet1_b_weight);
        t1->SetBranchAddress("jet2_b_weight", &jet2_b_weight);
        t1->SetBranchAddress("pu_SF", &pu_SF);
        if(flag == FLAG_MUONS){

            t1->SetBranchAddress("bcdef_HLT_SF", &bcdef_HLT_SF);
            t1->SetBranchAddress("bcdef_iso_SF", &bcdef_iso_SF);
            t1->SetBranchAddress("bcdef_id_SF", &bcdef_id_SF);
            t1->SetBranchAddress("gh_HLT_SF", &gh_HLT_SF);
            t1->SetBranchAddress("gh_iso_SF", &gh_iso_SF);
            t1->SetBranchAddress("gh_id_SF", &gh_id_SF);

            TH2D *h_bcdef = (TH2D *)h->Clone("h_back_bcdef");
            TH2D *h_gh = (TH2D *)h->Clone("h_back_gh");

            for (int i=0; i<nEntries; i++) {
                t1->GetEntry(i);
                bool no_bjets = has_no_bjets(nJets, jet1_pt, jet2_pt, jet1_cmva, jet2_cmva);
                if(m >= m_low && m <= m_high && met_pt < 50. 
                        && no_bjets){

                    Double_t bcdef_weight = gen_weight *pu_SF * bcdef_HLT_SF * bcdef_iso_SF * bcdef_id_SF;
                    Double_t gh_weight = gen_weight *pu_SF * gh_HLT_SF * gh_iso_SF * gh_id_SF;
                    if (nJets >= 1){
                        bcdef_weight *= jet1_b_weight;
                        gh_weight *= jet1_b_weight;
                    }
                    if (nJets >= 2){
                        bcdef_weight *= jet2_b_weight;
                        gh_weight *= jet2_b_weight;
                    }
                    h_bcdef->Fill(xF, cost, bcdef_weight);
                    h_gh->Fill(xF, cost, gh_weight);
                }
            }
            h_bcdef->Scale(bcdef_lumi*1000);
            h_gh->Scale(gh_lumi*1000);
            h->Add(h_bcdef);
            h->Add(h_gh);
        }
        else if(flag == FLAG_ELECTRONS) {
            t1->SetBranchAddress("el_id_SF", &el_id_SF);
            t1->SetBranchAddress("el_reco_SF", &el_reco_SF);
            for (int i=0; i<nEntries; i++) {
                t1->GetEntry(i);
                bool no_bjets = has_no_bjets(nJets, jet1_pt, jet2_pt, jet1_cmva, jet2_cmva);
                if(m >= m_low && m <= m_high && met_pt < 50. 
                        && no_bjets){


                    Double_t evt_weight = gen_weight * pu_SF * el_id_SF * el_reco_SF;
                    if (nJets >= 1){
                        evt_weight *= jet1_b_weight;
                    }
                    if (nJets >= 2){
                        evt_weight *= jet2_b_weight;
                    }
                    h->Fill(xF, cost, evt_weight);
                }
            }
            //h->Scale(1000*(bcdef_lumi + gh_lumi));
        }

        t1->ResetBranchAddresses();
    }
    h->Scale(1./h->Integral());
    return 0;
}

void make_emu_m_hist(TTree *t1, TH1F *h_m, bool is_data = false, int flag1 = FLAG_MUONS){
    Long64_t size  =  t1->GetEntries();
    Double_t m, xF, cost, mu1_pt, mu2_pt, jet1_cmva, jet2_cmva, gen_weight;
    Double_t bcdef_HLT_SF, bcdef_iso_SF, bcdef_id_SF;
    Double_t gh_HLT_SF, gh_iso_SF, gh_id_SF, el_id_SF, el_reco_SF;
    Double_t jet1_pt, jet2_pt, jet1_b_weight, jet2_b_weight, pu_SF;
    Float_t met_pt;
    Int_t nJets;
    nJets = 2;
    pu_SF=1;
    TH2D *h_m_bcdef = (TH2D *)h_m->Clone("h_m_bcdef");
    TH2D *h_m_gh = (TH2D *)h_m->Clone("h_m_gh");
    TLorentzVector *el=0;
    TLorentzVector *mu=0;

    t1->SetBranchAddress("el", &el);
    t1->SetBranchAddress("mu", &mu);
    t1->SetBranchAddress("met_pt", &met_pt);
    t1->SetBranchAddress("jet2_CMVA", &jet2_cmva);
    t1->SetBranchAddress("jet1_CMVA", &jet1_cmva);
    t1->SetBranchAddress("jet1_pt", &jet1_pt);
    t1->SetBranchAddress("jet2_pt", &jet2_pt);
    if(!is_data){
        t1->SetBranchAddress("nJets", &nJets);
        t1->SetBranchAddress("gen_weight", &gen_weight);
        t1->SetBranchAddress("jet1_b_weight", &jet1_b_weight);
        t1->SetBranchAddress("jet2_b_weight", &jet2_b_weight);
        t1->SetBranchAddress("pu_SF", &pu_SF);
        t1->SetBranchAddress("bcdef_HLT_SF", &bcdef_HLT_SF);
        t1->SetBranchAddress("bcdef_iso_SF", &bcdef_iso_SF);
        t1->SetBranchAddress("bcdef_id_SF", &bcdef_id_SF);
        t1->SetBranchAddress("gh_HLT_SF", &gh_HLT_SF);
        t1->SetBranchAddress("gh_iso_SF", &gh_iso_SF);
        t1->SetBranchAddress("gh_id_SF", &gh_id_SF);
        t1->SetBranchAddress("el_id_SF", &el_id_SF);
        t1->SetBranchAddress("el_reco_SF", &el_reco_SF);     
    }

    for (int i=0; i<size; i++) {
        t1->GetEntry(i);
        bool no_bjets = has_no_bjets(nJets, jet1_pt, jet2_pt, jet1_cmva, jet2_cmva);
        TLorentzVector cm;
        cm = *el + *mu;
        m = cm.M();

        if(m >= 150. && met_pt < 50. && no_bjets){
            if(is_data){
                h_m->Fill(m);
            }
            else{
                Double_t evt_weight = gen_weight * pu_SF * el_id_SF * el_reco_SF;
                Double_t bcdef_weight = bcdef_iso_SF * bcdef_id_SF;
                Double_t gh_weight = gh_iso_SF * gh_id_SF;
                if(flag1 == FLAG_MUONS){
                    bcdef_weight *= bcdef_HLT_SF;
                    gh_weight *= gh_HLT_SF;
                }
                if (nJets >= 1){
                    evt_weight *= jet1_b_weight;
                }
                if (nJets >= 2){
                    evt_weight *= jet2_b_weight;
                }
                //printf(" %.2e %.2e %.2e \n", evt_weight, evt_weight *bcdef_weight, evt_weight *gh_weight);
                h_m_bcdef->Fill(m,evt_weight * bcdef_weight);
                h_m_gh->Fill(m, evt_weight * gh_weight);
            }


            
        }
    }
    if(!is_data){
        h_m_bcdef ->Scale(bcdef_lumi * 1000);
        h_m_gh ->Scale(gh_lumi * 1000);
        //Printf("%.1f %.1f", h_m_bcdef->Integral(), h_m_gh->Integral());
        h_m->Add(h_m_bcdef, h_m_gh);
    }
}


void make_m_cost_hist(TTree *t1, TH1F *h_m, TH1F *h_cost, bool is_data=false, int flag1 = FLAG_MUONS, int flag2 = FLAG_NORMAL){
    //read event data
    Long64_t size  =  t1->GetEntries();
    Double_t m, xF, cost, mu1_pt, mu2_pt, jet1_cmva, jet2_cmva, gen_weight;
    Double_t bcdef_HLT_SF, bcdef_iso_SF, bcdef_id_SF;
    Double_t gh_HLT_SF, gh_iso_SF, gh_id_SF, el_id_SF, el_reco_SF;
    Double_t jet1_pt, jet2_pt, jet1_b_weight, jet2_b_weight, pu_SF;
    Float_t met_pt;
    Int_t nJets;
    nJets = 2;
    pu_SF=1;
    TH2D *h_m_bcdef = (TH2D *)h_m->Clone("h_m_bcdef");
    TH2D *h_m_gh = (TH2D *)h_m->Clone("h_m_gh");
    TH2D *h_cost_bcdef = (TH2D *)h_cost->Clone("h_cost_bcdef");
    TH2D *h_cost_gh = (TH2D *)h_cost->Clone("h_cost_gh");
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
        t1->SetBranchAddress("jet1_b_weight", &jet1_b_weight);
        t1->SetBranchAddress("jet2_b_weight", &jet2_b_weight);
        t1->SetBranchAddress("pu_SF", &pu_SF);
    }
    if(flag1 == FLAG_MUONS){
        if(flag2 == FLAG_NORMAL){
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
                    if(is_data){
                        h_m->Fill(m);
                        h_cost->Fill(cost);
                    }
                    else{
                        Double_t bcdef_weight = gen_weight *pu_SF * bcdef_HLT_SF * bcdef_iso_SF * bcdef_id_SF;
                        Double_t gh_weight = gen_weight *pu_SF * gh_HLT_SF * gh_iso_SF * gh_id_SF;
                        if (nJets >= 1){
                            bcdef_weight *= jet1_b_weight;
                            gh_weight *= jet1_b_weight;
                        }
                        if (nJets >= 2){
                            bcdef_weight *= jet2_b_weight;
                            gh_weight *= jet2_b_weight;
                        }
                        //Double_t weight = gen_weight;
                        h_m_bcdef->Fill(m,bcdef_weight);
                        h_m_gh->Fill(m,gh_weight);
                        h_cost_bcdef->Fill(cost,bcdef_weight);
                        h_cost_gh->Fill(cost,gh_weight);
                    }


                }
            }
            if(!is_data){
                h_m_bcdef ->Scale(bcdef_lumi * 1000);
                h_cost_bcdef ->Scale(bcdef_lumi * 1000);
                h_m_gh ->Scale(gh_lumi * 1000);
                h_cost_gh ->Scale(gh_lumi * 1000);
                h_m->Add(h_m_bcdef, h_m_gh);
                h_cost->Add(h_cost_bcdef, h_cost_gh);
            }
        }
        else if(flag2 == FLAG_QCD){
            Double_t evt_fakerate, mu1_fakerate, mu2_fakerate;
            Bool_t double_muon_trig;
            //t1->SetBranchAddress("evt_fakerate", &evt_fakerate);
            t1->SetBranchAddress("mu1_fakerate", &mu1_fakerate);
            t1->SetBranchAddress("mu2_fakerate", &mu2_fakerate);
            //t1->SetBranchAddress("double_muon_tig", &double_muon_trig);
            for (int i=0; i<size; i++) {
                t1->GetEntry(i);
                bool no_bjets = has_no_bjets(nJets, jet1_pt, jet2_pt, jet1_cmva, jet2_cmva);

                if(m >= 150. && met_pt < 50. && no_bjets){
                    //mu1_fakerate = std::min(mu1_fakerate, 0.94);
                    //mu1_fakerate = 0.13;
                    //mu2_fakerate = 0.13;
                    /*
                       if(double_muon_trig){
                       mu1_fakerate = 0.91;
                       mu2_fakerate = 0.13;
                       }
                       */
                    //mu2_fakerate = std::min(mu2_fakerate, 0.9);
                    evt_fakerate = mu1_fakerate*mu2_fakerate/((1-mu1_fakerate)*(1-mu2_fakerate));
                    //printf("Evt %.2f %.2f %.2f \n",mu1_fakerate, mu2_fakerate, evt_fakerate);
                    h_m->Fill(m, evt_fakerate);
                    h_cost->Fill(cost, evt_fakerate);
                }
            }
            printf("Tot QCD is %.0f \n", h_m->Integral());
        }
        else if(flag2 == FLAG_WJETS){
            Double_t evt_fakerate, mu_fakerate, frac_WJet;
            Bool_t double_muon_trig;
            t1->SetBranchAddress("evt_fakerate", &evt_fakerate);
            t1->SetBranchAddress("mu_fakerate", &mu_fakerate);
            //t1->SetBranchAddress("double_muon_trig", &double_muon_trig);
            for (int i=0; i<size; i++) {
                t1->GetEntry(i);
                bool no_bjets = has_no_bjets(nJets, jet1_pt, jet2_pt, jet1_cmva, jet2_cmva);
                
                double_t WJets_MC = 7.38e-02 * 1000* tot_lumi;
                double_t background_MC = 22.3e-02 *1000*tot_lumi;
                double_t frac_WJet = WJets_MC/(WJets_MC + background_MC);
                /*
                   if(double_muon_trig){
                   frac_WJet = 4.06/(4.06 + 17.88);
                   mu_fakerate = 0.13;
                   }
                   else{
                   frac_WJet = 3.33/(3.33 + 2.53);
                   mu_fakerate = 0.13;
                   }
                   */

                if(m >= 150. && met_pt < 50. && no_bjets){
                    //mu_fakerate = std::min(mu_fakerate, 0.9);
                    evt_fakerate =  frac_WJet * mu_fakerate/(1-mu_fakerate);
                    //printf("Evt %.2f %.2f \n", mu_fakerate, evt_fakerate);
                    h_m->Fill(m, evt_fakerate);
                    h_cost->Fill(cost, evt_fakerate);
                }
            }
            printf("Tot WJet is %.0f \n", h_m->Integral());
        }
    }
    else if(flag1==FLAG_ELECTRONS) {
        if(flag2 == FLAG_NORMAL){
            if(!is_data) t1->SetBranchAddress("el_id_SF", &el_id_SF);
            if(!is_data) t1->SetBranchAddress("el_reco_SF", &el_reco_SF);
            for (int i=0; i<size; i++) {
                t1->GetEntry(i);
                bool no_bjets = has_no_bjets(nJets, jet1_pt, jet2_pt, jet1_cmva, jet2_cmva);
                if(m >= 150 && met_pt < 50.  && no_bjets){
                    if(is_data){
                        h_m->Fill(m);
                        h_cost->Fill(cost);
                    }
                    else{

                        Double_t evt_weight = gen_weight *pu_SF * el_id_SF * el_reco_SF;
                        if (nJets >= 1){
                            evt_weight *= jet1_b_weight;
                        }
                        if (nJets >= 2){
                            evt_weight *= jet2_b_weight;
                        }
                        h_m->Fill(m, evt_weight);
                        h_cost->Fill(cost, evt_weight);
                    }
                }
            }
            if(!is_data){
                Double_t el_lumi = 1000*tot_lumi;
                h_m->Scale(el_lumi);
                h_cost->Scale(el_lumi);
            }
        }
        if(flag2 == FLAG_WJETS){
            Double_t evt_fakerate, el_fakerate;
            t1->SetBranchAddress("evt_fakerate", &evt_fakerate);
            t1->SetBranchAddress("el_fakerate", &el_fakerate);
            double_t WJets_MC = 1.19e-01 * 1000* tot_lumi;
            double_t background_MC = 1.9e-01 *1000*tot_lumi;
            double_t frac_WJet = WJets_MC/(WJets_MC + background_MC);
            printf("N events, WJets_MC, background_MC: %i, %.0f, %.0f \n", size, WJets_MC, background_MC);
            for (int i=0; i<size; i++) {
                t1->GetEntry(i);
                bool no_bjets = has_no_bjets(nJets, jet1_pt, jet2_pt, jet1_cmva, jet2_cmva);


                if(m >= 150. && met_pt < 50. && no_bjets){
                    //el_fakerate = std::min(el_fakerate, 0.9);
                    evt_fakerate =  frac_WJet * el_fakerate/(1-el_fakerate);
                    //printf("Evt %.2f %.2f \n", el_fakerate, evt_fakerate);
                    h_m->Fill(m, evt_fakerate);
                    h_cost->Fill(cost, evt_fakerate);
                }
            }
            printf("Tot WJet is %.0f from %i events \n", h_m->Integral(), size);
        }
        else if(flag2 == FLAG_QCD){
            Double_t evt_fakerate, el1_fakerate, el2_fakerate;
            //t1->SetBranchAddress("evt_fakerate", &evt_fakerate);
            t1->SetBranchAddress("el1_fakerate", &el1_fakerate);
            t1->SetBranchAddress("el2_fakerate", &el2_fakerate);
            for (int i=0; i<size; i++) {
                t1->GetEntry(i);
                bool no_bjets = has_no_bjets(nJets, jet1_pt, jet2_pt, jet1_cmva, jet2_cmva);
                double_t WJets_MC = 1.19e-01 * 1000* tot_lumi;
                double_t background_MC = 1.9e-01 *1000*tot_lumi;
                double_t frac_Wjets = WJets_MC/(WJets_MC + background_MC);
                double_t frac_QCD = 1;

                if(m >= 150. && met_pt < 50. && no_bjets){
                    //el1_fakerate = std::min(el1_fakerate, 0.94);
                    //el2_fakerate = std::min(el2_fakerate, 0.94);
                    el1_fakerate = 0.44;
                    el2_fakerate = 0.44;
                    evt_fakerate = frac_QCD* el1_fakerate*el2_fakerate/((1-el1_fakerate)*(1-el2_fakerate));
                    //printf("Evt %.2f %.2f %.2f \n",el1_fakerate, el2_fakerate, evt_fakerate);
                    h_m->Fill(m, evt_fakerate);
                    h_cost->Fill(cost, evt_fakerate);
                }
            }
            printf("Tot QCD is %.0f from %i events\n", h_m->Integral(), size);
        }

    }


    t1->ResetBranchAddresses();
}

typedef struct {
    TH2D *h;
} FakeRate;
//static type means functions scope is only this file, to avoid conflicts
static void setup_new_el_fakerate(FakeRate *FR){
    TFile *f0 = TFile::Open("../analyze/FakeRate/root_files/SingleElectron_data_fake_rate_v2_corrected_nov29.root");
    TH2D *h1 = (TH2D *) gDirectory->Get("h_rate_new")->Clone();
    h1->SetDirectory(0);
    FR->h = h1;
    f0->Close();
}
static void setup_new_mu_fakerate(FakeRate *FR){
    TFile *f0 = TFile::Open("../analyze/FakeRate/root_files/SingleMuon_data_fake_rate_v2_corrected_nov29.root");
    f0->ls();
    TDirectory *subdir = gDirectory;
    TH2D *h1 = (TH2D *) subdir->Get("h_rate_new")->Clone();
    h1->SetDirectory(0);
    h1->Print();
    FR->h = h1;
    f0->Close();
}


static Double_t get_new_fakerate_prob(Double_t pt, Double_t eta, TH2D *h){
    //pt=35;
    if (pt >= 80) pt =80 ;


    TAxis* x_ax =  h->GetXaxis();
    TAxis *y_ax =  h->GetYaxis();
    
    int xbin = x_ax->FindBin(std::abs(eta));
    int ybin = y_ax->FindBin(pt);

    Double_t prob = h->GetBinContent(xbin, ybin);
    //printf("prob: %.2f \n", prob);
    if(prob < 0.001 || prob >= 0.99){
        printf("Warning: %.2f Rate for pt %.0f, eta %1.1f! \n", prob, pt, eta);
        ybin -=1;
        prob = h->GetBinContent(xbin, ybin);
        if(prob < 0.001) printf("Tried 1 lower pt bin and still 0 fakerate \n");
    }

    prob = min(prob, 0.98);
    prob = max(prob, 0.07);
    //printf("Efficiency is %f \n", eff);
    return prob;
}

void Fakerate_est_mu(TTree *t_WJets, TTree *t_QCD, TTree *t_MC, TH1F *h_m, TH1F *h_cost){
    FakeRate FR;
    //TH2D *FR;
    setup_new_mu_fakerate(&FR);
    FR.h->Print();
    for (int l=0; l<=2; l++){
        printf("l=%i\n", l);
        TTree *t;
        if (l==0) t = t_WJets;
        if (l==1) t = t_QCD;
        if (l==2) t = t_MC;
        Double_t m, xF, cost, jet1_cmva, jet2_cmva, gen_weight;
        Double_t bcdef_HLT_SF, bcdef_iso_SF, bcdef_id_SF;
        Double_t gh_HLT_SF, gh_iso_SF, gh_id_SF;
        Double_t jet1_pt, jet2_pt, jet1_b_weight, jet2_b_weight, pu_SF;
        Double_t evt_fakerate, mu1_fakerate, mu2_fakerate, mu1_eta, mu1_pt, mu2_eta, mu2_pt;
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
        if(l==0){
            t->SetBranchAddress("iso_muon", &iso_mu);
        }
        if(l==2){
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
                Double_t mc_weight = (bcdef_weight *bcdef_lumi + gh_weight * gh_lumi)/(bcdef_lumi + gh_lumi);
                if(iso_mu ==0) mu1_fakerate = get_new_fakerate_prob(mu1_pt, mu1_eta, FR.h);
                if(iso_mu ==1) mu1_fakerate = get_new_fakerate_prob(mu2_pt, mu2_eta, FR.h);
                evt_fakerate = -(mu1_fakerate * mc_weight)/(1-mu1_fakerate);
            }



            h_m->Fill(m, evt_fakerate);
            h_cost->Fill(cost, evt_fakerate);
        }

    }
    printf("Total fakerate est is %.0f \n", h_m->Integral());
}

void Fakerate_est_el(TTree *t_WJets, TTree *t_QCD, TTree *t_MC, TH1F *h_m, TH1F *h_cost){
    FakeRate FR;
    //TH2D *FR;
    setup_new_el_fakerate(&FR);
    FR.h->Print();
    for (int l=0; l<=2; l++){
        TTree *t;
        if (l==0) t = t_WJets;
        if (l==1) t = t_QCD;
        if (l==2) t = t_MC;
        Double_t m, xF, cost, jet1_cmva, jet2_cmva, gen_weight;
        Double_t jet1_pt, jet2_pt, jet1_b_weight, jet2_b_weight, pu_SF;
        Double_t el_id_SF, el_reco_SF;
        Double_t evt_fakerate, el1_fakerate, el2_fakerate, el1_eta, el1_pt, el2_eta, el2_pt;
        Int_t iso_el;
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
        //t1->SetBranchAddress("el_fakerate", &el1_fakerate);
        t->SetBranchAddress("el1_pt", &el1_pt);
        t->SetBranchAddress("el2_pt", &el2_pt);
        t->SetBranchAddress("el1_eta", &el1_eta);
        t->SetBranchAddress("el2_eta", &el2_eta);
        t->SetBranchAddress("nJets", &nJets);
        if(l==0){
            t->SetBranchAddress("iso_el", &iso_el);
        }
        if(l==2){
            t->SetBranchAddress("el_id_SF", &el_id_SF);
            t->SetBranchAddress("el_reco_SF", &el_reco_SF);
            t->SetBranchAddress("gen_weight", &gen_weight);
            t->SetBranchAddress("jet1_b_weight", &jet1_b_weight);
            t->SetBranchAddress("jet2_b_weight", &jet2_b_weight);
            t->SetBranchAddress("pu_SF", &pu_SF);
        }

        Long64_t size  =  t->GetEntries();
        for (int i=0; i<size; i++) {
            t->GetEntry(i);
            bool no_bjets = has_no_bjets(nJets, jet1_pt, jet2_pt, jet1_cmva, jet2_cmva);
            if(l==0){
                if(iso_el ==0) el1_fakerate = get_new_fakerate_prob(el1_pt, el1_eta, FR.h);
                if(iso_el ==1) el1_fakerate = get_new_fakerate_prob(el2_pt, el2_eta, FR.h);
                evt_fakerate = el1_fakerate/(1-el1_fakerate);
            }
            if(l==1){
                el1_fakerate = get_new_fakerate_prob(el1_pt, el1_eta, FR.h);
                el2_fakerate = get_new_fakerate_prob(el2_pt, el2_eta, FR.h);
                evt_fakerate = -(el1_fakerate/(1-el1_fakerate)) * (el2_fakerate/(1-el2_fakerate));
            }
            if(l==2){
                Double_t mc_weight = gen_weight * el_id_SF * el_reco_SF * jet1_b_weight * jet2_b_weight;
                if(iso_el ==0) el1_fakerate = get_new_fakerate_prob(el1_pt, el1_eta, FR.h);
                if(iso_el ==1) el1_fakerate = get_new_fakerate_prob(el2_pt, el2_eta, FR.h);
                evt_fakerate = -(el1_fakerate * mc_weight)/(1-el1_fakerate);
            }



            h_m->Fill(m, evt_fakerate);
            h_cost->Fill(cost, evt_fakerate);
        }

    }
    printf("Total fakerate est is %.0f \n", h_m->Integral());
}

