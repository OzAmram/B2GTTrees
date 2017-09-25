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

Double_t bcdef_lumi = 5.746 + 2.572 + 4.242 + 4.024 + 3.104;
// adding new Hv2 data set to get full 2016 luminosity
Double_t gh_lumi =  7.573 + 0.215 + 8.434;
Double_t tot_lumi = 35.9;
Double_t alpha = 0.0180;

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

    printf("Fit will be run on %i events in mass range [%1.0f, %1.0f] \n", 
            nEvents, m_low, m_high);
    h->Scale(1./h->Integral());
    t1->ResetBranchAddresses();
    return nEvents;
}


int gen_mc_template(TTree *t1, TH2F* h_sym, TH2F *h_asym, TH2F *h_count, 
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
        Double_t gh_HLT_SF, gh_iso_SF, gh_id_SF, el_id_S, el_reco_SF;
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


void make_m_cost_hist(TTree *t1, TH1F *h_m, TH1F *h_cost, bool is_data=false, int flag = FLAG_MUONS){
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
    if(flag == FLAG_MUONS){
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
    else if(flag == FLAG_QCD){
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
                //mu1_fakerate = std::min(mu1_fakerate, 0.91);
                mu1_fakerate = 0.91;
                mu2_fakerate = 0.13;
                /*
                if(double_muon_trig){
                    mu1_fakerate = 0.91;
                    mu2_fakerate = 0.13;
                }
                */
                //mu2_fakerate = std::min(mu2_fakerate, 0.9);
                evt_fakerate = mu1_fakerate*mu2_fakerate/((1-mu1_fakerate)*(1-mu2_fakerate));
                printf("Evt %.2f %.2f %.2f \n",mu1_fakerate, mu2_fakerate, evt_fakerate);
                h_m->Fill(m, evt_fakerate);
                h_cost->Fill(cost, evt_fakerate);
            }
        }
        printf("Tot QCD is %.0f \n", h_m->Integral());
    }
    else if(flag == FLAG_WJETS){
        Double_t evt_fakerate, mu_fakerate, frac_WJet;
        Bool_t double_muon_trig;
        t1->SetBranchAddress("evt_fakerate", &evt_fakerate);
        t1->SetBranchAddress("mu_fakerate", &mu_fakerate);
        //t1->SetBranchAddress("double_muon_trig", &double_muon_trig);
        for (int i=0; i<size; i++) {
            t1->GetEntry(i);
            bool no_bjets = has_no_bjets(nJets, jet1_pt, jet2_pt, jet1_cmva, jet2_cmva);
            frac_WJet = 7.38/(22.3+7.38);
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
    else if(flag==FLAG_ELECTRONS) {
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

    t1->ResetBranchAddresses();
}

