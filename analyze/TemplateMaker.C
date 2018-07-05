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
#include "TTree.h"
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
#include "BTagUtils.C"
#include "HistMaker.C"

void cleanup_template(TH2F *h){
    int n_xf_bins = 5;
    int n_cost_bins = 10;
    for(int i=1; i<= n_xf_bins; i++){
        for(int j=1; j<=n_cost_bins; j++){
            float val = h->GetBinContent(i,j);
            if(val<0.) h->SetBinContent(i,j,0.);
        }
    }
}

int gen_data_template(TTree *t1, TH2F* h, vector<double> *v_xF, vector<double> *v_cost, Double_t var_low, Double_t var_high, 
        int flag1 = FLAG_MUONS, int flag2 = FLAG_M_BINS, bool turn_on_RC = false){
    Long64_t nEntries  =  t1->GetEntries();
    //printf("size is %i \n", nEntries);
    Double_t m, xF, cost, gen_weight, jet1_cmva, jet2_cmva,
             jet1_pt, jet2_pt;
    Double_t lep1_pt_corr, lep2_pt_corr;
    Double_t pt;
    h->Sumw2();

    Float_t met_pt;
    Int_t nJets;
    TLorentzVector *lep_p = 0;
    TLorentzVector *lep_m = 0;
    TLorentzVector cm;
    t1->SetBranchAddress("m", &m);
    t1->SetBranchAddress("xF", &xF);
    t1->SetBranchAddress("cost", &cost);
    t1->SetBranchAddress("jet1_CMVA", &jet1_cmva);
    t1->SetBranchAddress("jet2_CMVA", &jet2_cmva);
    t1->SetBranchAddress("met_pt", &met_pt);
    t1->SetBranchAddress("jet1_pt", &jet1_pt);
    t1->SetBranchAddress("jet2_pt", &jet2_pt);
    t1->SetBranchAddress("nJets", &nJets);
    if(flag1 == FLAG_MUONS){
        t1->SetBranchAddress("mu_p", &lep_p);
        t1->SetBranchAddress("mu_m", &lep_m);
        t1->SetBranchAddress("mu1_pt_corr", &lep1_pt_corr);
        t1->SetBranchAddress("mu2_pt_corr", &lep2_pt_corr);
    }
    else{
        t1->SetBranchAddress("el_p", &lep_p);
        t1->SetBranchAddress("el_m", &lep_m);
    }
    //t1->SetBranchAddress("nJets", &nJets);
    nJets =2;
    int nEvents = 0;
    int n=0;
    for (int i=0; i<nEntries; i++) {
        t1->GetEntry(i);
        cost = get_cost_v2(*lep_p, *lep_m);
        bool no_bjets = has_no_bjets(nJets, jet1_pt, jet2_pt, jet1_cmva, jet2_cmva);
        if(turn_on_RC && flag1 == FLAG_MUONS){
            TLorentzVector mu_p_new, mu_m_new;
            if((lep_p->Pt() > lep_m->Pt() && lep1_pt_corr > lep2_pt_corr) ||
                    (lep_p->Pt() < lep_m->Pt() && lep1_pt_corr < lep2_pt_corr)){
                mu_p_new.SetPtEtaPhiE(lep1_pt_corr, lep_p->Eta(), lep_p->Phi(), lep_p->E());
                mu_m_new.SetPtEtaPhiE(lep2_pt_corr, lep_m->Eta(), lep_m->Phi(), lep_m->E());
            }
            else{
                mu_p_new.SetPtEtaPhiE(lep2_pt_corr, lep_p->Eta(), lep_p->Phi(), lep_p->E());
                mu_m_new.SetPtEtaPhiE(lep1_pt_corr, lep_m->Eta(), lep_m->Phi(), lep_m->E());
            }
            double new_cost = get_cost_v2(mu_p_new, mu_m_new);
            cm = mu_p_new + mu_m_new;
            //if(cm.M() < 150.) continue;
            cost = new_cost;
            xF = abs(2.*cm.Pz()/13000.); 
            m = cm.M();
            pt = cm.Pt();
               
        }
        if(flag2 == FLAG_M_BINS){
            if(m >= var_low && m <= var_high && met_pt < 50. && no_bjets){
                //if(m> 2000.) printf("m %.0f \n", m);
                n++;
                h->Fill(xF, cost, 1); 
                v_xF->push_back(xF);
                v_cost->push_back(cost);
                nEvents++;
            }
        }
        else{
            cm = *lep_p + *lep_m;
            pt = cm.Pt();
            if(m>= 150. &&  pt >= var_low && pt <= var_high && met_pt < 50. && no_bjets){
                n++;
                h->Fill(xF, cost, 1); 
                v_xF->push_back(xF);
                v_cost->push_back(cost);
                nEvents++;
            }
        }

    }

    h->Scale(1./h->Integral());
    t1->ResetBranchAddresses();
    return nEvents;
}



int gen_mc_template(TTree *t1, Double_t alpha, TH2F* h_sym, TH2F *h_asym, TH2F *h_count, 
        Double_t var_low, Double_t var_high, int flag1 = FLAG_MUONS, int flag2 = FLAG_M_BINS,
        bool turn_on_RC = false){
    Long64_t nEntries  =  t1->GetEntries();
    //printf("size is %i \n", nEntries);

    h_sym->Sumw2();
    h_asym->Sumw2();

    Double_t m, xF, cost, gen_weight, reweight, jet1_cmva, jet2_cmva, cost_st;
    Double_t bcdef_HLT_SF, bcdef_iso_SF, bcdef_id_SF;
    Double_t gh_HLT_SF, gh_iso_SF, gh_id_SF;
    Double_t gh_trk_SF, bcdef_trk_SF;
    Double_t el_id_SF, el_reco_SF, pu_SF, el_HLT_SF;
    Double_t jet1_pt, jet2_pt, jet1_eta, jet2_eta, jet1_b_weight, jet2_b_weight;
    Double_t mu_R_up, mu_R_down, mu_F_up, mu_F_down, 
             mu_RF_up, mu_RF_down, pdf_up, pdf_down;
    Double_t mu1_pt_corr, mu2_pt_corr;
    Float_t cost_pt, met_pt;
    TLorentzVector *lep_p=0;
    TLorentzVector *lep_m=0;
    Double_t pt;
    Int_t nJets, jet1_flavour, jet2_flavour;
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
    t1->SetBranchAddress("pu_SF", &pu_SF);
    t1->SetBranchAddress("jet1_eta", &jet1_eta);
    t1->SetBranchAddress("jet2_eta", &jet2_eta);
    t1->SetBranchAddress("jet1_flavour", &jet1_flavour);
    t1->SetBranchAddress("jet2_flavour", &jet2_flavour);
    //t1->SetBranchAddress("mu_R_up", &mu_R_up);
    //t1->SetBranchAddress("mu_R_down", &mu_R_down);
    //t1->SetBranchAddress("mu_F_up", &mu_F_up);
    //t1->SetBranchAddress("mu_F_down", &mu_F_down);
    //t1->SetBranchAddress("mu_RF_up", &mu_RF_up);
    //t1->SetBranchAddress("mu_RF_down", &mu_RF_down);
    //t1->SetBranchAddress("pdf_up", &pdf_up);
    //t1->SetBranchAddress("pdf_down", &pdf_down);
    int n = 0;
    Double_t one = 1.0;
    Double_t *systematic = &one;
    BTag_readers b_reader;
    BTag_effs btag_effs;

    bool do_bweight = true;
    jet1_b_weight = 1.;
    jet2_b_weight = 1.;

    setup_btag_SFs(&b_reader, &btag_effs);

    if(flag1 == FLAG_MUONS){
        t1->SetBranchAddress("mu_p", &lep_p);
        t1->SetBranchAddress("mu_m", &lep_m);
        t1->SetBranchAddress("mu1_pt_corr", &mu1_pt_corr);
        t1->SetBranchAddress("mu2_pt_corr", &mu2_pt_corr);
        t1->SetBranchAddress("bcdef_HLT_SF", &bcdef_HLT_SF);
        t1->SetBranchAddress("bcdef_iso_SF", &bcdef_iso_SF);
        t1->SetBranchAddress("bcdef_id_SF", &bcdef_id_SF);
        t1->SetBranchAddress("bcdef_trk_SF", &bcdef_trk_SF);
        t1->SetBranchAddress("gh_HLT_SF", &gh_HLT_SF);
        t1->SetBranchAddress("gh_iso_SF", &gh_iso_SF);
        t1->SetBranchAddress("gh_id_SF", &gh_id_SF);
        t1->SetBranchAddress("gh_trk_SF", &gh_trk_SF);
        for (int i=0; i<nEntries; i++) {
            t1->GetEntry(i);
            bool no_bjets = has_no_bjets(nJets, jet1_pt, jet2_pt, jet1_cmva, jet2_cmva);
            cost = get_cost_v2(*lep_p, *lep_m);
            if(cost_st>0.) cost_st = fabs(cost);
            else cost_st = -fabs(cost);
            if(flag2 == FLAG_PT_BINS){
                TLorentzVector cm = *lep_p + *lep_m;
                pt = cm.Pt();
            }
            if(turn_on_RC){
                TLorentzVector mu_p_new, mu_m_new, cm;
                if((lep_p->Pt() > lep_m->Pt() && mu1_pt_corr > mu2_pt_corr) ||
                        (lep_p->Pt() < lep_m->Pt() && mu1_pt_corr < mu2_pt_corr)){
                    mu_p_new.SetPtEtaPhiE(mu1_pt_corr, lep_p->Eta(), lep_p->Phi(), lep_p->E());
                    mu_m_new.SetPtEtaPhiE(mu2_pt_corr, lep_m->Eta(), lep_m->Phi(), lep_m->E());
                }
                else{
                    mu_p_new.SetPtEtaPhiE(mu2_pt_corr, lep_p->Eta(), lep_p->Phi(), lep_p->E());
                    mu_m_new.SetPtEtaPhiE(mu1_pt_corr, lep_m->Eta(), lep_m->Phi(), lep_m->E());
                }
                double new_cost = get_cost_v2(mu_p_new, mu_m_new);
                cm = mu_p_new + mu_m_new;
                //if(cm.M() < 150.) continue;
                cost = new_cost;
                m = cm.M();
                pt = cm.Pt();
                xF = abs(2.*cm.Pz()/13000.); 
                if(cost_st>0.) cost_st = fabs(cost);
                else cost_st = -fabs(cost);
                   
            }
            bool pass = ((flag2 == FLAG_M_BINS && m >= var_low && m <= var_high) ||
                    (flag2 == FLAG_PT_BINS && m >= 150. && pt >= var_low && pt <= var_high))
                && met_pt < 50.  && no_bjets;
            if(pass){
                reweight = (4./3.)*cost_st*(2. + alpha)/
                    (1. + cost_st*cost_st + alpha*(1.- cost_st*cost_st));
                n++;
                gen_weight *= *systematic;

                if(do_bweight){
                    jet1_b_weight = get_btag_weight(jet1_pt, jet1_eta,(Float_t) jet1_flavour , btag_effs, b_reader, 0);
                    jet2_b_weight = get_btag_weight(jet2_pt, jet2_eta,(Float_t) jet2_flavour , btag_effs, b_reader, 0);
                }

                Double_t bcdef_weight = gen_weight * pu_SF * bcdef_HLT_SF * bcdef_iso_SF * bcdef_id_SF * bcdef_trk_SF;
                Double_t gh_weight = gen_weight * pu_SF * gh_HLT_SF * gh_iso_SF * gh_id_SF * gh_trk_SF;
                
                if (nJets >= 1){
                    bcdef_weight *= jet1_b_weight;
                    gh_weight *= jet1_b_weight;
                }
                if (nJets >= 2){
                    bcdef_weight *= jet2_b_weight;
                    gh_weight *= jet2_b_weight;
                }
                Double_t final_weight = (bcdef_weight*bcdef_lumi + gh_weight*gh_lumi);


                h_sym->Fill(xF, cost, final_weight); 
                h_sym->Fill(xF, -cost, final_weight); 

                h_asym->Fill(xF, cost, reweight * final_weight);
                h_asym->Fill(xF, -cost, -reweight * final_weight);

                h_count->Fill(xF, cost, 1);
                h_count->Fill(xF, -cost, 1);
            }
        }


    }
    else if (flag1 == FLAG_ELECTRONS) {
        t1->SetBranchAddress("el_p", &lep_p);
        t1->SetBranchAddress("el_m", &lep_m);
        t1->SetBranchAddress("el_id_SF", &el_id_SF);
        t1->SetBranchAddress("el_HLT_SF", &el_HLT_SF);
        t1->SetBranchAddress("el_reco_SF", &el_reco_SF);
        for (int i=0; i<nEntries; i++) {
            t1->GetEntry(i);
            bool no_bjets = has_no_bjets(nJets, jet1_pt, jet2_pt, jet1_cmva, jet2_cmva);
            if(flag2 == FLAG_PT_BINS){
                TLorentzVector cm = *lep_p + *lep_m;
                pt = cm.Pt();
            }
            bool pass = ((flag2 == FLAG_M_BINS && m >= var_low && m <= var_high) ||
                    (flag2 == FLAG_PT_BINS && m >= 150. && pt >= var_low && pt <= var_high))
                && met_pt < 50.  && no_bjets;
            if(pass){
                cost = get_cost_v2(*lep_p, *lep_m);
                if(cost_st>0.) cost_st = fabs(cost);
                else cost_st = -fabs(cost);
                reweight = (4./3.)*cost_st*(2. + alpha)/
                    (1. + cost_st*cost_st + alpha*(1.- cost_st*cost_st));
                n++;
                gen_weight *= *systematic;

                if(do_bweight){
                    jet1_b_weight = get_btag_weight(jet1_pt, jet1_eta,(Float_t) jet1_flavour , btag_effs, b_reader, 0);
                    jet2_b_weight = get_btag_weight(jet2_pt, jet2_eta,(Float_t) jet2_flavour , btag_effs, b_reader, 0);
                }


                Double_t evt_weight = gen_weight * el_id_SF *el_reco_SF * pu_SF * el_HLT_SF;
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




void gen_fakes_template(TTree *t_WJets, TTree *t_QCD, TTree *t_WJets_contam, 
        TTree* t_QCD_contam, TH2F *h, Double_t var_low, Double_t var_high, 
        int flag1 = FLAG_MUONS, int flag2 = FLAG_M_BINS){
    h->Sumw2();
    TLorentzVector *lep_p=0;
    TLorentzVector *lep_m=0;
    Double_t pt;
    if(flag1 == FLAG_MUONS){
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
            t->SetBranchAddress("mu_p", &lep_p);
            t->SetBranchAddress("mu_m", &lep_m);
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
                bool pass = ((flag2 == FLAG_M_BINS && m >= var_low && m <= var_high) ||
                        (flag2 == FLAG_PT_BINS && m >= 150. && pt >= var_low && pt <= var_high))
                    && met_pt < 50.  && no_bjets;

                if(pass){
                    cost = get_cost_v2(*lep_p, *lep_m);
                    TLorentzVector cm = *lep_p + *lep_m;
                    pt = cm.Pt();
                    xF = abs(2.*cm.Pz()/13000.); 
                    //if(l==3) printf("Evt fr %.2e \n", evt_fakerate);
                    //if(l==3) printf("cost, fr %.2f %.2e \n", cost, evt_fakerate);
                    h->Fill(xF, cost, evt_fakerate);
                }


            }
            printf("After iter %i current fakerate est is %.0f \n", l, h->Integral());
        }
    }
    else{

        FakeRate FR;
        //TH2D *FR;
        setup_new_el_fakerate(&FR);
        FR.h->Print();
        for (int l=0; l<=3; l++){
            TTree *t;
            if (l==0) t = t_WJets;
            if (l==1) t = t_QCD;
            if (l==2) t = t_WJets_contam;
            if (l==3) t = t_QCD_contam;
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
            t->SetBranchAddress("el_p", &lep_p);
            t->SetBranchAddress("el_m", &lep_m);
            if(l==0 || l==2 ){
                t->SetBranchAddress("iso_el", &iso_el);
            }
            if(l==2 || l==3){
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

                    Double_t mc_weight = gen_weight * el_id_SF * el_reco_SF  * 1000. * tot_lumi;
                    if(iso_el ==0) el1_fakerate = get_new_fakerate_prob(el1_pt, el1_eta, FR.h);
                    if(iso_el ==1) el1_fakerate = get_new_fakerate_prob(el2_pt, el2_eta, FR.h);
                    evt_fakerate = -(el1_fakerate * mc_weight)/(1-el1_fakerate);
                }
                if(l==3){
                    Double_t mc_weight = gen_weight * el_id_SF * el_reco_SF * 1000. * tot_lumi;

                    el1_fakerate = get_new_fakerate_prob(el1_pt, el1_eta, FR.h);
                    el2_fakerate = get_new_fakerate_prob(el2_pt, el2_eta, FR.h);
                    evt_fakerate = mc_weight * (el1_fakerate/(1-el1_fakerate)) * (el2_fakerate/(1-el2_fakerate));
                }



                if(flag2 == FLAG_PT_BINS){
                    TLorentzVector cm = *lep_p + *lep_m;
                    pt = cm.Pt();
                }
                bool pass = ((flag2 == FLAG_M_BINS && m >= var_low && m <= var_high) ||
                        (flag2 == FLAG_PT_BINS && m >= 150. && pt >= var_low && pt <= var_high))
                    && met_pt < 50.  && no_bjets;
                if(pass){
                    cost = get_cost_v2(*lep_p, *lep_m);
                    TLorentzVector cm = *lep_p + *lep_m;
                    pt = cm.Pt();
                    xF = abs(2.*cm.Pz()/13000.); 
                    //if(l==3) printf("Evt fr %.2e \n", evt_fakerate);
                    //if(l==3) printf("cost, fr %.2f %.2e \n", cost, evt_fakerate);
                    h->Fill(xF, cost, evt_fakerate);
                }
            }

            printf("After iter %i current fakerate est is %.0f \n", l, h->Integral());
        }
    }
    printf("Performing fakes cleanup (removing neg. bins) \n");
    cleanup_template(h);
    printf("Total Fakerate weight Weight is %.2f \n", h->Integral());
    if(h->Integral() < 0.){
        h->Scale(0.);
        printf("zeroing Fakes template \n");
    }
    return;
}

int gen_combined_background_template(int nTrees, TTree **ts, TH2F* h,  
        Double_t var_low, Double_t var_high, int flag1 = FLAG_MUONS, int flag2 = FLAG_M_BINS,
        bool turn_on_RC = false){
    h->Sumw2();

    BTag_readers b_reader;
    BTag_effs btag_effs;

    setup_btag_SFs(&b_reader, &btag_effs);
    for(int i=0; i<nTrees; i++){
        TTree *t1 = ts[i];
        Long64_t nEntries  =  t1->GetEntries();

        Double_t m, xF, cost, mu1_pt, mu2_pt, jet1_cmva, jet2_cmva, gen_weight;
        Double_t bcdef_HLT_SF, bcdef_iso_SF, bcdef_id_SF;
        Double_t gh_HLT_SF, gh_iso_SF, gh_id_SF, el_id_SF, el_reco_SF, el_HLT_SF;
        Double_t gh_trk_SF, bcdef_trk_SF;
        Double_t jet1_pt, jet2_pt, jet1_b_weight, jet2_b_weight, pu_SF, jet1_eta, jet2_eta;
        Int_t jet1_flavour, jet2_flavour;
        Double_t mu1_pt_corr, mu2_pt_corr;
        Float_t cost_pt, met_pt;
        Int_t nJets;
        TLorentzVector *lep_p=0;
        TLorentzVector *lep_m=0;
        Double_t pt;

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
        t1->SetBranchAddress("jet1_eta", &jet1_eta);
        t1->SetBranchAddress("jet2_eta", &jet2_eta);
        t1->SetBranchAddress("jet1_flavour", &jet1_flavour);
        t1->SetBranchAddress("jet2_flavour", &jet2_flavour);
        t1->SetBranchAddress("gen_weight", &gen_weight);
        t1->SetBranchAddress("pu_SF", &pu_SF);

        bool do_bweight = true;
        jet1_b_weight = 1.;
        jet2_b_weight = 1.;
        if(flag1 == FLAG_MUONS){

            t1->SetBranchAddress("bcdef_HLT_SF", &bcdef_HLT_SF);
            t1->SetBranchAddress("bcdef_iso_SF", &bcdef_iso_SF);
            t1->SetBranchAddress("bcdef_id_SF", &bcdef_id_SF);
            t1->SetBranchAddress("gh_HLT_SF", &gh_HLT_SF);
            t1->SetBranchAddress("gh_iso_SF", &gh_iso_SF);
            t1->SetBranchAddress("gh_id_SF", &gh_id_SF);
            t1->SetBranchAddress("gh_trk_SF", &gh_trk_SF);
            t1->SetBranchAddress("bcdef_trk_SF", &bcdef_trk_SF);
            t1->SetBranchAddress("mu_p", &lep_p);
            t1->SetBranchAddress("mu_m", &lep_m);
            t1->SetBranchAddress("mu1_pt_corr", &mu1_pt_corr);
            t1->SetBranchAddress("mu2_pt_corr", &mu2_pt_corr);


            for (int i=0; i<nEntries; i++) {
                t1->GetEntry(i);
                bool no_bjets = has_no_bjets(nJets, jet1_pt, jet2_pt, jet1_cmva, jet2_cmva);
                if(flag2 == FLAG_PT_BINS){
                    TLorentzVector cm = *lep_p + *lep_m;
                    pt = cm.Pt();
                }
                if(turn_on_RC){
                    TLorentzVector mu_p_new, mu_m_new, cm;
                    if((lep_p->Pt() > lep_m->Pt() && mu1_pt_corr > mu2_pt_corr) ||
                            (lep_p->Pt() < lep_m->Pt() && mu1_pt_corr < mu2_pt_corr)){
                        mu_p_new.SetPtEtaPhiE(mu1_pt_corr, lep_p->Eta(), lep_p->Phi(), lep_p->E());
                        mu_m_new.SetPtEtaPhiE(mu2_pt_corr, lep_m->Eta(), lep_m->Phi(), lep_m->E());
                    }
                    else{
                        mu_p_new.SetPtEtaPhiE(mu2_pt_corr, lep_p->Eta(), lep_p->Phi(), lep_p->E());
                        mu_m_new.SetPtEtaPhiE(mu1_pt_corr, lep_m->Eta(), lep_m->Phi(), lep_m->E());
                    }
                    double new_cost = get_cost_v2(mu_p_new, mu_m_new);
                    cm = mu_p_new + mu_m_new;
                    //if(cm.M() < 150.) continue;
                    cost = new_cost;
                    xF = abs(2.*cm.Pz()/13000.); 
                    m = cm.M();
                    pt = cm.Pt();
                       
                }
                bool pass = ((flag2 == FLAG_M_BINS && m >= var_low && m <= var_high) ||
                        (flag2 == FLAG_PT_BINS && m >= 150. && pt >= var_low && pt <= var_high))
                    && met_pt < 50.  && no_bjets;

                if(pass){
                    cost = get_cost_v2(*lep_p, *lep_m);

                    if(do_bweight){
                        jet1_b_weight = get_btag_weight(jet1_pt, jet1_eta,(Float_t) jet1_flavour , btag_effs, b_reader, 0);
                        jet2_b_weight = get_btag_weight(jet2_pt, jet2_eta,(Float_t) jet2_flavour , btag_effs, b_reader, 0);
                    }
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
                    Double_t evt_weight = gh_weight * gh_lumi + bcdef_weight * bcdef_lumi * 1000 *emu_scaling;
                    h->Fill(xF, cost, evt_weight);
                }
            }
        }
        else if(flag1 == FLAG_ELECTRONS) {
            t1->SetBranchAddress("el_p", &lep_p);
            t1->SetBranchAddress("el_m", &lep_m);
            t1->SetBranchAddress("el_id_SF", &el_id_SF);
            t1->SetBranchAddress("el_HLT_SF", &el_HLT_SF);
            t1->SetBranchAddress("el_reco_SF", &el_reco_SF);
            for (int i=0; i<nEntries; i++) {
                t1->GetEntry(i);
                bool no_bjets = has_no_bjets(nJets, jet1_pt, jet2_pt, jet1_cmva, jet2_cmva);
                if(flag2 == FLAG_PT_BINS){
                    TLorentzVector cm = *lep_p + *lep_m;
                    pt = cm.Pt();
                }
                bool pass = ((flag2 == FLAG_M_BINS && m >= var_low && m <= var_high) ||
                        (flag2 == FLAG_PT_BINS && m >= 150. && pt >= var_low && pt <= var_high))
                    && met_pt < 50.  && no_bjets;
                if(pass){

                    cost = get_cost_v2(*lep_p, *lep_m);
                    if(do_bweight){
                        jet1_b_weight = get_btag_weight(jet1_pt, jet1_eta,(Float_t) jet1_flavour , btag_effs, b_reader, 0);
                        jet2_b_weight = get_btag_weight(jet2_pt, jet2_eta,(Float_t) jet2_flavour , btag_effs, b_reader, 0);
                    }

                    Double_t evt_weight = gen_weight * pu_SF * el_id_SF * el_reco_SF * el_HLT_SF;
                    if (nJets >= 1){
                        evt_weight *= jet1_b_weight;
                    }
                    if (nJets >= 2){
                        evt_weight *= jet2_b_weight;
                    }
                    evt_weight *= 1000*tot_lumi*emu_scaling;
                    h->Fill(xF, cost, evt_weight);
                }
            }
        }

        t1->ResetBranchAddresses();
    }
    printf("Performing templ. cleanup (removing neg. bins) \n");
    cleanup_template(h);
    printf("Tot Weight is %.2f \n", h->Integral());
    h->Scale(1./h->Integral());
    return 0;
}





