//construct templates from reconstructed events
//
//
//

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <cmath>
#include <cfloat>
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
#include "TH2.h"
#include "TH3.h"
#include "TTree.h"
#include "TFitter.h"
#include "TSystem.h"
#include "Math/Functor.h"
#include "ScaleFactors.C"
#include "BTagUtils.C"
#include "HistMaker.C"


using namespace std;

el_SFs el_SF;
mu_SFs runs_bcdef, runs_gh;
pileup_SFs pu_SFs;
pileup_systematics pu_sys;

BTag_readers b_reader;
BTag_effs btag_effs;
bool btag_setup = false;

void setup_all_SFs(){
    setup_btag_SFs(&b_reader, &btag_effs);
    setup_el_SF(&el_SF);
    setup_SFs(&runs_bcdef, &runs_gh, &pu_SFs);
    setup_pileup_systematic(&pu_sys); 
}


void cleanup_template(TH2F *h){
    //printf("%i %i \n", h->GetNbinsX(), h->GetNbinsY());
    for(int i=0; i<= h->GetNbinsX()+1; i++){
        for(int j=0; j<= h->GetNbinsY()+1; j++){
            //printf("%i %i \n", i,j);
            float val = h->GetBinContent(i,j);
            if(val< 0.) h->SetBinContent(i,j,0.);
        }
    }
}

void print_hist(TH2 *h){
    printf("\n");
    for(int i=1; i<= h->GetNbinsX(); i++){
        for(int j=1; j<= h->GetNbinsY(); j++){
            printf("%.2e ",   (float) h->GetBinContent(i,j));
        }
        printf("\n");
    }

}
//
//static type means functions scope is only this file, to avoid conflicts

int gen_data_template(TTree *t1, TH2F* h,  Double_t var_low, Double_t var_high, 
        int flag1 = FLAG_MUONS, int flag2 = FLAG_M_BINS, bool turn_on_RC = true, bool ss = false){
    h->Sumw2();
    Long64_t nEntries  =  t1->GetEntries();
    //printf("size is %i \n", nEntries);
    Double_t m, xF, cost, gen_weight, jet1_cmva, jet2_cmva,
             jet1_pt, jet2_pt, lep1_pt, lep2_pt, lep1_pt_corr, lep2_pt_corr,
             lep1_pt_alt, lep2_pt_alt, pt;
    Double_t mu_p_SF, mu_m_SF, mu_p_SF_alt, mu_m_SF_alt;
    Float_t met_pt;
    Int_t nJets;
    TLorentzVector *lep_p = 0;
    TLorentzVector *lep_m = 0;
    TLorentzVector cm;
    t1->SetBranchAddress("m", &m);
    t1->SetBranchAddress("xF", &xF);
    t1->SetBranchAddress("cost", &cost);
    t1->SetBranchAddress("nJets", &nJets);
    t1->SetBranchAddress("jet1_CMVA", &jet1_cmva);
    t1->SetBranchAddress("jet2_CMVA", &jet2_cmva);
    t1->SetBranchAddress("met_pt", &met_pt);
    t1->SetBranchAddress("jet1_pt", &jet1_pt);
    t1->SetBranchAddress("jet2_pt", &jet2_pt);
    if(flag1 == FLAG_MUONS){
        t1->SetBranchAddress("mu_p", &lep_p);
        t1->SetBranchAddress("mu_m", &lep_m);
        t1->SetBranchAddress("mu1_pt", &lep1_pt);
        t1->SetBranchAddress("mu2_pt", &lep2_pt);
        if(turn_on_RC){
            t1->SetBranchAddress("mu_p_SF", &mu_p_SF);
            t1->SetBranchAddress("mu_m_SF", &mu_m_SF);
            t1->SetBranchAddress("mu_p_SF_alt", &mu_p_SF_alt);
            t1->SetBranchAddress("mu_m_SF_alt", &mu_m_SF_alt);
        }
    }
    else{
        t1->SetBranchAddress("el_p", &lep_p);
        t1->SetBranchAddress("el_m", &lep_m);
        t1->SetBranchAddress("el1_pt", &lep1_pt);
        t1->SetBranchAddress("el2_pt", &lep2_pt);
    }
    nJets =2;
    int nEvents = 0;
    int n=0;
    for (int i=0; i<nEntries; i++) {
        t1->GetEntry(i);
        cost = get_cost(*lep_p, *lep_m);
        bool no_bjets = has_no_bjets(nJets, jet1_pt, jet2_pt, jet1_cmva, jet2_cmva);
        bool not_cosmic = notCosmic(*lep_p, *lep_m);
        cm = *lep_p + *lep_m;
        pt = cm.Pt();

        if(flag1 == FLAG_MUONS && turn_on_RC){
            TLorentzVector mu_p_new, mu_m_new;
            mu_p_new.SetPtEtaPhiE(lep_p->Pt() * mu_p_SF, lep_p->Eta(), lep_p->Phi(), mu_p_SF * lep_p->E());
            mu_m_new.SetPtEtaPhiE(lep_m->Pt() * mu_m_SF, lep_m->Eta(), lep_m->Phi(), mu_m_SF * lep_m->E());

            double new_cost = get_cost(mu_p_new, mu_m_new);
            cm = mu_p_new + mu_m_new;
            //if(cm.M() < 150.) continue;
            cost = new_cost;
            m = cm.M();
            pt = cm.Pt();
            xF = compute_xF(cm); 

        }
        if(flag2 == FLAG_M_BINS){

            if(m >= var_low && m <= var_high && met_pt < 50. && no_bjets && not_cosmic){
                n++;
                if(!ss) h->Fill(xF, cost, 1); 
                else{
                    h->Fill(xF, -abs(cost), 1);
                }

                //printf("size %i \n", (int) v_cost->size());
                //printf("xf\n");
                //printf("end\n");
                nEvents++;
            }
        }
        else{
            if(m>= 150. &&  pt >= var_low && pt <= var_high && met_pt < 50. && no_bjets && not_cosmic){
                n++;
                if(!ss) h->Fill(xF, cost, 1); 
                else{
                    h->Fill(xF, -abs(cost), 1);
                }
                nEvents++;
            }
        }

    }

    t1->ResetBranchAddresses();
    return nEvents;
}



int gen_mc_template(TTree *t1, Double_t alpha, TH2F* h_sym, TH2F *h_asym, 
        Double_t var_low, Double_t var_high, int flag1 = FLAG_MUONS, int flag2 = FLAG_M_BINS, bool turn_on_RC = true,
        const string &sys_label = "" ){
    Long64_t nEntries  =  t1->GetEntries();

    h_sym->Sumw2();
    h_asym->Sumw2();

    //TH2F* h_reweights = (TH2F *) h_sym->Clone("h_rw");

    Double_t m, xF, cost, gen_weight, reweight, jet1_cmva, jet2_cmva, cost_st;
    Double_t bcdef_HLT_SF, bcdef_iso_SF, bcdef_id_SF;
    Double_t gh_HLT_SF, gh_iso_SF, gh_id_SF;
    Double_t gh_trk_SF, bcdef_trk_SF;
    Double_t el_id_SF, el_reco_SF, pu_SF, el_HLT_SF;
    Double_t jet1_pt, jet2_pt, jet1_b_weight, jet2_b_weight, jet1_eta, jet2_eta;
    Double_t mu1_pt, mu1_eta, mu2_pt, mu2_eta;
    Double_t el1_pt, el1_eta, el2_pt, el2_eta;
    Double_t mu_R_up, mu_R_down, mu_F_up, mu_F_down, 
             mu_RF_up, mu_RF_down, pdf_up, pdf_down;
    Double_t mu_p_SF, mu_m_SF, mu_p_SF_up, mu_m_SF_up, mu_p_SF_down, mu_m_SF_down, alphaS_up, alphaS_down;
    Float_t cost_pt, met_pt;
    Float_t pdf_weights[60];
    TLorentzVector *lep_p=0;
    TLorentzVector *lep_m=0;
    Double_t pt;
    Int_t nJets, pu_NtrueInt, jet1_flavour, jet2_flavour;
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
    t1->SetBranchAddress("mu_R_up", &mu_R_up);
    t1->SetBranchAddress("mu_R_down", &mu_R_down);
    t1->SetBranchAddress("mu_F_up", &mu_F_up);
    t1->SetBranchAddress("mu_F_down", &mu_F_down);
    t1->SetBranchAddress("mu_RF_up", &mu_RF_up);
    t1->SetBranchAddress("mu_RF_down", &mu_RF_down);
    t1->SetBranchAddress("pu_NtrueInt", &pu_NtrueInt);
    t1->SetBranchAddress("jet1_eta", &jet1_eta);
    t1->SetBranchAddress("jet2_eta", &jet2_eta);
    t1->SetBranchAddress("jet1_flavour", &jet1_flavour);
    t1->SetBranchAddress("jet2_flavour", &jet2_flavour);
    t1->SetBranchAddress("pdf_weights", &pdf_weights);
    t1->SetBranchAddress("alpha_up", &alphaS_up);
    t1->SetBranchAddress("alpha_down", &alphaS_down);
    int n = 0;

    //SYSTEMATICS 
    Double_t one = 1.0;
    Double_t *systematic = &one;
    int shift = 0;
    if(!sys_label.empty()){
        //doing some systematic
        if(sys_label.find("Up") != string::npos){
            shift = 1;
        }
        else if(sys_label.find("Down") != string::npos){
            shift = -1;
        }

        else{
            printf("systematic label not empty, but doesn't have UP or DOWN \n");
        }
    }
    // do_sys: 0 = nominal, 1= var up, -1 = var down
    int do_pdf_sys = 0;
    int do_btag_sys = 0;
    int do_pileup_sys = 0;

    int do_muHLT_sys = 0;
    int do_muID_sys = 0;
    int do_muISO_sys = 0;
    int do_muTRK_sys = 0;
    int do_muRC_sys = 0;

    int do_elID_sys = 0;
    int do_elHLT_sys = 0;
    int do_elRECO_sys = 0;
    int do_elScale_sys = 0;
    int do_elSmear_sys = 0;
    float elp_rescale, elm_rescale;



    if(shift !=0){
        if(sys_label.find("BTAG") != string::npos) do_btag_sys = shift;
        else if(sys_label.find("Pu") != string::npos) do_pileup_sys = shift;
        else if(sys_label.find("muHLT") != string::npos) do_muHLT_sys = shift;
        else if(sys_label.find("muID") != string::npos) do_muID_sys = shift;
        else if(sys_label.find("muISO") != string::npos) do_muISO_sys = shift;
        else if(sys_label.find("muTRK") != string::npos) do_muTRK_sys = shift;
        else if(sys_label.find("muRC") != string::npos) do_muRC_sys = shift;

        else if(sys_label.find("elID") != string::npos) do_elID_sys = shift;
        else if(sys_label.find("elHLT") != string::npos) do_elHLT_sys = shift;
        else if(sys_label.find("elRECO") != string::npos) do_elRECO_sys = shift;
        else if(sys_label.find("elScale") != string::npos) do_elScale_sys = shift;
        else if(sys_label.find("elSmear") != string::npos) do_elSmear_sys = shift;


        else if(sys_label.find("RENORM") != string::npos && shift > 0) systematic = &mu_R_up;
        else if(sys_label.find("RENORM") != string::npos && shift < 0) systematic = &mu_R_down;
        else if(sys_label.find("FAC") != string::npos && shift > 0) systematic = &mu_F_up;
        else if(sys_label.find("FAC") != string::npos && shift < 0) systematic = &mu_F_down;
        else if(sys_label.find("alphaS") != string::npos && shift < 0) systematic = &alphaS_down;
        else if(sys_label.find("alphaS") != string::npos && shift > 0) systematic = &alphaS_up;
        else if(sys_label.find("alpha") != string::npos && shift < 0) systematic = &one;

        else if(sys_label.find("pdf") != string::npos){
            if(shift > 0) sscanf(sys_label.c_str(), "_pdf%iUp", &do_pdf_sys);
            else sscanf(sys_label.c_str(), "_pdf%iDown", &do_pdf_sys);
            printf("Doing pdf sys %i \n", do_pdf_sys);
        }
        
        else printf("COULDN'T PARSE SYSTEMATIC %s !!! \n \n", sys_label.c_str());
    }







    if(flag1 == FLAG_MUONS){
        t1->SetBranchAddress("mu_p", &lep_p);
        t1->SetBranchAddress("mu_m", &lep_m);
        t1->SetBranchAddress("mu1_pt", &mu1_pt);
        t1->SetBranchAddress("mu1_eta", &mu1_eta);
        t1->SetBranchAddress("mu2_pt", &mu2_pt);
        t1->SetBranchAddress("mu2_eta", &mu2_eta);
        t1->SetBranchAddress("bcdef_HLT_SF", &bcdef_HLT_SF);
        t1->SetBranchAddress("bcdef_iso_SF", &bcdef_iso_SF);
        t1->SetBranchAddress("bcdef_id_SF", &bcdef_id_SF);
        t1->SetBranchAddress("bcdef_trk_SF", &bcdef_trk_SF);
        t1->SetBranchAddress("gh_HLT_SF", &gh_HLT_SF);
        t1->SetBranchAddress("gh_iso_SF", &gh_iso_SF);
        t1->SetBranchAddress("gh_id_SF", &gh_id_SF);
        t1->SetBranchAddress("gh_trk_SF", &gh_trk_SF);
        if(turn_on_RC){
            t1->SetBranchAddress("mu_p_SF", &mu_p_SF);
            t1->SetBranchAddress("mu_m_SF", &mu_m_SF);
            t1->SetBranchAddress("mu_p_SF_up", &mu_p_SF_up);
            t1->SetBranchAddress("mu_m_SF_up", &mu_m_SF_up);
            t1->SetBranchAddress("mu_p_SF_down", &mu_p_SF_down);
            t1->SetBranchAddress("mu_m_SF_down", &mu_m_SF_down);
        }


        //separate SFs for runs BCDEF and GH


        for (int i=0; i<nEntries; i++) {
            t1->GetEntry(i);
            bool no_bjets = has_no_bjets(nJets, jet1_pt, jet2_pt, jet1_cmva, jet2_cmva);
            bool not_cosmic = notCosmic(*lep_p, *lep_m);
            if(flag2 == FLAG_PT_BINS){
                TLorentzVector cm = *lep_p + *lep_m;
                pt = cm.Pt();
            }
            if(turn_on_RC){
                TLorentzVector mu_p_new, mu_m_new, cm;

                if(do_muRC_sys > 0){
                    mu_p_SF = mu_p_SF_up;
                    mu_m_SF = mu_m_SF_up;
                }
                else if(do_muRC_sys < 0){
                    mu_p_SF = mu_p_SF_down;
                    mu_m_SF = mu_m_SF_down;
                }
                mu_p_new.SetPtEtaPhiE(lep_p->Pt() * mu_p_SF, lep_p->Eta(), lep_p->Phi(), mu_p_SF * lep_p->E());
                mu_m_new.SetPtEtaPhiE(lep_m->Pt() * mu_m_SF, lep_m->Eta(), lep_m->Phi(), mu_m_SF * lep_m->E());

                cost = get_cost(mu_p_new, mu_m_new);
                cm = mu_p_new + mu_m_new;
                m = cm.M();
                pt = cm.Pt();
                xF = compute_xF(cm); 
                //if(cm.M() < 150.) continue;

            }
            else cost = get_cost(*lep_p, *lep_m);
            if(cost_st>0.) cost_st = fabs(cost);
            else cost_st = -fabs(cost);
            bool pass = ((flag2 == FLAG_M_BINS && m >= var_low && m <= var_high) ||
                    (flag2 == FLAG_PT_BINS && m >= 150. && pt >= var_low && pt <= var_high))
                && met_pt < 50.  && no_bjets && not_cosmic;
            if(pass){
                reweight = (4./3.)*cost_st*(2. + alpha)/
                    (1. + cost_st*cost_st + alpha*(1.- cost_st*cost_st));
                n++;
                Double_t pu_SF_sys = 1.;
                if(do_pileup_sys == -1) pu_SF_sys = get_pileup_SF(pu_NtrueInt, pu_sys.pileup_down);
                if(do_pileup_sys == 1) pu_SF_sys = get_pileup_SF(pu_NtrueInt, pu_sys.pileup_up);
                if(do_pdf_sys != 0){
                    if(shift > 0) *systematic = pdf_weights[do_pdf_sys-1];
                    if(shift < 0) *systematic = 2. - pdf_weights[do_pdf_sys-1];
                }
                gen_weight *= *systematic * pu_SF * pu_SF_sys;


                if(do_muHLT_sys)   bcdef_HLT_SF = get_HLT_SF(mu1_pt, mu1_eta, mu2_pt, mu2_eta, runs_bcdef.HLT_SF, runs_bcdef.HLT_MC_EFF, do_muHLT_sys);
                if(do_muHLT_sys)   gh_HLT_SF = get_HLT_SF(mu1_pt, mu1_eta, mu2_pt, mu2_eta, runs_gh.HLT_SF, runs_gh.HLT_MC_EFF, do_muHLT_sys);

                if(do_muISO_sys)   bcdef_iso_SF = get_SF(mu1_pt, mu1_eta, runs_bcdef.ISO_SF, do_muISO_sys) * get_SF(mu2_pt, mu2_eta, runs_bcdef.ISO_SF, do_muISO_sys);
                if(do_muISO_sys)   gh_iso_SF = get_SF(mu1_pt, mu1_eta, runs_gh.ISO_SF, do_muISO_sys) * get_SF(mu2_pt, mu2_eta, runs_gh.ISO_SF, do_muISO_sys);

                if(do_muID_sys)    bcdef_id_SF = get_SF(mu1_pt, mu1_eta, runs_bcdef.ID_SF, do_muID_sys) * get_SF(mu2_pt, mu2_eta, runs_bcdef.ID_SF, do_muID_sys);
                if(do_muID_sys)    gh_id_SF = get_SF(mu1_pt, mu1_eta, runs_gh.ID_SF, do_muID_sys) * get_SF(mu2_pt, mu2_eta, runs_gh.ID_SF, do_muID_sys);

                if(do_muTRK_sys)   bcdef_trk_SF = get_Mu_trk_SF(abs(mu1_eta), runs_bcdef.TRK_SF, do_muTRK_sys) * get_Mu_trk_SF(abs(mu2_eta), runs_bcdef.TRK_SF, do_muTRK_sys);
                if(do_muTRK_sys)   gh_trk_SF = get_Mu_trk_SF(abs(mu1_eta), runs_gh.TRK_SF, do_muTRK_sys) * get_Mu_trk_SF(abs(mu2_eta), runs_gh.TRK_SF, do_muTRK_sys);
                //bcdef_HLT_SF = gh_HLT_SF = bcdef_iso_SF = gh_iso_SF = bcdef_id_SF = gh_id_SF = bcdef_trk_SF = gh_trk_SF = 1.0;
                jet1_b_weight = get_btag_weight(jet1_pt, jet1_eta, (Float_t) jet1_flavour , btag_effs, b_reader, do_btag_sys);
                jet2_b_weight = get_btag_weight(jet2_pt, jet2_eta, (Float_t) jet2_flavour , btag_effs, b_reader, do_btag_sys);


                Double_t bcdef_weight = gen_weight * bcdef_HLT_SF * bcdef_iso_SF * bcdef_id_SF * bcdef_trk_SF;
                Double_t gh_weight = gen_weight * gh_HLT_SF * gh_iso_SF * gh_id_SF * gh_trk_SF;
                if (nJets >= 1){
                    bcdef_weight *= jet1_b_weight;
                    gh_weight *= jet1_b_weight;
                }
                if (nJets >= 2){
                    bcdef_weight *= jet2_b_weight;
                    gh_weight *= jet2_b_weight;
                }


                Double_t final_weight = 1000*(bcdef_weight*bcdef_lumi + gh_weight*gh_lumi);
                h_sym->Fill(xF, cost, final_weight); 
                h_sym->Fill(xF, -cost, final_weight); 

                h_asym->Fill(xF, cost, reweight * final_weight);
                h_asym->Fill(xF, -cost, -reweight * final_weight);
                //h_reweights->Fill(xF, fabs(cost), fabs(reweight));

            }
        }

    }
    else if (flag1 == FLAG_ELECTRONS) {
        t1->SetBranchAddress("el_p", &lep_p);
        t1->SetBranchAddress("el_m", &lep_m);
        t1->SetBranchAddress("el1_pt", &el1_pt);
        t1->SetBranchAddress("el1_eta", &el1_eta);
        t1->SetBranchAddress("el2_pt", &el2_pt);
        t1->SetBranchAddress("el2_eta", &el2_eta);
        t1->SetBranchAddress("el_id_SF", &el_id_SF);
        t1->SetBranchAddress("el_HLT_SF", &el_HLT_SF);
        t1->SetBranchAddress("el_reco_SF", &el_reco_SF);
        if(do_elScale_sys || do_elSmear_sys){
            if(do_elScale_sys >0){
                t1->SetBranchAddress("elp_scale_up", &elp_rescale);
                t1->SetBranchAddress("elm_scale_up", &elm_rescale);
            }
            else if(do_elScale_sys < 0){
                t1->SetBranchAddress("elp_scale_down", &elp_rescale);
                t1->SetBranchAddress("elm_scale_down", &elm_rescale);
            }
            else if(do_elSmear_sys < 0){
                t1->SetBranchAddress("elp_smear_down", &elp_rescale);
                t1->SetBranchAddress("elm_smear_down", &elm_rescale);
            }
            else if(do_elSmear_sys > 0){
                t1->SetBranchAddress("elp_smear_up", &elp_rescale);
                t1->SetBranchAddress("elm_smear_up", &elm_rescale);
            }

        }


        //separate SFs for runs BCDEF and GH
        for (int i=0; i<nEntries; i++) {
            t1->GetEntry(i);
            bool no_bjets = has_no_bjets(nJets, jet1_pt, jet2_pt, jet1_cmva, jet2_cmva);
            bool not_cosmic = notCosmic(*lep_p, *lep_m);
            if(flag2 == FLAG_PT_BINS){
                TLorentzVector cm = *lep_p + *lep_m;
                pt = cm.Pt();
            }
            if(do_elScale_sys || do_elSmear_sys){
                lep_p->SetPtEtaPhiE(lep_p->Pt() * (elp_rescale), lep_p->Eta(), lep_p->Phi(), lep_p->E() *(elp_rescale));
                lep_m->SetPtEtaPhiE(lep_m->Pt() * (elm_rescale), lep_m->Eta(), lep_m->Phi(), lep_m->E() *(elm_rescale));
                TLorentzVector cm = *lep_p + *lep_m;
                m = cm.M();
                xF = compute_xF(cm);

            }
            bool pass = ((flag2 == FLAG_M_BINS && m >= var_low && m <= var_high) ||
                    (flag2 == FLAG_PT_BINS && m >= 150. && pt >= var_low && pt <= var_high))
                && met_pt < 50.  && no_bjets && not_cosmic;
            if(pass){
                cost = get_cost(*lep_p, *lep_m);
                if(cost_st>0.) cost_st = fabs(cost);
                else cost_st = -fabs(cost);
                reweight = (4./3.)*cost_st*(2. + alpha)/
                    (1. + cost_st*cost_st + alpha*(1.- cost_st*cost_st));
                n++;
                Double_t pu_SF_sys = 1.;
                if(do_pileup_sys == -1) pu_SF_sys = get_pileup_SF(pu_NtrueInt, pu_sys.pileup_down);
                if(do_pileup_sys == 1) pu_SF_sys = get_pileup_SF(pu_NtrueInt, pu_sys.pileup_up);

                if(do_pdf_sys != 0){
                    if(shift > 0) *systematic = pdf_weights[do_pdf_sys-1];
                    if(shift < 0) *systematic = 2. - pdf_weights[do_pdf_sys-1];
                }
                gen_weight *= *systematic * pu_SF_sys * pu_SF;


                if(do_elID_sys) el_id_SF = get_el_SF(el1_pt, el1_eta, el_SF.ID_SF, do_elID_sys) * get_el_SF(el2_pt, el2_eta, el_SF.ID_SF, do_elID_sys);
                if(do_elRECO_sys) el_reco_SF = get_el_SF(el1_pt, el1_eta, el_SF.RECO_SF, do_elRECO_sys) * get_el_SF(el2_pt, el2_eta, el_SF.RECO_SF, do_elRECO_sys);
                if(do_elHLT_sys) el_HLT_SF = get_el_HLT_SF(el1_pt, el1_eta, el2_pt, el2_eta, el_SF.HLT_SF, el_SF.HLT_MC_EFF, do_elHLT_sys);
                jet1_b_weight = get_btag_weight(jet1_pt, jet1_eta,(Float_t) jet1_flavour , btag_effs, b_reader, do_btag_sys);
                jet2_b_weight = get_btag_weight(jet2_pt, jet2_eta,(Float_t) jet2_flavour , btag_effs, b_reader, do_btag_sys);


                Double_t evt_weight = gen_weight * el_id_SF *el_reco_SF * el_HLT_SF * 1000. * el_lumi;
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

            }
        }

    }

    //float norm = h_sym -> Integral();
    //h_sym->Scale(1./norm);
    //h_asym->Scale(1./norm);
    t1->ResetBranchAddresses();
    printf("MC templates generated from %i events \n \n", n);
    return 0;
}




float gen_fakes_template(TTree *t_WJets, TTree *t_QCD, TTree *t_WJets_contam, 
        TTree* t_QCD_contam, TH2F *h, Double_t var_low, Double_t var_high, 
        int flag1 = FLAG_MUONS, int flag2 = FLAG_M_BINS, bool ss = false){
    h->Sumw2();
    TH2D *h_err;
    TLorentzVector *lep_p=0;
    TLorentzVector *lep_m=0;
    Double_t pt;
    FakeRate FR;
    float tot_weight_os = 0.;
    float tot_weight_ss = 0.;
    if(flag1 == FLAG_MUONS){
        //TH2D *FR;
        setup_new_mu_fakerate(&FR);
        h_err = (TH2D *) FR.h->Clone("h_err");
        h_err->Reset();
        for (int l=0; l<=3; l++){
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
            Float_t mu1_charge, mu2_charge;
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
            t->SetBranchAddress("mu1_charge", &mu1_charge);
            t->SetBranchAddress("mu2_charge", &mu2_charge);
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
                bool not_cosmic = notCosmic(*lep_p, *lep_m);
                if(l==0){
                    if(iso_mu ==1) mu1_fakerate = get_new_fakerate_prob(mu1_pt, mu1_eta, FR.h);
                    if(iso_mu ==0) mu1_fakerate = get_new_fakerate_prob(mu2_pt, mu2_eta, FR.h);
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
                    if(iso_mu ==1) mu1_fakerate = get_new_fakerate_prob(mu1_pt, mu1_eta, FR.h);
                    if(iso_mu ==0) mu1_fakerate = get_new_fakerate_prob(mu2_pt, mu2_eta, FR.h);
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
                cost = get_cost(*lep_p, *lep_m);
                TLorentzVector cm = *lep_p + *lep_m;
                pt = cm.Pt();
                xF = compute_xF(cm); 
                bool pass = ((flag2 == FLAG_M_BINS && m >= var_low && m <= var_high) ||
                        (flag2 == FLAG_PT_BINS && m >= 150. && pt >= var_low && pt <= var_high))
                    && met_pt < 50.  && no_bjets && not_cosmic;

                bool opp_sign = ((abs(mu1_charge - mu2_charge)) > 0.01);
                if(!ss) pass = pass && opp_sign;
                if(pass){
                    if(l==0 && iso_mu ==1) h_err->Fill(min(abs(mu1_eta), 2.3), min(mu1_pt, 150.));
                    if(l==0 && iso_mu ==0) h_err->Fill(min(abs(mu2_eta), 2.3), min(mu2_pt, 150.));
                    cost = get_cost(*lep_p, *lep_m);
                    //if(l==3) printf("Evt fr %.2e \n", evt_fakerate);
                    //if(l==3) printf("cost, fr %.2f %.2e \n", cost, evt_fakerate);
                    if(!ss) h->Fill(xF, cost, evt_fakerate);
                    else{
                        h->Fill(xF, -abs(cost), evt_fakerate);
                    }
                    if(opp_sign) tot_weight_os += evt_fakerate;
                    else tot_weight_ss += evt_fakerate;
                }


            }
            printf("After iter %i current fakerate est is %.0f \n", l, h->Integral());
        }
    }
    else{

        //TH2D *FR;
        setup_new_el_fakerate(&FR);
        h_err = (TH2D *) FR.h->Clone("h_err");
        h_err->Reset();
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
            Float_t el1_charge, el2_charge;
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
            t->SetBranchAddress("el1_charge", &el1_charge);
            t->SetBranchAddress("el2_charge", &el2_charge);
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
                bool not_cosmic = notCosmic(*lep_p, *lep_m);
                if(l==0){
                    if(iso_el ==1) el1_fakerate = get_new_fakerate_prob(el1_pt, el1_eta, FR.h);
                    if(iso_el ==0) el1_fakerate = get_new_fakerate_prob(el2_pt, el2_eta, FR.h);
                    evt_fakerate = el1_fakerate/(1-el1_fakerate);
                }
                if(l==1){
                    el1_fakerate = get_new_fakerate_prob(el1_pt, el1_eta, FR.h);
                    el2_fakerate = get_new_fakerate_prob(el2_pt, el2_eta, FR.h);
                    evt_fakerate = -(el1_fakerate/(1-el1_fakerate)) * (el2_fakerate/(1-el2_fakerate));
                }
                if(l==2){

                    Double_t mc_weight = gen_weight * el_id_SF * el_reco_SF  * 1000. * el_lumi;
                    if(iso_el ==1) el1_fakerate = get_new_fakerate_prob(el1_pt, el1_eta, FR.h);
                    if(iso_el ==0) el1_fakerate = get_new_fakerate_prob(el2_pt, el2_eta, FR.h);
                    evt_fakerate = -(el1_fakerate * mc_weight)/(1-el1_fakerate);
                }
                if(l==3){
                    Double_t mc_weight = gen_weight * el_id_SF * el_reco_SF * 1000. * el_lumi;

                    el1_fakerate = get_new_fakerate_prob(el1_pt, el1_eta, FR.h);
                    el2_fakerate = get_new_fakerate_prob(el2_pt, el2_eta, FR.h);
                    evt_fakerate = mc_weight * (el1_fakerate/(1-el1_fakerate)) * (el2_fakerate/(1-el2_fakerate));
                }


                cost = get_cost(*lep_p, *lep_m);
                TLorentzVector cm = *lep_p + *lep_m;
                pt = cm.Pt();
                xF = compute_xF(cm); 
                bool pass = ((flag2 == FLAG_M_BINS && m >= var_low && m <= var_high) ||
                        (flag2 == FLAG_PT_BINS && m >= 150. && pt >= var_low && pt <= var_high))
                    && met_pt < 50.  && no_bjets && not_cosmic;
                bool opp_sign = ((abs(el1_charge - el2_charge)) > 0.01);
                if(!ss) pass = pass && opp_sign;
                if(pass){
                    if(l==0 && iso_el ==1) h_err->Fill(min(abs(el1_eta), 2.3), min(el1_pt, 150.));
                    if(l==0 && iso_el ==0) h_err->Fill(min(abs(el2_eta), 2.3), min(el2_pt, 150.));
                    //if(l==3) printf("Evt fr %.2e \n", evt_fakerate);
                    //if(l==3) printf("cost, fr %.2f %.2e \n", cost, evt_fakerate);
                    if(!ss) h->Fill(xF, cost, evt_fakerate);
                    else{
                        h->Fill(xF, -abs(cost), evt_fakerate);
                    }
                    if(opp_sign) tot_weight_os += evt_fakerate;
                    else tot_weight_ss += evt_fakerate;
                }
            }

            printf("After iter %i current fakerate est is %.0f \n", l, h->Integral());
        }
    }
    printf("Performing fakes cleanup (removing neg. bins) \n");
    cleanup_template(h);
    
    set_fakerate_errors(h_err, FR.h, h);
    //scale ss electrons by 0.5, muons get a more complicated norm
    printf("Total Fakerate weight Weight is %.2f \n", h->Integral());
    if(h->Integral() < 0.){
        h->Scale(0.);
        printf("zeroing Fakes template \n");
    }
    float scaling = tot_weight_os / (tot_weight_ss + tot_weight_os);
    return scaling;
}

void gen_emu_fakes_template(TTree *t_WJets, TTree *t_QCD, TTree *t_WJets_MC, TH2F *h, 
        float m_low = 150., float m_high = 10000.){
    FakeRate el_FR, mu_FR;
    //TH2D *FR;
    setup_new_el_fakerate(&el_FR);
    setup_new_mu_fakerate(&mu_FR);
    //FR.h->Print();
    for (int l=0; l<=2; l++){
        TTree *t;
        if (l==0) t = t_WJets;
        if (l==1) t = t_QCD;
        if (l==2) t = t_WJets_MC;
        Double_t m, xF, cost, jet1_cmva, jet2_cmva, gen_weight;
        Double_t jet1_pt, jet2_pt, jet1_b_weight, jet2_b_weight, pu_SF;

        jet1_b_weight = jet2_b_weight =1.;
        Double_t el_id_SF, el_reco_SF;
        Double_t bcdef_iso_SF, bcdef_id_SF, bcdef_trk_SF;
        Double_t gh_iso_SF, gh_id_SF, gh_trk_SF;
        Double_t evt_fakerate, lep1_fakerate, lep2_fakerate, el1_eta, el1_pt, mu1_eta, mu1_pt;
        TLorentzVector *el = 0;
        TLorentzVector *mu = 0;
        Int_t iso_lep;
        Float_t met_pt;
        Int_t nJets;
        nJets = 2;
        pu_SF=1;
        t->SetBranchAddress("met_pt", &met_pt);
        t->SetBranchAddress("jet2_CMVA", &jet2_cmva);
        t->SetBranchAddress("jet1_CMVA", &jet1_cmva);
        t->SetBranchAddress("jet1_pt", &jet1_pt);
        t->SetBranchAddress("jet2_pt", &jet2_pt);
        //t1->SetBranchAddress("evt_fakerate", &evt_fakerate);
        //t1->SetBranchAddress("el_fakerate", &el1_fakerate);
        t->SetBranchAddress("el1_pt", &el1_pt);
        t->SetBranchAddress("mu1_pt", &mu1_pt);
        t->SetBranchAddress("el1_eta", &el1_eta);
        t->SetBranchAddress("mu1_eta", &mu1_eta);
        t->SetBranchAddress("nJets", &nJets);
        t->SetBranchAddress("el", &el);
        t->SetBranchAddress("mu", &mu);

        if(l==0 || l==2 ){
            t->SetBranchAddress("iso_lep", &iso_lep);
        }
        if(l==2){
            t->SetBranchAddress("el_id_SF", &el_id_SF);
            t->SetBranchAddress("el_reco_SF", &el_reco_SF);
            t->SetBranchAddress("gen_weight", &gen_weight);
            t->SetBranchAddress("pu_SF", &pu_SF);
            t->SetBranchAddress("bcdef_trk_SF", &bcdef_trk_SF);
            t->SetBranchAddress("bcdef_id_SF", &bcdef_id_SF);
            t->SetBranchAddress("gh_trk_SF", &gh_trk_SF);
            t->SetBranchAddress("gh_id_SF", &gh_id_SF);
        }

        Long64_t size  =  t->GetEntries();
        for (int i=0; i<size; i++) {
            t->GetEntry(i);
            bool no_bjets = has_no_bjets(nJets, jet1_pt, jet2_pt, jet1_cmva, jet2_cmva);
            if(l==0){
                //iso_lep: 0 = muons, 1 electrons
                if(iso_lep ==1) lep1_fakerate = get_new_fakerate_prob(mu1_pt, mu1_eta, mu_FR.h);
                if(iso_lep ==0) lep1_fakerate = get_new_fakerate_prob(el1_pt, el1_eta, el_FR.h);
                evt_fakerate = lep1_fakerate/(1-lep1_fakerate);
            }
            if(l==1){
                lep1_fakerate = get_new_fakerate_prob(el1_pt, el1_eta, el_FR.h);
                lep2_fakerate = get_new_fakerate_prob(mu1_pt, mu1_eta, mu_FR.h);
                evt_fakerate = -(lep1_fakerate/(1-lep1_fakerate)) * (lep2_fakerate/(1-lep2_fakerate));
            }
            if(l==2){
                Double_t bcdef_SF = gen_weight * bcdef_trk_SF *  bcdef_id_SF;
                Double_t gh_SF = gen_weight * gh_trk_SF * gh_id_SF;
                Double_t mu_SF = (bcdef_SF *bcdef_lumi + gh_SF * gh_lumi) / 
                    (bcdef_lumi + gh_lumi);

                Double_t mc_weight = gen_weight * el_id_SF * el_reco_SF *pu_SF *
                    mu_SF  * 1000. * mu_lumi;
                if(iso_lep ==1) lep1_fakerate = get_new_fakerate_prob(mu1_pt, mu1_eta, mu_FR.h);
                if(iso_lep ==0) lep1_fakerate = get_new_fakerate_prob(el1_pt, el1_eta, el_FR.h);
                evt_fakerate = -(lep1_fakerate * mc_weight)/(1-lep1_fakerate);
            }

            TLorentzVector cm = *el + *mu;
            cost = get_cost(*el, *mu);
            float pt = cm.Pt();
            float xF = compute_xF(cm); 
            float m = cm.M();
            bool pass = m>= m_low && m <= m_high && met_pt < 50.  && no_bjets && 
                mu1_pt > 27.;
            if(pass){
                h->Fill(xF, cost, 0.5*evt_fakerate);
                h->Fill(xF, -cost, 0.5*evt_fakerate);
            }
        }

        printf("After iter %i current fakerate est is %.0f \n", l, h->Integral());
    }
    cleanup_template(h);
    printf("Total fakerate est is %.0f \n", h->Integral());
}

int gen_combined_background_template(int nTrees, TTree **ts, TH2F* h,  
        Double_t var_low, Double_t var_high, int flag1 = FLAG_MUONS, int flag2 = FLAG_M_BINS,
        bool turn_on_RC = true, bool ss =false, const string &sys_label = ""){
    h->Sumw2();


    for(int i=0; i<nTrees; i++){
        TTree *t1 = ts[i];
        Long64_t nEntries  =  t1->GetEntries();

        Double_t m, xF, cost, gen_weight, reweight, jet1_cmva, jet2_cmva, cost_st;
        Double_t bcdef_HLT_SF, bcdef_iso_SF, bcdef_id_SF;
        Double_t gh_HLT_SF, gh_iso_SF, gh_id_SF;
        Double_t gh_trk_SF, bcdef_trk_SF;
        Double_t el_id_SF, el_reco_SF, pu_SF, el_HLT_SF;
        Double_t jet1_pt, jet2_pt, jet1_b_weight, jet2_b_weight, jet1_eta, jet2_eta;
        Double_t mu1_pt, mu1_eta, mu2_pt, mu2_eta;
        Double_t el1_pt, el1_eta, el2_pt, el2_eta;
        Double_t mu_R_up, mu_R_down, mu_F_up, mu_F_down, 
                 mu_RF_up, mu_RF_down, pdf_up, pdf_down;
        Double_t mu_p_SF, mu_m_SF, mu_p_SF_up, mu_m_SF_up, mu_p_SF_down, mu_m_SF_down, alphaS_up, alphaS_down;
        Float_t cost_pt, met_pt;
        Float_t pdf_weights[60];
        TLorentzVector *lep_p=0;
        TLorentzVector *lep_m=0;
        Double_t pt;
        Int_t nJets, pu_NtrueInt, jet1_flavour, jet2_flavour;
        t1->SetBranchAddress("m", &m);
        t1->SetBranchAddress("xF", &xF);
        t1->SetBranchAddress("cost", &cost);
        t1->SetBranchAddress("jet1_CMVA", &jet1_cmva);
        t1->SetBranchAddress("jet2_CMVA", &jet2_cmva);
        t1->SetBranchAddress("met_pt", &met_pt);
        t1->SetBranchAddress("jet1_pt", &jet1_pt);
        t1->SetBranchAddress("jet2_pt", &jet2_pt);
        t1->SetBranchAddress("nJets", &nJets);
        t1->SetBranchAddress("gen_weight", &gen_weight);
        t1->SetBranchAddress("pu_SF", &pu_SF);
        t1->SetBranchAddress("mu_R_up", &mu_R_up);
        t1->SetBranchAddress("mu_R_down", &mu_R_down);
        t1->SetBranchAddress("mu_F_up", &mu_F_up);
        t1->SetBranchAddress("mu_F_down", &mu_F_down);
        t1->SetBranchAddress("mu_RF_up", &mu_RF_up);
        t1->SetBranchAddress("mu_RF_down", &mu_RF_down);
        t1->SetBranchAddress("pu_NtrueInt", &pu_NtrueInt);
        t1->SetBranchAddress("jet1_eta", &jet1_eta);
        t1->SetBranchAddress("jet2_eta", &jet2_eta);
        t1->SetBranchAddress("jet1_flavour", &jet1_flavour);
        t1->SetBranchAddress("jet2_flavour", &jet2_flavour);
        t1->SetBranchAddress("pdf_weights", &pdf_weights);
        t1->SetBranchAddress("alpha_up", &alphaS_up);
        t1->SetBranchAddress("alpha_down", &alphaS_down);
        int n = 0;

        //SYSTEMATICS 
        Double_t one = 1.0;
        Double_t *systematic = &one;
        int shift = 0;
        if(!sys_label.empty()){
            //doing some systematic
            if(sys_label.find("Up") != string::npos){
                shift = 1;
            }
            else if(sys_label.find("Down") != string::npos){
                shift = -1;
            }

            else{
                printf("systematic label not empty, but doesn't have Up or Down \n");
            }
        }
        // do_sys: 0 = nominal, 1= var up, -1 = var down
        int do_pdf_sys = 0;
        int do_btag_sys = 0;
        int do_pileup_sys = 0;

        int do_muHLT_sys = 0;
        int do_muID_sys = 0;
        int do_muISO_sys = 0;
        int do_muTRK_sys = 0;
        int do_muRC_sys = 0;

        int do_elID_sys = 0;
        int do_elHLT_sys = 0;
        int do_elRECO_sys = 0;
        int do_elScale_sys = 0;
        int do_elSmear_sys = 0;
        float elp_rescale, elm_rescale;

        if(shift !=0){
            if(sys_label.find("BTAG") != string::npos) do_btag_sys = shift;
            else if(sys_label.find("Pu") != string::npos) do_pileup_sys = shift;
            else if(sys_label.find("muHLT") != string::npos) do_muHLT_sys = shift;
            else if(sys_label.find("muID") != string::npos) do_muID_sys = shift;
            else if(sys_label.find("muISO") != string::npos) do_muISO_sys = shift;
            else if(sys_label.find("muTRK") != string::npos) do_muTRK_sys = shift;
            else if(sys_label.find("muRC") != string::npos) do_muRC_sys = shift;

            else if(sys_label.find("elID") != string::npos) do_elID_sys = shift;
            else if(sys_label.find("elHLT") != string::npos) do_elHLT_sys = shift;
            else if(sys_label.find("elRECO") != string::npos) do_elRECO_sys = shift;
            else if(sys_label.find("elScale") != string::npos) do_elScale_sys = shift;
            else if(sys_label.find("elSmear") != string::npos) do_elSmear_sys = shift;


            else if(sys_label.find("RENORM") != string::npos && shift > 0) systematic = &mu_R_up;
            else if(sys_label.find("RENORM") != string::npos && shift < 0) systematic = &mu_R_down;
            else if(sys_label.find("FAC") != string::npos && shift > 0) systematic = &mu_F_up;
            else if(sys_label.find("FAC") != string::npos && shift < 0) systematic = &mu_F_down;
            else if(sys_label.find("alphaS") != string::npos && shift < 0) systematic = &alphaS_down;
            else if(sys_label.find("alphaS") != string::npos && shift > 0) systematic = &alphaS_up;
            else if(sys_label.find("alpha") != string::npos) systematic = &one;
            else if(sys_label.find("pdf") != string::npos){
                if(shift > 0) sscanf(sys_label.c_str(), "_pdf%iUp", &do_pdf_sys);
                else sscanf(sys_label.c_str(), "_pdf%iDown", &do_pdf_sys);
                printf("Doing pdf sys %i \n", do_pdf_sys);
            }

            else printf("COULDN'T PARSE SYSTEMATIC %s !!! \n \n", sys_label.c_str());

        }





        if(flag1 == FLAG_MUONS){

            t1->SetBranchAddress("mu_p", &lep_p);
            t1->SetBranchAddress("mu_m", &lep_m);
            t1->SetBranchAddress("mu1_pt", &mu1_pt);
            t1->SetBranchAddress("mu1_eta", &mu1_eta);
            t1->SetBranchAddress("mu2_pt", &mu2_pt);
            t1->SetBranchAddress("mu2_eta", &mu2_eta);
            t1->SetBranchAddress("bcdef_HLT_SF", &bcdef_HLT_SF);
            t1->SetBranchAddress("bcdef_iso_SF", &bcdef_iso_SF);
            t1->SetBranchAddress("bcdef_id_SF", &bcdef_id_SF);
            t1->SetBranchAddress("bcdef_trk_SF", &bcdef_trk_SF);
            t1->SetBranchAddress("gh_HLT_SF", &gh_HLT_SF);
            t1->SetBranchAddress("gh_iso_SF", &gh_iso_SF);
            t1->SetBranchAddress("gh_id_SF", &gh_id_SF);
            t1->SetBranchAddress("gh_trk_SF", &gh_trk_SF);
            if(turn_on_RC){
                t1->SetBranchAddress("mu_p_SF", &mu_p_SF);
                t1->SetBranchAddress("mu_m_SF", &mu_m_SF);
                t1->SetBranchAddress("mu_p_SF_up", &mu_p_SF_up);
                t1->SetBranchAddress("mu_m_SF_up", &mu_m_SF_up);
                t1->SetBranchAddress("mu_p_SF_down", &mu_p_SF_down);
                t1->SetBranchAddress("mu_m_SF_down", &mu_m_SF_down);
            }

            //separate SFs for runs BCDEF and GH

            for (int i=0; i<nEntries; i++) {
                t1->GetEntry(i);
                if(abs(*systematic) > 10.0 || abs(*systematic) < 0.01 || std::isnan(*systematic)){
                    //printf("sys is %.4f  \n", *systematic);
                    *systematic = 1.;
                }
                bool no_bjets = has_no_bjets(nJets, jet1_pt, jet2_pt, jet1_cmva, jet2_cmva);
                bool not_cosmic = notCosmic(*lep_p, *lep_m);
                if(flag2 == FLAG_PT_BINS){
                    TLorentzVector cm = *lep_p + *lep_m;
                    pt = cm.Pt();
                }
                cost = get_cost(*lep_p, *lep_m);
                if(turn_on_RC){
                    TLorentzVector mu_p_new, mu_m_new, cm;
                    if(do_muRC_sys > 0){
                        mu_p_SF = mu_p_SF_up;
                        mu_m_SF = mu_m_SF_up;
                    }
                    else if(do_muRC_sys < 0){
                        mu_p_SF = mu_p_SF_down;
                        mu_m_SF = mu_m_SF_down;
                    }
                    mu_p_new.SetPtEtaPhiE(lep_p->Pt() * mu_p_SF, lep_p->Eta(), lep_p->Phi(), mu_p_SF * lep_p->E());
                    mu_m_new.SetPtEtaPhiE(lep_m->Pt() * mu_m_SF, lep_m->Eta(), lep_m->Phi(), mu_m_SF * lep_m->E());
                    double new_cost = get_cost(mu_p_new, mu_m_new);
                    cm = mu_p_new + mu_m_new;
                    //if(cm.M() < 150.) continue;
                    xF = compute_xF(cm); 
                    cost = new_cost;
                    m = cm.M();
                    pt = cm.Pt();

                }
                bool pass = ((flag2 == FLAG_M_BINS && m >= var_low && m <= var_high) ||
                        (flag2 == FLAG_PT_BINS && m >= 150. && pt >= var_low && pt <= var_high))
                    && met_pt < 50.  && no_bjets && not_cosmic;

                if(pass){

                    Double_t pu_SF_sys = 1.;
                    if(do_pileup_sys == -1) pu_SF_sys = get_pileup_SF(pu_NtrueInt, pu_sys.pileup_down);
                    if(do_pileup_sys == 1) pu_SF_sys = get_pileup_SF(pu_NtrueInt, pu_sys.pileup_up);
                    if(do_pdf_sys != 0){
                        if(shift > 0) *systematic = pdf_weights[do_pdf_sys-1];
                        if(shift < 0) *systematic = 2. - pdf_weights[do_pdf_sys-1];
                        //printf("mumu sys is  %.6f pdf weight is %.4f \n", *systematic, pdf_weights[do_pdf_sys-1]);
                    }
                    gen_weight *= (*systematic) * pu_SF * pu_SF_sys;

                    if(do_muHLT_sys) bcdef_HLT_SF = get_HLT_SF(mu1_pt, mu1_eta, mu2_pt, mu2_eta, runs_bcdef.HLT_SF, runs_bcdef.HLT_MC_EFF, do_muHLT_sys);
                    if(do_muHLT_sys) gh_HLT_SF = get_HLT_SF(mu1_pt, mu1_eta, mu2_pt, mu2_eta, runs_gh.HLT_SF, runs_gh.HLT_MC_EFF, do_muHLT_sys);

                    if(do_muISO_sys)   bcdef_iso_SF = get_SF(mu1_pt, mu1_eta, runs_bcdef.ISO_SF, do_muISO_sys) * get_SF(mu2_pt, mu2_eta, runs_bcdef.ISO_SF, do_muISO_sys);
                    if(do_muISO_sys)   gh_iso_SF = get_SF(mu1_pt, mu1_eta, runs_gh.ISO_SF, do_muISO_sys) * get_SF(mu2_pt, mu2_eta, runs_gh.ISO_SF, do_muISO_sys);

                    if(do_muID_sys)  bcdef_id_SF = get_SF(mu1_pt, mu1_eta, runs_bcdef.ID_SF, do_muID_sys) * get_SF(mu2_pt, mu2_eta, runs_bcdef.ID_SF, do_muID_sys);
                    if(do_muID_sys)  gh_id_SF = get_SF(mu1_pt, mu1_eta, runs_gh.ID_SF, do_muID_sys) * get_SF(mu2_pt, mu2_eta, runs_gh.ID_SF, do_muID_sys);

                    if(do_muTRK_sys) bcdef_trk_SF = get_Mu_trk_SF(abs(mu1_eta), runs_bcdef.TRK_SF, do_muTRK_sys) * get_Mu_trk_SF(abs(mu2_eta), runs_bcdef.TRK_SF, do_muTRK_sys);
                    if(do_muTRK_sys) gh_trk_SF = get_Mu_trk_SF(abs(mu1_eta), runs_gh.TRK_SF, do_muTRK_sys) * get_Mu_trk_SF(abs(mu2_eta), runs_gh.TRK_SF, do_muTRK_sys);
                    jet1_b_weight = get_btag_weight(jet1_pt, jet1_eta, (Float_t) jet1_flavour , btag_effs, b_reader, do_btag_sys);
                    jet2_b_weight = get_btag_weight(jet2_pt, jet2_eta, (Float_t) jet2_flavour , btag_effs, b_reader, do_btag_sys);


                    Double_t bcdef_weight = gen_weight * bcdef_HLT_SF * bcdef_iso_SF * bcdef_id_SF * bcdef_trk_SF;
                    Double_t gh_weight = gen_weight * gh_HLT_SF * gh_iso_SF * gh_id_SF * gh_trk_SF;
                    if (nJets >= 1){
                        bcdef_weight *= jet1_b_weight;
                        gh_weight *= jet1_b_weight;
                    }
                    if (nJets >= 2){
                        bcdef_weight *= jet2_b_weight;
                        gh_weight *= jet2_b_weight;
                    }


                    Double_t evt_weight = 1000*(bcdef_lumi * bcdef_weight + gh_lumi *gh_weight);
                    if(!ss) h->Fill(xF, cost, evt_weight);
                    else{
                        h->Fill(xF, -abs(cost), evt_weight);
                    }
                }
            }
        }
        else if(flag1 == FLAG_ELECTRONS) {
            t1->SetBranchAddress("el_p", &lep_p);
            t1->SetBranchAddress("el_m", &lep_m);
            t1->SetBranchAddress("el1_pt", &el1_pt);
            t1->SetBranchAddress("el1_eta", &el1_eta);
            t1->SetBranchAddress("el2_pt", &el2_pt);
            t1->SetBranchAddress("el2_eta", &el2_eta);
            t1->SetBranchAddress("el_id_SF", &el_id_SF);
            t1->SetBranchAddress("el_HLT_SF", &el_HLT_SF);
            t1->SetBranchAddress("el_reco_SF", &el_reco_SF);
            if(do_elScale_sys || do_elSmear_sys){
                if(do_elScale_sys >0){
                    t1->SetBranchAddress("elp_scale_up", &elp_rescale);
                    t1->SetBranchAddress("elm_scale_up", &elm_rescale);
                }
                else if(do_elScale_sys < 0){
                    t1->SetBranchAddress("elp_scale_down", &elp_rescale);
                    t1->SetBranchAddress("elm_scale_down", &elm_rescale);
                }
                else if(do_elSmear_sys < 0){
                    t1->SetBranchAddress("elp_smear_down", &elp_rescale);
                    t1->SetBranchAddress("elm_smear_down", &elm_rescale);
                }
                else if(do_elSmear_sys > 0){
                    t1->SetBranchAddress("elp_smear_up", &elp_rescale);
                    t1->SetBranchAddress("elm_smear_up", &elm_rescale);
                }

            }
            for (int i=0; i<nEntries; i++) {
                t1->GetEntry(i);
                if(abs(*systematic) > 10.0 || abs(*systematic) < 0.01 || std::isnan(*systematic)){
                    //printf("sys is %.4f  \n", *systematic);
                    *systematic = 1.;
                }
                bool no_bjets = has_no_bjets(nJets, jet1_pt, jet2_pt, jet1_cmva, jet2_cmva);
                bool not_cosmic = notCosmic(*lep_p, *lep_m);
                if(flag2 == FLAG_PT_BINS){
                    TLorentzVector cm = *lep_p + *lep_m;
                    pt = cm.Pt();
                }
                if(do_elScale_sys || do_elSmear_sys){
                    lep_p->SetPtEtaPhiE(lep_p->Pt() * (elp_rescale), lep_p->Eta(), lep_p->Phi(), lep_p->E() *(elp_rescale));
                    lep_m->SetPtEtaPhiE(lep_m->Pt() * (elm_rescale), lep_m->Eta(), lep_m->Phi(), lep_m->E() *(elm_rescale));
                    TLorentzVector cm = *lep_p + *lep_m;
                    m = cm.M();
                    xF = compute_xF(cm);

                }
                bool pass = ((flag2 == FLAG_M_BINS && m >= var_low && m <= var_high) ||
                        (flag2 == FLAG_PT_BINS && m >= 150. && pt >= var_low && pt <= var_high))
                    && met_pt < 50.  && no_bjets && not_cosmic;
                if(pass){

                    cost = get_cost(*lep_p, *lep_m);

                    Double_t pu_SF_sys = 1.;
                    if(do_pileup_sys == -1) pu_SF_sys = get_pileup_SF(pu_NtrueInt, pu_sys.pileup_down);
                    if(do_pileup_sys == 1) pu_SF_sys = get_pileup_SF(pu_NtrueInt, pu_sys.pileup_up);
                    if(do_pdf_sys != 0){
                        if(shift > 0) *systematic = pdf_weights[do_pdf_sys-1];
                        if(shift < 0) *systematic = 2. - pdf_weights[do_pdf_sys-1];
                        //printf("elel sys is  %.6f pdf weight is %.4f \n", *systematic, pdf_weights[do_pdf_sys-1]);
                    }
                    gen_weight *= (*systematic) * pu_SF * pu_SF_sys;

                    if(do_elID_sys) el_id_SF = get_el_SF(el1_pt, el1_eta, el_SF.ID_SF, do_elID_sys) * get_el_SF(el2_pt, el2_eta, el_SF.ID_SF, do_elID_sys);
                    if(do_elRECO_sys) el_reco_SF = get_el_SF(el1_pt, el1_eta, el_SF.RECO_SF, do_elRECO_sys) * get_el_SF(el2_pt, el2_eta, el_SF.RECO_SF, do_elRECO_sys);
                    if(do_elHLT_sys) el_HLT_SF = get_el_HLT_SF(el1_pt, el1_eta, el2_pt, el2_eta, el_SF.HLT_SF, el_SF.HLT_MC_EFF, do_elHLT_sys);
                    jet1_b_weight = get_btag_weight(jet1_pt, jet1_eta,(Float_t) jet1_flavour , btag_effs, b_reader, do_btag_sys);
                    jet2_b_weight = get_btag_weight(jet2_pt, jet2_eta,(Float_t) jet2_flavour , btag_effs, b_reader, do_btag_sys);

                    Double_t evt_weight = gen_weight * el_id_SF * el_reco_SF * el_HLT_SF * 1000. * el_lumi;
                    if (nJets >= 1){
                        evt_weight *= jet1_b_weight;
                    }
                    if (nJets >= 2){
                        evt_weight *= jet2_b_weight;
                    }
                    if(!ss) h->Fill(xF, cost, evt_weight);
                    else{
                        h->Fill(xF, -abs(cost), evt_weight);
                    }
                }
            }
        }

        t1->ResetBranchAddresses();
    }
    printf("Performing templ. cleanup (removing neg. bins) \n");
    cleanup_template(h);
    printf("Tot Weight is %.2f \n", h->Integral());
    return 0;
}


void gen_emu_template(TTree *t1, TH2F *h, bool is_data = false, 
        float m_low = 150., float m_high = 999999.){
    Long64_t size  =  t1->GetEntries();
    Double_t m, xF, cost, mu1_pt, mu2_pt, jet1_cmva, jet2_cmva, gen_weight;
    Double_t bcdef_HLT_SF, bcdef_iso_SF, bcdef_id_SF;
    Double_t gh_HLT_SF, gh_iso_SF, gh_id_SF, el_id_SF, el_reco_SF, el_HLT_SF;
    Double_t bcdef_trk_SF, gh_trk_SF;
    Double_t jet1_pt, jet2_pt, jet1_b_weight, jet2_b_weight, pu_SF;
    jet1_b_weight = jet2_b_weight =1.;
    Float_t met_pt;
    Int_t nJets;
    nJets = 2;
    pu_SF=1;
    TLorentzVector *el=0;
    TLorentzVector *mu=0;
    int nEvents =0;

    t1->SetBranchAddress("el", &el);
    t1->SetBranchAddress("mu", &mu);
    t1->SetBranchAddress("met_pt", &met_pt);
    t1->SetBranchAddress("jet2_CMVA", &jet2_cmva);
    t1->SetBranchAddress("jet1_CMVA", &jet1_cmva);
    t1->SetBranchAddress("jet1_pt", &jet1_pt);
    t1->SetBranchAddress("jet2_pt", &jet2_pt);
    t1->SetBranchAddress("nJets", &nJets);
    if(!is_data){
        t1->SetBranchAddress("gen_weight", &gen_weight);
        t1->SetBranchAddress("pu_SF", &pu_SF);
        t1->SetBranchAddress("bcdef_HLT_SF", &bcdef_HLT_SF);
        t1->SetBranchAddress("bcdef_iso_SF", &bcdef_iso_SF);
        t1->SetBranchAddress("bcdef_id_SF", &bcdef_id_SF);
        t1->SetBranchAddress("gh_HLT_SF", &gh_HLT_SF);
        t1->SetBranchAddress("gh_iso_SF", &gh_iso_SF);
        t1->SetBranchAddress("gh_id_SF", &gh_id_SF);
        t1->SetBranchAddress("gh_trk_SF", &gh_trk_SF);
        t1->SetBranchAddress("bcdef_trk_SF", &bcdef_trk_SF);
        t1->SetBranchAddress("el_id_SF", &el_id_SF);
        t1->SetBranchAddress("el_reco_SF", &el_reco_SF);     
    }

    for (int i=0; i<size; i++) {
        t1->GetEntry(i);
        bool no_bjets = has_no_bjets(nJets, jet1_pt, jet2_pt, jet1_cmva, jet2_cmva);
        TLorentzVector cm;
        cm = *el + *mu;
        m = cm.M();
        double cost = get_cost(*el, *mu);
        xF  = compute_xF(cm); 

        if(m >= m_low && m <= m_high && met_pt < 50. && no_bjets){
            if(is_data){
                h->Fill(xF, cost, 0.5);
                h->Fill(xF, -cost, 0.5);
            }
            else{
                nEvents++;
                Double_t evt_weight = gen_weight * pu_SF * el_id_SF * el_reco_SF;
                Double_t bcdef_weight = bcdef_iso_SF * bcdef_id_SF * bcdef_trk_SF;
                Double_t gh_weight = gh_iso_SF * gh_id_SF * gh_trk_SF;
                bcdef_weight *= bcdef_HLT_SF;
                gh_weight *= gh_HLT_SF;

                if (nJets >= 1){
                    evt_weight *= jet1_b_weight;
                }
                if (nJets >= 2){
                    evt_weight *= jet2_b_weight;
                }
                //printf(" %.2e %.2e %.2e \n", evt_weight, evt_weight *bcdef_weight, evt_weight *gh_weight);
                double total_weight = 1000 * (evt_weight * gh_weight * gh_lumi + evt_weight * bcdef_weight * bcdef_lumi);
                h->Fill(xF, cost, 0.5 *total_weight);
                h->Fill(xF, -cost, 0.5 *total_weight);
            }



        }
    }
    cleanup_template(h);
}



