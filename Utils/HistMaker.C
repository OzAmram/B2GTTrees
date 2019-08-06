
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
#include "TH2D.h"
#include "Math/Functor.h"

using namespace std;

#define FLAG_MUONS  0
#define FLAG_ELECTRONS  1
#define FLAG_EMU  2

#define FLAG_M_BINS 0
#define FLAG_PT_BINS 1

Double_t bcdef_lumi16 = 3.119 + 4.035 + 4.270 +  2.615 + 5.809;
Double_t gh_lumi16 =  8.754 + 7.655;
Double_t mu_lumi16 = bcdef_lumi16 + gh_lumi16;
Double_t el_lumi16 = 35.9;

Double_t mu_lumi17 = 42.13;
Double_t el_lumi17 = 41.52;

Double_t mu_lumi18 = 60.43;
Double_t el_lumi18 = 58.83;


void cleanup_hist(TH1 *h){
    int nBins = h->GetNbinsX();
    for(int i=1; i<= nBins; i++){
        float val = h->GetBinContent(i);
        if(val<0.) h->SetBinContent(i,0.);
    }
}

void setHistError(TH1 *h, float e){
    int nBins = h->GetNbinsX();
    for(int i=1; i<= nBins; i++){
        float val = h->GetBinContent(i);
        h->SetBinError(i, val*e);
    }
}

float computeChi2(TH1 *h){
    // only use on ratio plots, with expected value of 1
    float sum = 0.;
    int nBins = h->GetNbinsX();
    for(int i=1; i<= nBins; i++){
        float val = h->GetBinContent(i);
        float err = h->GetBinError(i);
        if (val > 0. && err > 0.){
            sum += std::pow((val-1)/err,2);
        }
    }
    return sum;
}


void symmetrize1d(TH1F *h){
    int n_cost_bins = h->GetNbinsX();

        for(int j=1; j<= n_cost_bins/2; j++){
            float content = h->GetBinContent(j);
            float error = h->GetBinError(j);

            int opp_j = (n_cost_bins + 1) -j;
            float content2 = h->GetBinContent(opp_j);
            float error2 = h->GetBinError(opp_j);

            float new_content = (content + content2)/2.0;
            float new_error = pow((error*error + error2*error2), 0.5)/2.0;
            h->SetBinContent(j, new_content);
            h->SetBinContent(opp_j, new_content);

            h->SetBinError(j, new_error);
            h->SetBinError(opp_j, new_error);


        }
    
}


void set_fakerate_errors(TH2D *h_errs, TH2D *h_fr, TH1F *h){
    float err_sum = 0.;
    for(int i=1; i<= h_errs->GetNbinsX(); i++){
        for(int j=1; j<= h_errs->GetNbinsY(); j++){
            float err = h_fr->GetBinError(i,j);
            float num = h_errs->GetBinContent(i,j);
            err_sum += err*err*num*num;
        }
    }
    //increase error sum by 50% to approximate subtraction of 2 fail template
    err_sum *= 1.50;
    float n_err_events = h_errs->Integral();
    for(int i=1; i<= h->GetNbinsX(); i++){
        float bin_num = h->GetBinContent(i);
        float scaling = pow(bin_num/n_err_events, 2);
        float num_err = pow(h->GetBinError(i),2);
        float weight_err = scaling * err_sum;
        float new_err = sqrt(num_err + weight_err);
        //printf("Err was %.1f, is now %.1f \n", sqrt(num_err), new_err);

        h->SetBinError(i, new_err);
    }
}

void set_fakerate_errors(TH2D *h_errs, TH2D *h_fr, TH2F *h){
    float err_sum = 0.;
    for(int i=1; i<= h_errs->GetNbinsX(); i++){
        for(int j=1; j<= h_errs->GetNbinsY(); j++){
            float err = h_fr->GetBinError(i,j);
            float num = h_errs->GetBinContent(i,j);
            err_sum += err*err*num*num;
        }
    }
    //increase error sum by 50% to approximate subtraction of 2 fail template
    err_sum *= 1.50;
    float n_err_events = h_errs->Integral();
    for(int i=1; i<= h->GetNbinsX(); i++){
        for(int j=1; j<= h->GetNbinsY(); j++){
            float bin_num = h->GetBinContent(i, j);
            float scaling = pow(bin_num/n_err_events, 2);
            float num_err = pow(h->GetBinError(i,j),2);
            float weight_err = scaling * err_sum;
            float new_err = sqrt(num_err + weight_err);
            //printf("Err was %.1f, is now %.1f \n", sqrt(num_err), new_err);

            h->SetBinError(i,j, new_err);
        }
    }
}
    



double get_cost(TLorentzVector lep_p, TLorentzVector lep_m, bool do_reco_flip = true){

    TLorentzVector cm = lep_p + lep_m;
    double root2 = sqrt(2);
    double lep_p_pls = (lep_p.E()+lep_p.Pz())/root2;
    double lep_p_min = (lep_p.E()-lep_p.Pz())/root2;
    double lep_m_pls = (lep_m.E()+lep_m.Pz())/root2;
    double lep_m_min = (lep_m.E()-lep_m.Pz())/root2;
    double qt2 = cm.Px()*cm.Px()+cm.Py()*cm.Py();
    double cm_m2 = cm.M2();
    //cost_p = cos(theta)_r (reconstructed collins soper angle, sign
    //may be 'wrong' if lepton pair direction is not the same as inital
    //quark direction)
    double cost = 2*(lep_m_pls*lep_p_min - lep_m_min*lep_p_pls)/sqrt(cm_m2*(cm_m2 + qt2));
    if(cm.Pz() <0. && do_reco_flip) cost = -cost;
    return cost;
}

bool notCosmic(TLorentzVector v1, TLorentzVector v2){
    Double_t ang = v1.Angle(v2.Vect());
    return ang < (TMath::Pi() - 0.005);
}
Double_t compute_xF(TLorentzVector cm){
    return  abs(2.*cm.Pz()/13000.); 
}

bool goodElEta(Float_t el1_eta){
    el1_eta = abs(el1_eta);
    return (el1_eta < 1.4442) || (el1_eta > 1.566 && el1_eta < 2.5); 
}

double get_cost_v2(TLorentzVector lep_p, TLorentzVector lep_m){

    double Ebeam = 6500.;
    double Pbeam = sqrt(Ebeam*Ebeam - 0.938*0.938);
    TLorentzVector cm = lep_p + lep_m;
    TLorentzVector p1(0., 0., Pbeam, Ebeam);
    TLorentzVector p2(0., 0., -Pbeam, Ebeam);

    if(cm.Pz() < 0. ){
        TLorentzVector p = p1;
        p1 = p2;
        p2 = p;
    }

    TVector3 beta = -cm.BoostVector();
    lep_m.Boost(beta);
    lep_p.Boost(beta);
    p1.Boost(beta);
    p2.Boost(beta);

    // Now calculate the direction of the new z azis

    TVector3 p1u = p1.Vect();
    p1u.SetMag(1.0);
    TVector3 p2u = p2.Vect();
    p2u.SetMag(1.0);
    TVector3 bisec = p1u - p2u;
    bisec.SetMag(1.0);
    TVector3 plep = lep_m.Vect();
    double cost = TMath::Cos(plep.Angle(bisec));
    return cost;
}

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

typedef struct {
    TH2D *h;
} FakeRate;
//
//static type means functions scope is only this file, to avoid conflicts
void setup_new_el_fakerate(FakeRate *FR, int year){
    TFile *f0;
    if (year == 2016) f0 = TFile::Open("../analyze/FakeRate/root_files/2016/SingleElectron_data_fakerate_corrected_june25.root");
    if (year == 2017) f0 = TFile::Open("../analyze/FakeRate/root_files/2016/SingleElectron_data_fakerate_corrected_june25.root");
    if (year == 2018) f0 = TFile::Open("../analyze/FakeRate/root_files/2016/SingleElectron_data_fakerate_corrected_june25.root");
    TH2D *h1 = (TH2D *) gDirectory->Get("h_rate_new")->Clone();
    h1->SetDirectory(0);
    FR->h = h1;
    f0->Close();
}
void setup_new_mu_fakerate(FakeRate *FR, int year){
    TFile *f0;
    if (year == 2016) f0 = TFile::Open("../analyze/FakeRate/root_files/2016/SingleMuon_data_fakerate_corrected_june25.root");
    if (year == 2017) f0 = TFile::Open("../analyze/FakeRate/root_files/2016/SingleMuon_data_fakerate_corrected_june25.root");
    if (year == 2018) f0 = TFile::Open("../analyze/FakeRate/root_files/2016/SingleMuon_data_fakerate_corrected_june25.root");
    TDirectory *subdir = gDirectory;
    TH2D *h1 = (TH2D *) subdir->Get("h_rate_new")->Clone();
    h1->SetDirectory(0);
    FR->h = h1;
    f0->Close();
}


Double_t get_new_fakerate_prob(Double_t pt, Double_t eta, TH2D *h){
    //pt=35;
    if (pt >= 150) pt =150 ;
    if(abs(eta) > 2.4) eta = 2.3;

    Double_t variation = 0.;
    Double_t err;


    TAxis* x_ax =  h->GetXaxis();
    TAxis *y_ax =  h->GetYaxis();

    int xbin = x_ax->FindBin(std::abs(eta));
    int ybin = y_ax->FindBin(pt);

    Double_t prob = h->GetBinContent(xbin, ybin);
    err = h->GetBinError(xbin, ybin);
    //printf("prob err %.2f %.2f \n", prob, err);
    prob += err * variation;
    //printf("prob: %.2f \n", prob);
    if(prob < 0.001 || prob >= 0.99){
        printf("Warning: %.2f Rate for pt %.0f, eta %1.1f! \n"
                "err %.2f varr %.2f \n",
                prob, pt, eta, err, variation);
        ybin -=1;
        prob = h->GetBinContent(xbin, ybin);
        err = h->GetBinError(xbin, ybin);
        prob += err * variation;
        if(prob < 0.001) printf("Tried 1 lower pt bin and still 0 fakerate \n");
    }

    prob = min(prob, 0.98);
    prob = max(prob, 0.02);
    //printf("Efficiency is %f \n", eff);
    return prob;
}
Double_t get_new_fakerate_unc(Double_t pt, Double_t eta, TH2D *h){
    //pt=35;
    if (pt >= 150) pt =150 ;
    if(abs(eta) > 2.4) eta = 2.3;

    Double_t err;


    TAxis* x_ax =  h->GetXaxis();
    TAxis *y_ax =  h->GetYaxis();

    int xbin = x_ax->FindBin(std::abs(eta));
    int ybin = y_ax->FindBin(pt);

    err = h->GetBinError(xbin, ybin);
    return err;
}


void make_emu_m_cost_pt_xf_hist(TTree *t1, TH1F *h_m, TH1F *h_cost,  TH1F *h_pt, TH1F *h_xf, bool is_data = false, 
        int year=2016, float m_low = 150., float m_high = 999999., bool ss = false){
    Long64_t size  =  t1->GetEntries();
    Double_t m, xF, cost, mu1_pt, mu2_pt, jet1_cmva, jet2_cmva, gen_weight;
    Double_t era1_HLT_SF, era1_iso_SF, era1_id_SF;
    Double_t era2_HLT_SF, era2_iso_SF, era2_id_SF, el_id_SF, el_reco_SF, el_HLT_SF;
    Double_t jet1_pt, jet2_pt, pu_SF;
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
        t1->SetBranchAddress("era1_HLT_SF", &era1_HLT_SF);
        t1->SetBranchAddress("era1_iso_SF", &era1_iso_SF);
        t1->SetBranchAddress("era1_id_SF", &era1_id_SF);
        t1->SetBranchAddress("era2_HLT_SF", &era2_HLT_SF);
        t1->SetBranchAddress("era2_iso_SF", &era2_iso_SF);
        t1->SetBranchAddress("era2_id_SF", &era2_id_SF);
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
        double xf = compute_xF(cm); 

        if(m >= m_low && m <= m_high && met_pt < 50. && no_bjets){
            if(is_data){
                h_m->Fill(m);
                if(ss) h_cost->Fill(abs(cost));
                else h_cost->Fill(cost);
                h_xf->Fill(xf);
                h_pt->Fill(cm.Pt());
            }
            else{
                nEvents++;
                Double_t evt_weight = gen_weight * pu_SF * el_id_SF * el_reco_SF;
                Double_t era1_weight = era1_iso_SF * era1_id_SF;
                Double_t era2_weight = era2_iso_SF * era2_id_SF;
                era1_weight *= era1_HLT_SF;
                era2_weight *= era2_HLT_SF;

                //printf(" %.2e %.2e %.2e \n", tot_weight, tot_weight *era1_weight, tot_weight *era2_weight);
                Double_t tot_weight;
                if(year ==2016) tot_weight = 1000 * (evt_weight * era2_weight * gh_lumi16 + evt_weight * era1_weight * bcdef_lumi16);
                if(year ==2017) tot_weight = 1000 * evt_weight * era1_weight * mu_lumi17;
                if(year ==2018) tot_weight = 1000 * evt_weight * era1_weight * mu_lumi18;
                h_m->Fill(m, tot_weight);
                if(ss) h_cost->Fill(abs(cost), tot_weight);
                else h_cost->Fill(cost, tot_weight);
                h_xf->Fill(xf, tot_weight);
                h_pt->Fill(cm.Pt(), tot_weight);
            }



        }
    }
    if(!is_data){
        //Printf("%.1f %.1f", h_m_era1->Integral(), h_m_gh->Integral());
        printf("%i MC events. Poission unc. is %.2f \n", nEvents, 1./(sqrt(nEvents)));
    }
}
void make_qcd_from_emu_m_cost_pt_xf_hist(TTree *t_data, TTree *t_ttbar, TTree *t_diboson, TTree *t_dy, 
        TH1F *h_m, TH1F *h_cost, TH1F *h_pt, TH1F *h_xf,  float m_low = 150., float m_high = 99999.){
    TH1F *data_m = (TH1F*)h_m->Clone("data_m");
    TH1F *ttbar_m = (TH1F*)h_m->Clone("ttbar_m");
    TH1F *diboson_m = (TH1F*)h_m->Clone("diboson_m");
    TH1F *dy_m = (TH1F*)h_m->Clone("dy_m");
    TH1F *data_cost = (TH1F*)h_cost->Clone("data_cost");
    TH1F *ttbar_cost = (TH1F*)h_cost->Clone("ttbar_cost");
    TH1F *diboson_cost = (TH1F*)h_cost->Clone("diboson_cost");
    TH1F *dy_cost = (TH1F*)h_cost->Clone("dy_cost");
    TH1F *data_xf = (TH1F*)h_xf->Clone("data_xf");
    TH1F *ttbar_xf = (TH1F*)h_xf->Clone("ttbar_xf");
    TH1F *diboson_xf = (TH1F*)h_xf->Clone("diboson_xf");
    TH1F *dy_xf = (TH1F*)h_xf->Clone("dy_xf");
    TH1F *data_pt = (TH1F*)h_pt->Clone("data_pt");
    TH1F *ttbar_pt = (TH1F*)h_pt->Clone("ttbar_pt");
    TH1F *diboson_pt = (TH1F*)h_pt->Clone("diboson_pt");
    TH1F *dy_pt = (TH1F*)h_pt->Clone("dy_pt");

    bool ss = true;
    make_emu_m_cost_pt_xf_hist(t_data, data_m, data_cost, data_pt, data_xf,  true, m_low, m_high, ss);
    make_emu_m_cost_pt_xf_hist(t_ttbar, ttbar_m, ttbar_cost,  ttbar_pt,ttbar_xf, false, m_low,  m_high, ss);
    make_emu_m_cost_pt_xf_hist(t_diboson, diboson_m, diboson_cost, diboson_pt, diboson_xf, false,  m_low, m_high, ss);
    make_emu_m_cost_pt_xf_hist(t_dy, dy_m, dy_cost, dy_pt, dy_xf,  false,  m_low, m_high, ss);


    //h_m = data_m - ttbar_m  - dy_m - diboson_m;
    //h_xf = data_xf - ttbar_xf  - dy_xf - diboson_xf;
    //h_cost = data_cost - ttbar_cost  - dy_cost - diboson_cost;
    //
    TH1F *temp1, *temp2;

    temp1 = (TH1F*) data_m->Clone("temp1");
    temp2 = (TH1F*) data_m->Clone("temp2");
    temp1->Add(data_m, ttbar_m ,1, -1);
    temp2->Add(dy_m, diboson_m ,-1, -1);
    h_m->Add(temp1, temp2);

    temp1 = (TH1F*) data_cost->Clone("temp1");
    temp2 = (TH1F*) data_cost->Clone("temp2");
    temp1->Add(data_cost, ttbar_cost ,1, -1);
    temp2->Add(dy_cost, diboson_cost ,-1, -1);
    h_cost->Add(temp1, temp2);

    temp1 = (TH1F*) data_xf->Clone("temp1");
    temp2 = (TH1F*) data_xf->Clone("temp2");
    temp1->Add(data_xf, ttbar_xf ,1, -1);
    temp2->Add(dy_xf, diboson_xf ,-1, -1);
    h_xf->Add(temp1, temp2);

    temp1 = (TH1F*) data_pt->Clone("temp1");
    temp2 = (TH1F*) data_pt->Clone("temp2");
    temp1->Add(data_pt, ttbar_pt ,1, -1);
    temp2->Add(dy_pt, diboson_pt ,-1, -1);
    h_pt->Add(temp1, temp2);
    cleanup_hist(h_m);
    cleanup_hist(h_cost);
    cleanup_hist(h_xf);
    cleanup_hist(h_pt);

}


void make_m_cost_pt_xf_hist(TTree *t1, TH1F *h_m, TH1F *h_cost, TH1F *h_pt, TH1F *h_xf, bool is_data=false, int flag1 = FLAG_MUONS, bool turn_on_RC = true,
        int year = 2016, Double_t m_low = 150., Double_t m_high = 9999999., bool ss = false){
    //read event data
    Long64_t size  =  t1->GetEntries();
    Double_t m, xF, cost, mu1_pt, mu2_pt, jet1_cmva, jet2_cmva, gen_weight;
    Double_t era1_HLT_SF, era1_iso_SF, era1_id_SF;
    Double_t era2_HLT_SF, era2_iso_SF, era2_id_SF, el_id_SF, el_reco_SF, el_HLT_SF;
    Double_t jet1_pt, jet2_pt, pu_SF;
    Double_t mu_p_SF, mu_m_SF, mu_p_SF_alt, mu_m_SF_alt;
    TLorentzVector *mu_p = 0;
    TLorentzVector *mu_m = 0;
    TLorentzVector *el_p = 0;
    TLorentzVector *el_m = 0;
    TLorentzVector cm;
    Float_t met_pt;
    Int_t nJets;
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
    t1->SetBranchAddress("nJets", &nJets);
    if(!is_data){
        t1->SetBranchAddress("gen_weight", &gen_weight);
        t1->SetBranchAddress("pu_SF", &pu_SF);
    }
    if(flag1 == FLAG_MUONS){
        t1->SetBranchAddress("mu_p", &mu_p);
        t1->SetBranchAddress("mu_m", &mu_m);
        t1->SetBranchAddress("mu1_pt", &mu1_pt);
        t1->SetBranchAddress("mu2_pt", &mu2_pt);
        if(!is_data){
            t1->SetBranchAddress("era1_HLT_SF", &era1_HLT_SF);
            t1->SetBranchAddress("era1_iso_SF", &era1_iso_SF);
            t1->SetBranchAddress("era1_id_SF", &era1_id_SF);
            t1->SetBranchAddress("era2_HLT_SF", &era2_HLT_SF);
            t1->SetBranchAddress("era2_iso_SF", &era2_iso_SF);
            t1->SetBranchAddress("era2_id_SF", &era2_id_SF);
        }
        if(turn_on_RC){
            t1->SetBranchAddress("mu_p_SF", &mu_p_SF);
            t1->SetBranchAddress("mu_m_SF", &mu_m_SF);
            t1->SetBranchAddress("mu_p_SF_alt", &mu_p_SF_alt);
            t1->SetBranchAddress("mu_m_SF_alt", &mu_m_SF_alt);
        }
        for (int i=0; i<size; i++) {
            t1->GetEntry(i);
            bool no_bjets = has_no_bjets(nJets, jet1_pt, jet2_pt, jet1_cmva, jet2_cmva);
            bool not_cosmic = notCosmic(*mu_p, *mu_m);

            cm = *mu_p + *mu_m;
            Double_t pt = cm.Pt();
            cost = get_cost(*mu_p, *mu_m);
            if(turn_on_RC){
                TLorentzVector mu_p_new, mu_m_new;
                mu_p_new.SetPtEtaPhiE(mu_p->Pt() * mu_p_SF, mu_p->Eta(), mu_p->Phi(), mu_p_SF * mu_p->E());
                mu_m_new.SetPtEtaPhiE(mu_m->Pt() * mu_m_SF, mu_m->Eta(), mu_m->Phi(), mu_m_SF * mu_m->E());
                //double new_cost = get_cost(mu_p_new, mu_m_new);
                double new_cost = get_cost(mu_p_new, mu_m_new);
                cm = mu_p_new + mu_m_new;
                //if(cm.M() < 150.) continue;
                cost = new_cost;
                m = cm.M();
                pt = cm.Pt();
                xF = compute_xF(cm);
                //printf("old mass %.3f new mass %.3f \n", mu_p->M(), mu_p_new.M());
                   

            }
            if(m >= m_low && m <= m_high && met_pt < 50. && no_bjets && not_cosmic){

                if(is_data){
                    h_m->Fill(m);
                    if(ss){
                        h_cost->Fill(cost, 0.5);
                        h_cost->Fill(-cost, 0.5);
                    }
                    else h_cost->Fill(cost);
                    h_pt->Fill(pt);
                    h_xf->Fill(xF);
                }
                else{

                    Double_t era1_weight = gen_weight *pu_SF * era1_HLT_SF * era1_iso_SF * era1_id_SF;
                    Double_t era2_weight = gen_weight *pu_SF * era2_HLT_SF * era2_iso_SF * era2_id_SF;
                    Double_t tot_weight;
                    if(year == 2016) tot_weight = 1000*(bcdef_lumi16 * era1_weight + gh_lumi16 *era2_weight);
                    if(year == 2017) tot_weight = 1000*(mu_lumi17 * era1_weight);
                    if(year == 2018) tot_weight = 1000*(mu_lumi18 * era1_weight);
                    //Double_t weight = gen_weight;
                    h_m->Fill(m,tot_weight);
                    if(ss){
                        h_cost->Fill(-cost, 0.5* tot_weight);
                        h_cost->Fill(cost, 0.5* tot_weight);
                    }
                    else h_cost->Fill(cost, tot_weight);
                    h_pt->Fill(pt,tot_weight);
                    h_xf->Fill(xF, tot_weight);
                }


            }
        }
    }
    else if(flag1==FLAG_ELECTRONS) {
        t1->SetBranchAddress("el_p", &el_p);
        t1->SetBranchAddress("el_m", &el_m);
        if(!is_data) t1->SetBranchAddress("el_id_SF", &el_id_SF);
        if(!is_data) t1->SetBranchAddress("el_reco_SF", &el_reco_SF);
        if(!is_data) t1->SetBranchAddress("el_HLT_SF", &el_HLT_SF);
        for (int i=0; i<size; i++) {
            t1->GetEntry(i);
            bool no_bjets = has_no_bjets(nJets, jet1_pt, jet2_pt, jet1_cmva, jet2_cmva);
            bool not_cosmic = notCosmic(*el_p, *el_m);
            if(m >= m_low && m <= m_high && met_pt < 50.  && no_bjets && not_cosmic){
                TLorentzVector cm = *el_p + *el_m;
                Double_t pt = cm.Pt();
                cost = get_cost(*el_p, *el_m);
                if(is_data){
                    h_m->Fill(m);
                    if(ss){
                        h_cost->Fill(cost, 0.5);
                        h_cost->Fill(-cost, 0.5);
                    }
                    else h_cost->Fill(cost);
                    h_pt->Fill(pt);
                    h_xf->Fill(xF);
                }
                else{

                    Double_t el_lumi;
                    if(year == 2016) el_lumi = el_lumi16;
                    if(year == 2017) el_lumi = el_lumi17;
                    if(year == 2018) el_lumi = el_lumi18;

                    Double_t tot_weight = gen_weight *pu_SF * el_id_SF * el_reco_SF * el_HLT_SF * 1000. * el_lumi;
                    h_m->Fill(m, tot_weight);
                    if(ss){
                        h_cost->Fill(-cost, 0.5* tot_weight);
                        h_cost->Fill(cost, 0.5* tot_weight);
                    }
                    else h_cost->Fill(cost, tot_weight);
                    h_pt->Fill(pt, tot_weight);
                    h_xf->Fill(xF, tot_weight);
                }
            }
        }
    }


    t1->ResetBranchAddresses();
}


void Fakerate_est_mu(TTree *t_WJets, TTree *t_QCD, TTree *t_WJets_contam, TTree *t_QCD_contam, TH1F *h_m, TH1F *h_cost, TH1F *h_pt, TH1F *h_xf, 
        int year = 2016, float m_low=150., float m_high = 99999., bool ss = false, bool in_os_region=true){
    
    FakeRate FR;
    //TH2D *FR;
    setup_new_mu_fakerate(&FR, year);
    TH2D *h_err = (TH2D *) FR.h->Clone("h_err");
    h_err->Reset();
    float tot_weight_os = 0.;
    float tot_weight_ss = 0.;
    //FR.h->Print();
    for (int l=0; l<=3; l++){
        printf("l=%i\n", l);
        TTree *t;
        if (l==0) t = t_WJets;
        if (l==1) t = t_QCD;
        if (l==2) t = t_WJets_contam;
        if (l==3) t = t_QCD_contam;
        Double_t m, xF, cost, jet1_cmva, jet2_cmva, gen_weight;
        Double_t era1_HLT_SF, era1_iso_SF, era1_id_SF;
        Double_t era2_HLT_SF, era2_iso_SF, era2_id_SF;
        Double_t jet1_pt, jet2_pt,  pu_SF;
        Float_t mu1_charge, mu2_charge;
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
        t->SetBranchAddress("mu1_pt", &mu1_pt);
        t->SetBranchAddress("mu2_pt", &mu2_pt);
        t->SetBranchAddress("mu1_eta", &mu1_eta);
        t->SetBranchAddress("mu2_eta", &mu2_eta);
        t->SetBranchAddress("mu1_charge", &mu1_charge);
        t->SetBranchAddress("mu2_charge", &mu2_charge);
        t->SetBranchAddress("nJets", &nJets);
        t->SetBranchAddress("mu_p", &mu_p);
        t->SetBranchAddress("mu_m", &mu_m);
        if(l==0 || l==2 ){
            t->SetBranchAddress("iso_muon", &iso_mu);
        }
        if(l==2 || l==3){
            t->SetBranchAddress("gen_weight", &gen_weight);
            t->SetBranchAddress("pu_SF", &pu_SF);
            t->SetBranchAddress("era1_HLT_SF", &era1_HLT_SF);
            t->SetBranchAddress("era1_iso_SF", &era1_iso_SF);
            t->SetBranchAddress("era1_id_SF", &era1_id_SF);
            t->SetBranchAddress("era2_HLT_SF", &era2_HLT_SF);
            t->SetBranchAddress("era2_iso_SF", &era2_iso_SF);
            t->SetBranchAddress("era2_id_SF", &era2_id_SF);
        }

        Long64_t size  =  t->GetEntries();
        for (int i=0; i<size; i++) {
            t->GetEntry(i);
            bool no_bjets = has_no_bjets(nJets, jet1_pt, jet2_pt, jet1_cmva, jet2_cmva);
            bool not_cosmic = notCosmic(*mu_p, *mu_m);
            if(l==0){
                if(iso_mu ==1) mu1_fakerate = get_new_fakerate_prob(mu1_pt, mu1_eta, FR.h);
                else           mu1_fakerate = get_new_fakerate_prob(mu2_pt, mu2_eta, FR.h);
                evt_fakerate = mu1_fakerate/(1-mu1_fakerate);
            }
            if(l==1){
                mu1_fakerate = get_new_fakerate_prob(mu1_pt, mu1_eta, FR.h);
                mu2_fakerate = get_new_fakerate_prob(mu2_pt, mu2_eta, FR.h);
                evt_fakerate = -(mu1_fakerate/(1-mu1_fakerate)) * (mu2_fakerate/(1-mu2_fakerate));
            }
            if(l==2){
                Double_t era1_weight = gen_weight * era1_HLT_SF *  era1_id_SF * era1_iso_SF;
                Double_t era2_weight = gen_weight * era2_HLT_SF * era2_id_SF * era2_iso_SF;

                Double_t mc_weight;
                if(year ==2016) mc_weight = 1000 * (era2_weight * gh_lumi16 + era1_weight * bcdef_lumi16);
                if(year ==2017) mc_weight = 1000 * era1_weight * mu_lumi17;
                if(year ==2018) mc_weight = 1000 * era1_weight * mu_lumi18;


                if(iso_mu ==1) mu1_fakerate = get_new_fakerate_prob(mu1_pt, mu1_eta, FR.h);
                else           mu1_fakerate = get_new_fakerate_prob(mu2_pt, mu2_eta, FR.h);
                evt_fakerate = -(mu1_fakerate * mc_weight)/(1-mu1_fakerate);
            }
            if(l==3){
                Double_t era1_weight = gen_weight * era1_HLT_SF *  era1_id_SF;
                Double_t era2_weight = gen_weight * era2_HLT_SF * era2_id_SF;
                Double_t mc_weight;
                if(year ==2016) mc_weight = 1000 * (era2_weight * gh_lumi16 + era1_weight * bcdef_lumi16);
                if(year ==2017) mc_weight = 1000 * era1_weight * mu_lumi17;
                if(year ==2018) mc_weight = 1000 * era1_weight * mu_lumi18;

                mu1_fakerate = get_new_fakerate_prob(mu1_pt, mu1_eta, FR.h);
                mu2_fakerate = get_new_fakerate_prob(mu2_pt, mu2_eta, FR.h);
                evt_fakerate = mc_weight * (mu1_fakerate/(1-mu1_fakerate)) * (mu2_fakerate/(1-mu2_fakerate));
            }


            bool pass = m >= m_low && m <= m_high && met_pt < 50.  && no_bjets && not_cosmic;
            bool opp_sign = ((abs(mu1_charge - mu2_charge)) > 0.01);
            if(!ss) pass = pass && opp_sign;

            if(pass){
                //if(l==3) printf("Evt rate %.2e \n", evt_fakerate);
                TLorentzVector cm = *mu_p + *mu_m;

                if(l==0 && iso_mu ==1) h_err->Fill(min(abs(mu1_eta), 2.3), min(mu1_pt, 150.));
                if(l==0 && iso_mu ==0) h_err->Fill(min(abs(mu2_eta), 2.3), min(mu2_pt, 150.));
                h_m->Fill(m, evt_fakerate);
                if(ss){
                    h_cost->Fill(cost, 0.5*evt_fakerate);
                    h_cost->Fill(-cost, 0.5*evt_fakerate);
                }
                else h_cost->Fill(cost, evt_fakerate);
                h_pt->Fill(cm.Pt(), evt_fakerate);
                h_xf->Fill(compute_xF(cm), evt_fakerate);
                if(opp_sign) tot_weight_os += evt_fakerate;
                else tot_weight_ss += evt_fakerate;
            }
        }
        printf("After iter %i current fakerate est is %.0f \n", l, h_cost->Integral());
    }
    cleanup_hist(h_m);
    cleanup_hist(h_pt);
    cleanup_hist(h_xf);
    cleanup_hist(h_cost);
    set_fakerate_errors(h_err, FR.h, h_cost);
    if(ss){
        float scaling = 1.0;
        if(in_os_region) scaling = tot_weight_os / (tot_weight_ss + tot_weight_os);
        else scaling = tot_weight_ss / (tot_weight_ss + tot_weight_os);
        printf("Scaling is %.2f \n", scaling);
        h_m->Scale(scaling);
        h_pt->Scale(scaling);
        h_xf->Scale(scaling);
        h_cost->Scale(scaling);
    }
    Double_t err;
    printf("Total fakerate est is %.0f +/- %.0f \n", h_cost->IntegralAndError(1, h_cost->GetNbinsX(), err), err);
}

void Fakerate_est_el(TTree *t_WJets, TTree *t_QCD, TTree *t_WJets_MC, TTree *t_QCD_MC, TH1F *h_m, TH1F *h_cost, TH1F *h_pt, TH1F *h_xf,
        int year=2016, float m_low = 150., float m_high = 99999., bool ss = false, bool in_os_region = true){
    FakeRate FR;
    //TH2D *FR;
    setup_new_el_fakerate(&FR, year);
    TH2D *h_err = (TH2D *) FR.h->Clone("h_err");
    h_err->Reset();
    TH1F *h_rw;
    //FR.h->Print();
    float tot_weight_os = 0.;
    float tot_weight_ss = 0.;
    for (int l=0; l<=3; l++){
        TTree *t;
        if (l==0) t = t_WJets;
        if (l==1) t = t_QCD;
        if (l==2) t = t_WJets_MC;
        if (l==3) t = t_QCD_MC;
        Double_t m, xF, cost, jet1_cmva, jet2_cmva, gen_weight;
        Double_t jet1_pt, jet2_pt, pu_SF;
        Double_t el_id_SF, el_reco_SF;
        Double_t evt_fakerate, el1_fakerate, el2_fakerate, el1_eta, el1_pt, el2_eta, el2_pt;
        Float_t el1_charge, el2_charge;
        TLorentzVector *el_p = 0;
        TLorentzVector *el_m = 0;
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
        t->SetBranchAddress("el_p", &el_p);
        t->SetBranchAddress("el_m", &el_m);

        if(l==0 || l==2 ){
            t->SetBranchAddress("iso_el", &iso_el);
        }
        if(l==2 || l==3){
            t->SetBranchAddress("el_id_SF", &el_id_SF);
            t->SetBranchAddress("el_reco_SF", &el_reco_SF);
            t->SetBranchAddress("gen_weight", &gen_weight);
            t->SetBranchAddress("pu_SF", &pu_SF);
        }
        Double_t el_lumi;
        if(year == 2016) el_lumi = el_lumi16;
        if(year == 2017) el_lumi = el_lumi17;
        if(year == 2018) el_lumi = el_lumi18;

        Long64_t size  =  t->GetEntries();
        for (int i=0; i<size; i++) {
            t->GetEntry(i);
            bool no_bjets = has_no_bjets(nJets, jet1_pt, jet2_pt, jet1_cmva, jet2_cmva);
            bool not_cosmic = notCosmic(*el_p, *el_m);
            if(l==0){
                double err;
                if(iso_el ==1) el1_fakerate = get_new_fakerate_prob(el1_pt, el1_eta, FR.h); 
                else           el1_fakerate = get_new_fakerate_prob(el2_pt, el2_eta, FR.h); 
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
                else           el1_fakerate = get_new_fakerate_prob(el2_pt, el2_eta, FR.h);
                evt_fakerate = -(el1_fakerate * mc_weight)/(1-el1_fakerate);
            }
            if(l==3){
                Double_t mc_weight = gen_weight * el_id_SF * el_reco_SF * 1000. * el_lumi;

                el1_fakerate = get_new_fakerate_prob(el1_pt, el1_eta, FR.h);
                el2_fakerate = get_new_fakerate_prob(el2_pt, el2_eta, FR.h);
                evt_fakerate = mc_weight * (el1_fakerate/(1-el1_fakerate)) * (el2_fakerate/(1-el2_fakerate));
            }



            bool pass = m >= m_low && m <= m_high && met_pt < 50.  && no_bjets && not_cosmic;
            bool opp_sign = ((abs(el1_charge - el2_charge)) > 0.01);
            if(!ss) pass = pass && opp_sign;

            if(pass){
                //if(l==3) printf("Evt fr %.2e \n", evt_fakerate);
                //if(l==3) printf("cost, fr %.2f %.2e \n", cost, evt_fakerate);
                TLorentzVector cm = *el_p + *el_m;
                cost = get_cost(*el_p, *el_m);
                float pt = cm.Pt();

                if(l==0 && iso_el ==1) h_err->Fill(min(abs(el1_eta), 2.3), min(el1_pt, 150.));
                if(l==0 && iso_el ==0) h_err->Fill(min(abs(el2_eta), 2.3), min(el2_pt, 150.));
                h_m->Fill(m, evt_fakerate);
                if(ss){
                    h_cost->Fill(cost, 0.5*evt_fakerate);
                    h_cost->Fill(-cost, 0.5*evt_fakerate);
                }
                else h_cost->Fill(cost, evt_fakerate);
                h_pt->Fill(cm.Pt(), evt_fakerate);
                h_xf->Fill(compute_xF(cm), evt_fakerate);
                if(opp_sign) tot_weight_os += evt_fakerate;
                else tot_weight_ss += evt_fakerate;
            }
        }

        printf("After iter %i current fakerate est is %.0f \n", l, h_cost->Integral());
    }
    cleanup_hist(h_m);
    cleanup_hist(h_pt);
    cleanup_hist(h_xf);
    cleanup_hist(h_cost);
    set_fakerate_errors(h_err, FR.h, h_cost);
    if(ss){
        float scaling = 1.0;
        if(in_os_region) scaling = tot_weight_os / (tot_weight_ss + tot_weight_os);
        else scaling = tot_weight_ss / (tot_weight_ss + tot_weight_os);
        printf("Scaling is %.2f \n", scaling);
        h_m->Scale(scaling);
        h_pt->Scale(scaling);
        h_xf->Scale(scaling);
        h_cost->Scale(scaling);
    }
    Double_t err;
    printf("Total fakerate est is %.0f +/- %.0f \n", h_cost->IntegralAndError(1, h_cost->GetNbinsX(), err), err);
}


void Fakerate_est_emu(TTree *t_WJets, TTree *t_QCD, TTree *t_WJets_MC, TH1F *h_m, int flag1 = FLAG_MUONS, 
        int year=2016, float m_low = 150., float m_high = 10000.){
    FakeRate el_FR, mu_FR;
    //TH2D *FR;
    setup_new_el_fakerate(&el_FR, year);
    setup_new_mu_fakerate(&mu_FR, year);
    //FR.h->Print();
    for (int l=0; l<=2; l++){
        TTree *t;
        if (l==0) t = t_WJets;
        if (l==1) t = t_QCD;
        if (l==2) t = t_WJets_MC;
        Double_t m, xF, cost, jet1_cmva, jet2_cmva, gen_weight;
        Double_t jet1_pt, jet2_pt, pu_SF;

        Double_t el_id_SF, el_reco_SF;
        Double_t era1_iso_SF, era1_id_SF;
        Double_t era2_iso_SF, era2_id_SF;
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
            t->SetBranchAddress("era1_id_SF", &era1_id_SF);
            t->SetBranchAddress("era2_id_SF", &era2_id_SF);
        }

        Long64_t size  =  t->GetEntries();
        for (int i=0; i<size; i++) {
            t->GetEntry(i);
            bool no_bjets = has_no_bjets(nJets, jet1_pt, jet2_pt, jet1_cmva, jet2_cmva);
            if(l==0){
                //iso_lep: 0 = muons, 1 electrons
                if(iso_lep ==1) lep1_fakerate = get_new_fakerate_prob(mu1_pt, mu1_eta, mu_FR.h);
                else  lep1_fakerate = get_new_fakerate_prob(el1_pt, el1_eta, el_FR.h);
                evt_fakerate = lep1_fakerate/(1-lep1_fakerate);
            }
            if(l==1){
                lep1_fakerate = get_new_fakerate_prob(el1_pt, el1_eta, el_FR.h);
                lep2_fakerate = get_new_fakerate_prob(mu1_pt, mu1_eta, mu_FR.h);
                evt_fakerate = -(lep1_fakerate/(1-lep1_fakerate)) * (lep2_fakerate/(1-lep2_fakerate));
            }
            if(l==2){
                Double_t era1_SF = era1_id_SF;
                Double_t era2_SF = era2_id_SF;

                Double_t mu_SF, mu_lumi;
                if(year ==2016){
                    mu_SF = (era1_SF *bcdef_lumi16 + era2_SF * gh_lumi16) / 
                    (bcdef_lumi16 + gh_lumi16);
                    mu_lumi=mu_lumi16;
                }
                if(year==2017){
                    mu_SF = era1_SF;
                    mu_lumi=mu_lumi17;
                }
                if(year==2018){
                    mu_SF = era1_SF;
                    mu_lumi=mu_lumi18;
                }

                Double_t mc_weight = gen_weight * el_id_SF * el_reco_SF *pu_SF *
                    mu_SF  * 1000. * mu_lumi;
                if(iso_lep ==1) lep1_fakerate = get_new_fakerate_prob(mu1_pt, mu1_eta, mu_FR.h);
                else            lep1_fakerate = get_new_fakerate_prob(el1_pt, el1_eta, el_FR.h);
                evt_fakerate = -(lep1_fakerate * mc_weight)/(1-lep1_fakerate);
            }

            TLorentzVector cm = *el + *mu;
            float m = cm.M();
            bool pass = m>= m_low && m <= m_high && met_pt < 50.  && no_bjets && 
                ((flag1 == FLAG_MUONS && mu1_pt > 27.) || (flag1 == FLAG_ELECTRONS && el1_pt > 29.));
            if(pass){
                //if(l==3) printf("Evt fr %.2e \n", evt_fakerate);
                h_m->Fill(m, evt_fakerate);
            }
        }

        printf("After iter %i current fakerate est is %.0f \n", l, h_m->Integral());
    }
    cleanup_hist(h_m);
    printf("Total fakerate est is %.0f \n", h_m->Integral());
}
