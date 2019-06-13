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

int make_gen_cost(TTree *t1, TH1F *h_cost_st, TH1F *h_cost_r, TH1F* h_pt,  TH1F *h_xf, float m_low = 150., float m_high = 100000., bool phot_ind=false){
    //read event data
    h_cost_st->Sumw2();
    Long64_t size  =  t1->GetEntries();
    double norm = 1.0;


    int nSelected=0;
    double root2 = sqrt(2);
    
    Int_t lep1_id, lep2_id, q1_id, q2_id;
    float cost_st, cost_r, gen_weight;
    TLorentzVector *lep_pls = 0;
    TLorentzVector *lep_mns = 0;
    TLorentzVector *q1 = 0;
    TLorentzVector *q2 = 0;
    TLorentzVector cm;
    t1->SetBranchAddress("lep_pls", &lep_pls);
    t1->SetBranchAddress("lep_mns", &lep_mns);
    t1->SetBranchAddress("lep1_id", &lep1_id);
    t1->SetBranchAddress("lep2_id", &lep2_id);
    t1->SetBranchAddress("q1_id", &q1_id);
    t1->SetBranchAddress("q2_id", &q2_id);
    t1->SetBranchAddress("q1", &q1);
    t1->SetBranchAddress("q2", &q2);
    t1->SetBranchAddress("gen_weight", &gen_weight);
    for (int i=0; i<size; i++) {
        t1->GetEntry(i);
        cm = *lep_pls + *lep_mns;
        //printf("M= %.2f \n", cm.M());
        if(cm.M() >= m_low && cm.M() < m_high && (lep1_id < 14) && (lep2_id < 14)
            && (!phot_ind || q1_id == 22 || q2_id ==22)
           //&& ((abs(lep_pls->Eta()) < 2.5) && (abs(lep_mns->Eta()) < 2.5))
                ){
            nSelected++;

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
            cost_r = cost_st = gen_cost;
            if((q1_id > 0 && q1_id < 6) || (q2_id < 0 && q2_id > -6)){
                if(q1->Pz() < 0) cost_st = -gen_cost;
            }
            else{
                if(q2->Pz() < 0) cost_st = -gen_cost;
            }
            if(cm.Pz() < 0) cost_r = -gen_cost;
            double xf = abs(2.*cm.Pz()/13000.);

            //if(gen_weight >0) gen_weight = 1.0;
            //else gen_weight = -1.0;
            h_cost_st->Fill(cost_st, gen_weight *norm);
            h_cost_r->Fill(cost_r, gen_weight*norm);
            h_pt->Fill(cm.Pt(), gen_weight*norm);
            h_xf->Fill(xf, gen_weight*norm);
        }
        



    }
    printf("Selected %i events \n", nSelected);


    t1->ResetBranchAddresses();
    return nSelected;
}
void fit_gen_cost(){
    gStyle->SetOptStat(0);
    TFile *f1= TFile::Open("../generator_stuff/root_files/madgraph_m200_may23.root");
    //TFile *f= TFile::Open("../generator_stuff/root_files/powheg_m150_may6.root");
    TTree *t_gen1 = (TTree *)f1->Get("T_lhe");


    TH1F *h_cost = new TH1F("h_mad_cost", "", 20, -1., 1.);
    TH1F *h_cost_r = new TH1F("h_mad_cost_r", "", 20, -1., 1.);
    TH1F *h_pt = new TH1F("h_pt", "", 20, 0., 300.);
    TH1F *h_xf = new TH1F("h_xf", "", 20, 0., 1.);
    float m_low = 200.;
    float m_high = 400.;

    int nEvents = make_gen_cost(t_gen1,  h_cost, h_cost_r, h_pt, h_xf, m_low, m_high);


    //TF1 *func = new TF1("func", "(1 + x*x + [1]*(1-x*x) + (4./3.)*(2. + [1])*[0]*x) /(8./3. + 4.*[1]/3.)", -1., 1.);
    TF1 *func = new TF1("func", "3./8.*(1 + x*x + ([1]/2.)*(1-3*x*x)) + [0]*x", -1., 1.);
    func->SetParameter(0,0.6);
    func->SetParameter(1,0.1);
    //func->SetParLimits(1, -1.0, 1.0);

    //func->FixParameter(1, 0.098);
    
    Double_t nB = h_cost->Integral(1,10);
    Double_t nF = h_cost->Integral(11,20);
    
    //bin size is 0.1, so 1/bin_size = 10.
    h_cost->Scale(10./h_cost->Integral());
    h_cost->Fit(func);

    Double_t AFB = ((nF - nB))/((nF+nB));
    Double_t dAFB = (1.-AFB*AFB)/sqrt((nEvents));

    printf("Mass range from %.0f to %.0f \n", m_low, m_high);
    printf("AFB: %.4f +/- %.4f \n", func->GetParameter(0), func->GetParError(0));
    printf("A0: %.3f +/- %.3f \n", func->GetParameter(1), func->GetParError(1));
    printf("Counting: NF %.0f NB %.0f \n", nF, nB);
    printf("Counting: AFB %.4f +/- %.4f \n", AFB, dAFB);
}
