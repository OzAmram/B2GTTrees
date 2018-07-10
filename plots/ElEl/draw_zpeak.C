

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

const int type = FLAG_ELECTRONS;


void make_m_zpeak_hist(TTree *t1, TH1F *h_m, bool is_data=false, int flag1 = FLAG_MUONS){
    //read event data
    Long64_t size  =  t1->GetEntries();
    Double_t m, xF, cost, mu1_pt, mu2_pt, jet1_cmva, jet2_cmva, gen_weight;
    Double_t bcdef_HLT_SF, bcdef_iso_SF, bcdef_id_SF;
    Double_t gh_HLT_SF, gh_iso_SF, gh_id_SF, el_id_SF, el_reco_SF, el_HLT_SF;
    Double_t bcdef_trk_SF, gh_trk_SF;
    Double_t jet1_pt, jet2_pt, jet1_b_weight, jet2_b_weight, pu_SF;
    jet1_b_weight = jet2_b_weight =1.;
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
    if(!is_data){
        t1->SetBranchAddress("nJets", &nJets);
        t1->SetBranchAddress("gen_weight", &gen_weight);
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
            t1->SetBranchAddress("gh_trk_SF", &gh_trk_SF);
            t1->SetBranchAddress("bcdef_trk_SF", &bcdef_trk_SF);
        }
        for (int i=0; i<size; i++) {
            t1->GetEntry(i);
            bool no_bjets = has_no_bjets(nJets, jet1_pt, jet2_pt, jet1_cmva, jet2_cmva);

            if(m >= 50. && met_pt < 50. && no_bjets){
                cm = *mu_p + *mu_m;
                Double_t pt = cm.Pt();
                if(is_data){
                    h_m->Fill(m);
                }
                else{
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
                    double final_weight = 1000.*(bcdef_weight*bcdef_lumi + gh_weight*gh_lumi);
                    //Double_t weight = gen_weight;
                    h_m->Fill(m,final_weight);
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
            if(m >= 50. && met_pt < 50.  && no_bjets){
                TLorentzVector cm = *el_p + *el_m;
                Double_t pt = cm.Pt();
                if(is_data){
                    h_m->Fill(m);
                }
                else{

                    Double_t evt_weight = gen_weight *pu_SF * el_id_SF * el_reco_SF * el_HLT_SF;
                    if (nJets >= 1){
                        evt_weight *= jet1_b_weight;
                    }
                    if (nJets >= 2){
                        evt_weight *= jet2_b_weight;
                    }
                    h_m->Fill(m, evt_weight);
                }
            }
        }
        if(!is_data){
            Double_t el_lumi = 1000*tot_lumi;
            h_m->Scale(el_lumi);
        }
    }


    t1->ResetBranchAddresses();
}


void Fakerate_est_zpeak_el(TTree *t_WJets, TTree *t_QCD, TTree *t_WJets_MC, TTree *t_QCD_MC, TH1F *h_m){
    FakeRate FR;
    //TH2D *FR;
    setup_new_el_fakerate(&FR);
    //FR.h->Print();
    for (int l=0; l<=3; l++){
        TTree *t;
        if (l==0) t = t_WJets;
        if (l==1) t = t_QCD;
        if (l==2) t = t_WJets_MC;
        if (l==3) t = t_QCD_MC;
        Double_t m, xF, cost, jet1_cmva, jet2_cmva, gen_weight;
        Double_t jet1_pt, jet2_pt, jet1_b_weight, jet2_b_weight, pu_SF;
        Double_t el_id_SF, el_reco_SF;
        Double_t evt_fakerate, el1_fakerate, el2_fakerate, el1_eta, el1_pt, el2_eta, el2_pt;
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



            if(m >= 105. && met_pt < 50.  && no_bjets){
                //if(l==3) printf("Evt fr %.2e \n", evt_fakerate);
                //if(l==3) printf("cost, fr %.2f %.2e \n", cost, evt_fakerate);
                TLorentzVector cm = *el_p + *el_m;
                h_m->Fill(m, evt_fakerate);
            }
        }

        printf("After iter %i current fakerate est is %.0f \n", l, h_m->Integral());
    }
    printf("Total fakerate est is %.0f \n", h_m->Integral());
}
void draw_zpeak(){
    init();
    setTDRStyle();

    f_data = TFile::Open("../analyze/output_files/SingleElectron_data_july6.root");
    t_data = (TTree *)f_data->Get("T_data");

    f_mc = TFile::Open("../analyze/output_files/ElEl_DY_unbinned_june20.root");
    t_mc = (TTree *)f_mc->Get("T_data");
    t_mc_nosig = (TTree *)f_mc->Get("T_back");



    
    int n_bins = 23;
    float bins[] = {50., 55., 60., 64., 68., 72., 76., 81., 86., 91., 96., 101., 106., 110., 115., 120., 126., 133., 141., 150., 160., 171., 185., 200};

    //int n_bins = 50;
    //double start = 50;
    //double end = 2000.;
    TH1F *data_m = new TH1F("data_m", "Data Dimuon Mass Distribution", n_bins, bins);
    TH1F *mc_m = new TH1F("mc_m", "MC Signal (qqbar, qglu, qbarglu)", n_bins, bins);
    TH1F *diboson_m = new TH1F("diboson_m", "DiBoson (WW, WZ, ZZ)", n_bins, bins);
    TH1F *QCD_m = new TH1F("QCD_m", "QCD", n_bins, bins);
    TH1F *WJets_m = new TH1F("WJets_m", "WJets", n_bins, bins);
    TH1F *wt_m = new TH1F("wt_m", "tw + #bar{t}w", n_bins, bins);
    TH1F *mc_nosig_m = new TH1F("mc_nosig_m", "MC no signal (qq, gluglu qbarqbar)", n_bins, bins);
    TH1F *ttbar_m = new TH1F("ttbar_m", "TTBar Background", n_bins, bins);
    mc_m->SetFillColor(kRed+1);
    mc_m->SetMarkerStyle(21);
    mc_m->SetMarkerColor(kRed+1);
    mc_nosig_m->SetFillColor(kMagenta);
    mc_nosig_m->SetMarkerStyle(21);
    mc_nosig_m->SetMarkerColor(kMagenta);
    ttbar_m->SetFillColor(kBlue);
    ttbar_m->SetMarkerStyle(21);
    ttbar_m->SetMarkerColor(kBlue);






    wt_m->SetFillColor(kOrange+7); 
    diboson_m->SetFillColor(kGreen+3);
    QCD_m->SetFillColor(kRed -7);

    make_m_zpeak_hist(t_data, data_m, true, type);
    make_m_zpeak_hist(t_mc, mc_m, false, type);
    make_m_zpeak_hist(t_mc_nosig, mc_nosig_m, false, type);
    make_m_zpeak_hist(t_ttbar, ttbar_m, false, type);
    make_m_zpeak_hist(t_wt, wt_m, false, type);
    make_m_zpeak_hist(t_diboson, diboson_m, false, type);

    //Fakerate_est_zpeak_el(t_WJets, t_QCD, t_WJets_mc, t_QCD_mc, QCD_m);




    Double_t EMu_ratio= 1.05;
    ttbar_m->Scale(EMu_ratio);
    diboson_m->Scale(EMu_ratio);
    wt_m->Scale(EMu_ratio);


    int nBins_x = QCD_m->GetXaxis()->GetNbins();
    int nBins_y = QCD_m->GetYaxis()->GetNbins();
    //printf("Get size %i \n", nBins);
    for (int i=1; i <= nBins_x; i++){
        for (int j=1; j <= nBins_y; j++){

            Double_t m_val = QCD_m->GetBinContent(i,j);

            QCD_m->SetBinError(i,j, 0.2*m_val);
        }
    }

    

    THStack *m_stack = new THStack("m_stack", "MuMu Mass Distribution: Data vs MC ; m_{#mu^{+}#mu^{-}} (GeV)");
    m_stack->Add(diboson_m);
    m_stack->Add(wt_m);
    m_stack->Add(QCD_m);
    m_stack->Add(ttbar_m);
    m_stack->Add(mc_nosig_m);
    m_stack->Add(mc_m);



    TCanvas *c_m = new TCanvas("c_m", "Histograms", 200, 10, 900, 700);
    TPad *pad1 = new TPad("pad1", "pad1", 0.,0.3,0.98,1.);
    pad1->SetBottomMargin(0);
    pad1->Draw();
    pad1->cd();
    pad1->SetLogy();
    m_stack->Draw("hist");
    //m_stack->SetMaximum(65000);
    m_stack->SetMinimum(0.1);
    gStyle->SetEndErrorSize(4);
    data_m->SetMarkerStyle(kFullCircle);
    data_m->SetMarkerColor(1);
    data_m->DrawCopy("P E same");


    gStyle->SetLegendBorderSize(0);
    TLegend *leg1 = new TLegend(0.5, 0.35, 0.75, 0.5);
    leg1->AddEntry(data_m, "data", "p");
    leg1->AddEntry(mc_m, "DY (q#bar{q}, qg #bar{q}g)", "f");
    leg1->AddEntry(mc_nosig_m, "DY no asymmety(gg, qq, #bar{q}#bar{q})", "f");
    leg1->AddEntry(diboson_m, "WW + WZ + ZZ", "f");
    leg1->AddEntry(wt_m, "tW + #bar{t}W", "f");
    leg1->AddEntry(QCD_m, "QCD + WJets", "f");
    leg1->AddEntry(ttbar_m, "t#bar{t}", "f");
    leg1->Draw();

    //gPad->BuildLegend();
    c_m->cd();
    TPad *pad2 = new TPad("pad2", "pad2", 0.,0,.98,0.3);
    //pad2->SetTopMargin(0);
    pad2->SetBottomMargin(0.2);
    pad2->SetGridy();
    pad2->Draw();
    pad2->cd();



    TList *stackHists = m_stack->GetHists();
    TH1* m_mc_sum = (TH1*)stackHists->At(0)->Clone();
    m_mc_sum->Reset();

    for (int i=0;i<stackHists->GetSize();++i) {
      m_mc_sum->Add((TH1*)stackHists->At(i));
    }
    auto ratio = (TH1F *) data_m->Clone("h_ratio");
    ratio->SetMinimum(0.75);
    ratio->SetMaximum(1.25);
    ratio->Sumw2();
    ratio->SetStats(0);
    ratio->Divide(m_mc_sum);
    ratio->SetMarkerStyle(21);
    ratio->Draw("ep");
    TLine *l1 = new TLine(50,1,2000,1);
    l1->SetLineStyle(2);
    l1->Draw();
    c_m->cd();

    ratio->SetTitle("");
    // Y axis ratio plot settings
   ratio->GetYaxis()->SetTitle("Data/MC");
   ratio->GetYaxis()->SetNdivisions(505);
   ratio->GetYaxis()->SetTitleSize(20);
   ratio->GetYaxis()->SetTitleFont(43);
   ratio->GetYaxis()->SetTitleOffset(1.2);
   ratio->GetYaxis()->SetLabelFont(43); // Absolute font size in pixel (precision 3)
   ratio->GetYaxis()->SetLabelSize(15);
   // X axis ratio plot settings
   ratio->GetXaxis()->SetTitle("M_{ee} (GeV)");
   ratio->GetXaxis()->SetTitleSize(20);
   ratio->GetXaxis()->SetTitleFont(43);
   ratio->GetXaxis()->SetTitleOffset(3.);
   ratio->GetXaxis()->SetLabelFont(43); // Absolute font size in pixel (precision 3)
   ratio->GetXaxis()->SetLabelSize(20);
 
    //lumi_sqrtS = "";       // used with iPeriod = 0, e.g. for simulation-only plots (default is an empty string)
    int iPeriod = 4; 
    CMS_lumi(pad1, iPeriod, 33 );
    c_m->Print("ElEl_zpeak.pdf");


 
}

    
    
