
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
#include "Math/Functor.h"
#include "../../analyze/HistMaker.C"
#include "../tdrstyle.C"
#include "../CMS_lumi.C"
#include "root_files.h"


typedef struct{
    TH1F *m_hist;
    TH1F *pt_hist;
    TH1F *eta_hist;
    TH1F *phi_hist;
} kin_hists;

void make_4vec_hists(TTree *t1, kin_hists *k, bool is_data=false){
    //read event data
    
    Long64_t size  =  t1->GetEntries();
    Double_t m, xF, cost, mu1_pt, mu2_pt, jet1_cmva, jet2_cmva, gen_weight;
    Double_t bcdef_HLT_SF, bcdef_iso_SF, bcdef_id_SF;
    Double_t gh_HLT_SF, gh_iso_SF, gh_id_SF, el_id_SF, el_reco_SF, el_HLT_SF;
    Double_t bcdef_trk_SF, gh_trk_SF;
    Double_t jet1_pt, jet2_pt, jet1_b_weight, jet2_b_weight, pu_SF;
    TLorentzVector *el_p = 0;
    TLorentzVector *el_m = 0;
    TLorentzVector obs;
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
        t1->SetBranchAddress("jet1_b_weight", &jet1_b_weight);
        t1->SetBranchAddress("jet2_b_weight", &jet2_b_weight);
        t1->SetBranchAddress("pu_SF", &pu_SF);
    }
        t1->SetBranchAddress("el_p", &el_p);
        t1->SetBranchAddress("el_m", &el_m);
        if(!is_data) t1->SetBranchAddress("el_id_SF", &el_id_SF);
        if(!is_data) t1->SetBranchAddress("el_reco_SF", &el_reco_SF);
        if(!is_data) t1->SetBranchAddress("el_HLT_SF", &el_HLT_SF);
        for (int i=0; i<size; i++) {
            t1->GetEntry(i);
            bool no_bjets = has_no_bjets(nJets, jet1_pt, jet2_pt, jet1_cmva, jet2_cmva);

            if(m >= 150. && met_pt < 50. && no_bjets){
                cm = *el_p + *el_m;

                //mu_p->Print();
                //change definitition of obs for different observables


                //if(mu_p->Pt() > mu_m->Pt()) obs = *mu_p;
                //else obs = *mu_m;
                //obs.SetPxPyPzE(mu_p->Px(), mu_p->Py(), mu_p->Pz(), mu_p->E());
                obs=cm;
                

                if(is_data){
                    //k->m_hist->Fill(obs.M());
                    //k->pt_hist->Fill(obs.Pt());
                    //k->eta_hist->Fill(obs.Eta());
                    //k->phi_hist->Fill(obs.Phi());

                    k->m_hist->Fill(el_p->E());
                    k->pt_hist->Fill(el_p->Pt());
                    k->eta_hist->Fill(el_p->Eta());
                    k->phi_hist->Fill(el_p->Phi());

                    k->m_hist->Fill(el_m->E());
                    k->pt_hist->Fill(el_m->Pt());
                    k->eta_hist->Fill(el_m->Eta());
                    k->phi_hist->Fill(el_m->Phi());
                }
                else{
                    Double_t evt_weight = gen_weight *pu_SF * el_id_SF * el_reco_SF * el_HLT_SF * 1000. * tot_lumi;
                    if (nJets >= 1){
                        evt_weight *= jet1_b_weight;
                    }
                    if (nJets >= 2){
                        evt_weight *= jet2_b_weight;
                    }
                    //k->m_hist->Fill(obs.M(), evt_weight);
                    //k->pt_hist->Fill(obs.Pt(), evt_weight);
                    //k->eta_hist->Fill(obs.Eta(), evt_weight);
                    //k->phi_hist->Fill(obs.Phi(), evt_weight);

                    k->m_hist->Fill(el_p->E(), evt_weight);
                    k->pt_hist->Fill(el_p->Pt(), evt_weight);
                    k->eta_hist->Fill(el_p->Eta(), evt_weight);
                    k->phi_hist->Fill(el_p->Phi(), evt_weight);

                    k->m_hist->Fill(el_m->E(), evt_weight);
                    k->pt_hist->Fill(el_m->Pt(), evt_weight);
                    k->eta_hist->Fill(el_m->Eta(), evt_weight);
                    k->phi_hist->Fill(el_m->Phi(), evt_weight);

                }


            }
    }


    t1->ResetBranchAddresses();
}


void setup_kin_hists(kin_hists *k1, char name[20]){

    Double_t m_min = 0.;
    Double_t m_max = 2000.;
    int m_bins = 50;

    Double_t pt_min = 0.;
    Double_t pt_max = 700.;
    int pt_bins = 35;
    
    Double_t eta_min = -3.;
    Double_t eta_max = 3.;
    int eta_bins = 30;

    Double_t phi_min = -4.0;
    Double_t phi_max = 4.0;
    int phi_bins = 30;
    char m_name[40], pt_name[40], eta_name[40], phi_name[40];

    sprintf(m_name, "%s_m_hist", name);
    sprintf(pt_name, "%s_pt_hist", name);
    sprintf(eta_name, "%s_eta_hist", name);
    sprintf(phi_name, "%s_phi_hist", name);

    TH1F *m_hist = new TH1F(m_name, "", m_bins, m_min, m_max);
    TH1F *pt_hist = new TH1F(pt_name, "", pt_bins, pt_min, pt_max);
    TH1F *eta_hist = new TH1F(eta_name, "", eta_bins, eta_min, eta_max);
    TH1F *phi_hist = new TH1F(phi_name, "", phi_bins, phi_min, phi_max);

    m_hist->SetDirectory(0);
    pt_hist->SetDirectory(0);
    eta_hist->SetDirectory(0);
    phi_hist->SetDirectory(0);

    k1->m_hist = m_hist;
    k1->pt_hist = pt_hist;
    k1->eta_hist = eta_hist;
    k1->phi_hist = phi_hist;

    //k1->m_hist->Print();
    return;
}

void setcolors(kin_hists *k, Int_t col){
    k->m_hist->SetFillColor(col);
    k->pt_hist->SetFillColor(col);
    k->eta_hist->SetFillColor(col);
    k->phi_hist->SetFillColor(col);

    k->m_hist->SetLineColor(kBlack);
    k->pt_hist->SetLineColor(kBlack);
    k->eta_hist->SetLineColor(kBlack);
    k->phi_hist->SetLineColor(kBlack);

    k->m_hist->SetMarkerColor(col);
    k->pt_hist->SetMarkerColor(col);
    k->eta_hist->SetMarkerColor(col);
    k->phi_hist->SetMarkerColor(col);

    k->m_hist->SetMarkerStyle(21);
    k->pt_hist->SetMarkerStyle(21);
    k->eta_hist->SetMarkerStyle(21);
    k->phi_hist->SetMarkerStyle(21);
}

void do_emu_scaling(kin_hists *k){
    Double_t emu_scaling = 0.97;
    k->m_hist->Scale(emu_scaling);
    k->pt_hist->Scale(emu_scaling);
    k->eta_hist->Scale(emu_scaling);
    k->phi_hist->Scale(emu_scaling);
}


void make_plots(char name[80], kin_hists *k_data, kin_hists *k_mc, kin_hists *k_mc_nosig, 
                                   kin_hists *k_ttbar, kin_hists *k_diboson, 
                                   kin_hists *k_wt){

    setcolors(k_data, kBlack);
    setcolors(k_mc, kRed+1);
    setcolors(k_mc_nosig, kMagenta);
    setcolors(k_ttbar, kBlue);
    setcolors(k_diboson, kGreen+3);
    setcolors(k_wt, kOrange+7);
    char m_label[80], pt_label[80], eta_label[80], phi_label[80];
    sprintf(m_label, "%s Energy (Gev)", name);
    sprintf(pt_label, "%s Pt (Gev)", name);
    sprintf(eta_label, "%s #eta", name);
    sprintf(phi_label, "%s #phi", name);
    THStack *m_stack = new THStack("m_stack", "");
    THStack *pt_stack = new THStack("pt_stack", "");
    THStack *eta_stack = new THStack("eta_stack", "");
    THStack *phi_stack = new THStack("phi_stack", "");

    char *labels[] = {m_label, pt_label, eta_label, phi_label};

    kin_hists *kins[] = {k_mc_nosig, k_wt, k_diboson, k_ttbar, k_mc};

    for(int i=0; i<5; i++){
        kins[i]->m_hist->Print();
        m_stack->Add(kins[i]->m_hist);
        pt_stack->Add(kins[i]->pt_hist);
        eta_stack->Add(kins[i]->eta_hist);
        phi_stack->Add(kins[i]->phi_hist);
    }

    printf("made stack \n");

    TLegend *leg1 = new TLegend(0.15, 0.75, 0.35, 0.9);
    leg1->AddEntry(k_data->m_hist, "data", "p");
    leg1->AddEntry(k_mc->m_hist, "DY (q#bar{q}, qg #bar{q}g)", "f");
    leg1->AddEntry(k_ttbar->m_hist, "t#bar{t}", "f");
    leg1->AddEntry(k_diboson->m_hist, "WW + WZ + ZZ", "f");
    leg1->AddEntry(k_wt->m_hist, "tW + #bar{t}W", "f");
    leg1->AddEntry(k_mc_nosig->m_hist, "DY no asymmety(gg, qq, #bar{q}#bar{q})", "f");

    THStack *h_stacks[4] = {m_stack, pt_stack, eta_stack, phi_stack};
    TH1F *data_hists[4] = {k_data->m_hist, k_data->pt_hist, k_data->eta_hist, k_data->phi_hist};
    
    //setup canvas for output
    TCanvas *c1 = new TCanvas("c", "", 800, 800);
    TPad *pad1 = new TPad("pad1", "pad1", 0,0.3,1,1.0);
    gStyle->SetEndErrorSize(4);
    pad1->SetBottomMargin(0);
    pad1->Draw();
    pad1->cd();
    pad1->SetLogy();
    gStyle->SetLegendBorderSize(0);
    c1->cd();
    TPad *pad2 = new TPad("pad2", "pad2", 0,0.05,1,0.3);
    //pad2->SetTopMargin(0);
    pad2->SetBottomMargin(0.2);
    pad2->Draw();
    pad2->SetGridy();

    char outfile[4][80];
    sprintf(outfile[0], "%s_M.pdf", name);
    sprintf(outfile[1], "%s_pt.pdf", name);
    sprintf(outfile[2], "%s_eta.pdf", name);
    sprintf(outfile[3], "%s_phi.pdf", name);
    for(int i=0; i<4; i++){
        printf("loop %i \n", i);
        pad1->cd();
        data_hists[i]->SetMarkerStyle(kFullCircle);
        data_hists[i]->SetMarkerColor(1);
        h_stacks[i]->SetMaximum(10*data_hists[i]->GetMaximum());
        h_stacks[i]->SetMinimum(1);
        h_stacks[i]->Draw("hist");
        data_hists[i]->DrawCopy("P E same");

        leg1->Draw();

        TList *stackHists = h_stacks[i]->GetHists();
        TH1* stack_sum = (TH1*)stackHists->At(0)->Clone();
        stack_sum->Reset();

        for (int i=0;i<stackHists->GetSize();++i) {
          stack_sum->Add((TH1*)stackHists->At(i));
        }
        pad2->cd();
        auto ratio = (TH1F *) data_hists[i]->Clone("h_ratio");
        ratio->SetMinimum(0.6);
        ratio->SetMaximum(1.4);
        ratio->Sumw2();
        ratio->SetStats(0);
        ratio->Divide(stack_sum);
        ratio->SetMarkerStyle(21);
        ratio->Draw("ep");
        //TLine *l1 = new TLine(150,1,2000,1);
        //l1->SetLineStyle(2);
        //l1->Draw();

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
       ratio->GetXaxis()->SetTitle(labels[i]);
       ratio->GetXaxis()->SetTitleSize(20);
       ratio->GetXaxis()->SetTitleFont(43);
       ratio->GetXaxis()->SetTitleOffset(3.);
       ratio->GetXaxis()->SetLabelFont(43); // Absolute font size in pixel (precision 3)
       ratio->GetXaxis()->SetLabelSize(20);
     
        //lumi_sqrtS = "";       // used with iPeriod = 0, e.g. for simulation-only plots (default is an empty string)
        int iPeriod = 4; 
        CMS_lumi(pad1, iPeriod, 33 );

        printf("outputing \n");

        c1->Print(outfile[i]);
        //c1->Print("out1.pdf");

    }

    return;
}
void draw_kinematics(){
    init(); 
    kin_hists k_data, k_mc, k_mc_nosig, k_ttbar, k_diboson, k_wt;

    setup_kin_hists(&k_data, "data");
    setup_kin_hists(&k_mc, "mc");
    setup_kin_hists(&k_mc_nosig, "mc_nosig");
    setup_kin_hists(&k_ttbar, "ttbar");
    setup_kin_hists(&k_diboson, "diboson");
    setup_kin_hists(&k_wt, "wt");



    make_4vec_hists(t_data, &k_data, true);
    make_4vec_hists(t_mc, &k_mc);
    make_4vec_hists(t_mc_nosig, &k_mc_nosig);
    make_4vec_hists(t_ttbar, &k_ttbar);
    make_4vec_hists(t_diboson, &k_diboson);
    make_4vec_hists(t_wt, &k_wt);

    printf("doing emu scaling \n");

    do_emu_scaling(&k_ttbar);
    do_emu_scaling(&k_diboson);
    do_emu_scaling(&k_wt);

    printf("making plots \n");
    make_plots("El", &k_data, &k_mc, &k_mc_nosig, &k_ttbar, &k_diboson, &k_wt);

    return;
}
