#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <iostream>
#include <cstring>
#include <algorithm>
#include "TFile.h"
#include "TH2D.h"
#include "TGraphAsymmErrors.h"





typedef struct {
    TH2D *HLT_SF;
    TH2D *HLT_MC_EFF;
    TH2D *ISO_SF;
    TH2D *ID_SF;
    TGraphAsymmErrors *TRK_SF;
} mu_SFs;


typedef struct {
    TH1D *data_pileup;
    TH1D *pileup_ratio;
} pileup_SFs;

typedef struct {
    TH1D *pileup_up;
    TH1D *pileup_down;
} pileup_systematics;

typedef struct{
    TH2D *ID_SF;
    TH2D *RECO_SF;
    TH2D *HLT_SF;
    TH2D *HLT_MC_EFF;
} el_SFs;

double get_var(Double_t vals[100]){
    float mean(0.), var(0.);
    int n_vars = 100;
    int n_entries = n_vars;
    for(int i =0; i< n_vars; i++){
        //printf("%.2f \n", vals[i]);
        if(std::isnan((float)vals[i])) n_entries--;
        else{
            //printf("val %.2f \n", vals[i]);
            mean += vals[i];
        }
        //printf("%.3e \n", vals[i]);
    }
    mean = mean / n_entries;
    //printf("mean %.3f n_entries %i\n", mean, n_entries);

    for(int  i=0; i< n_vars; i++){
        if(std::isnan((float)vals[i])) continue;
        else var += pow(vals[i] - mean, 2);
    }
    var = var/(n_entries -1);
    //printf("std %.3f \n\n\n", sqrt(var));
    return var;
}

Double_t get_pileup_SF(Int_t n_int, TH1D *h){

    TAxis* x_ax =  h->GetXaxis();
    int xbin = x_ax->FindBin(n_int);

    Double_t result = h->GetBinContent(xbin);
    //if(result < 0.0001) printf("0 pileup SF for %i vertices\n", n_int);
    return result;
}

Double_t get_Mu_trk_SF(Double_t eta, TGraphAsymmErrors *h, int systematic = 0){
    eta = abs(eta);
    Double_t result = h->Eval(eta);
    
    if(systematic !=0){
        TAxis* x_ax =  h->GetXaxis();
        int xbin = x_ax->FindBin(eta);
        Double_t err;
        if(systematic >0) err = h->GetErrorYhigh(xbin);
        else if(systematic <0) err = h->GetErrorYlow(xbin);
        result += err*systematic;
    }

    //if(result < 0.0001) printf("0 pileup SF for %i vertices\n", n_int);
    return result;
}

void get_pdf_avg_std_dev(Float_t pdf_Weights[100], Float_t *pdf_avg, Float_t *pdf_std_dev){
    Float_t sum = 0;
    for (int i=0; i<100; i++){
        sum+= pdf_Weights[i];
    }
    *pdf_avg = sum/100.;
    Float_t var = 0;
    for (int i=0; i<100; i++){
        var += pow(*pdf_avg - pdf_Weights[i], 2);
    }
    *pdf_std_dev = sqrt(var/99.);
    return;
}





Double_t get_SF(Double_t pt, Double_t eta, TH2D *h, int systematic = 0){
    //stay in range of histogram
    if (pt >= 90) pt = 90.;
    if (pt <= 22.5) pt = 22.5;
    eta = abs(eta);
    TAxis* x_ax =  h->GetXaxis();
    TAxis *y_ax =  h->GetYaxis();
    int xbin = x_ax->FindBin(eta);
    int ybin = y_ax->FindBin(pt);

    Double_t result = h->GetBinContent(xbin, ybin);
    if(systematic != 0){
        Double_t err = h->GetBinError(xbin, ybin);
        err = sqrt(err*err + 0.01*0.01);
        //printf("SF is %.3f +/- %.3f \n", result, err);
        result += (systematic * err );
    }
    if(result < 0.001){ 
        printf("0 muon SF for Pt %.1f, Eta %1.2f \n", pt, eta);
        result = 1;
    }
    return result;
}

Double_t get_HLT_SF_1mu(Double_t mu1_pt, Double_t mu1_eta, TH2D *h_SF){
    //get HLT SF for event with just 1 muon
    //stay in range of histogram
    if (mu1_pt >= 350.) mu1_pt = 350.;
    mu1_eta = abs(mu1_eta);
    TAxis *x_ax_SF =  h_SF->GetXaxis();
    TAxis *y_ax_SF =  h_SF->GetYaxis();
    int xbin1_SF = x_ax_SF->FindBin(mu1_pt);
    int ybin1_SF = y_ax_SF->FindBin(std::abs(mu1_eta));


    Double_t SF1 = h_SF->GetBinContent(xbin1_SF, ybin1_SF);

    Double_t result = SF1;
    if(result < 0.01) printf("0 HLT SF for Pt %.1f, Eta %1.2f \n", mu1_pt, mu1_eta);
    if(TMath::IsNaN(result)){ 
        printf("Nan HLT SF 1 mu for Pt %.1f, Eta %1.2f \n", mu1_pt, mu1_eta);
        result = 1;
    }
    //printf("Result, SF1 = (%0.3f, %0.3f) \n", result, SF1);
    return result;
}

Double_t get_HLT_SF_1el(Double_t el1_pt, Double_t el1_eta, TH2D *h_SF){
    //get HLT SF for event with just 1 elon
    //stay in range of histogram
    if (el1_pt >= 350.) el1_pt = 350.;
    el1_eta = abs(el1_eta);
    TAxis *x_ax_SF =  h_SF->GetXaxis();
    TAxis *y_ax_SF =  h_SF->GetYaxis();
    int xbin1_SF = x_ax_SF->FindBin(el1_eta);
    int ybin1_SF = y_ax_SF->FindBin(el1_pt);


    Double_t SF1 = h_SF->GetBinContent(xbin1_SF, ybin1_SF);

    Double_t result = SF1;
    if(result < 0.01) printf("0 HLT SF for Pt %.1f, Eta %1.2f \n", el1_pt, el1_eta);
    if(TMath::IsNaN(result)){ 
        printf("Nan HLT SF for Pt %.1f, Eta %1.2f \n", el1_pt, el1_eta);
        result = 1;
    }
    //printf("Result, SF1 = (%0.3f, %0.3f) \n", result, SF1);
    return result;
}

Double_t get_HLT_SF(Double_t mu1_pt, Double_t mu1_eta, Double_t mu2_pt, Double_t mu2_eta, TH2D *h_SF, TH2D *h_MC_EFF, int systematic = 0){
    //printf("Getting HLT for %.2f %.2f %.2f %.2f \n", mu1_pt, mu1_eta, mu2_pt, mu2_eta);
    //Get HLT SF for event with 2 Muons
    //stay in range of histogram
    if(mu1_pt < 26.) mu1_pt = 26.01;
    if (mu1_pt >= 350.) mu1_pt = 350.;
    if (mu2_pt >= 350.) mu2_pt = 350.;
    mu1_eta = abs(mu1_eta);
    mu2_eta = abs(mu2_eta);
    TAxis *x_ax_SF =  h_SF->GetXaxis();
    TAxis *y_ax_SF =  h_SF->GetYaxis();
    int xbin1_SF = x_ax_SF->FindBin(std::fabs(mu1_eta));
    int ybin1_SF = y_ax_SF->FindBin(mu1_pt);

    int xbin2_SF = x_ax_SF->FindBin(std::fabs(mu2_eta));
    int ybin2_SF = y_ax_SF->FindBin(mu2_pt);

    Double_t SF1 = h_SF->GetBinContent(xbin1_SF, ybin1_SF);
    Double_t SF2 = h_SF->GetBinContent(xbin2_SF, ybin2_SF);
    if(systematic != 0){
        Double_t SF1_err = h_SF->GetBinError(xbin1_SF, ybin1_SF);
        Double_t SF2_err = h_SF->GetBinError(xbin2_SF, ybin2_SF);
        SF1_err = sqrt(SF1_err*SF1_err + 0.005*0.005);
        SF2_err = sqrt(SF2_err*SF2_err + 0.005*0.005);
        //printf("SF is %.3f +/- %.3f \n", SF1, SF1_err);
        //printf("SF is %.3f +/- %.3f \n", SF2, SF2_err);
        SF1 += SF1_err * systematic;
        SF2 += SF2_err * systematic;
    }


    TAxis *x_ax_MC_EFF =  h_MC_EFF->GetXaxis();
    TAxis *y_ax_MC_EFF =  h_MC_EFF->GetYaxis();
    int xbin1_MC_EFF = x_ax_MC_EFF->FindBin(std::fabs(mu1_eta));
    int ybin1_MC_EFF = y_ax_MC_EFF->FindBin(mu1_pt);

    int xbin2_MC_EFF = x_ax_MC_EFF->FindBin(std::fabs(mu2_eta));
    int ybin2_MC_EFF = y_ax_MC_EFF->FindBin(mu2_pt);

    Double_t MC_EFF1 = h_MC_EFF->GetBinContent(xbin1_MC_EFF, ybin1_MC_EFF);
    Double_t MC_EFF2 = h_MC_EFF->GetBinContent(xbin2_MC_EFF, ybin2_MC_EFF);
    Double_t result = (1 - (1-MC_EFF1*SF1)*(1-MC_EFF2*SF2))/
                      (1 - (1-MC_EFF1)*(1-MC_EFF2));
    if(result < 0.01) printf("0 HLT SF for Pt %.1f, Eta %1.2f \n", mu1_pt, mu1_eta);
    if(TMath::IsNaN(result)){ 
        printf("Nan mu HLT SF for Pt1 %.2f, Eta1 %.2f Pt2 %.2f Eta2 %.2f  %.2f %.2f %.2f %.2f \n", mu1_pt, mu1_eta, mu2_pt, mu2_eta, MC_EFF1, MC_EFF2, SF1, SF2);
        result = 1;
    }
    //printf("Result, SF1 = (%0.3f, %0.3f) \n", result, SF1);
    return result;
}

Double_t get_el_HLT_SF(Double_t el1_pt, Double_t el1_eta, Double_t el2_pt, Double_t el2_eta, 
        TH2D *h_SF, TH2D *h_MC_EFF, int systematic = 0){
    //Get HLT SF for event with 2 elons
    //stay in range of histogram
    //restrict to < 200 pt
    if (el1_pt >= 350.) el1_pt = 350.;
    if (el2_pt >= 350.) el2_pt = 350.;
    TAxis *x_ax_SF =  h_SF->GetXaxis();
    TAxis *y_ax_SF =  h_SF->GetYaxis();
    int xbin1_SF = x_ax_SF->FindBin(el1_eta);
    int ybin1_SF = y_ax_SF->FindBin(el1_pt);

    int xbin2_SF = x_ax_SF->FindBin(el2_eta);
    int ybin2_SF = y_ax_SF->FindBin(el2_pt);

    Double_t SF1 = h_SF->GetBinContent(xbin1_SF, ybin1_SF);
    Double_t SF2 = h_SF->GetBinContent(xbin2_SF, ybin2_SF);

    if(systematic != 0){
        Double_t SF1_err = h_SF->GetBinError(xbin1_SF, ybin1_SF);
        Double_t SF2_err = h_SF->GetBinError(xbin2_SF, ybin2_SF);
        //printf("%.3f %.3f \n", SF1_err, SF2_err);
        //SF1_err = min(SF1_err, 0.001);
        //SF2_err = min(SF2_err, 0.001);
        SF1_err = std::max(SF1_err, 0.01);
        SF2_err = std::max(SF2_err, 0.01);
        SF1 += SF1_err * systematic;
        SF2 += SF2_err * systematic;
    }


    TAxis *x_ax_MC_EFF =  h_MC_EFF->GetXaxis();
    TAxis *y_ax_MC_EFF =  h_MC_EFF->GetYaxis();
    int xbin1_MC_EFF = x_ax_MC_EFF->FindBin(el1_eta);
    int ybin1_MC_EFF = y_ax_MC_EFF->FindBin(el1_pt);

    int xbin2_MC_EFF = x_ax_MC_EFF->FindBin(el2_eta);
    int ybin2_MC_EFF = y_ax_MC_EFF->FindBin(el2_pt);

    Double_t MC_EFF1 = h_MC_EFF->GetBinContent(xbin1_MC_EFF, ybin1_MC_EFF);
    Double_t MC_EFF2 = h_MC_EFF->GetBinContent(xbin2_MC_EFF, ybin2_MC_EFF);
    Double_t result = (1 - (1-MC_EFF1*SF1)*(1-MC_EFF2*SF2))/
                      (1 - (1-MC_EFF1)*(1-MC_EFF2));
    if(result < 0.01) printf("0 EL HLT SF for Pt %.1f, Eta %1.2f \n", el1_pt, el1_eta);
    
    if(TMath::IsNaN(result)){ 
        printf("Nan EL HLT SF for Pt %.1f, Eta %1.2f  %.2f %.2f %.2f %.2f \n", el1_pt, el1_eta, MC_EFF1, MC_EFF2, SF1, SF2);
        result = 1;
    }
    //if(abs(result -1.0) > 0.2) printf("Result, SF1 = (%0.3f, %0.3f) \n", result, SF1);
    return result;
}


Double_t get_el_SF(Double_t pt, Double_t eta, TH2D *h, int systematic = 0){
    if( pt <= 25.) pt = 25;
    if (pt >= 150.) pt = 149.;
    TAxis* x_ax =  h->GetXaxis();

    TAxis *y_ax =  h->GetYaxis();
    int xbin = x_ax->FindBin(eta);
    int ybin = y_ax->FindBin(pt);
    Double_t result = h->GetBinContent(xbin, ybin);
    if(systematic != 0){
        Double_t err = h->GetBinError(xbin, ybin);
        err =abs(err);
        //Systematics already included for ELectrons
        result += (systematic * err);
    }

    if(result < 0.01){
        printf("0 el SF for Pt %.1f, Eta %1.2f \n", pt, eta);
        result = 1;
    }
    return result;
}







void setup_SFs(mu_SFs *runs_BCDEF, mu_SFs *runs_GH, pileup_SFs *pu_SF){
    TH1::AddDirectory(kFALSE);
    TFile *f1 = TFile::Open("SFs/2016/Mu_BCDEF_HLT.root");
    f1->cd("IsoMu24_OR_IsoTkMu24_PtEtaBins");
    TDirectory *subdir1 = gDirectory;
    TH2D *HLT_1 = (TH2D *) subdir1->Get("abseta_pt_ratio")->Clone();
    HLT_1->SetDirectory(0);
    runs_BCDEF->HLT_SF = HLT_1;
    subdir1->cd("efficienciesMC");
    TDirectory *subdir12 = gDirectory;
    TH2D *MC_EFF1 = (TH2D *) subdir12->Get("abseta_pt_MC")->Clone();
    MC_EFF1->SetDirectory(0);
    runs_BCDEF->HLT_MC_EFF = MC_EFF1;
    f1->Close();


    TFile *f2 = TFile::Open("SFs/2016/Mu_BCDEF_ID.root");
    TH2D *ID_1 = (TH2D *) f2->Get("NUM_HighPtID_DEN_genTracks_eta_pair_newTuneP_probe_pt")->Clone();
    ID_1->SetDirectory(0);
    runs_BCDEF->ID_SF = ID_1;
    f2->Close();


    TFile *f3 = TFile::Open("SFs/2016/Mu_BCDEF_ISO.root");
    TH2D *ISO_1 = (TH2D *) f3->Get("NUM_LooseRelTkIso_DEN_HighPtIDandIPCut_eta_pair_newTuneP_probe_pt")->Clone();
    ISO_1->SetDirectory(0);
    runs_BCDEF->ISO_SF = ISO_1;
    f3->Close();


    TFile *f4 = TFile::Open("SFs/2016/Mu_GH_HLT.root");
    f4->cd("IsoMu24_OR_IsoTkMu24_PtEtaBins");
    TDirectory *subdir4 = gDirectory;
    TH2D *HLT_2 = (TH2D *) subdir4->Get("abseta_pt_ratio")->Clone();
    HLT_2->SetDirectory(0);
    runs_GH->HLT_SF = HLT_2;
    subdir4->cd("efficienciesMC");
    TDirectory *subdir42 = gDirectory;
    runs_GH->HLT_MC_EFF = (TH2D *) subdir42->Get("abseta_pt_MC")->Clone();
    runs_GH->HLT_MC_EFF ->SetDirectory(0);
    f4->Close();


    TFile *f5 = TFile::Open("SFs/2016/Mu_GH_ID.root");
    TH2D *ID_2 = (TH2D *) f5->Get("NUM_HighPtID_DEN_genTracks_eta_pair_newTuneP_probe_pt")->Clone();
    ID_2->SetDirectory(0);
    runs_GH->ID_SF = ID_2;
    f5->Close();


    TFile *f6 = TFile::Open("SFs/2016/Mu_GH_ISO.root");
    TH2D *ISO_2 = (TH2D *) f6->Get("NUM_LooseRelTkIso_DEN_HighPtIDandIPCut_eta_pair_newTuneP_probe_pt")->Clone();
    ISO_2->SetDirectory(0);
    runs_GH->ISO_SF = ISO_2;
    f6->Close();

    TFile *f6a = TFile::Open("SFs/2016/Muon_tracking_SF_BCDEF.root");
    TDirectory *subdir6a = gDirectory;
    TGraphAsymmErrors *TRK_1 = (TGraphAsymmErrors *) subdir6a->Get("ratio_eff_aeta_dr030e030_corr");
    runs_BCDEF->TRK_SF = TRK_1;
    f6a->Close();

    TFile *f6b = TFile::Open("SFs/2016/Muon_tracking_SF_GH.root");
    TDirectory *subdir6b = gDirectory;
    TGraphAsymmErrors *TRK_2 = (TGraphAsymmErrors *) subdir6b->Get("ratio_eff_aeta_dr030e030_corr");
    runs_GH->TRK_SF = TRK_2;
    f6b->Close();



    TFile *f7 = TFile::Open("SFs/2016/DataPileupHistogram_69200.root");
    TH1D *data_pileup = (TH1D *) f7->Get("pileup")->Clone();
    data_pileup->Scale(1./data_pileup->Integral());
    data_pileup->SetDirectory(0);
    pu_SF->data_pileup = data_pileup;
    pu_SF->pileup_ratio = (TH1D *) data_pileup->Clone("pileup_ratio");
    f7->Close();
}

void setup_el_SF(el_SFs *sf){
    //Setup electron SF's
    TFile *f7 = TFile::Open("SFs/2016/El_ID.root");
    TDirectory *subdir7 = gDirectory;
    TH2D *h1 = (TH2D *) subdir7->Get("EGamma_SF2D")->Clone();
    h1->SetDirectory(0);
    sf->ID_SF = h1;
    f7->Close();
    //el_SF->Print();
    //

    TFile *f8 = TFile::Open("SFs/2016/El_RECO.root");
    TDirectory *subdir8 = gDirectory;
    TH2D *h2 = (TH2D *) subdir8->Get("EGamma_SF2D")->Clone();
    h2->SetDirectory(0);
    sf->RECO_SF = h2;
    f8->Close();

    TFile *f9 = TFile::Open("SFs/2016/El_HLT.root");
    TDirectory *subdir9 = gDirectory;
    TH2D *h_hltsf = (TH2D *) subdir9->Get("EGamma_SF2D")->Clone();
    h_hltsf->SetDirectory(0);
    sf->HLT_SF = h_hltsf;
    TH2D *h_hltmc = (TH2D *) subdir9->Get("EGamma_EffMC2D")->Clone();
    h_hltmc->SetDirectory(0);
    sf->HLT_MC_EFF = h_hltmc;
    f9->Close();
    
}

void setup_pileup_systematic(pileup_systematics *pu_sys){

    TFile *f7 = TFile::Open("SFs/2016/DataPileupHistogram_69200.root");
    TH1D * pileup_data_nom = (TH1D *) f7->Get("pileup")->Clone();
    //pileup_data_nom->Scale(1./pileup_data_nom->Integral());
    pileup_data_nom->SetDirectory(0);

    TFile *f8 = TFile::Open("SFs/2016/DataPileupHistogram_66017.root");
    TH1D * pileup_data_up = (TH1D *) f8->Get("pileup")->Clone();
    //pileup_data_up->Scale(1./pileup_data_up->Integral());
    pileup_data_up->SetDirectory(0);

    TFile *f9 = TFile::Open("SFs/2016/DataPileupHistogram_72383.root");
    TH1D * pileup_data_down = (TH1D *) f9->Get("pileup")->Clone();
    //pileup_data_down->Scale(1./pileup_data_down->Integral());
    pileup_data_down->SetDirectory(0);

    //Printf(" Ints are %.4e %.4e %.4e \n", pileup_data_nom->Integral(), pileup_data_up->Integral(), pileup_data_down->Integral());
    pu_sys->pileup_up = (TH1D *) pileup_data_nom->Clone();
    pu_sys->pileup_down = (TH1D *) pileup_data_nom->Clone();
    pu_sys->pileup_up->Divide(pileup_data_up, pileup_data_nom);
    pu_sys->pileup_down->Divide(pileup_data_down, pileup_data_nom);
    pu_sys->pileup_up->SetDirectory(0);
    pu_sys->pileup_down->SetDirectory(0);
    f7->Close();
    f8->Close();
    f9->Close();
}

