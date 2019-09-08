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
} mu_SFs;


typedef struct {
    TH1D *data_pileup;
    TH1D *pileup_ratio;
} pileup_SFs;

typedef struct {
    TH1D *data_pileup_up;
    TH1D *data_pileup_down;
    TH1D *data_pileup_nom;

    TH1D *ratio_pileup_up;
    TH1D *ratio_pileup_down;
    TH1D *ratio_pileup_nom;
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
    if(n_int > 99) n_int = 99;

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





Double_t get_mu_SF(Double_t pt, Double_t eta, int year, TH2D *h, int systematic = 0){
    if(year <2016 || year > 2018) printf("Year is not from 2016-2018. This is bad!! \n");
    //stay in range of histogram
    float sys_unc = 0.;
    if (pt >= 115.){
        sys_unc = 0.01;
        pt = 90.;
    }
    if (pt <= 22.5) pt = 22.5;
    TAxis* x_ax =  h->GetXaxis();
    TAxis *y_ax =  h->GetYaxis();

    int xbin(0),ybin(0);
    if(year == 2016){
        //2016 hists are abs(eta) and pt
        xbin = x_ax->FindBin(eta);
        ybin = y_ax->FindBin(pt);
    }
    else{
        //2017 & 18 hists are pt and eta
        xbin = x_ax->FindBin(pt);
        ybin = y_ax->FindBin(abs(eta));
    }


    Double_t result = h->GetBinContent(xbin, ybin);
    if(systematic != 0){
        Double_t err = h->GetBinError(xbin, ybin);
        err = sqrt(err*err + sys_unc*sys_unc);
        //printf("SF is %.3f +/- %.3f \n", result, err);
        result += (systematic * err );
    }
    if(result < 0.001){ 
        printf("0 muon SF for Pt %.1f, Eta %1.2f \n", pt, eta);
        result = 1;
    }
    return result;
}

Double_t get_el_SF(Double_t pt, Double_t eta, TH2D *h, int systematic = 0){
    float sys_unc = 0.;
    if( pt <= 25.) pt = 25;
    if (pt >= 500.){
        pt = 200.;
        sys_unc = 0.01;
    }
    TAxis* x_ax =  h->GetXaxis();

    TAxis *y_ax =  h->GetYaxis();
    int xbin = x_ax->FindBin(eta);
    int ybin = y_ax->FindBin(pt);
    Double_t result = h->GetBinContent(xbin, ybin);
    if(systematic != 0){
        Double_t err = h->GetBinError(xbin, ybin);
        err =abs(err);
        err = sqrt(err*err + sys_unc*sys_unc);
        result += (systematic * err);
    }

    if(result < 0.01){
        printf("0 el SF for Pt %.1f, Eta %1.2f \n", pt, eta);
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
    int xbin1_SF = x_ax_SF->FindBin(std::abs(mu1_eta));
    int ybin1_SF = y_ax_SF->FindBin(mu1_pt);


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
        //SF1_err = std::max(SF1_err, 0.01);
        //SF2_err = std::max(SF2_err, 0.01);
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




void setup_pu_SFs(pileup_SFs *pu_SF, int year){
    TFile *f7;

    //all use 2016 for now
    if(year == 2016) f7 = TFile::Open("SFs/2016/Data16PileupHistogram_69200.root");
    else if(year == 2017) f7 = TFile::Open("SFs/2017/Data17PileupHistogram_69200.root");
    else if(year == 2018) f7 = TFile::Open("SFs/2018/Data18PileupHistogram_69200.root");

    TH1D *data_pileup = (TH1D *) f7->Get("pileup")->Clone();
    data_pileup->Scale(1./data_pileup->Integral());
    data_pileup->SetDirectory(0);
    pu_SF->data_pileup = data_pileup;
    pu_SF->pileup_ratio = (TH1D *) data_pileup->Clone("pileup_ratio");
    pu_SF->pileup_ratio->SetDirectory(0);
    f7->Close();

    if(pu_SF->data_pileup == NULL) printf("Something wrong getting Pileup SF!\n\n\n");
}




void setup_mu_SFs(mu_SFs *era1, mu_SFs *era2, int year){
    TH1::AddDirectory(kFALSE);

    TFile *f1, *f2, *f3, *f4, *f5, *f6;


    if(year == 2016){
        f1 = TFile::Open("SFs/2016/Mu_BCDEF_HLT.root");
        f1->cd("IsoMu24_OR_IsoTkMu24_PtEtaBins");
        era1->HLT_SF = (TH2D *) gDirectory->Get("abseta_pt_ratio")->Clone();
        era1->HLT_SF->SetDirectory(0);
        gDirectory->cd("efficienciesMC");
        era1->HLT_MC_EFF = (TH2D *) gDirectory->Get("abseta_pt_MC")->Clone();
        era1->HLT_MC_EFF->SetDirectory(0);
        f1->Close();

        f4 = TFile::Open("SFs/2016/Mu_GH_HLT.root");
        f4->cd("IsoMu24_OR_IsoTkMu24_PtEtaBins");
        era2->HLT_SF = (TH2D *) gDirectory->Get("abseta_pt_ratio")->Clone();
        era2->HLT_SF ->SetDirectory(0);
        gDirectory->cd("efficienciesMC");
        era2->HLT_MC_EFF = (TH2D *) gDirectory->Get("abseta_pt_MC")->Clone();
        era2->HLT_MC_EFF ->SetDirectory(0);
        f4->Close();

        f2 = TFile::Open("SFs/2016/Mu_BCDEF_ID.root");
        TH2D *ID_1 = (TH2D *) f2->Get("NUM_TightID_DEN_genTracks_eta_pt")->Clone();
        ID_1->SetDirectory(0);
        era1->ID_SF = ID_1;
        f2->Close();

        f3 = TFile::Open("SFs/2016/Mu_BCDEF_ISO.root");
        TH2D *ISO_1 = (TH2D *) f3->Get("NUM_TightRelIso_DEN_TightIDandIPCut_eta_pt")->Clone();
        ISO_1->SetDirectory(0);
        era1->ISO_SF = ISO_1;
        f3->Close();


        f5 = TFile::Open("SFs/2016/Mu_GH_ID.root");
        TH2D *ID_2 = (TH2D *) f5->Get("NUM_TightID_DEN_genTracks_eta_pt")->Clone();
        ID_2->SetDirectory(0);
        era2->ID_SF = ID_2;
        f5->Close();

        f6 = TFile::Open("SFs/2016/Mu_GH_ISO.root");
        TH2D *ISO_2 = (TH2D *) f6->Get("NUM_TightRelIso_DEN_TightIDandIPCut_eta_pt")->Clone();
        ISO_2->SetDirectory(0);
        era2->ISO_SF = ISO_2;
        f6->Close();

    }
    else if(year == 2017){
        printf("hlt\n");
        f1 = TFile::Open("SFs/2017/Mu_HLT.root");
        f1->cd("IsoMu27_PtEtaBins");
        era1->HLT_SF = (TH2D *) gDirectory->Get("abseta_pt_ratio")->Clone();
        era1->HLT_SF->SetDirectory(0);
        gDirectory->cd("efficienciesMC");
        era1->HLT_MC_EFF = (TH2D *) gDirectory->Get("abseta_pt_MC")->Clone();
        era1->HLT_MC_EFF->SetDirectory(0);
        f1->Close();
        era2->HLT_SF = era1->HLT_SF;
        era2->HLT_MC_EFF = era1->HLT_MC_EFF;

        printf("id\n");
        f2 = TFile::Open("SFs/2017/Mu_ID.root");
        TH2D *ID_1 = (TH2D *) f2->Get("NUM_TightID_DEN_genTracks_pt_abseta")->Clone();
        ID_1->SetDirectory(0);
        era1->ID_SF = ID_1;
        era2->ID_SF = ID_1;
        f2->Close();

        printf("iso\n");
        f3 = TFile::Open("SFs/2017/Mu_ISO.root");
        TH2D *ISO_1 = (TH2D *) f3->Get("NUM_TightRelIso_DEN_TightIDandIPCut_pt_abseta")->Clone();
        ISO_1->SetDirectory(0);
        era1->ISO_SF = ISO_1;
        era2->ISO_SF = ISO_1;
        f3->Close();


    }
    else if(year == 2018){
        f1 = TFile::Open("SFs/2018/Mu_per1_HLT.root");
        f1->cd("IsoMu24_PtEtaBins");
        era1->HLT_SF = (TH2D *) gDirectory->Get("abseta_pt_ratio")->Clone();
        era1->HLT_SF->SetDirectory(0);
        gDirectory->cd("efficienciesMC");
        era1->HLT_MC_EFF = (TH2D *) gDirectory->Get("abseta_pt_MC")->Clone();
        era1->HLT_MC_EFF->SetDirectory(0);
        f1->Close();

        f4 = TFile::Open("SFs/2018/Mu_per2_HLT.root");
        f4->cd("IsoMu24_PtEtaBins");
        era2->HLT_SF = (TH2D *) gDirectory->Get("abseta_pt_ratio")->Clone();
        era2->HLT_SF->SetDirectory(0);
        gDirectory->cd("efficienciesMC");
        era2->HLT_MC_EFF = (TH2D *) gDirectory->Get("abseta_pt_MC")->Clone();
        era2->HLT_MC_EFF->SetDirectory(0);
        f4->Close();

        f2 = TFile::Open("SFs/2018/Mu_ID.root");
        TH2D *ID_1 = (TH2D *) f2->Get("NUM_TightID_DEN_TrackerMuons_pt_abseta")->Clone();
        ID_1->SetDirectory(0);
        era1->ID_SF = ID_1;
        era2->ID_SF = ID_1;
        f2->Close();

        f3 = TFile::Open("SFs/2018/Mu_ISO.root");
        TH2D *ISO_1 = (TH2D *) f3->Get("NUM_TightRelIso_DEN_TightIDandIPCut_pt_abseta")->Clone();
        ISO_1->SetDirectory(0);
        era1->ISO_SF = ISO_1;
        era2->ISO_SF = ISO_1;
        f3->Close();


    }

    if(era1->HLT_SF == NULL || era1->ID_SF == NULL || era1->ISO_SF == NULL  ||
       era2->HLT_SF == NULL || era2->ID_SF == NULL || era2->ISO_SF == NULL) printf("Something wrong setup muon SF'S !!!! \n\n\n");
}

void setup_el_SF(el_SFs *sf, int year){
    //Setup electron SF's
    TFile *f_id, *f_reco, *f_hlt;
    if(year == 2016){
        f_hlt = TFile::Open("SFs/2016/El_HLT.root");
        f_id = TFile::Open("SFs/2016/El_ID.root");
        f_reco = TFile::Open("SFs/2016/El_RECO.root");
    }
    else if(year == 2017){
        f_hlt = TFile::Open("SFs/2016/El_HLT.root");
        f_id = TFile::Open("SFs/2017/El_ID.root");
        f_reco = TFile::Open("SFs/2017/El_RECO.root");
    }
    else if(year == 2018){
        f_hlt = TFile::Open("SFs/2016/El_HLT.root");
        f_id = TFile::Open("SFs/2018/El_ID.root");
        f_reco = TFile::Open("SFs/2018/El_RECO.root");
    }

    TH2D *h1 = (TH2D *) f_id->Get("EGamma_SF2D")->Clone();
    h1->SetDirectory(0);
    sf->ID_SF = h1;
    f_id->Close();

    TH2D *h2 = (TH2D *) f_reco->Get("EGamma_SF2D")->Clone();
    h2->SetDirectory(0);
    sf->RECO_SF = h2;
    f_reco->Close();

    TH2D *h_hltsf = (TH2D *) f_hlt->Get("EGamma_SF2D")->Clone();
    h_hltsf->SetDirectory(0);
    sf->HLT_SF = h_hltsf;
    TH2D *h_hltmc = (TH2D *) f_hlt->Get("EGamma_EffMC2D")->Clone();
    h_hltmc->SetDirectory(0);
    sf->HLT_MC_EFF = h_hltmc;
    f_hlt->Close();
    
    if(sf->HLT_SF == NULL || sf->ID_SF == NULL || sf->RECO_SF == NULL ) printf("Something wrong setup electron SF'S !!!! \n\n\n");
}

void setup_pileup_systematic(pileup_systematics *pu_sys, int year){


    TFile *f7, *f8, *f9;
    if(year == 2016){
        f7 = TFile::Open("SFs/2016/Data16PileupHistogram_69200.root");
        f8 = TFile::Open("SFs/2016/Data16PileupHistogram_66017.root");
        f9 = TFile::Open("SFs/2016/Data16PileupHistogram_72383.root");
    }
    if(year == 2017){
        f7 = TFile::Open("SFs/2017/Data17PileupHistogram_69200.root");
        f8 = TFile::Open("SFs/2017/Data17PileupHistogram_66017.root");
        f9 = TFile::Open("SFs/2017/Data17PileupHistogram_72383.root");
    }
    if(year == 2018){
        f7 = TFile::Open("SFs/2018/Data18PileupHistogram_69200.root");
        f8 = TFile::Open("SFs/2018/Data18PileupHistogram_66017.root");
        f9 = TFile::Open("SFs/2018/Data18PileupHistogram_72383.root");
    }


    TH1D * pileup_data_nom = (TH1D *) f7->Get("pileup")->Clone();
    pileup_data_nom->Scale(1./pileup_data_nom->Integral());
    pileup_data_nom->SetDirectory(0);

    TH1D * pileup_data_up = (TH1D *) f8->Get("pileup")->Clone();
    pileup_data_up->Scale(1./pileup_data_up->Integral());
    pileup_data_up->SetDirectory(0);

    TH1D * pileup_data_down = (TH1D *) f9->Get("pileup")->Clone();
    pileup_data_down->Scale(1./pileup_data_down->Integral());
    pileup_data_down->SetDirectory(0);

    //Printf(" Ints are %.4e %.4e %.4e \n", pileup_data_nom->Integral(), pileup_data_up->Integral(), pileup_data_down->Integral());
    pu_sys->data_pileup_up = (TH1D *) pileup_data_up->Clone("data_pu_up");
    pu_sys->data_pileup_down = (TH1D *) pileup_data_down->Clone("data_pu_down");
    pu_sys->data_pileup_nom = (TH1D *) pileup_data_nom->Clone("data_pu_nom");
    pu_sys->data_pileup_up->SetDirectory(0);
    pu_sys->data_pileup_down->SetDirectory(0);
    pu_sys->data_pileup_nom->SetDirectory(0);

    pu_sys->ratio_pileup_up = (TH1D *) pileup_data_up->Clone("ratio_pu_up");
    pu_sys->ratio_pileup_down = (TH1D *) pileup_data_down->Clone("ratio_pu_down");
    pu_sys->ratio_pileup_nom = (TH1D *) pileup_data_nom->Clone("ratio_pu_nom");
    pu_sys->ratio_pileup_up->SetDirectory(0);
    pu_sys->ratio_pileup_down->SetDirectory(0);
    pu_sys->ratio_pileup_nom->SetDirectory(0);
    f7->Close();
    f8->Close();
    f9->Close();
}

