#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <iostream>
#include <cstring>
#include <algorithm>
#include "TFile.h"

#include "CondFormats/BTauObjects/interface/BTagCalibration.h"
#include "CondTools/BTau/interface/BTagCalibrationReader.h"




typedef struct {
    TH2D *HLT_SF;
    TH2D *ISO_SF;
    TH2D *ID_SF;
} SFs;

typedef struct {
    BTagCalibrationReader b_reader;
    BTagCalibrationReader c_reader;
    BTagCalibrationReader udsg_reader;
} BTag_readers;


Double_t get_SF(Double_t pt, Double_t eta, TH2D *h){
    //stay in range of histogram
    if (pt >= 100) pt = 100.;
    if (pt <= 20) pt = 20.;
    TAxis* x_ax =  h->GetXaxis();
    TAxis *y_ax =  h->GetYaxis();
    int xbin = x_ax->FindBin(pt);
    int ybin = y_ax->FindBin(std::abs(eta));

    Double_t result = h->GetBinContent(xbin, ybin);
    if(result < 0.01) printf("0 SF for Pt %.1f, Eta %1.2f \n", pt, eta);
    return result;
}
Double_t get_HLT_SF(Double_t pt, Double_t eta, TH2D *h){
    //stay in range of histogram
    if (pt >= 350.) pt = 350.;
    TAxis* x_ax =  h->GetXaxis();
    TAxis *y_ax =  h->GetYaxis();
    int xbin = x_ax->FindBin(pt);
    int ybin = y_ax->FindBin(std::abs(eta));

    Double_t result = h->GetBinContent(xbin, ybin);
    if(result < 0.01) printf("0 SF for Pt %.1f, Eta %1.2f \n", pt, eta);
    return result;
}

Double_t get_btag_weight(Double_t pt, Double_t eta, Double_t SF, TH2D *mc_eff){
    //compute weighting from btagging scale factors
    return 1;
}
    


void setup_SFs(SFs *runs_BCDEF, SFs *runs_GH, BTag_readers *btag_r){
    BTagCalibration calib("csvv1", "SFs/cMVAv2_Moriond17_B_H.csv");
    btag_r->b_reader = BTagCalibrationReader(BTagEntry::OP_MEDIUM, "central");
    btag_r->b_reader.load(calib, BTagEntry::FLAV_B, "ttbar");
    btag_r->c_reader = BTagCalibrationReader(BTagEntry::OP_MEDIUM, "central");
    btag_r->c_reader.load(calib, BTagEntry::FLAV_C, "ttbar");
    btag_r->udsg_reader = BTagCalibrationReader(BTagEntry::OP_MEDIUM, "central");
    btag_r->udsg_reader.load(calib, BTagEntry::FLAV_UDSG, "incl");


    TFile *f1 = TFile::Open("SFs/EfficienciesAndSF_RunBtoF.root");
    f1->cd("IsoMu24_OR_IsoTkMu24_PtEtaBins");
    TDirectory *subdir1 = gDirectory;
    TH2D *HLT_1 = (TH2D *) subdir1->Get("pt_abseta_ratio")->Clone();
    HLT_1->SetDirectory(0);
    runs_BCDEF->HLT_SF = HLT_1;
    f1->Close();

    
    TFile *f2 = TFile::Open("SFs/EfficienciesAndSF_BCDEF_ID.root");
    f2->cd("MC_NUM_TightID_DEN_genTracks_PAR_pt_eta");
    TDirectory *subdir2 = gDirectory;
    TH2D *ID_1 = (TH2D *) subdir2->Get("pt_abseta_ratio")->Clone();
    ID_1->SetDirectory(0);
    runs_BCDEF->ID_SF = ID_1;
    f2->Close();


    TFile *f3 = TFile::Open("SFs/EfficienciesAndSF_BCDEF_ISO.root");
    f3->cd("TightISO_TightID_pt_eta");
    TDirectory *subdir3 = gDirectory;
    TH2D *ISO_1 = (TH2D *) subdir3->Get("pt_abseta_ratio")->Clone();
    ISO_1->SetDirectory(0);
    runs_BCDEF->ISO_SF = ISO_1;
    f3->Close();


    TFile *f4 = TFile::Open("SFs/EfficienciesAndSF_Period4.root");
    f4->cd("IsoMu24_OR_IsoTkMu24_PtEtaBins");
    TDirectory *subdir4 = gDirectory;
    TH2D *HLT_2 = (TH2D *) subdir4->Get("pt_abseta_ratio")->Clone();
    HLT_2->SetDirectory(0);
    runs_GH->HLT_SF = HLT_2;
    f4->Close();

    
    TFile *f5 = TFile::Open("SFs/EfficienciesAndSF_GH_ID.root");
    f5->cd("MC_NUM_TightID_DEN_genTracks_PAR_pt_eta");
    TDirectory *subdir5 = gDirectory;
    TH2D *ID_2 = (TH2D *) subdir5->Get("pt_abseta_ratio")->Clone();
    ID_2->SetDirectory(0);
    runs_GH->ID_SF = ID_2;
    f5->Close();


    TFile *f6 = TFile::Open("SFs/EfficienciesAndSF_GH_ISO.root");
    f6->cd("TightISO_TightID_pt_eta");
    TDirectory *subdir6 = gDirectory;
    TH2D *ISO_2 = (TH2D *) subdir6->Get("pt_abseta_ratio")->Clone();
    ISO_2->SetDirectory(0);
    runs_GH->ISO_SF = ISO_2;
    f6->Close();
}

