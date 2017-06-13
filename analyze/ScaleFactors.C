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

typedef struct {
    TH2D *b_eff;
    TH2D *c_eff;
    TH2D *udsg_eff;
} BTag_effs;

typedef struct{
    TH2D *h;
} el_SFs;


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
Double_t get_el_SF(Double_t pt, Double_t eta, TH2D *h){
    if (pt >= 150.) pt = 148.;
    TAxis* x_ax =  h->GetXaxis();

    TAxis *y_ax =  h->GetYaxis();
    int xbin = x_ax->FindBin(eta);
    int ybin = y_ax->FindBin(pt);

    Double_t result = h->GetBinContent(xbin, ybin);
    if(result < 0.01) printf("0 SF for Pt %.1f, Eta %1.2f \n", pt, eta);
    return result;
}

Double_t btag_weight_helper(Double_t pt, Double_t eta, Double_t SF, TH2D *mc_eff){
    if (pt >=500) pt = 450;

    TAxis* x_ax =  mc_eff->GetXaxis();
    TAxis *y_ax =  mc_eff->GetYaxis();
    int xbin = x_ax->FindBin(pt);
    int ybin = y_ax->FindBin(std::abs(eta));

    Double_t eff = mc_eff->GetBinContent(xbin, ybin);
    if(eff == 0) printf("Warning: 0 efficiency for pt %.0f, eta %1.1f \n!", pt, eta);
    //printf("Efficiency is %f \n", eff);
    Double_t weight = (1-SF*eff)/(1-eff);
    return weight;
}


Double_t btag_eff(Double_t pt, Double_t eta,TH2D *mc_eff){
    if (pt >=500) pt = 450;

    TAxis* x_ax =  mc_eff->GetXaxis();
    TAxis *y_ax =  mc_eff->GetYaxis();
    int xbin = x_ax->FindBin(pt);
    int ybin = y_ax->FindBin(std::abs(eta));

    Double_t eff = mc_eff->GetBinContent(xbin, ybin);
    if(eff == 0) printf("Warning: 0 efficiency for pt %.0f, eta %1.1f \n!", pt, eta);
    //printf("Efficiency is %f \n", eff);
    return eff;
}

Double_t get_btag_weight(Double_t pt, Double_t eta, Float_t flavour, BTag_effs btag_effs, BTag_readers b_readers){
    //compute weighting from btagging scale factors

    Double_t weight, bjet_SF;

    if(std::abs(flavour - 5.) < 0.01){ //bjet
        bjet_SF = b_readers.b_reader.eval_auto_bounds("central", BTagEntry::FLAV_B, eta, pt);
        weight = btag_weight_helper(pt, eta, bjet_SF, btag_effs.b_eff);
        //printf("B jet, SF is %0.3f, weight is %.4f \n", bjet_SF, weight);

    }
    else if (std::abs(flavour - 4.) < 0.01){ //cjet
        bjet_SF = b_readers.c_reader.eval_auto_bounds("central", BTagEntry::FLAV_C, eta, pt);
        weight = btag_weight_helper(pt, eta, bjet_SF, btag_effs.c_eff);
        //printf("C jet, SF is %0.3f, weight is %.4f \n", bjet_SF, weight);
    }
    else{ //udsg jet
        bjet_SF = b_readers.udsg_reader.eval_auto_bounds("central", BTagEntry::FLAV_UDSG, eta,pt);
        weight = btag_weight_helper(pt, eta, bjet_SF, btag_effs.udsg_eff);
        //printf("UDSG jet, SF is %0.3f, weight is %.4f \n", bjet_SF, weight);
    }
    if(bjet_SF == 0) printf("WARNING: Scale factor return 0 for Flavour %1.0f pt %.0f eta %1.1f \n!",
                            flavour, pt, eta);
    return weight;
}


Double_t get_emu_btag_weight(Double_t pt1, Double_t eta1, Float_t flavour1, Double_t pt2, Double_t eta2, Float_t flavour2, BTag_effs btag_effs, BTag_readers b_readers){
    //compute weighting from btagging scale factors

    Double_t bjet1_SF, bjet1_eff, bjet2_SF, bjet2_eff;

    if(std::abs(flavour1 - 5.) < 0.01){ //bjet
        bjet1_SF = b_readers.b_reader.eval_auto_bounds("central", BTagEntry::FLAV_B, eta1, pt1);
        bjet1_eff= btag_eff(pt1, eta1, btag_effs.b_eff);
        //printf("B jet, SF is %0.3f, weight is %.4f \n", bjet_SF, weight);

    }
    else if (std::abs(flavour1 - 4.) < 0.01){ //cjet
        bjet1_SF = b_readers.c_reader.eval_auto_bounds("central", BTagEntry::FLAV_C, eta1, pt1);
        bjet1_eff= btag_eff(pt1, eta1, btag_effs.c_eff);
        //printf("C jet, SF is %0.3f, weight is %.4f \n", bjet_SF, weight);
    }
    else{ //udsg jet
        bjet1_SF = b_readers.udsg_reader.eval_auto_bounds("central", BTagEntry::FLAV_UDSG, eta1,pt1);
        bjet1_eff= btag_eff(pt1, eta1, btag_effs.udsg_eff);
        //printf("UDSG jet, SF is %0.3f, weight is %.4f \n", bjet_SF, weight);
    }


    if(std::abs(flavour2 - 5.) < 0.01){ //bjet
        bjet2_SF = b_readers.b_reader.eval_auto_bounds("central", BTagEntry::FLAV_B, eta2, pt2);
        bjet2_eff= btag_eff(pt2, eta2, btag_effs.b_eff);
        //printf("B jet, SF is %0.3f, weight is %.4f \n", bjet_SF, weight);

    }
    else if (std::abs(flavour2 - 4.) < 0.01){ //cjet
        bjet2_SF = b_readers.c_reader.eval_auto_bounds("central", BTagEntry::FLAV_C, eta2, pt2);
        bjet2_eff= btag_eff(pt2, eta2, btag_effs.c_eff);
        //printf("C jet, SF is %0.3f, weight is %.4f \n", bjet_SF, weight);
    }
    else{ //udsg jet
        bjet2_SF = b_readers.udsg_reader.eval_auto_bounds("central", BTagEntry::FLAV_UDSG, eta2,pt2);
        bjet2_eff= btag_eff(pt2, eta2, btag_effs.udsg_eff);
        //printf("UDSG jet, SF is %0.3f, weight is %.4f \n", bjet_SF, weight);
    }
    
    Double_t P_mc = bjet1_eff + bjet2_eff - bjet1_eff*bjet2_eff;
    Double_t P_data = bjet1_SF*bjet1_eff + bjet2_SF*bjet2_eff - bjet1_SF*bjet2_SF*bjet1_eff*bjet2_eff;
    return P_data/P_mc;
}



void setup_SFs(SFs *runs_BCDEF, SFs *runs_GH, BTag_readers *btag_r, BTag_effs *b_effs){
    TH1::AddDirectory(kFALSE);
    BTagCalibration calib("csvv1", "SFs/cMVAv2_Moriond17_B_H.csv");
    btag_r->b_reader = BTagCalibrationReader(BTagEntry::OP_MEDIUM, "central");
    btag_r->b_reader.load(calib, BTagEntry::FLAV_B, "ttbar");
    btag_r->c_reader = BTagCalibrationReader(BTagEntry::OP_MEDIUM, "central");
    btag_r->c_reader.load(calib, BTagEntry::FLAV_C, "ttbar");
    btag_r->udsg_reader = BTagCalibrationReader(BTagEntry::OP_MEDIUM, "central");
    btag_r->udsg_reader.load(calib, BTagEntry::FLAV_UDSG, "incl");

    TFile *f0 = TFile::Open("SFs/BTag_efficiency_may24.root");
    TDirectory *subdir0 = gDirectory;
    TH2D *b_eff = (TH2D *) subdir0->Get("b_eff")->Clone();
    TH2D *c_eff = (TH2D *) subdir0->Get("c_eff")->Clone();
    TH2D *udsg_eff = (TH2D *) subdir0->Get("udsg_eff")->Clone();
    b_eff->SetDirectory(0);
    c_eff->SetDirectory(0);
    udsg_eff->SetDirectory(0);
    b_effs->b_eff = b_eff;
    b_effs->c_eff = c_eff;
    b_effs->udsg_eff = udsg_eff;



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

void setup_el_SF(el_SFs *sf){
    //Setup electron SF's
    TFile *f7 = TFile::Open("SFs/egammaEffi.txt_EGM2D.root");
    TDirectory *subdir7 = gDirectory;
    TH2D *h = (TH2D *) subdir7->Get("EGamma_SF2D")->Clone();
    h->SetDirectory(0);
    sf->h = h;
    f7->Close();
    //el_SF->Print();
}

