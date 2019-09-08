
#include "CondFormats/BTauObjects/interface/BTagCalibration.h"
#include "CondTools/BTau/interface/BTagCalibrationReader.h"
#include "TH2D.h"

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
Double_t get_btag_weight(Double_t pt, Double_t eta, Float_t flavour, BTag_effs btag_effs, BTag_readers b_readers, int systematic = 0){
    //compute weighting from btagging scale factors

    Double_t weight, bjet_SF;

    char const *sys;
    if (systematic == 0) sys = "central";
    if (systematic == 1) sys = "up";
    if (systematic == -1) sys = "down";

    if(std::abs(flavour - 5.) < 0.01){ //bjet
        bjet_SF = b_readers.b_reader.eval_auto_bounds(sys, BTagEntry::FLAV_B, eta, pt);
        weight = btag_weight_helper(pt, eta, bjet_SF, btag_effs.b_eff);
        //printf("B jet, SF is %0.3f, weight is %.4f \n", bjet_SF, weight);

    }
    else if (std::abs(flavour - 4.) < 0.01){ //cjet
        bjet_SF = b_readers.c_reader.eval_auto_bounds(sys, BTagEntry::FLAV_C, eta, pt);
        weight = btag_weight_helper(pt, eta, bjet_SF, btag_effs.c_eff);
        //printf("C jet, SF is %0.3f, weight is %.4f \n", bjet_SF, weight);
    }
    else{ //udsg jet
        bjet_SF = b_readers.udsg_reader.eval_auto_bounds(sys, BTagEntry::FLAV_UDSG, eta,pt);
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

void setup_btag_SFs(BTag_readers *btag_r, BTag_effs *b_effs, int year){
    TH1::AddDirectory(kFALSE);
    char file[80];
    // USE OLD FOR NOW, CHANGE WHEN REMAKE NTUPLES
    BTagCalibration calib;
    //if(year == 2016) sprintf(file, "%s", "SFs/2016/DeepCSV_2016LegacySF_V1.csv");
    if(year == 2016) calib = BTagCalibration("DeepCSV", "SFs/2016/DeepCSV_2016LegacySF_V1.csv");
    if(year == 2017) calib = BTagCalibration("DeepCSV", "SFs/2017/DeepCSV_94XSF_V4_B_F.csv");
    if(year == 2018) calib = BTagCalibration("DeepCSV", "SFs/2018/DeepCSV_102XSF_V1.csv");


    btag_r->b_reader = BTagCalibrationReader (BTagEntry::OP_MEDIUM, "central", {"up", "down"});
    btag_r->b_reader.load(calib, BTagEntry::FLAV_B, "comb");
    btag_r->c_reader = BTagCalibrationReader (BTagEntry::OP_MEDIUM, "central", {"up", "down"});
    btag_r->c_reader.load(calib, BTagEntry::FLAV_C, "comb");
    btag_r->udsg_reader = BTagCalibrationReader (BTagEntry::OP_MEDIUM, "central", {"up", "down"});
    btag_r->udsg_reader.load(calib, BTagEntry::FLAV_UDSG, "incl");

    /*
    btag_r->b_reader_up = BTagCalibrationReader(BTagEntry::OP_MEDIUM, "up");
    btag_r->b_reader_up.load(calib, BTagEntry::FLAV_B, "comb");
    btag_r->c_reader_up = BTagCalibrationReader(BTagEntry::OP_MEDIUM, "up");
    btag_r->c_reader_up.load(calib, BTagEntry::FLAV_C, "comb");
    btag_r->udsg_reader_up = BTagCalibrationReader(BTagEntry::OP_MEDIUM, "up");
    btag_r->udsg_reader_up.load(calib, BTagEntry::FLAV_UDSG, "incl");
    printf("3 \n");

    btag_r->b_reader_down = BTagCalibrationReader(BTagEntry::OP_MEDIUM, "down");
    btag_r->b_reader_down.load(calib, BTagEntry::FLAV_B, "comb");
    btag_r->c_reader_down = BTagCalibrationReader(BTagEntry::OP_MEDIUM, "down");
    btag_r->c_reader_down.load(calib, BTagEntry::FLAV_C, "comb");
    btag_r->udsg_reader_down = BTagCalibrationReader(BTagEntry::OP_MEDIUM, "down");
    btag_r->udsg_reader_down.load(calib, BTagEntry::FLAV_UDSG, "incl");
    printf("4 \n");
    */

    TFile *f0 = TFile::Open("SFs/2016/BTag_efficiency_may24.root");
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
    printf("Btag SF's set up \n");
}
