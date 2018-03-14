#include "TROOT.h"
#include "TFile.h"

TFile *f_elel_mc, *f_elel_back, *f_elel_data, *f_elel_QCD, *f_elel_WJets, *f_elel_WJets_contam, *f_elel_QCD_contam;
TTree *t_elel_mc, *t_elel_back, *t_elel_data, *t_elel_QCD, *t_elel_WJets, *t_elel_WJets_contam, *t_elel_QCD_contam, *t_elel_nosig;

TFile *f_mumu_mc, *f_mumu_back, *f_mumu_data, *f_mumu_QCD, *f_mumu_WJets, *f_mumu_WJets_contam, *f_mumu_QCD_contam;
TTree *t_mumu_mc, *t_mumu_back, *t_mumu_data, *t_mumu_QCD, *t_mumu_WJets, *t_mumu_WJets_contam, *t_mumu_QCD_contam, *t_mumu_nosig;

void init(){
    //MC templates
    printf("init \n");
    f_elel_mc = (TFile*) TFile::Open("output_files/ElEl_DY_mar8.root");
    t_elel_mc = (TTree *) f_elel_mc ->Get("T_data");
    t_elel_nosig = (TTree *) f_elel_mc ->Get("T_back");
    f_elel_back = (TFile*) TFile::Open("output_files/ElEl_combined_back_mar8.root");
    t_elel_back = (TTree *) f_elel_back ->Get("T_data");

    f_elel_data = TFile::Open("output_files/SingleElectron_data_jan22.root");
    t_elel_data = (TTree *)f_elel_data->Get("T_data"); 

    f_elel_QCD = TFile::Open("../analyze/output_files/ElEl_QCD_est_nov2.root");
    t_elel_QCD = (TTree *)f_elel_QCD->Get("T_data");

    f_elel_WJets = TFile::Open("../analyze/output_files/ElEl_WJets_est_nov2.root");
    t_elel_WJets = (TTree *)f_elel_WJets->Get("T_data");

    f_elel_WJets_contam = TFile::Open("../analyze/FakeRate/root_files/ElEl_fakerate_WJets_MC_dec4.root");
    t_elel_WJets_contam = (TTree *)f_elel_WJets_contam->Get("T_data");

    f_elel_QCD_contam = TFile::Open("../analyze/FakeRate/root_files/ElEl_fakerate_QCD_MC_dec4.root");
    t_elel_QCD_contam = (TTree *)f_elel_QCD_contam->Get("T_data");
    ////////////////////////////////////////

    f_mumu_mc = (TFile*) TFile::Open("output_files/MuMu_DY_jan16.root");
    t_mumu_mc = (TTree *) f_mumu_mc ->Get("T_data");
    t_mumu_nosig = (TTree *) f_mumu_mc ->Get("T_back");
    f_mumu_back = (TFile*) TFile::Open("output_files/MuMu_combined_back_jan22.root");
    t_mumu_back = (TTree *) f_mumu_back ->Get("T_data");

    f_mumu_data = TFile::Open("output_files/SingleMuon_data_jan22.root");
    t_mumu_data = (TTree *)f_mumu_data->Get("T_data"); 

    f_mumu_QCD = TFile::Open("../analyze/output_files/MuMu_QCD_est_mar8.root");
    t_mumu_QCD = (TTree *)f_mumu_QCD->Get("T_data");

    f_mumu_WJets = TFile::Open("../analyze/output_files/MuMu_WJets_est_mar13.root");
    t_mumu_WJets = (TTree *)f_mumu_WJets->Get("T_data");

    f_mumu_WJets_contam = TFile::Open("../analyze/FakeRate/root_files/MuMu_fakerate_Wjets_MC_mar8.root");
    t_mumu_WJets_contam = (TTree *)f_mumu_WJets_contam->Get("T_data");

    f_mumu_QCD_contam = TFile::Open("../analyze/FakeRate/root_files/MuMu_fakerate_QCD_MC_mar8.root");
    t_mumu_QCD_contam = (TTree *)f_mumu_QCD_contam->Get("T_data");

    return;
}
