#include "TROOT.h"
#include "TFile.h"

TFile *f_data, *f_mc, *f_mc_nosig, *f_ttbar, *f_QCD, *f_WJets, *f_WJets_mc, *f_QCD_mc, *f_diboson, *f_wt;
TTree *t_data, *t_mc, *t_mc_nosig, *t_ttbar, *t_QCD, *t_WJets, *t_WJets_mc, *t_QCD_mc, *t_diboson, *t_wt;


void init(){
    f_data = TFile::Open("../analyze/output_files/SingleMuon_data_slim_nov26.root");
    t_data = (TTree *)f_data->Get("T_data");
    f_mc = TFile::Open("../analyze/output_files/MuMu_DY_slim_nov2.root");
    t_mc = (TTree *)f_mc->Get("T_data");
    t_mc_nosig = (TTree *)f_mc->Get("T_back");
    f_ttbar = TFile::Open("../analyze/output_files/MuMu_ttbar_sep4.root");
    t_ttbar = (TTree *)f_ttbar->Get("T_data");

    f_QCD = TFile::Open("../analyze/output_files/MuMu_QCD_est_nov26.root");
    t_QCD = (TTree *)f_QCD->Get("T_data");

    f_WJets = TFile::Open("../analyze/output_files/MuMu_WJets_est_nov26.root");
    t_WJets = (TTree *)f_WJets->Get("T_data");
    f_WJets_mc = TFile::Open("../analyze/output_files/MuMu_WJets_MC_contam_nov5.root");
    t_WJets_mc = (TTree *)f_WJets_mc->Get("T_data");

    f_QCD_mc = TFile::Open("../analyze/output_files/MuMu_QCD_MC_contam_nov5.root");
    t_QCD_mc = (TTree *)f_QCD_mc->Get("T_data");

    f_diboson = TFile::Open("../analyze/output_files/MuMu_diboson_sep4.root");
    t_diboson = (TTree *)f_diboson->Get("T_data");

    f_wt = TFile::Open("../analyze/output_files/MuMu_WT_sep4.root");
    t_wt = (TTree *)f_wt->Get("T_data");
}
