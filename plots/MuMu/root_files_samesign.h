
#include "TROOT.h"
#include "TFile.h"

TFile *f_data, *f_mc, *f_mc_nosig, *f_ttbar, *f_QCD, *f_WJets, *f_WJets_mc, *f_QCD_mc, *f_diboson, *f_wt;
TTree *t_data, *t_mc, *t_mc_nosig, *t_ttbar, *t_QCD, *t_WJets, *t_WJets_mc, *t_QCD_mc, *t_diboson, *t_wt;


void init(){
    f_data = TFile::Open("../analyze/output_files/MuMu_samesign_data_dec3.root");
    t_data = (TTree *)f_data->Get("T_data");

    f_QCD = TFile::Open("../analyze/output_files/MuMu_samesign_fakerate_qcd_est_dec3.root");
    t_QCD = (TTree *)f_QCD->Get("T_data");

    f_WJets = TFile::Open("../analyze/output_files/MuMu_samesign_fakerate_wjets_est_dec3.root");
    t_WJets = (TTree *)f_WJets->Get("T_data");
    f_WJets_mc = TFile::Open("../analyze/output_files/MuMu_samesign_fakerate_wjets_MC_dec3.root");
    t_WJets_mc = (TTree *)f_WJets_mc->Get("T_data");

    //dummy tree
    t_QCD_mc = new TTree();

    f_diboson = TFile::Open("../analyze/output_files/MuMu_samesign_background_dec3.root");
    t_diboson = (TTree *)f_diboson->Get("T_data");

}
