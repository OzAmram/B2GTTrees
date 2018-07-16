#include "TROOT.h"
#include "TTree.h"
#include "TFile.h"

TFile *f_elel_mc, *f_elel_back, *f_elel_data, *f_elel_QCD, *f_elel_WJets, *f_elel_WJets_contam, *f_elel_QCD_contam;
TTree *t_elel_mc, *t_elel_back, *t_elel_data, *t_elel_QCD, *t_elel_WJets, *t_elel_WJets_contam, *t_elel_QCD_contam, *t_elel_nosig;

TFile *f_mumu_mc, *f_mumu_back, *f_mumu_data, *f_mumu_QCD, *f_mumu_WJets, *f_mumu_WJets_contam, *f_mumu_QCD_contam;
TTree *t_mumu_mc, *t_mumu_back, *t_mumu_data, *t_mumu_QCD, *t_mumu_WJets, *t_mumu_WJets_contam, *t_mumu_QCD_contam, *t_mumu_nosig;

int n_xf_bins = 5;
Float_t xf_bins[] = {0., 0.02, 0.04, 0.07, 0.10, 1.0};
//int n_cost_bins = 6;
//Float_t cost_bins[] = {-1.0, -.667, -.333, 0., 0.33, 0.667,  1.0};
//int n_cost_bins = 8;
//Float_t cost_bins[] = {-1.0, -.75, -.5, -.25, 0., 0.25, 0.5,  0.75, 1.0};
int n_cost_bins = 10;
Float_t cost_bins[] = {-1.0, -.8, -.6, -.4, -.2, 0., 0.2, 0.4, 0.6, 0.8, 1.0};
//int n_cost_bins = 12;
//Float_t cost_bins[] = {-1.0, -0.8333, -0.6667, -0.5, -0.3333, -0.1667, 0., 0.1667, 0.3333, 0.5, 0.6667, 0.8333, 1.0};
//int n_cost_bins = 14;
//Float_t cost_bins[] = {-1.0, -.857, -.714, -.571, -.429, -0.286, -.143,  0., 0.143, .286, 0.429, 0.571, 0.714, 0.857, 1.0};
int n_m_bins = 6;
int n_pt_bins = 6;
Double_t m_bins[] = {150,200,   250,    350,    500,    700, 100000};
Double_t alphas[6] = {0.109, 0.078, 0.0762, 0.112, 0.065, 0.06};
Double_t alpha_unc[6] = {0.015, 0.015, 0.02, 0.03,   0.02, 0.02};
Double_t m_alphas[6] = {0.109, 0.078, 0.0762, 0.112, 0.065, 0.06};
Double_t m_alpha_unc[6] = {0.015, 0.015, 0.02, 0.03,   0.02, 0.02};
Double_t pt_bins[] =        {0.,25.,  50., 80.,   120.,   200., 10000.};
Double_t pt_alphas[6] =    {0.007, 0.136, 0.337, 0.546, 0.776, 0.945};
Double_t pt_alpha_unc[6] = {0.006, 0.015, 0.035, 0.05, 0.08, .15};
Double_t alpha;

void init(){
    //MC templates
    printf("init \n");
    f_elel_mc = (TFile*) TFile::Open("output_files/ElEl_DY_slim_june25.root");
    //f_elel_mc = (TFile*) TFile::Open("output_files/ElEl_DY_mar8.root");
    t_elel_mc = (TTree *) f_elel_mc ->Get("T_data");
    t_elel_nosig = (TTree *) f_elel_mc ->Get("T_back");
    //f_elel_back = (TFile*) TFile::Open("output_files/ElEl_combined_back_mar8.root");
    f_elel_back = (TFile*) TFile::Open("output_files/ElEl_combined_back_july05.root");
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

    f_mumu_mc = (TFile*) TFile::Open("output_files/MuMu_DY_slim_july10.root");
    //f_mumu_mc = (TFile*) TFile::Open("output_files/MuMu_DY_jan16.root");
    t_mumu_mc = (TTree *) f_mumu_mc ->Get("T_data");
    t_mumu_nosig = (TTree *) f_mumu_mc ->Get("T_back");
    //f_mumu_back = (TFile*) TFile::Open("output_files/MuMu_combined_back_jan22.root");
    f_mumu_back = (TFile*) TFile::Open("output_files/MuMu_combined_back_july10.root");
    t_mumu_back = (TTree *) f_mumu_back ->Get("T_data");

    f_mumu_data = TFile::Open("output_files/SingleMuon_data_slim_july10.root");
    //f_mumu_data = TFile::Open("output_files/SingleMuon_data_jan22.root");
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
