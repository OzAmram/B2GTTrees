#include "TROOT.h"
#include "TTree.h"
#include "TFile.h"

TFile *f_elel_mc, *f_elel_back, *f_elel_data, *f_elel_QCD, *f_elel_WJets, *f_elel_WJets_contam, *f_elel_QCD_contam;
TTree *t_elel_mc, *t_elel_back, *t_elel_data, *t_elel_QCD, *t_elel_WJets, *t_elel_WJets_contam, *t_elel_QCD_contam, *t_elel_nosig;

TFile *f_mumu_mc, *f_mumu_back, *f_mumu_data, *f_mumu_QCD, *f_mumu_WJets, *f_mumu_WJets_contam, *f_mumu_QCD_contam;
TTree *t_mumu_mc, *t_mumu_back, *t_mumu_data, *t_mumu_QCD, *t_mumu_WJets, *t_mumu_WJets_contam, *t_mumu_QCD_contam, *t_mumu_nosig;

TFile *f_emu_dy, *f_emu_back, *f_emu_data, *f_emu_QCD, *f_emu_WJets, *f_emu_WJets_contam;
TTree *t_emu_dy, *t_emu_back, *t_emu_data, *t_emu_QCD, *t_emu_WJets, *t_emu_WJets_contam;

TFile *f_elel_ss_dy, *f_elel_ss_back, *f_elel_ss_data, *f_elel_ss_QCD, *f_elel_ss_WJets, *f_elel_ss_WJets_contam, *f_elel_ss_QCD_contam;
TTree *t_elel_ss_dy, *t_elel_ss_back, *t_elel_ss_data, *t_elel_ss_QCD, *t_elel_ss_WJets, *t_elel_ss_WJets_contam, *t_elel_ss_QCD_contam;

TFile *f_mumu_ss_dy, *f_mumu_ss_back, *f_mumu_ss_data, *f_mumu_ss_QCD, *f_mumu_ss_WJets, *f_mumu_ss_WJets_contam, *f_mumu_ss_QCD_contam;
TTree *t_mumu_ss_dy, *t_mumu_ss_back, *t_mumu_ss_data, *t_mumu_ss_QCD, *t_mumu_ss_WJets, *t_mumu_ss_WJets_contam, *t_mumu_ss_QCD_contam;

TFile *f_emu_ss_data, *f_emu_ss_ttbar, *f_emu_ss_dy, *f_emu_ss_diboson;
TTree *t_emu_ss_data, *t_emu_ss_ttbar, *t_emu_ss_dy, *t_emu_ss_diboson;

TFile *f_mumu_gamgam, *f_elel_gamgam;
TTree *t_mumu_gamgam, *t_elel_gamgam;

int n_xf_bins = 4;
Float_t xf_bins[] = {0., 0.04, 0.07, 0.10, 1.0};
int n_cost_bins = 10;
Float_t cost_bins[] = {-1.0, -.8, -.6, -.4, -.2, 0., 0.2, 0.4, 0.6, 0.8, 1.0};
int n_cost_ss_bins = n_cost_bins/2;
Float_t cost_ss_bins[] = {-1.0, -0.8, -0.6, -0.4, -0.2, 0.0};
int n_m_bins = 6;
Double_t m_bins[] = {150, 200,   250,    350,    500,    700, 100000};
Double_t alphas_num[6] = {0.125, 0.115, 0.07, 0.065, 0.05, 0.04};
Double_t alpha_num_unc[6] = {0.04, 0.04, 0.04, 0.04,   0.03, 0.03};
Double_t alphas_denom[6] = {0.15, 0.13, 0.08, 0.065, 0.07, 0.045};
Double_t alpha_denom_unc[6] = {0.05, 0.03, 0.03, 0.04,   0.04, 0.04};

void init(int year){
    //MC templates
    printf("init \n");
    if(year == 2016){
        f_elel_data = TFile::Open("../analyze/output_files/2016/ElEl_data16_june23.root");
        t_elel_data = (TTree *)f_elel_data->Get("T_data"); 

        f_elel_mc = (TFile*) TFile::Open("../analyze/output_files/2016/ElEl_dy_july1.root");
        t_elel_mc = (TTree *) f_elel_mc ->Get("T_data");
        t_elel_nosig = (TTree *) f_elel_mc ->Get("T_back");

        f_elel_back = (TFile*) TFile::Open("../analyze/output_files/2016/ElEl_comb_back_july1.root");
        t_elel_back = (TTree *) f_elel_back ->Get("T_data");

        f_elel_gamgam = TFile::Open("../analyze/output_files/2016/ElEl_gamgam_back_june25.root");
        t_elel_gamgam = (TTree *)f_elel_gamgam->Get("T_data");

        f_elel_QCD = TFile::Open("../analyze/output_files/2016/ElEl_QCD_est_june25.root");
        t_elel_QCD = (TTree *)f_elel_QCD->Get("T_data");

        f_elel_WJets = TFile::Open("../analyze/output_files/2016/ElEl_WJets_est_june25.root");
        t_elel_WJets = (TTree *)f_elel_WJets->Get("T_data");

        f_elel_WJets_contam = TFile::Open("../analyze/output_files/2016/ElEl_WJets_MC_july8.root");
        t_elel_WJets_contam = (TTree *)f_elel_WJets_contam->Get("T_data");

        f_elel_QCD_contam = TFile::Open("../analyze/output_files/2016/ElEl_QCD_MC_july8.root");
        t_elel_QCD_contam = (TTree *)f_elel_QCD_contam->Get("T_data");
        ////////////////////////////////////////

        f_mumu_data = TFile::Open("../analyze/output_files/2016/MuMu_data_july12.root");
        t_mumu_data = (TTree *)f_mumu_data->Get("T_data"); 

        f_mumu_mc = (TFile*) TFile::Open("../analyze/output_files/2016/MuMu_dy_july15.root");
        t_mumu_mc = (TTree *) f_mumu_mc ->Get("T_data");
        t_mumu_nosig = (TTree *) f_mumu_mc ->Get("T_back");
        f_mumu_back = (TFile*) TFile::Open("../analyze/output_files/2016/MuMu_comb_back_july15.root");
        t_mumu_back = (TTree *) f_mumu_back ->Get("T_data");


        f_mumu_gamgam = TFile::Open("../analyze/output_files/2016/MuMu_gamgam_july15.root");
        t_mumu_gamgam = (TTree *)f_mumu_gamgam->Get("T_data");

        f_mumu_QCD = TFile::Open("../analyze/output_files/2016/MuMu_QCD_est_june25.root");
        t_mumu_QCD = (TTree *)f_mumu_QCD->Get("T_data");

        f_mumu_WJets = TFile::Open("../analyze/output_files/2016/MuMu_WJets_est_june25.root");
        t_mumu_WJets = (TTree *)f_mumu_WJets->Get("T_data");

        f_mumu_WJets_contam = TFile::Open("../analyze/output_files/2016/MuMu_WJets_MC_july8.root");
        t_mumu_WJets_contam = (TTree *)f_mumu_WJets_contam->Get("T_data");

        f_mumu_QCD_contam = TFile::Open("../analyze/output_files/2016/MuMu_QCD_MC_july8.root");
        t_mumu_QCD_contam = (TTree *)f_mumu_QCD_contam->Get("T_data");

        return;
    }
    if (year == 2017){
        f_mumu_mc = (TFile*) TFile::Open("../analyze/output_files/2017/MuMu_dy17_sep3.root");
        t_mumu_mc = (TTree *)f_mumu_mc->Get("T_sig");

        f_elel_mc = (TFile*) TFile::Open("../analyze/output_files/2017/ElEl_dy17_sep3.root");
        t_elel_mc = (TTree *)f_elel_mc->Get("T_sig");

        return;
    }

}
void init_emu(int year){


    f_emu_data = TFile::Open("../analyze/output_files/2016/EMu_data_july1.root");
    t_emu_data = (TTree *)f_emu_data->Get("T_data");

                         
    f_emu_back = TFile::Open("../analyze/output_files/2016/EMu_comb_back_july2.root");
    t_emu_back = (TTree *)f_emu_back->Get("T_data");

    f_emu_dy = TFile::Open("../analyze/output_files/2016/EMu_dy_july2.root");
    t_emu_dy = (TTree *)f_emu_dy->Get("T_data");

    f_emu_QCD = TFile::Open("../analyze/output_files/2016/EMu_QCD_est_july1.root");
    t_emu_QCD = (TTree *)f_emu_QCD->Get("T_data");

    f_emu_WJets = TFile::Open("../analyze/output_files/2016/EMu_WJets_est_july1.root");
    t_emu_WJets = (TTree *)f_emu_WJets->Get("T_data");

    f_emu_WJets_contam = TFile::Open("../analyze/output_files/2016/EMu_WJets_MC_july12.root");
    t_emu_WJets_contam = (TTree *)f_emu_WJets_contam->Get("T_data");

}



void init_ss(int year){
    f_elel_ss_data = TFile::Open("../analyze/output_files/2016/ElEl_data_ss_july3.root");
    t_elel_ss_data = (TTree *)f_elel_ss_data->Get("T_data");

    f_elel_ss_dy = TFile::Open("../analyze/output_files/2016/ElEl_dy_ss_july3.root");
    t_elel_ss_dy = (TTree *)f_elel_ss_dy->Get("T_data");

    f_elel_ss_back = TFile::Open("../analyze/output_files/2016/ElEl_comb_back_ss_july8.root");
    t_elel_ss_back = (TTree *)f_elel_ss_back->Get("T_data");

//////////////////////////////////////////////
    f_mumu_ss_data = TFile::Open("../analyze/output_files/2016/MuMu_data_ss_july12.root");
    t_mumu_ss_data = (TTree *)f_mumu_ss_data->Get("T_data");

    f_mumu_ss_dy = TFile::Open("../analyze/output_files/2016/MuMu_dy_ss_july3.root");
    t_mumu_ss_dy = (TTree *)f_mumu_ss_dy->Get("T_data");

    f_mumu_ss_back = TFile::Open("../analyze/output_files/2016/MuMu_comb_back_ss_july8.root");
    t_mumu_ss_back = (TTree *)f_mumu_ss_back->Get("T_data");
}

