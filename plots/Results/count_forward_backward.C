#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <iostream>
#include <cstring>
#include <algorithm>
#include <vector>

#include "TROOT.h"
#include "TFile.h"
#include "TObject.h"
#include "TRatioPlot.h"
#include "TH1F.h"
#include "TProfile.h"
#include "TStyle.h"
#include "TCanvas.h"
#include "TPostScript.h"
#include "TLorentzVector.h"
#include "TVector3.h"
#include "TFitter.h"
#include "TSystem.h"
#include "Math/Functor.h"
#include "../../Utils/HistMaker.C"
#include "../tdrstyle.C"
#include "../CMS_lumi.C"
#include "../../Utils/root_files.h"






void count_forward_backward(){
    init();
    init_gamgam();

    int flag_type = FLAG_M_BINS;



    
    Float_t m_bins[] = {150,200,   250,    350,    500,    700, 100000};
    for(int i =0; i<6; i++){

        double var_low, var_high;

            var_low = m_bins[i];
            var_high = m_bins[i+1];

        TH1F *mumu_data_cost = new TH1F("mumu_data_cost", "Data", 2, -1.,1.);
        TH1F *mumu_mc_back_cost = new TH1F("mumu_mc_cost", "DiBoson (WW, WZ,ZZ)", 2, -1.,1);
        TH1F *mumu_dy_cost = new TH1F("mumu_dy_cost", "DiBoson (WW, WZ,ZZ)", 2, -1.,1);
        TH1F *mumu_QCD_cost = new TH1F("mumu_QCD_cost", "QCD", 2, -1.,1);
        TH1F *mumu_gam_cost = new TH1F("mumu_gamgam_cost", "gamma", 2, -1.,1);

        TH1F *elel_data_cost = new TH1F("elel_data_cost", "Data", 2, -1.,1.);
        TH1F *elel_mc_back_cost = new TH1F("elel_mc_cost", "DiBoson (WW, WZ,ZZ)", 2, -1.,1);
        TH1F *elel_dy_cost = new TH1F("elel_dy_cost", "DiBoson (WW, WZ,ZZ)", 2, -1.,1);
        TH1F *elel_QCD_cost = new TH1F("elel_QCD_cost", "QCD", 2, -1.,1);
        TH1F *elel_gam_cost = new TH1F("elel_gamgam_cost", "gamma", 2, -1.,1);
        
        TH1F *h_dummy = new TH1F("dummy", "", 100, 0, 100.);
        bool do_RC = true;
        int type = FLAG_MUONS;
        make_m_cost_pt_xf_hist(t_mumu_data, h_dummy, mumu_data_cost, h_dummy, h_dummy, true, type, do_RC, var_low, var_high);
        make_m_cost_pt_xf_hist(t_mumu_back, h_dummy, mumu_mc_back_cost, h_dummy, h_dummy, false, type, do_RC, var_low, var_high);
        make_m_cost_pt_xf_hist(t_mumu_gamgam, h_dummy, mumu_gam_cost, h_dummy, h_dummy, false, type, do_RC, var_low, var_high);
        make_m_cost_pt_xf_hist(t_mumu_nosig, h_dummy, mumu_mc_back_cost, h_dummy, h_dummy, false, type, do_RC, var_low, var_high);
        make_m_cost_pt_xf_hist(t_mumu_mc, h_dummy, mumu_dy_cost, h_dummy, h_dummy, false, type, do_RC, var_low, var_high);
        bool ss_qcd = true;
        bool in_os_region = true;
        Fakerate_est_mu(t_mumu_WJets, t_mumu_QCD, t_mumu_WJets_contam, t_mumu_QCD_contam, h_dummy, mumu_QCD_cost, h_dummy, h_dummy, var_low, var_high, ss_qcd, in_os_region);


        type = FLAG_ELECTRONS;
        make_m_cost_pt_xf_hist(t_elel_data, h_dummy, elel_data_cost, h_dummy, h_dummy, true, type, do_RC, var_low, var_high);
        make_m_cost_pt_xf_hist(t_elel_back, h_dummy, elel_mc_back_cost, h_dummy, h_dummy, false, type, do_RC, var_low, var_high);
        make_m_cost_pt_xf_hist(t_elel_gamgam, h_dummy, elel_gam_cost, h_dummy, h_dummy, false, type, do_RC, var_low, var_high);
        make_m_cost_pt_xf_hist(t_elel_nosig, h_dummy, elel_mc_back_cost, h_dummy, h_dummy, false, type, do_RC, var_low, var_high);
        make_m_cost_pt_xf_hist(t_elel_mc, h_dummy, elel_dy_cost, h_dummy, h_dummy, false, type, do_RC, var_low, var_high);

        Fakerate_est_el(t_elel_WJets, t_elel_QCD, t_elel_WJets_contam, t_elel_QCD_contam, h_dummy, elel_QCD_cost, h_dummy, h_dummy, var_low, var_high, ss_qcd);

        double n_mumu_f_data = mumu_data_cost->GetBinContent(2);
        double n_mumu_b_data = mumu_data_cost->GetBinContent(1);

        double n_mumu_f_dy = mumu_dy_cost->GetBinContent(2);
        double n_mumu_b_dy = mumu_dy_cost->GetBinContent(1);
        
        double n_mumu_f_back = mumu_mc_back_cost->GetBinContent(2);
        double n_mumu_b_back = mumu_mc_back_cost->GetBinContent(1);

        double n_mumu_f_gam = mumu_gam_cost->GetBinContent(2);
        double n_mumu_b_gam = mumu_gam_cost->GetBinContent(1);

        double n_mumu_f_qcd = mumu_QCD_cost->GetBinContent(2);
        double n_mumu_b_qcd = mumu_QCD_cost->GetBinContent(1);


        printf("MUONS \n");
        printf("Var from %.0f to %.0f  \n", var_low, var_high);
        printf("total data events %.0f \n", mumu_data_cost->Integral());
        printf("Data backward %0.f,forward %.0f \n", n_mumu_b_data, n_mumu_f_data);
        printf("DY backward %0.f,forward %.0f \n", n_mumu_b_dy, n_mumu_f_dy);
        printf("MC Backgrounds backward %0.f,forward %.0f \n", n_mumu_b_back, n_mumu_f_back);
        printf("Photon Induced Backgrounds backward %0.f,forward %.0f \n", n_mumu_b_gam, n_mumu_f_gam);
        printf("Fakes Backgrounds backward %0.f,forward %.0f \n", n_mumu_b_qcd, n_mumu_f_qcd);

        printf("%.0f-%.0f & %.0f & %.0f &\n", var_low, var_high, n_mumu_b_data, n_mumu_f_data);
        printf("%.0f & %.0f &\n", var_low, var_high, n_mumu_b_dy, n_mumu_f_dy);
        printf("%.0f & %.0f &\n", var_low, var_high, n_mumu_b_back, n_mumu_f_back);
        printf("%.0f & %.0f &\n", var_low, var_high, n_mumu_b_gam, n_mumu_f_gam);
        printf("%.0f & %.0f \\\\ \n", var_low, var_high, n_mumu_b_qcd, n_mumu_f_qcd);


        double n_elel_f_data = elel_data_cost->GetBinContent(2);
        double n_elel_b_data = elel_data_cost->GetBinContent(1);

        double n_elel_f_dy = elel_dy_cost->GetBinContent(2);
        double n_elel_b_dy = elel_dy_cost->GetBinContent(1);
        
        double n_elel_f_back = elel_mc_back_cost->GetBinContent(2);
        double n_elel_b_back = elel_mc_back_cost->GetBinContent(1);

        double n_elel_f_gam = elel_gam_cost->GetBinContent(2);
        double n_elel_b_gam = elel_gam_cost->GetBinContent(1);

        double n_elel_f_qcd = elel_QCD_cost->GetBinContent(2);
        double n_elel_b_qcd = elel_QCD_cost->GetBinContent(1);


        printf("ELECTRONS \n");
        printf("Var from %.0f to %.0f  \n", var_low, var_high);
        printf("total data events %.0f \n", elel_data_cost->Integral());
        printf("Data backward %0.f,forward %.0f \n", n_elel_b_data, n_elel_f_data);
        printf("DY backward %0.f,forward %.0f \n", n_elel_b_dy, n_elel_f_dy);
        printf("MC Backgrounds backward %0.f,forward %.0f \n", n_elel_b_back, n_elel_f_back);
        printf("Photon Induced Backgrounds backward %0.f,forward %.0f \n", n_elel_b_gam, n_elel_f_gam);
        printf("Fakes Backgrounds backward %0.f,forward %.0f \n", n_elel_b_qcd, n_elel_f_qcd);


        printf("%.0f-%.0f & %.0f & %.0f & \n", var_low, var_high, n_elel_b_data, n_elel_f_data);
        printf("%.0f & %.0f & \n", var_low, var_high, n_elel_b_dy, n_elel_f_dy);
        printf("%.0f & %.0f & \n", var_low, var_high, n_elel_b_back, n_elel_f_back);
        printf("%.0f & %.0f & \n", var_low, var_high, n_elel_b_gam, n_elel_f_gam);
        printf("%.0f & %.0f \\\\ \n", var_low, var_high, n_elel_b_qcd, n_elel_f_qcd);

    }

}


