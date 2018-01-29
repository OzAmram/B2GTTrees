
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
#include "../../analyze/TemplateMaker.C"





void get_ratio(){
    /*
    int type = FLAG_MUONS;
    TFile *f_data = TFile::Open("../analyze/output_files/EMu_SingleMuon_data_nov3.root");
    TTree *t_data = (TTree *)f_data->Get("T_data");

    TFile *f_mc = TFile::Open("../analyze/output_files/EMu_combined_back_Mu_jan22.root");
    TTree *t_mc = (TTree *)f_mc->Get("T_data");
    */
                                
    int type = FLAG_ELECTRONS;
    TFile *f_data = TFile::Open("../analyze/output_files/EMu_SingleElectron_data_nov3.root");
    TTree *t_data = (TTree *)f_data->Get("T_data");

    TFile *f_mc = TFile::Open("../analyze/output_files/EMu_combined_back_El_jan22.root");
    TTree *t_mc = (TTree *)f_mc->Get("T_data");

    TFile *f_QCD = TFile::Open("../analyze/output_files/EMu_QCD_fakerate_est_jan24.root");
    TTree *t_QCD = (TTree *)f_QCD->Get("T_data");

    TFile *f_WJets = TFile::Open("../analyze/output_files/EMu_WJets_fakerate_est_jan24.root");
    TTree *t_WJets = (TTree *)f_WJets->Get("T_data");

    TFile *f_WJets_mc = TFile::Open("../analyze/output_files/EMu_WJets_MC_jan24.root");
    TTree *t_WJets_mc = (TTree *)f_WJets_mc->Get("T_data");
     
    TH1F *data_m = new TH1F("data_m", "WJets", 30, 150, 2000);


    TH1F *mc_m = new TH1F("mc_m", "WJets", 30, 150, 2000);
    TH1F *qcd_m = new TH1F("qcd_m", "WJets", 30, 150, 2000);

    make_emu_m_hist(t_data, data_m, true, type);
    make_emu_m_hist(t_mc, mc_m, false, type);
    Fakerate_est_emu(t_WJets, t_QCD,t_WJets_mc, qcd_m, type);

    Double_t data_count = data_m->Integral();
    Double_t fake_count = qcd_m->Integral();
    Double_t mc_count = mc_m->Integral();

    Double_t fake_unc = 0.3 * fake_count;
    Double_t mc_unc = sqrt(mc_count);

    printf("Data count %.0f \n", data_count);
    printf("MC count %.0f \n", mc_count);
    printf("Fake count %.0f \n", fake_count);
    Double_t ratio = data_count / (mc_count + fake_count);
    Double_t unc = sqrt( (data_count/(mc_count + fake_count)/(mc_count + fake_count)) +  
                          pow(data_count/(mc_count + fake_count)/(mc_count + fake_count), 2) * (fake_unc*fake_unc + mc_unc*mc_unc));
    printf("Ratio is %1.2f +/- %1.2f \n", ratio, unc);
    return;
}









