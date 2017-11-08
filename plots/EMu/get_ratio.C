
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
    int type = FLAG_ELECTRONS;
    /*
    TFile *f_data = TFile::Open("../analyze/output_files/EMu_SingleMuon_data_nov3.root");
    TTree *t_data = (TTree *)f_data->Get("T_data");

    TFile *f_mc = TFile::Open("../analyze/output_files/EMu_background_mu_nov3.root");
    TTree *t_mc = (TTree *)f_mc->Get("T_data");
    */
                                
    TFile *f_data = TFile::Open("../analyze/output_files/EMu_SingleElectron_data_nov3.root");
    TTree *t_data = (TTree *)f_data->Get("T_data");

    TFile *f_mc = TFile::Open("../analyze/output_files/EMu_background_El_nov3.root");
    TTree *t_mc = (TTree *)f_mc->Get("T_data");

     
    TH1F *data_m = new TH1F("data_m", "WJets", 30, 150, 2000);


    TH1F *mc_m = new TH1F("mc_m", "WJets", 30, 150, 2000);

    make_emu_m_hist(t_data, data_m, true, type);
    make_emu_m_hist(t_mc, mc_m, false, type);

    Double_t data_count = data_m->Integral();
    Double_t mc_count = mc_m->Integral();

    printf("Data count %.0f \n", data_count);
    printf("MC count %.0f \n", mc_count);
    Double_t ratio = data_count / mc_count;
    Double_t unc = sqrt(data_count/mc_count/mc_count + pow((data_count/mc_count/mc_count),2) * mc_count);
    printf("Ratio is %1.2f +/- %1.2f \n", ratio, unc);
    return;
}









