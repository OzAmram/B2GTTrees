
#include "TROOT.h"
#include "TFile.h"

TFile *f_ttbar, *f_diboson, *f_wt;
TTree *t_ttbar, *t_diboson, *t_wt;


void ee_init(){
    f_ttbar = TFile::Open("../analyze/output_files/ElEl_ttbar_sep4.root");
    t_ttbar = (TTree *)f_ttbar->Get("T_data");


    f_diboson = TFile::Open("../analyze/output_files/ElEl_diboson_sep4.root");
    t_diboson = (TTree *)f_diboson->Get("T_data");

    f_wt = TFile::Open("../analyze/output_files/ElEl_WT_sep4.root");
    t_wt = (TTree *)f_wt->Get("T_data");
}
