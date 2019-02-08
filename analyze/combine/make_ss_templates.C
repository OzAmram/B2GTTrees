
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <iostream>
#include <cstring>
#include <algorithm>
#include <vector>

#include "TObject.h"
#include "TH1F.h"
#include "TProfile.h"
#include "TStyle.h"
#include "TPostScript.h"
#include "TLorentzVector.h"
#include "TVector3.h"
#include "TFitter.h"
#include "TSystem.h"
#include "HiggsAnalysis/CombinedLimit/interface/RooParametricHist.h"
#include "RooWorkspace.h"
//#include"Minuit2/Minuit2Minimizer.h"
#include "Math/Functor.h"
#include "../TemplateMaker_systematics.C"
//#include "../TemplateMaker.C"
#include "../AFB_fit/FitUtils.C"

Double_t m_low;
Double_t m_high;

bool print = true;
const TString fout_name("combine/templates/feb7_qcd_fakerate_param.root");
TFile * fout;
char dirname[40];




RooWorkspace *w;
RooRealVar *var = new RooRealVar("var", "var", 0.,n_cost_bins * n_xf_bins);
RooRealVar *qcd_norm = new RooRealVar("Rqcd", "QCD normalization", 1, 0., 10.);



TH1F * dummy1 = new TH1F("dummy", "", 1, 0, 1);
bool do_emu_scale = false;
bool ss = true;
bool do_RC = true;


void write_roo_hist(TH1F *h){
    RooDataHist r(h->GetName(), h->GetName(), *var, h);
    w->import(r);
}


void convert_ss_qcd_to_param_hist(TH1F *h, FILE *f_log){
    //convert a hist to a parametric hist 
    RooArgList *bin_list = new RooArgList();

    char h_name[40];
    sprintf(h_name, "%s_param", h->GetName());
    for(int j=1; j <= h->GetNbinsX(); j++){

        float content = h->GetBinContent(j);
        if(content<0) printf("Bin %i Content is %.0f \n", j, content);
        float error = h->GetBinError(j);
        printf("Bin %.1f error %.1f \n", content,error);
        char bin_name[40];
        char form_name[40];
        sprintf(bin_name, "%s_bin%i",h_name, j); 
        sprintf(form_name, "%s_form%i",h_name, j); 
        RooRealVar *bin = new RooRealVar(bin_name, bin_name, content);
        fprintf(f_log, "%s param %.2f %.2f \n", bin_name, content, error);
        RooFormulaVar *form = new RooFormulaVar(form_name, form_name, "@0*@1", RooArgList(*bin, *qcd_norm));
        //bin->Print();
        bin_list->add(*form);
    }
    bin_list->Print();

    RooParametricHist *p= new RooParametricHist (h_name, h_name, *var, *bin_list, *h);
    char norm_name[40];
    sprintf(norm_name, "%s_norm", h_name);
    RooAddition *n = new RooAddition(norm_name, norm_name, *bin_list);

    
    w->import(*p);
    w->import(*n,RooFit::RecycleConflictNodes());
}


void make_ss_data_templates(){
    auto h_elel_data = new TH2F("ee_ss_data_obs", "Data template of (x_f, cost_r)",
            n_xf_bins, xf_bins, n_cost_bins, cost_bins);
    h_elel_data->SetDirectory(0);
    auto h_mumu_data = new TH2F("mumu_ss_data_obs", "Data template of (x_f, cost_r)",
            n_xf_bins, xf_bins, n_cost_bins, cost_bins);
    h_mumu_data->SetDirectory(0);

    nElEl_DataEvents = gen_data_template(t_elel_ss_data, h_elel_data,  m_low, m_high, FLAG_ELECTRONS, FLAG_M_BINS, do_RC, ss);
    nMuMu_DataEvents = gen_data_template(t_mumu_ss_data, h_mumu_data,  m_low, m_high, FLAG_MUONS, FLAG_M_BINS, do_RC, ss);
    auto h1_elel_data = convert2d(h_elel_data);
    auto h1_mumu_data = convert2d(h_mumu_data);


    printf("Integral of data templates were %.2f %.2f \n", h_elel_data->Integral(), h_mumu_data->Integral()); 
    fout->cd();
    gDirectory->cd(dirname);
    //h_elel_data->Write();
    //h_mumu_data->Write();
    write_roo_hist(h1_elel_data);
    write_roo_hist(h1_mumu_data);
    printf("Made data templates \n");
}

void make_ss_qcd_templates(FILE *f_log){
    auto h_elel_qcd = new TH2F((string("ee_ss_qcd") ).c_str(), "Combined background template",
            n_xf_bins, xf_bins, n_cost_bins, cost_bins);
    h_elel_qcd->SetDirectory(0);
    auto h_mumu_qcd = new TH2F((string("mumu_ss_qcd") ).c_str(), "Combined background template",
            n_xf_bins, xf_bins, n_cost_bins, cost_bins);
    h_mumu_qcd->SetDirectory(0);


    gen_fakes_template(t_elel_ss_WJets, t_elel_ss_QCD, t_elel_ss_WJets_contam, t_elel_ss_QCD_contam, h_elel_qcd, m_low, m_high, FLAG_ELECTRONS, FLAG_M_BINS, ss);
    gen_fakes_template(t_mumu_ss_WJets, t_mumu_ss_QCD, t_mumu_ss_WJets_contam, t_mumu_ss_QCD_contam, h_mumu_qcd, m_low, m_high, FLAG_MUONS, FLAG_M_BINS, ss);
    printf("Integral of qcd templates are %.2f %.2f \n", h_elel_qcd->Integral(), h_mumu_qcd->Integral()); 
    auto h1_elel_qcd = convert2d(h_elel_qcd);
    auto h1_mumu_qcd = convert2d(h_mumu_qcd);

    RooParametricHist *p_elel_qcd, *p_mumu_qcd;
    RooAddition *n_elel_qcd, *n_mumu_qcd;
    convert_ss_qcd_to_param_hist(h1_elel_qcd, f_log);
    convert_ss_qcd_to_param_hist(h1_mumu_qcd, f_log);
    printf("Made qcd templates \n");
}

void make_ss_mc_templates(){
    auto h_elel_dy = new TH2F((string("ee_ss_dy") ).c_str(), "ss DY",
            n_xf_bins, xf_bins, n_cost_bins, cost_bins);
    h_elel_dy->SetDirectory(0);
    auto h_mumu_dy = new TH2F((string("mumu_ss_dy") ).c_str(), "ss DY",
            n_xf_bins, xf_bins, n_cost_bins, cost_bins);
    auto h_elel_bk = new TH2F((string("ee_ss_bk") ).c_str(), "Combined background template",
            n_xf_bins, xf_bins, n_cost_bins, cost_bins);
    h_elel_bk->SetDirectory(0);
    auto h_mumu_bk = new TH2F((string("mumu_ss_bk") ).c_str(), "Combined background template",
            n_xf_bins, xf_bins, n_cost_bins, cost_bins);



    TTree *mumu_ts[1] = {t_mumu_ss_back};
    gen_combined_background_template(1, mumu_ts, h_mumu_bk, m_low, m_high, FLAG_MUONS, FLAG_M_BINS, do_RC, ss);
    mumu_ts[0] = t_mumu_ss_dy;
    gen_combined_background_template(1, mumu_ts, h_mumu_dy, m_low, m_high, FLAG_MUONS, FLAG_M_BINS, do_RC, ss);

    TTree *elel_ts[1] = {t_elel_ss_back};
    gen_combined_background_template(1, elel_ts, h_elel_bk, m_low, m_high, FLAG_ELECTRONS, FLAG_M_BINS, do_RC, ss);
    elel_ts[0] = t_elel_ss_dy;
    gen_combined_background_template(1, elel_ts, h_elel_dy, m_low, m_high, FLAG_ELECTRONS, FLAG_M_BINS, do_RC, ss);

    auto h1_elel_bk = convert2d(h_elel_bk);
    auto h1_mumu_bk = convert2d(h_mumu_bk);

    auto h1_elel_dy = convert2d(h_elel_dy);
    auto h1_mumu_dy = convert2d(h_mumu_dy);




    printf("Integral of dy templates are %.2f %.2f \n", h_elel_dy->Integral(), h_mumu_dy->Integral()); 
    printf("Integral of bkg templates are %.2f %.2f \n", h_elel_bk->Integral(), h_mumu_bk->Integral()); 
    write_roo_hist(h1_elel_dy);
    write_roo_hist(h1_mumu_dy);
    printf("Made dy templates \n");

    write_roo_hist(h1_elel_bk);
    write_roo_hist(h1_mumu_bk);
    printf("Made bk templates \n");


}





void make_ss_templates(){

    init_ss();
    init_emu_ss();


    fout = TFile::Open(fout_name, "RECREATE");
    FILE *f_log;
    char f_log_name[80];


    for(int i=0; i<n_m_bins; i++){
        fout->cd();
        snprintf(dirname, 10, "w%i", i);
        gDirectory->mkdir(dirname);
        gDirectory->cd(dirname);
        w = new RooWorkspace("w", "w");

        m_low = m_bins[i];
        m_high = m_bins[i+1];

        sprintf(f_log_name, "combine/SameSign_fits/cards/mbin%i_bins.txt", i);
        f_log = fopen(f_log_name, "w");



        make_ss_data_templates();
        make_ss_qcd_templates(f_log);
        make_ss_mc_templates();
        fout->cd();
        gDirectory->cd(dirname);
        w->Write();
        fclose(f_log);
    }
    printf("Templates written to %s \n", fout_name.Data());
}





