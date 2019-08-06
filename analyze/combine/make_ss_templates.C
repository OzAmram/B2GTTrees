
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
//#include"Minuit2/Minuit2Minimizer.h"
#include "Math/Functor.h"
#include "TemplateUtils.h"


const TString ss_fout_name("combine/templates/feb12_ss_qcd_fakerate_param.root");
TFile * ss_fout;




void convert_ss_qcd_to_param_hist(TH1F *h, FILE *f_log, int flag){
    //convert a hist to a parametric hist 
    RooArgList *bin_list = new RooArgList();

    char h_name[40];
    sprintf(h_name, "%s_param", h->GetName());
    RooRealVar *Rqcd;
    if(flag == FLAG_MUONS) Rqcd = Rqcd_mumu_ss;
    else Rqcd = Rqcd_ee_ss;
    for(int j=1; j <= h->GetNbinsX(); j++){

        double content = h->GetBinContent(j);
        if(content<0) printf("Bin %i Content is %.0f \n", j, content);
        double error = h->GetBinError(j);
        //printf("Bin %.1f error %.1f \n", content,error);
        //prevent underflowing
        content = max(content, error);
        char bin_name[40];
        char form_name[40];
        sprintf(bin_name, "%s_bin%i",h_name, j); 
        sprintf(form_name, "%s_form%i",h_name, j); 
        RooRealVar *bin = new RooRealVar(bin_name, bin_name, content, 0., 10000.);
        fprintf(f_log, "%s param %.2f %.2f \n", bin_name, content, error);
        //bin->Print();
        bin_list->add(*bin);
        //RooFormulaVar *form = new RooFormulaVar(form_name, form_name, "@0*@1", RooArgList(*bin, *Rqcd));
        //bin_list->add(*form);
    }
    bin_list->Print();

    RooParametricHist *p= new RooParametricHist (h_name, h_name, *var, *bin_list, *h);
    char norm_name[40];
    sprintf(norm_name, "%s_norm", h_name);
    RooAddition *n = new RooAddition(norm_name, norm_name, *bin_list);

    
    w->import(*p);
    w->import(*n,RooFit::RecycleConflictNodes());
}


void make_ss_data_templates(int year){

    char title[100];
    sprintf(title, "ee%i_ss_data_obs", year %2000);
    auto h_elel_data = new TH2F(title, "Data template of (x_f, cost_r)",
            n_xf_bins, xf_bins, n_cost_ss_bins, cost_ss_bins);
    h_elel_data->SetDirectory(0);
    sprintf(title, "mumu%i_ss_data_obs", year %2000);
    auto h_mumu_data = new TH2F("mumu_ss_data_obs", "Data template of (x_f, cost_r)",
            n_xf_bins, xf_bins, n_cost_ss_bins, cost_ss_bins);
    h_mumu_data->SetDirectory(0);
    bool ss = true;

    gen_data_template(t_elel_ss_data, h_elel_data,  year, m_low, m_high, FLAG_ELECTRONS, do_RC, ss);
    gen_data_template(t_mumu_ss_data, h_mumu_data,  year, m_low, m_high, FLAG_MUONS,  do_RC, ss);
    auto h1_elel_data = convert2d(h_elel_data);
    auto h1_mumu_data = convert2d(h_mumu_data);


    printf("Integral of data templates were %.2f %.2f \n", h_elel_data->Integral(), h_mumu_data->Integral()); 
    //h_elel_data->Write();
    //h_mumu_data->Write();
    write_roo_hist(h1_elel_data, var_ss);
    write_roo_hist(h1_mumu_data, var_ss);
    printf("Made ss data templates \n");
}

void make_ss_qcd_templates(int year, FILE *f_log){
    char title[100];
    sprintf(title, "ee%i_ss_qcd", year %2000);
    auto h_elel_qcd = new TH2F(title, "Combined background template",
            n_xf_bins, xf_bins, n_cost_ss_bins, cost_ss_bins);
    h_elel_qcd->SetDirectory(0);
    sprintf(title, "mumu%i_ss_qcd", year %2000);
    auto h_mumu_qcd = new TH2F(title, "Combined background template",
            n_xf_bins, xf_bins, n_cost_ss_bins, cost_ss_bins);
    h_mumu_qcd->SetDirectory(0);


    bool ss = true;
    gen_fakes_template(t_elel_WJets, t_elel_QCD, t_elel_WJets_contam, t_elel_QCD_contam, h_elel_qcd, year, m_low, m_high, FLAG_ELECTRONS,  ss);
    gen_fakes_template(t_mumu_WJets, t_mumu_QCD, t_mumu_WJets_contam, t_mumu_QCD_contam, h_mumu_qcd, year, m_low, m_high, FLAG_MUONS, ss);
    printf("Integral of qcd templates are %.2f %.2f \n", h_elel_qcd->Integral(), h_mumu_qcd->Integral()); 
    auto h1_elel_qcd = convert2d(h_elel_qcd);
    auto h1_mumu_qcd = convert2d(h_mumu_qcd);

    convert_ss_qcd_to_param_hist(h1_elel_qcd, f_log, FLAG_ELECTRONS);
    convert_ss_qcd_to_param_hist(h1_mumu_qcd, f_log, FLAG_MUONS);
    printf("Made ss qcd templates \n");
}

void make_ss_mc_templates(int year){
    char title[100];
    sprintf(title, "ee%i_ss_dy", year %2000);
    auto h_elel_dy = new TH2F(title, "ss DY",
            n_xf_bins, xf_bins, n_cost_ss_bins, cost_ss_bins);
    h_elel_dy->SetDirectory(0);
    sprintf(title, "mumu%i_ss_dy", year %2000);
    auto h_mumu_dy = new TH2F(title, "ss DY",
            n_xf_bins, xf_bins, n_cost_ss_bins, cost_ss_bins);

    sprintf(title, "ee%i_ss_bk", year %2000);
    auto h_elel_bk = new TH2F(title, "Combined background template",
            n_xf_bins, xf_bins, n_cost_ss_bins, cost_ss_bins);
    h_elel_bk->SetDirectory(0);
    sprintf(title, "mumu%i_ss_bk", year %2000);
    auto h_mumu_bk = new TH2F((string("mumu_ss_bk") ).c_str(), "Combined background template",
            n_xf_bins, xf_bins, n_cost_ss_bins, cost_ss_bins);

    bool ss = true;


    TTree *mumu_ts[1] = {t_mumu_ss_back};
    gen_combined_background_template(1, mumu_ts, h_mumu_bk, year, m_low, m_high, FLAG_MUONS,  do_RC, ss);
    mumu_ts[0] = t_mumu_ss_dy;
    gen_combined_background_template(1, mumu_ts, h_mumu_dy, year, m_low, m_high, FLAG_MUONS,  do_RC, ss);

    TTree *elel_ts[1] = {t_elel_ss_back};
    gen_combined_background_template(1, elel_ts, h_elel_bk, year, m_low, m_high, FLAG_ELECTRONS,  do_RC, ss);
    elel_ts[0] = t_elel_ss_dy;
    gen_combined_background_template(1, elel_ts, h_elel_dy, year, m_low, m_high, FLAG_ELECTRONS,  do_RC, ss);

    auto h1_elel_bk = convert2d(h_elel_bk);
    auto h1_mumu_bk = convert2d(h_mumu_bk);

    auto h1_elel_dy = convert2d(h_elel_dy);
    auto h1_mumu_dy = convert2d(h_mumu_dy);




    printf("Integral of dy templates are %.2f %.2f \n", h_elel_dy->Integral(), h_mumu_dy->Integral()); 
    printf("Integral of bkg templates are %.2f %.2f \n", h_elel_bk->Integral(), h_mumu_bk->Integral()); 
    write_roo_hist(h1_elel_dy, var_ss);
    write_roo_hist(h1_mumu_dy, var_ss);
    printf("Made ss dy templates \n");

    write_roo_hist(h1_elel_bk, var_ss);
    write_roo_hist(h1_mumu_bk, var_ss);
    printf("Made ss bk templates \n");


}





void make_ss_templates(int year=2016){

    init(year);
    init_ss(year);


    ss_fout = TFile::Open(ss_fout_name, "RECREATE");
    FILE *f_log;
    char f_log_name[80];
    char dirname[40];


    for(int i=0; i<n_m_bins; i++){
        ss_fout->cd();
        snprintf(dirname, 10, "w%i", i);
        gDirectory->mkdir(dirname);
        gDirectory->cd(dirname);
        w = new RooWorkspace("w", "w");

        m_low = m_bins[i];
        m_high = m_bins[i+1];

        sprintf(f_log_name, "combine/SameSign_fits/cards/y%i_mbin%i_bins.txt", year, i);
        f_log = fopen(f_log_name, "w");



        make_ss_data_templates(year);
        make_ss_qcd_templates(year, f_log);
        make_ss_mc_templates(year);
        ss_fout->cd();
        gDirectory->cd(dirname);
        w->Write();
        fclose(f_log);
    }
    printf("Templates written to %s \n", ss_fout_name.Data());
}





