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
#include "../TemplateMaker_systematics.C"
//#include "../TemplateMaker.C"
#include "../AFB_fit/FitUtils.C"

Double_t m_low;
Double_t m_high;

bool print = true;
const TString fout_name("combine/templates/jan22_test.root");
TFile * fout;
char dirname[10];


Double_t med_btag = 0.4432;

TH2F *h_elel_asym, *h_elel_sym, *h_elel_back,  *h_elel_data, *h_elel_mc, *h_elel_qcd;
TH1F *h1_elel_mn, *h1_elel_pl, *h1_elel_back,  *h1_elel_data, *h1_elel_mc, *h1_elel_qcd;
TH2F *h_elel_mc_count, *h_elel_sym_count;
TH2F *h_mumu_asym, *h_mumu_sym, *h_mumu_back,  *h_mumu_data, *h_mumu_mc, *h_mumu_qcd;
TH1F *h1_mumu_mn, *h1_mumu_pl, *h1_mumu_back,  *h1_mumu_data, *h1_mumu_mc, *h1_mumu_qcd;
TH2F *h_mumu_mc_count, *h_mumu_sym_count;
//void *x_axp, *y_axp, *z_axp;

//define axis globally for convinience 
TAxis *x_ax, *y_ax, *z_ax;



vector<double> v_elel_xF;
vector<double> v_elel_cost;
vector<double> v_mumu_xF;
vector<double> v_mumu_cost;
unsigned int nElEl_DataEvents;
unsigned int nMuMu_DataEvents;



void make_data_templates(){
    h_elel_data = new TH2F("ee_data_obs", "Data template of (x_f, cost_r)",
            n_xf_bins, xf_bins, n_cost_bins, cost_bins);
    h_elel_data->SetDirectory(0);
    h_mumu_data = new TH2F("mumu_data_obs", "Data template of (x_f, cost_r)",
            n_xf_bins, xf_bins, n_cost_bins, cost_bins);
    h_mumu_data->SetDirectory(0);

    nElEl_DataEvents = gen_data_template(t_elel_data, h_elel_data,  &v_elel_xF, &v_elel_cost, m_low, m_high, FLAG_ELECTRONS);
    nMuMu_DataEvents = gen_data_template(t_mumu_data, h_mumu_data, &v_mumu_xF, &v_mumu_cost, m_low, m_high, FLAG_MUONS);
    h1_elel_data = convert2d(h_elel_data);
    h1_mumu_data = convert2d(h_mumu_data);

    printf("Integral of data templates are %.2f %.2f \n", h1_elel_data->Integral(), h1_mumu_data->Integral()); 
    fout->cd();
    gDirectory->cd(dirname);
    h1_elel_data->Write();
    h1_mumu_data->Write();
    printf("Made data templates \n");
}

void make_qcd_templates(){
    h_elel_qcd = new TH2F((string("ee_qcd") ).c_str(), "Combined background template",
            n_xf_bins, xf_bins, n_cost_bins, cost_bins);
    h_elel_qcd->SetDirectory(0);
    h_mumu_qcd = new TH2F((string("mumu_qcd") ).c_str(), "Combined background template",
            n_xf_bins, xf_bins, n_cost_bins, cost_bins);
    h_mumu_qcd->SetDirectory(0);
    gen_fakes_template(t_elel_WJets, t_elel_QCD, t_elel_WJets_contam, t_elel_QCD_contam, h_elel_qcd, m_low, m_high, FLAG_ELECTRONS);
    gen_fakes_template(t_mumu_WJets, t_mumu_QCD, t_mumu_WJets_contam, t_mumu_QCD_contam, h_mumu_qcd, m_low, m_high, FLAG_MUONS);

    h_elel_qcd->Scale(1./h_elel_qcd->Integral());
    h_mumu_qcd->Scale(1./h_mumu_qcd->Integral());

    h1_elel_qcd = convert2d(h_elel_qcd);
    h1_mumu_qcd = convert2d(h_mumu_qcd);
    printf("Integral of QCD templates are %.2f %.2f \n", h1_elel_qcd->Integral(), h1_mumu_qcd->Integral());


    fout->cd();
    gDirectory->cd(dirname);
    h1_elel_qcd->Write();
    h1_mumu_qcd->Write();
    printf("Made qcd templates \n");
}

void make_mc_templates(const string &sys_label){
    h_elel_mc_count = new TH2F("h_elel_mc_count", "Events in bins for MC templates",
            n_xf_bins, xf_bins, n_cost_bins, cost_bins);
    h_elel_mc_count->SetDirectory(0);
    h_elel_sym_count = new TH2F("h_elel_sym_count", "Events in bins for MC templates",
            n_xf_bins, xf_bins, n_cost_bins, cost_bins);
    h_elel_sym_count->SetDirectory(0);
    h_elel_sym = new TH2F((string("ee_sym") + sys_label).c_str(), "Symmetric template of mc",
            n_xf_bins, xf_bins, n_cost_bins, cost_bins);
    h_elel_sym->SetDirectory(0);
    h_elel_asym = new TH2F((string("ee_asym") + sys_label).c_str(), "Asymmetric template of mc",
            n_xf_bins, xf_bins, n_cost_bins, cost_bins);
    h_elel_asym->SetDirectory(0);
    h_elel_back = new TH2F((string("ee_bk") + sys_label).c_str(), "Combined background template",
            n_xf_bins, xf_bins, n_cost_bins, cost_bins);
    h_elel_back->SetDirectory(0);






    h_mumu_mc_count = new TH2F("h_mumu_mc_count", "Events in bins for MC templates",
            n_xf_bins, xf_bins, n_cost_bins, cost_bins);
    h_mumu_mc_count->SetDirectory(0);
    h_mumu_sym_count = new TH2F("h_mumu_sym_count", "Events in bins for MC templates",
            n_xf_bins, xf_bins, n_cost_bins, cost_bins);
    h_mumu_sym_count->SetDirectory(0);
    h_mumu_sym = new TH2F((string("mumu_sym") + sys_label).c_str(), "Symmetric template of mc",
            n_xf_bins, xf_bins, n_cost_bins, cost_bins);
    h_mumu_sym->SetDirectory(0);
    h_mumu_asym = new TH2F((string("mumu_asym") + sys_label).c_str(), "Asymmetric template of mc",
            n_xf_bins, xf_bins, n_cost_bins, cost_bins);
    h_mumu_asym->SetDirectory(0);
    h_mumu_back = new TH2F((string("mumu_bk") + sys_label).c_str(), "Combined background template",
            n_xf_bins, xf_bins, n_cost_bins, cost_bins);
    h_mumu_back->SetDirectory(0);


    gen_mc_template(t_mumu_mc, alpha, h_mumu_sym, h_mumu_asym, h_mumu_sym_count, m_low, m_high, FLAG_MUONS);
    TTree *mumu_ts[2] = {t_mumu_back, t_mumu_nosig};
    gen_combined_background_template(2, mumu_ts, h_mumu_back, m_low, m_high, FLAG_MUONS);


    gen_mc_template(t_elel_mc, alpha, h_elel_sym, h_elel_asym, h_elel_sym_count, m_low, m_high, FLAG_ELECTRONS);
    TTree *elel_ts[2] = {t_elel_back, t_elel_nosig};

    gen_combined_background_template(2, elel_ts, h_elel_back, m_low, m_high, FLAG_ELECTRONS);


}



void convert_mc_templates(const string &sys_label){
    h1_mumu_back = convert2d(h_mumu_back);
    auto h_mumu_pl = *h_mumu_sym + *h_mumu_asym;
    auto h_mumu_mn = *h_mumu_sym - *h_mumu_asym;
    h_mumu_pl.Scale(0.5);
    h_mumu_mn.Scale(0.5);
    h1_mumu_pl = convert2d(&h_mumu_pl);
    h1_mumu_mn = convert2d(&h_mumu_mn);
    h1_mumu_pl->SetName((string("mumu_fpl") + sys_label).c_str());
    h1_mumu_mn->SetName((string("mumu_fmn") + sys_label).c_str());

    h1_elel_back = convert2d(h_elel_back);
    auto h_elel_pl = *h_elel_sym + *h_elel_asym;
    auto h_elel_mn = *h_elel_sym - *h_elel_asym;
    h_elel_pl.Scale(0.5);
    h_elel_mn.Scale(0.5);
    h1_elel_pl = convert2d(&h_elel_pl);
    h1_elel_mn = convert2d(&h_elel_mn);
    h1_elel_pl->SetName((string("ee_fpl") + sys_label).c_str());
    h1_elel_mn->SetName((string("ee_fmn") + sys_label).c_str());
}

void write_out_mc_templates(){

    fout->cd();
    gDirectory->cd(dirname);
    h1_mumu_back->Write();
    h1_mumu_pl->Write();
    h1_mumu_mn->Write();

    h1_elel_back->Write();
    h1_elel_pl->Write();
    h1_elel_mn->Write();

}



void make_templates(){

    init();

    vector<string> sys_labels {""};

    fout = TFile::Open(fout_name, "RECREATE");



    for(int i=0; i<n_m_bins; i++){
    //for(int i=0; i<1; i++){
        fout->cd();
        snprintf(dirname, 10, "%i", i);
        gDirectory->mkdir(dirname);
        gDirectory->cd(dirname);

        m_low = m_bins[i];
        m_high = m_bins[i+1];
        alpha = alphas[i];
        alpha = alphas[i] + alpha_unc[i];
        printf("Start making templates for mass bin %.0f-%.0f \n", m_low, m_high);

        make_data_templates();
        make_qcd_templates();
        for(auto iter = sys_labels.begin(); iter !=sys_labels.end(); iter++){
            printf("Making MC templates for sys %s \n", (*iter).c_str());

            make_mc_templates(*iter);
            convert_mc_templates(*iter);
            write_out_mc_templates();
        }
    }


    printf("Templates written to %s \n", fout_name.Data());

}

