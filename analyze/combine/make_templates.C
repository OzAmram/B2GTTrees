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
#include "make_ss_templates.C"
#include "make_emu_templates.C"
#include "TemplateUtils.h"




TH2F *h_elel_asym, *h_elel_sym, *h_elel_back,  *h_elel_dy_gg, *h_elel_data, *h_elel_mc, *h_elel_qcd;
TH1F *h1_elel_mn, *h1_elel_pl, *h1_elel_back,  *h1_elel_dy_gg, *h1_elel_data, *h1_elel_mc, *h1_elel_qcd;
TH2F *h_elel_mc_count, *h_elel_sym_count;
TH2F *h_mumu_asym, *h_mumu_sym, *h_mumu_back,  *h_mumu_dy_gg, *h_mumu_data, *h_mumu_mc, *h_mumu_qcd;
TH1F *h1_mumu_mn, *h1_mumu_pl, *h1_mumu_back,  *h1_mumu_dy_gg, *h1_mumu_data, *h1_mumu_mc, *h1_mumu_qcd;




vector<double> v_elel_xF;
vector<double> v_elel_cost;
vector<double> v_mumu_xF;
vector<double> v_mumu_cost;
unsigned int nElEl_DataEvents;
unsigned int nMuMu_DataEvents;

void convert_qcd_to_param_hist(TH2F *h, FILE *f_log, float sign_scaling, int flag){
    //convert a hist to a parametric hist 
    RooArgList *bin_list = new RooArgList();
    RooArgList *bin_list_os = new RooArgList();
    RooArgList *bin_list_ss = new RooArgList();

    TH1F *h1 = convert2d(h);

    char h_name[40];
    char h_ss_name[40];
    char R_sign_param[40];
    sprintf(h_name, "%s_param", h->GetName());
    RooRealVar *R_qcd_sign_fraction;
    fprintf(f_log, "\n");
    if(flag == FLAG_MUONS){
        sprintf(h_ss_name, "%s", "mumu_ss_qcd_param");
        sprintf(R_sign_param, "%s", "R_mumu_os_fakes");
    }
    else{
        sprintf(h_ss_name, "%s", "ee_ss_qcd_param");
        sprintf(R_sign_param, "%s", "R_ee_os_fakes");
        R_qcd_sign_fraction = new RooRealVar(R_sign_param, "Fraction of os fakes events", sign_scaling , 0., 1.);
        fprintf(f_log, "%s param %.4f 0.05 \n", R_sign_param, sign_scaling);
    }
    for(int i=1; i <= n_xf_bins; i++){
        for(int j=1; j <= n_cost_bins; j++){



            int g_idx = TwoDToOneDIdx(n_cost_bins, i, j);
            int sym1_idx, sym2_idx;
            TwoDToSymIdxs(n_cost_bins, i,j, sym1_idx, sym2_idx);
            printf("i,j: %i %i ", i,j);
            printf("g_idx, sym1, sym2: %i %i %i  \n", g_idx, sym1_idx, sym2_idx);

            double content = h1->GetBinContent(g_idx);
            double error = h1->GetBinError(g_idx);
            if(content<0) printf("Bin %i Content is %.0f \n", j, content);

            //printf("Bin %.1f error %.1f \n", content,error);
            char bin_name[40];
            char form_name_ss[40], form_name1_os[40], form_name2_os[40];
            sprintf(bin_name, "%s_bin%i",h_name, g_idx); 
            sprintf(form_name_ss, "%s_form_%i",h_ss_name, g_idx); 
            sprintf(form_name1_os, "%s_form_%i",h_name, sym1_idx); 
            sprintf(form_name2_os, "%s_form_%i",h_name, sym2_idx); 
            //prevent underflowing by fixing super small bins
            content = max(content, 0.001);
            if (content < error){
                content = error/2.;
                error = 0.1*content;
            }
            else if(content < 2.5 * error){
                error = 0.3*content;
            }
            if(j<=(n_cost_bins/2)){
                //printf("first fill \n");
                RooRealVar *bin = new RooRealVar(bin_name, bin_name, content, 0., 10000.);
                fprintf(f_log, "%s param %.4f %.4f \n", bin_name, content, error);
                RooFormulaVar *form1 = new RooFormulaVar(form_name1_os, form_name1_os, "0.5*@0*@1", RooArgList(*bin, *R_qcd_sign_fraction));
                RooFormulaVar *form_ss = new RooFormulaVar(form_name_ss, form_name_ss, "@0*(1.0 - @1)", RooArgList(*bin, *R_qcd_sign_fraction));
                //form1->Print();
                //form_ss->Print();
                bin_list->add(*bin);
                bin_list_ss->add(*form_ss);
                bin_list_os->add(*form1);
            }

            else{
                //printf("2nd fill \n");
                int old_j = sym2_idx % n_cost_bins;
                int old_g_idx = TwoDToOneDIdx(n_cost_bins, i, old_j);
                sprintf(bin_name, "%s_bin%i",h_name, old_g_idx); 
                printf("Looking for bin %s \n", bin_name);
                RooRealVar *bin = (RooRealVar *) bin_list->find(bin_name);
                if(bin==nullptr) printf("NULL lookup of %s from bin list \n", bin_name);
                RooFormulaVar *form1 = new RooFormulaVar(form_name1_os, form_name1_os, "0.5*@0*@1", RooArgList(*bin, *R_qcd_sign_fraction));
                //form1->Print();
                bin_list_os->add(*form1);
            }
        }
    
    }
    bin_list_ss->Print();
    bin_list_os->Print();
    char norm_ss_name[40], norm_name[40];
    sprintf(norm_ss_name, "%s_norm", h_ss_name);
    sprintf(norm_name, "%s_norm", h_name);
    RooAddition *norm_ss = new RooAddition(norm_ss_name, norm_ss_name, *bin_list_ss);

    RooAddition *norm = new RooAddition(norm_name, norm_name, *bin_list_os);

    RooParametricHist *p= new RooParametricHist (h_name, h_name, *var, *bin_list_os, *h_dummy);
    RooParametricHist *p_ss= new RooParametricHist (h_ss_name, h_ss_name, *var_ss, *bin_list_ss, *h1);

    
    w->import(*p_ss);
    w->import(*p, RooFit::RecycleConflictNodes());
    w->import(*norm,RooFit::RecycleConflictNodes());
    w->import(*norm_ss,RooFit::RecycleConflictNodes());
}


void make_data_templates(){
    h_elel_data = new TH2F("ee_data_obs", "Data template of (x_f, cost_r)",
            n_xf_bins, xf_bins, n_cost_bins, cost_bins);
    h_elel_data->SetDirectory(0);
    h_mumu_data = new TH2F("mumu_data_obs", "Data template of (x_f, cost_r)",
            n_xf_bins, xf_bins, n_cost_bins, cost_bins);
    h_mumu_data->SetDirectory(0);

    nElEl_DataEvents = gen_data_template(t_elel_data, h_elel_data,  m_low, m_high, FLAG_ELECTRONS, FLAG_M_BINS, do_RC);
    nMuMu_DataEvents = gen_data_template(t_mumu_data, h_mumu_data,  m_low, m_high, FLAG_MUONS, FLAG_M_BINS, do_RC);
    h1_elel_data = convert2d(h_elel_data);
    h1_mumu_data = convert2d(h_mumu_data);
    

    printf("Integral of data templates are %.2f %.2f \n", h1_elel_data->Integral(), h1_mumu_data->Integral()); 
    write_roo_hist(h1_elel_data, var);
    write_roo_hist(h1_mumu_data, var);
    printf("Made data templates \n");
}

void make_qcd_templates(FILE* f_log){
    h_elel_qcd = new TH2F((string("ee_qcd") ).c_str(), "Combined background template",
            n_xf_bins, xf_bins, n_cost_ss_bins, cost_ss_bins);
    h_elel_qcd->SetDirectory(0);
    h_mumu_qcd = new TH2F((string("mumu_qcd") ).c_str(), "Combined background template",
            n_xf_bins, xf_bins, n_cost_ss_bins, cost_ss_bins);
    h_mumu_qcd->SetDirectory(0);
    bool ss = true;
    float elel_sign_scaling = gen_fakes_template(t_elel_WJets, t_elel_QCD, t_elel_WJets_contam, t_elel_QCD_contam, h_elel_qcd, m_low, m_high, FLAG_ELECTRONS, FLAG_M_BINS, ss);
    float mumu_sign_scaling = gen_fakes_template(t_mumu_WJets, t_mumu_QCD, t_mumu_WJets_contam, t_mumu_QCD_contam, h_mumu_qcd, m_low, m_high, FLAG_MUONS, FLAG_M_BINS,ss);

    //combined os and ss regions to estimate qcd, scale it to estimate amount
    //in os region 
    //float scaling = 1./(1. + R_mu_ss_os);
    //h_mumu_qcd->Scale(scaling);

    printf("Integral of QCD templates are %.2f %.2f \n", h_elel_qcd->Integral(), h_mumu_qcd->Integral());

    convert_qcd_to_param_hist(h_elel_qcd, f_log, elel_sign_scaling, FLAG_ELECTRONS);
    convert_qcd_to_param_hist(h_mumu_qcd, f_log, mumu_sign_scaling, FLAG_MUONS);

    printf("Made qcd templates \n");
}

void make_mc_templates(const string &sys_label){
    bool do_mu, do_el;
    if(sys_label.find("mu") != string::npos){
        printf("Doing mu only \n");
        do_mu = true;
        do_el = false;
    }
    if(sys_label.find("el") != string::npos){
        printf("Doing el only \n");
        do_mu = false;
        do_el = true;
    }
    else{
        do_mu = true;
        do_el = true;
    }
    bool ss= false;
    
    if(do_mu){
        printf("making muon mc templates \n");
        h_mumu_sym = new TH2F((string("mumu_sym") + sys_label).c_str(), "Symmetric template of mc",
                n_xf_bins, xf_bins, n_cost_bins, cost_bins);
        h_mumu_sym->SetDirectory(0);
        h_mumu_asym = new TH2F((string("mumu_asym") + sys_label).c_str(), "Asymmetric template of mc",
                n_xf_bins, xf_bins, n_cost_bins, cost_bins);
        h_mumu_asym->SetDirectory(0);
        h_mumu_back = new TH2F((string("mumu_bk") + sys_label).c_str(), "Combined background template",
                n_xf_bins, xf_bins, n_cost_bins, cost_bins);
        h_mumu_back->SetDirectory(0);
        h_mumu_dy_gg = new TH2F((string("mumu_dy_gg") + sys_label).c_str(), "Combined background template",
                n_xf_bins, xf_bins, n_cost_bins, cost_bins);
        h_mumu_dy_gg->SetDirectory(0);

        gen_mc_template(t_mumu_mc, alpha, h_mumu_sym, h_mumu_asym, m_low, m_high, FLAG_MUONS, FLAG_M_BINS, do_RC, sys_label );
        TTree *mumu_ts[1] = {t_mumu_back};
        gen_combined_background_template(1, mumu_ts, h_mumu_back, m_low, m_high, FLAG_MUONS, FLAG_M_BINS, do_RC, ss, sys_label);
        mumu_ts[0] = t_mumu_nosig;
        gen_combined_background_template(1, mumu_ts, h_mumu_dy_gg, m_low, m_high, FLAG_MUONS, FLAG_M_BINS, do_RC, ss, sys_label);
    }
    if(do_el){
        printf("making electron mc templates \n");
        h_elel_sym = new TH2F((string("ee_sym") + sys_label).c_str(), "Symmetric template of mc",
                n_xf_bins, xf_bins, n_cost_bins, cost_bins);
        h_elel_sym->SetDirectory(0);
        h_elel_asym = new TH2F((string("ee_asym") + sys_label).c_str(), "Asymmetric template of mc",
                n_xf_bins, xf_bins, n_cost_bins, cost_bins);
        h_elel_asym->SetDirectory(0);
        h_elel_back = new TH2F((string("ee_bk") + sys_label).c_str(), "Combined background template",
                n_xf_bins, xf_bins, n_cost_bins, cost_bins);
        h_elel_back->SetDirectory(0);
        h_elel_dy_gg = new TH2F((string("ee_dy_gg") + sys_label).c_str(), "Combined background template",
                n_xf_bins, xf_bins, n_cost_bins, cost_bins);
        h_elel_dy_gg->SetDirectory(0);

        gen_mc_template(t_elel_mc, alpha, h_elel_sym, h_elel_asym,  m_low, m_high, FLAG_ELECTRONS, FLAG_M_BINS, do_RC, sys_label);
        TTree *elel_ts[1] = {t_elel_back};
        gen_combined_background_template(1, elel_ts, h_elel_back, m_low, m_high, FLAG_ELECTRONS, FLAG_M_BINS,do_RC,ss, sys_label);
        elel_ts[0] = t_elel_nosig;
        gen_combined_background_template(1, elel_ts, h_elel_dy_gg, m_low, m_high, FLAG_ELECTRONS, FLAG_M_BINS,do_RC,ss, sys_label);
    }

}



void convert_mc_templates(const string &sys_label){
    bool do_mu, do_el;
    if(sys_label.find("mu") != string::npos){
        do_mu = true;
        do_el = false;
    }
    else if(sys_label.find("el") != string::npos){
        do_mu = false;
        do_el = true;
    }
    else{
        do_mu = true;
        do_el = true;
    }
    if(do_mu){
        h1_mumu_back = convert2d(h_mumu_back);
        h1_mumu_dy_gg = convert2d(h_mumu_dy_gg);
        auto h_mumu_pl = *h_mumu_sym + *h_mumu_asym;
        auto h_mumu_mn = *h_mumu_sym - *h_mumu_asym;
        h_mumu_pl.Scale(0.5);
        h_mumu_mn.Scale(0.5);
        h1_mumu_pl = convert2d(&h_mumu_pl);
        h1_mumu_mn = convert2d(&h_mumu_mn);
        h1_mumu_pl->SetName((string("mumu_fpl") + sys_label).c_str());
        h1_mumu_mn->SetName((string("mumu_fmn") + sys_label).c_str());

        write_roo_hist(h1_mumu_back, var);
        write_roo_hist(h1_mumu_dy_gg, var);
        write_roo_hist(h1_mumu_pl, var);
        write_roo_hist(h1_mumu_mn, var);
    }

    if(do_el){
        h1_elel_back = convert2d(h_elel_back);
        h1_elel_dy_gg = convert2d(h_elel_dy_gg);
        auto h_elel_pl = *h_elel_sym + *h_elel_asym;
        auto h_elel_mn = *h_elel_sym - *h_elel_asym;
        h_elel_pl.Scale(0.5);
        h_elel_mn.Scale(0.5);
        h1_elel_pl = convert2d(&h_elel_pl);
        h1_elel_mn = convert2d(&h_elel_mn);
        h1_elel_pl->SetName((string("ee_fpl") + sys_label).c_str());
        h1_elel_mn->SetName((string("ee_fmn") + sys_label).c_str());


        write_roo_hist(h1_elel_back, var);
        write_roo_hist(h1_elel_dy_gg, var);
        write_roo_hist(h1_elel_pl, var);
        write_roo_hist(h1_elel_mn, var);
    }
}
void write_groups(FILE *f_log){

    char label[4][40], intro[4][40];
    int sizes[4] = {40,20,20,60};
    sprintf(intro[0], "emu_fake_shape group = ");
    sprintf(intro[1], "ee_fake_shape group = ");
    sprintf(intro[2], "mumu_fake_shape group = ");
    sprintf(intro[3], "pdfs group = ");

    sprintf(label[0], "emu_qcd_param_bin");
    sprintf(label[1], "ee_qcd_param_bin");
    sprintf(label[2], "mumu_qcd_param_bin");
    sprintf(label[3], "pdf");

    for(int i=0; i<4; i++){
        fprintf(f_log, "\n %s", intro[i]);
        for(int j=1; j<=sizes[i]; j++){

            fprintf(f_log, " %s%i", label[i], j);
        }
    }
        fprintf(f_log, "\n");
}




void make_templates(int nJobs = 6, int iJob =-1){
    const TString fout_name("combine/templates/april22_no_sys.root");
    TFile * fout;

    init();
    init_ss();
    init_emu();
    printf("Setting up SFs... ");
    setup_all_SFs();
    printf("   done \n");

    vector<string> sys_labels {""};
    //vector<string> sys_labels {"_FACUp", "_pdfDown", "_RENORMDown"};
    
    
    /*
    vector<string> sys_labels {"_elScaleUp", "_elScaleDown", "_elSmearUp", "_elSmearDown", 
        "_muRCUp", "_muRCDown", "_muHLTUp", "_muHLTDown", "_muIDUp", "_muIDDown", "_muISOUp", "_muISODown", "_muTRKUp", "_muTRKDown",  
        "_elHLTUp", "_elHLTDown", "_elIDUp", "_elIDDown", "_elRECOUp", "_elRECODown", 
        "_RENORMUp", "_RENORMDown", "_FACUp", "_FACDown",
        "_PuUp", "_PuDown", "_BTAGUp", "_BTAGDown", "_alphaUp", "_alphaDown", "_alphaSUp", "_alphaSDown" };
        */
        
        

    fout = TFile::Open(fout_name, "RECREATE");
    FILE *f_log;
    char f_log_name[80];

    char dirname[40];

    int i_start=0;
    int i_max = n_m_bins;
    if(iJob >0){
        i_start =iJob;
        i_max = iJob +1;
    }


    for(int i=i_start; i<i_max; i++){
    //for(int i=0; i<1; i++){
        fout->cd();
        snprintf(dirname, 10, "w%i", i);
        gDirectory->mkdir(dirname);
        gDirectory->cd(dirname);
        w = new RooWorkspace("w", "w");
        sprintf(f_log_name, "combine/AFB_fits/cards/mbin%i_bins.txt", i);
        f_log = fopen(f_log_name, "w");

        m_low = m_bins[i];
        m_high = m_bins[i+1];
        printf("\n \n Start making templates for mass bin %.0f-%.0f \n", m_low, m_high);

        make_data_templates();
        make_ss_data_templates();
        make_ss_mc_templates();
        make_emu_data_templates();
        make_emu_qcd_templates(f_log);
        make_emu_mc_templates();
        make_qcd_templates(f_log);
        for(auto iter = sys_labels.begin(); iter !=sys_labels.end(); iter++){
            printf("Making MC templates for sys %s \n", (*iter).c_str());
            alpha = alphas[i];
            if(iter->find("alphaUp") != string::npos) alpha = alphas[i] + alpha_unc[i];
            if(iter->find("alphaDown") != string::npos) alpha = alphas[i] - alpha_unc[i];

            make_mc_templates(*iter);
            convert_mc_templates(*iter);
        }
        fout->cd();
        gDirectory->cd(dirname);
        w->Write();
        write_groups(f_log);
        fclose(f_log);
    }


    fout->Close();
    printf("Templates written to %s \n", fout_name.Data());

}

