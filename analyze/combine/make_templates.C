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




TH2F *h_elel_asym, *h_elel_sym, *h_elel_alpha, *h_elel_back,  *h_elel_dy_gg, *h_elel_data, *h_elel_mc, *h_elel_qcd, *h_elel_gam;
TH2F *h_elel_mc_count, *h_elel_sym_count;
TH2F *h_mumu_asym, *h_mumu_sym, *h_mumu_alpha, *h_mumu_back,  *h_mumu_dy_gg, *h_mumu_data, *h_mumu_mc, *h_mumu_qcd, *h_mumu_gam;





void convert_qcd_to_param_hist(TH2F *h, FILE *f_log, float sign_scaling, int flag){
    //convert a hist to a parametric hist 
    RooArgList *bin_list = new RooArgList();
    RooArgList *bin_list_os = new RooArgList();
    RooArgList *bin_list_ss = new RooArgList();

    TH1F *h1 = convert2d(h);

    char h_name[40];
    char h_ss_name[40];
    char R_sign_param[40];
    sprintf(h_name, "%s_qcd_param", h->GetName());
    RooRealVar *R_qcd_sign_fraction;
    fprintf(f_log, "\n");
    sprintf(h_ss_name, "%s_ss_qcd_param", h->GetName());
    sprintf(R_sign_param, "R_%s_os_fakes", h->GetName());
    if(flag == FLAG_ELECTRONS){
        R_qcd_sign_fraction = new RooRealVar(R_sign_param, "Fraction of os fakes events", sign_scaling , 0., 1.);
        fprintf(f_log, "%s param %.4f 0.05 \n", R_sign_param, sign_scaling);
    }
    for(int i=1; i <= n_xf_bins; i++){
        for(int j=1; j <= n_cost_bins; j++){



            int g_idx = TwoDToOneDIdx(n_cost_bins, i, j);
            int sym1_idx, sym2_idx;
            TwoDToSymIdxs(n_cost_bins, i,j, sym1_idx, sym2_idx);
            //printf("i,j: %i %i ", i,j);
            //printf("g_idx, sym1, sym2: %i %i %i  \n", g_idx, sym1_idx, sym2_idx);

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
                //printf("Looking for bin %s \n", bin_name);
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


void make_data_templates(int year){
    char title[100];
    sprintf(title, "ee%i_data_obs", year %2000);
    h_elel_data = new TH2F(title, "Data template of (x_f, cost_r)",
            n_xf_bins, xf_bins, n_cost_bins, cost_bins);
    h_elel_data->SetDirectory(0);
    sprintf(title, "mumu%i_data_obs", year %2000);
    h_mumu_data = new TH2F(title, "Data template of (x_f, cost_r)",
            n_xf_bins, xf_bins, n_cost_bins, cost_bins);
    h_mumu_data->SetDirectory(0);

    int nElEl_DataEvents = gen_data_template(t_elel_data, h_elel_data,  year, m_low, m_high, FLAG_ELECTRONS,  do_RC);
    int nMuMu_DataEvents = gen_data_template(t_mumu_data, h_mumu_data,  year, m_low, m_high, FLAG_MUONS, do_RC);
    auto h1_elel_data = convert2d(h_elel_data);
    auto h1_mumu_data = convert2d(h_mumu_data);
    

    printf("Integral of data templates are %.2f %.2f \n", h1_elel_data->Integral(), h1_mumu_data->Integral()); 
    write_roo_hist(h1_elel_data, var);
    write_roo_hist(h1_mumu_data, var);
    printf("Made data templates \n");
}

void make_qcd_templates(int year, FILE* f_log){
    char title[100];
    //titles taken care of in conversion

    sprintf(title, "ee%i", year %2000);
    h_elel_qcd = new TH2F(title, "Combined background template",
            n_xf_bins, xf_bins, n_cost_ss_bins, cost_ss_bins);
    h_elel_qcd->SetDirectory(0);
    sprintf(title, "mumu%i", year %2000);
    h_mumu_qcd = new TH2F(title, "Combined background template",
            n_xf_bins, xf_bins, n_cost_ss_bins, cost_ss_bins);
    h_mumu_qcd->SetDirectory(0);
    bool ss = true;
    float elel_sign_scaling = gen_fakes_template(t_elel_WJets, t_elel_QCD, t_elel_WJets_contam, t_elel_QCD_contam, h_elel_qcd, year, m_low, m_high, FLAG_ELECTRONS, ss);
    float mumu_sign_scaling = gen_fakes_template(t_mumu_WJets, t_mumu_QCD, t_mumu_WJets_contam, t_mumu_QCD_contam, h_mumu_qcd, year, m_low, m_high, FLAG_MUONS, ss);

    //combined os and ss regions to estimate qcd, scale it to estimate amount
    //in os region 
    //float scaling = 1./(1. + R_mu_ss_os);
    //h_mumu_qcd->Scale(scaling);

    printf("Integral of QCD templates are %.2f %.2f \n", h_elel_qcd->Integral(), h_mumu_qcd->Integral());

    convert_qcd_to_param_hist(h_elel_qcd, f_log, elel_sign_scaling, FLAG_ELECTRONS);
    convert_qcd_to_param_hist(h_mumu_qcd, f_log, mumu_sign_scaling, FLAG_MUONS);

    printf("Made qcd templates \n");
}

void cleanup_mc_templates(){
    delete h_elel_back; 
    delete h_elel_dy_gg; 
    delete h_elel_gam; 
    delete h_elel_sym; 
    delete h_elel_asym; 
    delete h_elel_alpha;
    delete h_mumu_back; 
    delete h_mumu_dy_gg; 
    delete h_mumu_gam; 
    delete h_mumu_sym; 
    delete h_mumu_asym; 
    delete h_mumu_alpha;
}

void make_mc_templates(int year, Double_t alpha_denom, const string &sys_label){
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
    
    char title[100];
    if(do_mu){
        printf("making muon mc templates \n");
        sprintf(title, "mumu%i_sym%s", year %2000, sys_label.c_str());
        h_mumu_sym = new TH2F(title, "Symmetric template of mc",
                n_xf_bins, xf_bins, n_cost_bins, cost_bins);
        h_mumu_sym->SetDirectory(0);
        sprintf(title, "mumu%i_alpha%s", year %2000, sys_label.c_str());
        h_mumu_alpha = new TH2F(title, "Gauge boson polarization template of mc",
                n_xf_bins, xf_bins, n_cost_bins, cost_bins);
        h_mumu_alpha->SetDirectory(0);
        sprintf(title, "mumu%i_asym%s", year %2000, sys_label.c_str());
        h_mumu_asym = new TH2F(title, "Asymmetric template of mc",
                n_xf_bins, xf_bins, n_cost_bins, cost_bins);
        h_mumu_asym->SetDirectory(0);
        sprintf(title, "mumu%i_bk%s", year %2000, sys_label.c_str());
        h_mumu_back = new TH2F(title, "Combined background template",
                n_xf_bins, xf_bins, n_cost_bins, cost_bins);
        h_mumu_back->SetDirectory(0);
        sprintf(title, "mumu%i_dy_gg%s", year %2000, sys_label.c_str());
        h_mumu_dy_gg = new TH2F(title, "Combined background template",
                n_xf_bins, xf_bins, n_cost_bins, cost_bins);
        h_mumu_dy_gg->SetDirectory(0);
        sprintf(title, "mumu%i_gam%s", year %2000, sys_label.c_str());
        h_mumu_gam = new TH2F(title, "Combined background template",
                n_xf_bins, xf_bins, n_cost_bins, cost_bins);
        h_mumu_gam->SetDirectory(0);

        printf("Making mumu mc \n");
        gen_mc_template(t_mumu_mc, alpha_denom, h_mumu_sym, h_mumu_asym, h_mumu_alpha, year, m_low, m_high, FLAG_MUONS, do_RC, sys_label );
        TTree *mumu_ts[1] = {t_mumu_back};
        printf("Making mumu back \n");
        gen_combined_background_template(1, mumu_ts, h_mumu_back, year, m_low, m_high, FLAG_MUONS,  do_RC, ss, sys_label);
        mumu_ts[0] = t_mumu_nosig;
        printf("Making mumu nosig \n");
        gen_combined_background_template(1, mumu_ts, h_mumu_dy_gg, year, m_low, m_high, FLAG_MUONS,  do_RC, ss, sys_label);

        mumu_ts[0] = t_mumu_gamgam;
        printf("Making mumu gamgam \n");
        gen_combined_background_template(1, mumu_ts, h_mumu_gam, year, m_low, m_high, FLAG_MUONS,  do_RC, ss, sys_label);
    }
    if(do_el){
        printf("making electron mc templates \n");
        sprintf(title, "ee%i_sym%s", year %2000, sys_label.c_str());
        h_elel_sym = new TH2F(title, "Symmetric template of mc",
                n_xf_bins, xf_bins, n_cost_bins, cost_bins);
        h_elel_sym->SetDirectory(0);
        sprintf(title, "ee%i_alpha%s", year %2000, sys_label.c_str());
        h_elel_alpha = new TH2F(title, "Gauge boson polarization template of mc",
                n_xf_bins, xf_bins, n_cost_bins, cost_bins);
        h_elel_alpha->SetDirectory(0);
        sprintf(title, "ee%i_asym%s", year %2000, sys_label.c_str());
        h_elel_asym = new TH2F(title, "Asymmetric template of mc",
                n_xf_bins, xf_bins, n_cost_bins, cost_bins);
        h_elel_asym->SetDirectory(0);
        sprintf(title, "ee%i_bk%s", year %2000, sys_label.c_str());
        h_elel_back = new TH2F(title, "Combined background template",
                n_xf_bins, xf_bins, n_cost_bins, cost_bins);
        h_elel_back->SetDirectory(0);
        sprintf(title, "ee%i_dy_gg%s", year %2000, sys_label.c_str());
        h_elel_dy_gg = new TH2F(title, "Combined background template",
                n_xf_bins, xf_bins, n_cost_bins, cost_bins);
        h_elel_dy_gg->SetDirectory(0);
        sprintf(title, "ee%i_gam%s", year %2000, sys_label.c_str());
        h_elel_gam = new TH2F(title, "Combined background template",
                n_xf_bins, xf_bins, n_cost_bins, cost_bins);
        h_elel_gam->SetDirectory(0);

        printf("starting elel dy \n");
        gen_mc_template(t_elel_mc, alpha_denom, h_elel_sym, h_elel_asym, h_elel_alpha, year, m_low, m_high, FLAG_ELECTRONS,  do_RC, sys_label);
        TTree *elel_ts[1] = {t_elel_back};
        gen_combined_background_template(1, elel_ts, h_elel_back, year, m_low, m_high, FLAG_ELECTRONS, do_RC,ss, sys_label);
        elel_ts[0] = t_elel_nosig;
        gen_combined_background_template(1, elel_ts, h_elel_dy_gg, year, m_low, m_high, FLAG_ELECTRONS, do_RC,ss, sys_label);

        elel_ts[0] = t_elel_gamgam;
        printf("Making ee gamgam \n");
        gen_combined_background_template(1, elel_ts, h_elel_gam, year, m_low, m_high, FLAG_ELECTRONS, do_RC,ss, sys_label);
    }

}



void convert_mc_templates(int year, const string &sys_label){
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
        symmetrize2d(h_mumu_gam);
        auto h1_mumu_back = convert2d(h_mumu_back);
        auto h1_mumu_dy_gg = convert2d(h_mumu_dy_gg);
        auto h1_mumu_gam = convert2d(h_mumu_gam);
        auto h1_mumu_alpha = convert2d(h_mumu_alpha);

        auto h_mumu_pl = *h_mumu_sym + *h_mumu_asym;
        auto h_mumu_mn = *h_mumu_sym - *h_mumu_asym;
        h_mumu_pl.Scale(0.5);
        h_mumu_mn.Scale(0.5);
        auto h1_mumu_pl = convert2d(&h_mumu_pl);
        auto h1_mumu_mn = convert2d(&h_mumu_mn);
        char title[100];
        sprintf(title, "mumu%i_fpl%s", year%2000, sys_label.c_str());
        h1_mumu_pl->SetName(title);
        sprintf(title, "mumu%i_fmn%s", year%2000, sys_label.c_str());
        h1_mumu_mn->SetName(title);

        write_roo_hist(h1_mumu_back, var);
        write_roo_hist(h1_mumu_dy_gg, var);
        write_roo_hist(h1_mumu_gam, var);
        write_roo_hist(h1_mumu_alpha, var);
        write_roo_hist(h1_mumu_pl, var);
        write_roo_hist(h1_mumu_mn, var);
    }

    if(do_el){
        symmetrize2d(h_elel_gam);
        auto h1_elel_back = convert2d(h_elel_back);
        auto h1_elel_dy_gg = convert2d(h_elel_dy_gg);
        auto h1_elel_gam = convert2d(h_elel_gam);
        auto h1_elel_alpha = convert2d(h_elel_alpha);

        auto h_elel_pl = *h_elel_sym + *h_elel_asym;
        auto h_elel_mn = *h_elel_sym - *h_elel_asym;
        h_elel_pl.Scale(0.5);
        h_elel_mn.Scale(0.5);
        auto h1_elel_pl = convert2d(&h_elel_pl);
        auto h1_elel_mn = convert2d(&h_elel_mn);

        char title[100];
        sprintf(title, "ee%i_fpl%s", year%2000, sys_label.c_str());
        h1_elel_pl->SetName(title);
        sprintf(title, "ee%i_fmn%s", year%2000, sys_label.c_str());
        h1_elel_mn->SetName(title);


        write_roo_hist(h1_elel_back, var);
        write_roo_hist(h1_elel_dy_gg, var);
        write_roo_hist(h1_elel_gam, var);
        write_roo_hist(h1_elel_alpha, var);
        write_roo_hist(h1_elel_pl, var);
        write_roo_hist(h1_elel_mn, var);
    }
}

void write_groups(int year, FILE *f_log){

    char label[4][40], intro[4][40];
    int sizes[4] = {40,20,20,60};
    sprintf(intro[0], "emu%i_fake_shape group = ", year%2000);
    sprintf(intro[1], "ee%i_fake_shape group = ", year%2000);
    sprintf(intro[2], "mumu%i_fake_shape group = ", year%2000);
    sprintf(intro[3], "pdfs group = ");

    sprintf(label[0], "emu%i_qcd_param_bin", year%2000);
    sprintf(label[1], "ee%i_qcd_param_bin", year%2000);
    sprintf(label[2], "mumu%i_qcd_param_bin", year%2000);
    sprintf(label[3], "pdf");

    for(int i=0; i<4; i++){
        fprintf(f_log, "\n %s", intro[i]);
        for(int j=1; j<=sizes[i]; j++){

            fprintf(f_log, " %s%i", label[i], j);
        }
    }
        fprintf(f_log, "\n");
}




void make_templates(int year = 2016, int nJobs = 6, int iJob =-1){
    const TString fout_name("combine/templates/sep6_test.root");
    year = 2016;



    printf("Initializing files \n");
    init(year);
    init_ss(year);
    init_emu(year);
    printf("Setting up SFs... ");
    setup_all_SFs(year);
    printf("   done \n");

    vector<string> sys_labels {""};
        

    TFile * fout;
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
        sprintf(f_log_name, "combine/AFB_fits/cards/y%i_mbin%i_bins.txt", year, i);
        f_log = fopen(f_log_name, "w");

        m_low = m_bins[i];
        m_high = m_bins[i+1];
        printf("\n \n Start making templates for mass bin %.0f-%.0f \n", m_low, m_high);

        make_data_templates(year);
        make_ss_data_templates(year);
        make_ss_mc_templates(year);
        make_emu_data_templates(year);
        make_emu_qcd_templates(year, f_log);
        make_emu_mc_templates(year);
        make_qcd_templates(year,f_log);
        for(auto iter = sys_labels.begin(); iter !=sys_labels.end(); iter++){
            printf("Making MC templates for sys %s \n", (*iter).c_str());
            Double_t alpha_denom = alphas_denom[i];

            make_mc_templates(year, alpha_denom, *iter);
            convert_mc_templates(year, *iter);
        }
        fout->cd();
        gDirectory->cd(dirname);
        w->Write();
        write_groups(year, f_log);
        fclose(f_log);
        cleanup_mc_templates();
    }


    fout->Close();
    printf("Templates written to %s \n", fout_name.Data());

}

