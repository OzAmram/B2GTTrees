

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <iostream>
#include <cstring>
#include <algorithm>
#include <vector>
#include <string>

#include "TROOT.h"
#include "TFile.h"
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
#include "../TemplateMaker.C"

const int n_bins = 6;
const int n_sys =6;
vector<double> v_AFB_nom, v_AFB_stat_unc;

double get_var(Float_t vals[100]){
    float mean, var;
    int n_vars = 100;
    int n_entries = n_vars;
    for(int i =0; i< n_vars; i++){
        //printf("%.2f \n", vals[i]);
        if(vals[i] < 0. || vals[i] > 1. || std::isnan((float)vals[i])) n_entries--;
        else{
            //printf("val %.2f \n", vals[i]);
            mean += vals[i];
        }
        //printf("%.3e \n", vals[i]);
    }
    mean = mean / n_entries;
    printf("mean %.3f n_entries %i\n", mean, n_entries);

    for(int  i=0; i< n_vars; i++){
        if(vals[i] < 0. || vals[i] > 1. || std::isnan((float)vals[i])) continue;
        else var += pow(vals[i] - mean, 2);
    }
    var = var/(n_entries -1);
    printf("std %.3f \n\n\n", sqrt(var));
    return var;
}

void eval_systematic(TTree *t, double *afb_var_up, double *afb_var_down){
    Double_t afb;
    t->SetBranchAddress("AFB", &afb); 
    for(int i=0; i < n_bins; i++){
        t->GetEntry(i);
        double afb_diff = afb - v_AFB_nom[i];
        if(std::isnan(afb_diff)) printf("hi \n");
        //if (abs(afb_diff) > 0.01) printf ("%.3f \n", afb_diff);
        afb_var_up[i] += pow((std::max(0., afb_diff)),2.);
        afb_var_down[i] += pow(std::max(0., -afb_diff),2.);
        if(afb_var_up[i] < 0. || afb_var_down[i] < 0.) printf("hey \n");

    }

}

void eval_scale_systematic(TTree *t, double *afb_var_up, double *afb_var_down){
    Double_t afb;
    t->SetBranchAddress("AFB", &afb); 
    for(int i=0; i < n_bins; i++){
        t->GetEntry(i);
        double afb_diff = afb - v_AFB_nom[i];
        if(afb_diff > 0) afb_var_up[i] = max(afb_diff*afb_diff, afb_var_up[i]);
        if(afb_diff < 0) afb_var_down[i] = max(afb_diff*afb_diff, afb_var_down[i]);

    }

}

void eval_finite_mc_systematic(TTree *t_mc_fixed, TTree *t_mc, double *afb_var_up, double *afb_var_down){
    Double_t afb_fixed_err;
    Double_t afb_err;
    t_mc_fixed->SetBranchAddress("AFB_err", &afb_fixed_err); 
    t_mc->SetBranchAddress("AFB_err", &afb_err); 
    for(int i=0; i < n_bins; i++){
        t_mc_fixed->GetEntry(i);
        t_mc->GetEntry(i);
        double var_diff = abs(afb_fixed_err*afb_fixed_err - afb_err*afb_err);
        printf("Fixed err = %.4f , err = %.4f \n", afb_fixed_err, afb_err);

        afb_var_up[i] = var_diff;
        afb_var_down[i] = var_diff;

    }
    return;
}


void eval_pdf_systematic(TTree *t, double *afb_var_up, double *afb_var_down){
    Float_t afb[100];
    float var_low, var_high;
    t->SetBranchAddress("AFB", &afb); 
    for(int i=0; i < n_bins; i++){
        t->GetEntry(i);
        float afb_err = get_var(afb);
        //printf("AFB err %.2f \n", afb_err);

        afb_var_up[i] += afb_err;
        afb_var_down[i] += afb_err;

    }

}

void combine_sys(int n_sys, double *afb_up, double * afb_down, double afb_var_up[n_sys][n_bins], double afb_var_down[n_sys][n_bins]){
    for(int i =0; i< n_bins; i++){
        for(int j=0; j<n_sys; j++){
            if(afb_var_up[j][i] < 0. || afb_var_down[j][i] < 0.) printf("hey bin %i sys %i is nan \n", i,j);
            afb_up[i] += afb_var_up[j][i];
            afb_down[i] += afb_var_down[j][i];
        }
        afb_up[i] = sqrt(afb_up[i]);
        afb_down[i] = sqrt(afb_down[i]);
    }
    return;
}



void combined_systematics(){
    TTree *t_nom;
    bool output_file = true;
    FILE * fout;
    if(output_file) fout = fopen("AFB_fit/systematics/Combined_systematics_july16.txt", "w");


    string f_nominal("AFB_fit/fit_results/m_bins/combined_nominal_july16.root");
    TFile *f1 = TFile ::Open(f_nominal.c_str());
    t_nom  = (TTree *) f1->Get("T_fit_res");

    string f_pdf_str = string("AFB_fit/fit_results/m_bins/combined_m_fit_mar22_pdf_sys.root"); 
    TFile *f_pdf = TFile::Open(f_pdf_str.c_str());
    TTree * t_pdf = (TTree *) f_pdf->Get("T_fit_res");

    string f_bin_up_str = string("AFB_fit/fit_results/m_bins/combined_bin_up_july16.root"); 
    TFile *f_bin_up = TFile::Open(f_bin_up_str.c_str());
    TTree * t_bin_up = (TTree *) f_bin_up->Get("T_fit_res");
    string f_bin_down_str = string("AFB_fit/fit_results/m_bins/combined_bin_down_july16.root"); 
    TFile *f_bin_down = TFile::Open(f_bin_down_str.c_str());
    TTree * t_bin_down = (TTree *) f_bin_down->Get("T_fit_res");


    string f_alpha_up_str = string("AFB_fit/fit_results/m_bins/combined_alpha_up_july16.root"); 
    TFile *f_alpha_up = TFile::Open(f_alpha_up_str.c_str());
    TTree * t_alpha_up = (TTree *) f_alpha_up->Get("T_fit_res");
    string f_alpha_down_str = string("AFB_fit/fit_results/m_bins/combined_alpha_down_july16.root"); 
    TFile *f_alpha_down = TFile::Open(f_alpha_down_str.c_str());
    TTree * t_alpha_down = (TTree *) f_alpha_down->Get("T_fit_res");


    /*
    string f_emu_up_str = string("AFB_fit/fit_results/m_bins/combined_emu_up_july16.root"); 
    TFile *f_emu_up = TFile::Open(f_emu_up_str.c_str());
    TTree * t_emu_up = (TTree *) f_emu_up->Get("T_fit_res");
    string f_emu_down_str = string("AFB_fit/fit_results/m_bins/combined_emu_down_july16.root"); 
    TFile *f_emu_down = TFile::Open(f_emu_down_str.c_str());
    TTree * t_emu_down = (TTree *) f_emu_down->Get("T_fit_res");
    */


    string f_el_fakes_up_str = string("AFB_fit/fit_results/m_bins/combined_el_fake_up_july16.root"); 
    TFile *f_el_fakes_up = TFile::Open(f_el_fakes_up_str.c_str());
    TTree * t_el_fakes_up = (TTree *) f_el_fakes_up->Get("T_fit_res");
    string f_el_fakes_down_str = string("AFB_fit/fit_results/m_bins/combined_el_fake_down_july16.root"); 
    TFile *f_el_fakes_down = TFile::Open(f_el_fakes_down_str.c_str());
    TTree * t_el_fakes_down = (TTree *) f_el_fakes_down->Get("T_fit_res");

    string f_mu_fakes_up_str = string("AFB_fit/fit_results/m_bins/combined_mu_fake_up_july16.root"); 
    TFile *f_mu_fakes_up = TFile::Open(f_mu_fakes_up_str.c_str());
    TTree * t_mu_fakes_up = (TTree *) f_mu_fakes_up->Get("T_fit_res");
    string f_mu_fakes_down_str = string("AFB_fit/fit_results/m_bins/combined_mu_fake_down_july16.root"); 
    TFile *f_mu_fakes_down = TFile::Open(f_mu_fakes_down_str.c_str());
    TTree * t_mu_fakes_down = (TTree *) f_mu_fakes_down->Get("T_fit_res");

    string f_rc_alt_str = string("AFB_fit/fit_results/m_bins/combined_RC_alt_july16.root"); 
    TFile *f_rc_alt = TFile::Open(f_rc_alt_str.c_str());
    TTree * t_rc_alt = (TTree *) f_rc_alt->Get("T_fit_res");

    /*
    string f_btag_up_str = string("AFB_fit/fit_results/m_bins/combined_btag_up_july16.root"); 
    TFile *f_btag_up = TFile::Open(f_btag_up_str.c_str());
    TTree * t_btag_up = (TTree *) f_btag_up->Get("T_fit_res");
    string f_btag_down_str = string("AFB_fit/fit_results/m_bins/combined_btag_down_july16.root"); 
    TFile *f_btag_down = TFile::Open(f_btag_down_str.c_str());
    TTree * t_btag_down = (TTree *) f_btag_down->Get("T_fit_res");
    */

    string f_pileup_up_str = string("AFB_fit/fit_results/m_bins/combined_pu_up_july16.root"); 
    TFile *f_pileup_up = TFile::Open(f_pileup_up_str.c_str());
    TTree * t_pileup_up = (TTree *) f_pileup_up->Get("T_fit_res");
    string f_pileup_down_str = string("AFB_fit/fit_results/m_bins/combined_pu_down_july16.root"); 
    TFile *f_pileup_down = TFile::Open(f_pileup_down_str.c_str());
    TTree * t_pileup_down = (TTree *) f_pileup_down->Get("T_fit_res");

    /*
    const int n_vars =6;
    TTree *t_scale[n_vars];
    string f_scale[n_vars];
    f_scale[0] = string("AFB_fit/fit_results/m_bins/combined_mu_RF_up_july16.root");
    f_scale[1] = string("AFB_fit/fit_results/m_bins/combined_mu_RF_down_july16.root");
    f_scale[2] = string("AFB_fit/fit_results/m_bins/combined_mu_F_up_july16.root");
    f_scale[3] = string("AFB_fit/fit_results/m_bins/combined_mu_F_down_july16.root");
    f_scale[4] = string("AFB_fit/fit_results/m_bins/combined_mu_R_up_july16.root");
    f_scale[5] = string("AFB_fit/fit_results/m_bins/combined_mu_R_down_july16.root");
    */

    const int n_SFs =2;
    TTree *t_SF[n_SFs];
    string f_SF[n_SFs];
    f_SF[0] = string("AFB_fit/fit_results/m_bins/combined_el_SF_off_july16.root");
    f_SF[1] = string("AFB_fit/fit_results/m_bins/combined_mu_SF_off_july16.root");
    /*
    const int n_SFs =14;
    TTree *t_SF[n_SFs];
    string f_SF[n_SFs];
    f_SF[0] = string("AFB_fit/fit_results/m_bins/combined_fit_el_HLT_down_mar24.root");
    f_SF[1] = string("AFB_fit/fit_results/m_bins/combined_fit_el_HLT_up_july9.root");
    f_SF[2] = string("AFB_fit/fit_results/m_bins/combined_fit_el_ID_up_mar24.root");
    f_SF[3] = string("AFB_fit/fit_results/m_bins/combined_fit_el_ID_down_mar24.root");
    f_SF[4] = string("AFB_fit/fit_results/m_bins/combined_fit_el_reco_down_mar24.root");
    f_SF[5] = string("AFB_fit/fit_results/m_bins/combined_fit_el_reco_up_mar24.root");
    f_SF[6] = string("AFB_fit/fit_results/m_bins/combined_fit_mu_HLT_down_mar24.root");
    f_SF[7] = string("AFB_fit/fit_results/m_bins/combined_fit_mu_HLT_up_mar24.root");
    f_SF[8] = string("AFB_fit/fit_results/m_bins/combined_fit_mu_ID_up_mar24.root");
    f_SF[9] = string("AFB_fit/fit_results/m_bins/combined_fit_mu_ID_down_mar24.root");
    f_SF[10] = string("AFB_fit/fit_results/m_bins/combined_fit_mu_reco_down_mar24.root");
    f_SF[11] = string("AFB_fit/fit_results/m_bins/combined_fit_mu_reco_up_mar24.root");
    f_SF[12] = string("AFB_fit/fit_results/m_bins/combined_fit_mu_iso_down_mar24.root");
    f_SF[13] = string("AFB_fit/fit_results/m_bins/combined_fit_mu_iso_up_mar24.root");
    */

    const int n_finite_mc=2;
    TTree *t_finite_mc[n_finite_mc];
    string f_finite_mc[n_finite_mc];
    f_finite_mc[0] = string("AFB_fit/fit_results/m_bins/Combined_fit_finite_mc_stat_fixed_july9.root");
    f_finite_mc[1] = string("AFB_fit/fit_results/m_bins/Combined_fit_finite_mc_stat_july9.root");
    TFile *f_mc_fixed = TFile::Open(f_finite_mc[0].c_str());
    TTree *t_mc_stat_fixed = (TTree *) f_mc_fixed->Get("T_fit_res");
    TFile *f_mc = TFile::Open(f_finite_mc[1].c_str());
    TTree *t_mc_stat = (TTree *) f_mc->Get("T_fit_res");

    /*
    for(int i=0; i<n_vars; i++){
        TFile *f = TFile::Open(f_scale[i].c_str());
        t_scale[i] = (TTree *) f->Get("T_fit_res");
    }
    */

    /*
    for(int i=0; i<n_SFs; i++){
        TFile *f = TFile::Open(f_SF[i].c_str());
        t_SF[i] = (TTree *) f->Get("T_fit_res");
    }
    */


    Double_t afb_nom, afb_err;
    t_nom->SetBranchAddress("AFB", &afb_nom); 
    t_nom->SetBranchAddress("AFB_err", &afb_err); 
    int nEntries = t_nom->GetEntries();
    for(int i=0; i<nEntries; i++){
        t_nom->GetEntry(i);
        v_AFB_nom.push_back(afb_nom);

        v_AFB_stat_unc.push_back(afb_err);
    }


    double afb_var_up[n_sys][n_bins];
    double afb_var_down[n_sys][n_bins];
    double afb_up[n_bins];
    double afb_down[n_bins];
    for(int i =0; i<n_sys; i++){
        for(int j=0; j<n_bins; j++){
            afb_var_up[i][j] =0;
            afb_var_down[i][j] =0;
            afb_up[j] =0;
            afb_down[j] =0;
        }
    }


    eval_pdf_systematic(t_pdf, afb_var_up[0], afb_var_down[0]);

    eval_systematic(t_bin_up, afb_var_up[1], afb_var_down[1]);
    eval_systematic(t_bin_down, afb_var_up[1], afb_var_down[1]);

    eval_systematic(t_el_fakes_up, afb_var_up[2], afb_var_down[2]);
    eval_systematic(t_el_fakes_down, afb_var_up[2], afb_var_down[2]);
    eval_systematic(t_mu_fakes_up, afb_var_up[2], afb_var_down[2]);
    eval_systematic(t_mu_fakes_down, afb_var_up[2], afb_var_down[2]);

    eval_systematic(t_alpha_up, afb_var_up[3], afb_var_down[3]);
    eval_systematic(t_alpha_down, afb_var_up[3], afb_var_down[3]);

    eval_systematic(t_rc_alt, afb_var_up[4], afb_var_down[4]);

    //eval_systematic(t_emu_up, afb_var_up[4], afb_var_down[4]);
    //eval_systematic(t_emu_down, afb_var_up[4], afb_var_down[4]);

    //eval_systematic(t_btag_up, afb_var_up[5], afb_var_down[5]);
    //eval_systematic(t_btag_down, afb_var_up[5], afb_var_down[5]);

    eval_systematic(t_pileup_up, afb_var_up[5], afb_var_down[5]);
    eval_systematic(t_pileup_down, afb_var_up[5], afb_var_down[5]);


    //for(int j=0; j<n_vars; j++){
        //eval_scale_systematic(t_scale[j], afb_var_up[7], afb_var_down[7]);
    //}

    for(int j=0; j<n_SFs; j++){
        //if(j==10) continue;
        //eval_systematic(t_SF[j], afb_var_up[8], afb_var_down[8]);
        //printf("stds are %.3f %.3f \n", sqrt(afb_var_up[8][3]), sqrt(afb_var_down[8][3]));
    }
    //eval_finite_mc_systematic(t_mc_stat_fixed, t_mc_stat, afb_var_up[9], afb_var_down[9]);

    combine_sys(n_sys, afb_up, afb_down, afb_var_up, afb_var_down);

    if(output_file){
        fprintf(fout, "Systematic Errors  from files: \n");
        fprintf(fout, "%s \n", f_nominal.c_str());
        fprintf(fout, "%s \n", f_pdf_str.c_str());
        fprintf(fout, "%s %s \n", f_bin_up_str.c_str(), f_bin_down_str.c_str());
        fprintf(fout, "%s %s \n", f_el_fakes_up_str.c_str(), f_el_fakes_down_str.c_str());
        fprintf(fout, "%s %s \n", f_mu_fakes_up_str.c_str(), f_mu_fakes_down_str.c_str());
        fprintf(fout, "%s \n", f_rc_alt_str.c_str());
        fprintf(fout, "%s %s \n", f_alpha_up_str.c_str(), f_alpha_down_str.c_str());
     //   fprintf(fout, "%s %s \n", f_emu_up_str.c_str(), f_emu_down_str.c_str());
        //fprintf(fout, "%s %s \n", f_btag_up_str.c_str(), f_btag_down_str.c_str());
        fprintf(fout, "%s %s \n", f_pileup_up_str.c_str(), f_pileup_down_str.c_str());
        //for(int j=0; j<n_vars; j++){
            //fprintf(fout, "%s \n", f_scale[j].c_str());
        //}
        for(int j=0; j<n_SFs; j++){
            //fprintf(fout, "%s \n", f_SF[j].c_str());
        }
        for(int j=0; j<n_finite_mc; j++){
            fprintf(fout, "%s \n", f_finite_mc[j].c_str());
        }
    }


    if(output_file){ 
        fprintf(fout, "\n\n\n");
        fprintf(fout, "Bin AFB    S_UNC  SYS_U   SYS_D  UNC_U  UNC_D  N_SYS SYS1_U SYS1_D  ... \n");
    }

    for(int i=0; i < n_bins; i++){
        //fix to print systematic first then different sources 
        printf("Bin %i: AFB = %.3f +- %.3f + %.3f - %.3f \n",
                i, v_AFB_nom[i], v_AFB_stat_unc[i], afb_up[i], afb_down[i]);
        if(output_file) {
            fprintf(fout, "%i   %.3f  %.3f  %.3f  %.3f  %.3f %.3f    %i   ",
                i, v_AFB_nom[i], v_AFB_stat_unc[i], afb_up[i], -afb_down[i], 
                sqrt(pow(v_AFB_stat_unc[i],2) + pow(afb_up[i],2)), sqrt(pow(v_AFB_stat_unc[i],2) + pow(afb_down[i],2)), n_sys);
            for(int j=0; j < n_sys; j++){
                fprintf(fout, "%.3f  %.3f  ", sqrt(afb_var_up[j][i]), -sqrt(afb_var_down[j][i]));
            }
            fprintf(fout, "\n");
        }
    }

    if(output_file) fclose(fout);
    return;


}
        
    
    


