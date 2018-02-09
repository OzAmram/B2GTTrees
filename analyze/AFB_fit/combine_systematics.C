

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
const int n_sys =8;
vector<double> v_AFB_nom, v_AFB_stat_unc;

double get_var(Float_t vals[100]){
    double mean, var;
    for(int i =0; i< 100; i++){
        //printf("%.2f \n", vals[i]);
        mean += vals[i];
    }
    mean = mean / 100.;

    for(int  i=0; i< 100; i++){
        var += pow(vals[i] - mean, 2);
    }
    return var/99.;
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

void eval_pdf_systematic(TTree *t, double *afb_var_up, double *afb_var_down){
    Float_t afb[100];
    float var_low, var_high;
    t->SetBranchAddress("AFB", &afb); 
    for(int i=0; i < n_bins; i++){
        t->GetEntry(i);
        double afb_err = get_var(afb);

        afb_var_up[i] += afb_err;
        afb_var_down[i] += afb_err;

    }

}

void combine_sys(int n_sys, double *afb_up, double * afb_down, double afb_var_up[n_sys][n_bins], double afb_var_down[n_sys][n_bins]){
    for(int i =0; i< n_bins; i++){
        for(int j=0; j<n_sys; j++){
            if(afb_var_up[j][i] < 0. || afb_var_down[j][i] < 0.) printf("hey %i %i \n", i,j);
            afb_up[i] += afb_var_up[j][i];
            afb_down[i] += afb_var_down[j][i];
        }
        afb_up[i] = sqrt(afb_up[i]);
        afb_down[i] = sqrt(afb_down[i]);
    }
    return;
}



void combine_systematics(){
    TTree *t_nom;
    bool output_file = true;
    FILE * fout;
    if(output_file) fout = fopen("AFB_fit/systematics/Muon_systematics_feb6.txt", "w");


    string f_nominal("AFB_fit/fit_results/m_bins/MuMu_fit_jan31_nominal.root");
    TFile *f1 = TFile ::Open(f_nominal.c_str());
    t_nom  = (TTree *) f1->Get("T_fit_res");

    string f_pdf_str = string("AFB_fit/fit_results/m_bins/MuMu_m_fit_feb2_pdf_sys.root"); 
    TFile *f_pdf = TFile::Open(f_pdf_str.c_str());
    TTree * t_pdf = (TTree *) f_pdf->Get("T_fit_res");

    string f_bin_up_str = string("AFB_fit/fit_results/m_bins/MuMu_fit_feb6_bin_up.root"); 
    TFile *f_bin_up = TFile::Open(f_bin_up_str.c_str());
    TTree * t_bin_up = (TTree *) f_bin_up->Get("T_fit_res");
    string f_bin_down_str = string("AFB_fit/fit_results/m_bins/MuMu_fit_feb6_bin_down.root"); 
    TFile *f_bin_down = TFile::Open(f_bin_down_str.c_str());
    TTree * t_bin_down = (TTree *) f_bin_down->Get("T_fit_res");


    string f_alpha_up_str = string("AFB_fit/fit_results/m_bins/MuMu_fit_feb6_alpha_up.root"); 
    TFile *f_alpha_up = TFile::Open(f_alpha_up_str.c_str());
    TTree * t_alpha_up = (TTree *) f_alpha_up->Get("T_fit_res");
    string f_alpha_down_str = string("AFB_fit/fit_results/m_bins/MuMu_fit_feb6_alpha_down.root"); 
    TFile *f_alpha_down = TFile::Open(f_alpha_down_str.c_str());
    TTree * t_alpha_down = (TTree *) f_alpha_down->Get("T_fit_res");


    string f_emu_up_str = string("AFB_fit/fit_results/m_bins/MuMu_fit_feb6_emu_up.root"); 
    TFile *f_emu_up = TFile::Open(f_emu_up_str.c_str());
    TTree * t_emu_up = (TTree *) f_emu_up->Get("T_fit_res");
    string f_emu_down_str = string("AFB_fit/fit_results/m_bins/MuMu_fit_feb6_emu_down.root"); 
    TFile *f_emu_down = TFile::Open(f_emu_down_str.c_str());
    TTree * t_emu_down = (TTree *) f_emu_down->Get("T_fit_res");


    string f_fakes_up_str = string("AFB_fit/fit_results/m_bins/MuMu_fit_feb6_fakes_up.root"); 
    TFile *f_fakes_up = TFile::Open(f_fakes_up_str.c_str());
    TTree * t_fakes_up = (TTree *) f_fakes_up->Get("T_fit_res");
    string f_fakes_down_str = string("AFB_fit/fit_results/m_bins/MuMu_fit_feb6_fakes_down.root"); 
    TFile *f_fakes_down = TFile::Open(f_fakes_down_str.c_str());
    TTree * t_fakes_down = (TTree *) f_fakes_down->Get("T_fit_res");

    string f_btag_up_str = string("AFB_fit/fit_results/m_bins/MuMu_fit_feb6_btag_up.root"); 
    TFile *f_btag_up = TFile::Open(f_btag_up_str.c_str());
    TTree * t_btag_up = (TTree *) f_btag_up->Get("T_fit_res");
    string f_btag_down_str = string("AFB_fit/fit_results/m_bins/MuMu_fit_feb6_btag_down.root"); 
    TFile *f_btag_down = TFile::Open(f_btag_down_str.c_str());
    TTree * t_btag_down = (TTree *) f_btag_down->Get("T_fit_res");

    string f_pileup_up_str = string("AFB_fit/fit_results/m_bins/MuMu_fit_feb6_pileup_up.root"); 
    TFile *f_pileup_up = TFile::Open(f_pileup_up_str.c_str());
    TTree * t_pileup_up = (TTree *) f_pileup_up->Get("T_fit_res");
    string f_pileup_down_str = string("AFB_fit/fit_results/m_bins/MuMu_fit_feb6_pileup_down.root"); 
    TFile *f_pileup_down = TFile::Open(f_pileup_down_str.c_str());
    TTree * t_pileup_down = (TTree *) f_pileup_down->Get("T_fit_res");

    const int n_vars =6;
    TTree *t_scale[n_vars];
    string f_scale[n_vars];
    f_scale[0] = string("AFB_fit/fit_results/m_bins/MuMu_fit_feb6_mu_F_down.root");
    f_scale[1] = string("AFB_fit/fit_results/m_bins/MuMu_fit_feb6_mu_F_up.root");
    f_scale[2] = string("AFB_fit/fit_results/m_bins/MuMu_fit_feb6_mu_R_down.root");
    f_scale[3] = string("AFB_fit/fit_results/m_bins/MuMu_fit_feb6_mu_R_up.root");
    f_scale[4] = string("AFB_fit/fit_results/m_bins/MuMu_fit_feb6_mu_RF_down.root");
    f_scale[5] = string("AFB_fit/fit_results/m_bins/MuMu_fit_feb6_mu_RF_up.root");

    for(int i=0; i<n_vars; i++){
        TFile *f = TFile::Open(f_scale[i].c_str());
        t_scale[i] = (TTree *) f->Get("T_fit_res");
    }


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

    eval_systematic(t_fakes_up, afb_var_up[2], afb_var_down[2]);
    eval_systematic(t_fakes_down, afb_var_up[2], afb_var_down[2]);

    eval_systematic(t_alpha_up, afb_var_up[3], afb_var_down[3]);
    eval_systematic(t_alpha_down, afb_var_up[3], afb_var_down[3]);

    eval_systematic(t_emu_up, afb_var_up[4], afb_var_down[4]);
    eval_systematic(t_emu_down, afb_var_up[4], afb_var_down[4]);

    eval_systematic(t_btag_up, afb_var_up[5], afb_var_down[5]);
    eval_systematic(t_btag_down, afb_var_up[5], afb_var_down[5]);

    eval_systematic(t_pileup_up, afb_var_up[6], afb_var_down[6]);
    eval_systematic(t_pileup_down, afb_var_up[6], afb_var_down[6]);


    for(int j=0; j<n_vars; j++){
        eval_scale_systematic(t_scale[j], afb_var_up[7], afb_var_down[7]);
    }


    combine_sys(n_sys, afb_up, afb_down, afb_var_up, afb_var_down);

    if(output_file){
        fprintf(fout, "Systematic Errors  from files: \n");
        fprintf(fout, "%s \n", f_nominal.c_str());
        fprintf(fout, "%s \n", f_pdf_str.c_str());
        fprintf(fout, "%s %s \n", f_bin_up_str.c_str(), f_bin_down_str.c_str());
        fprintf(fout, "%s %s \n", f_fakes_up_str.c_str(), f_fakes_down_str.c_str());
        fprintf(fout, "%s %s \n", f_alpha_up_str.c_str(), f_alpha_down_str.c_str());
        fprintf(fout, "%s %s \n", f_emu_up_str.c_str(), f_emu_down_str.c_str());
        fprintf(fout, "%s %s \n", f_btag_up_str.c_str(), f_btag_down_str.c_str());
        fprintf(fout, "%s %s \n", f_pileup_up_str.c_str(), f_pileup_down_str.c_str());
        for(int j=0; j<n_vars; j++){
            fprintf(fout, "%s \n", f_scale[j].c_str());
        }
    }


    if(output_file){ 
        fprintf(fout, "\n\n\n");
        fprintf(fout, "Bin AFB    S_UNC  SYS_U   SYS_D  N_SYS SYS1_U SYS1_D  ... \n");
    }

    for(int i=0; i < n_bins; i++){
        //fix to print systematic first then different sources 
        printf("Bin %i: AFB = %.3f +- %.3f + %.3f - %.3f \n",
                i, v_AFB_nom[i], v_AFB_stat_unc[i], afb_up[i], afb_down[i]);
        if(output_file) {
            fprintf(fout, "%i   %.3f  %.3f  %.3f  %.3f    %i   ",
                i, v_AFB_nom[i], v_AFB_stat_unc[i], afb_up[i], -afb_down[i], n_sys);
            for(int j=0; j < n_sys; j++){
                fprintf(fout, "%.3f  %.3f  ", sqrt(afb_var_up[j][i]), -sqrt(afb_var_down[j][i]));
            }
            fprintf(fout, "\n");
        }
    }

    if(output_file) fclose(fout);
    return;


}
        
    
    


