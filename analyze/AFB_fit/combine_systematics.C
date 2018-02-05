

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
const int n_sys =2;
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
        afb_var_up[i] += pow((std::max(0., afb_diff)),2.);
        afb_var_down[i] += pow(std::max(0., -afb_diff),2.);

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
    if(output_file) fout = fopen("AFB_fit/systematics/Electron_systematics_feb5_test.txt", "w");


    string f_nominal("AFB_fit/fit_results/m_bins/ElEl_fit_jan31_nominal.root");
    TFile *f1 = TFile ::Open(f_nominal.c_str());
    t_nom  = (TTree *) f1->Get("T_fit_res");

    string f_pdf_str = string("AFB_fit/fit_results/m_bins/ElEl_m_fit_feb2_pdf_sys.root"); 
    TFile *f_pdf = TFile::Open(f_pdf_str.c_str());
    TTree * t_pdf = (TTree *) f_pdf->Get("T_fit_res");


    const int n_vars =2;
    TTree *t_vars[n_vars];
    string f_vars[n_vars];
    f_vars[0] = string("AFB_fit/fit_results/m_bins/ElEl_fit_jan31_nominal.root");
    f_vars[1] = string("AFB_fit/fit_results/m_bins/ElEl_fit_jan31_nominal.root");


    for(int i=0; i<n_vars; i++){
        TFile *f = TFile::Open(f_vars[i].c_str());
        t_vars[i] = (TTree *) f->Get("T_fit_res");
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
    for(int j=0; j<n_vars; j++){
        eval_systematic(t_vars[j], afb_var_up[0], afb_var_down[0]);
    }
    eval_pdf_systematic(t_pdf, afb_var_up[1], afb_var_down[1]);

    double afb_up[n_bins];
    double afb_down[n_bins];

    combine_sys(n_sys, afb_up, afb_down, afb_var_up, afb_var_down);

    printf("Systematic Errors  from files: \n");
    printf("%s \n", f_nominal.c_str());
    printf("%s \n", f_pdf.c_str());
    if(output_file){
        fprintf(fout, "Systematic Errors  from files: \n");
        fprintf(fout, "%s \n", f_nominal.c_str());
    }

    for(int j=0; j<n_vars; j++){
        printf("%s \n", f_vars[j].c_str());
        if(output_file) fprintf(fout, "%s \n", f_vars[j].c_str());
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
        
    
    


