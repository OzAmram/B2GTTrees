
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




void combined_fit_systematics(){
    TTree *t_nom;
    const int n_vars =2;
    TTree *t_vars[n_vars];
    const int n_bins = 8;

    vector<double> v_AFB_nom;
    vector<double> v_elel_nom;
    vector<double> v_mumu_nom;


    string f_nominal("AFB_fit/fit_results/combined_fit_dec11_nominal.root");
    string f_vars[n_vars];


    f_vars[0] = string("AFB_fit/fit_results/combined_fit_dec11_pdf_up.root");
    f_vars[1] = string("AFB_fit/fit_results/combined_fit_dec11_pdf_down.root");

    TFile *f1 = TFile ::Open(f_nominal.c_str());
    t_nom  = (TTree *) f1->Get("T_fit_res");

    for(int i=0; i<n_vars; i++){
        TFile *f = TFile::Open(f_vars[i].c_str());
        t_vars[i] = (TTree *) f->Get("T_fit_res");
    }
    Double_t afb_nom, r_elel_nom, r_mumu_nom;
    t_nom->SetBranchAddress("AFB", &afb_nom); 
    t_nom->SetBranchAddress("r_elel_back", &r_elel_nom); 
    t_nom->SetBranchAddress("r_mumu_back", &r_mumu_nom); 
    int nEntries = t_nom->GetEntries();
    for(int i=0; i<nEntries; i++){
        t_nom->GetEntry(i);
        v_AFB_nom.push_back(afb_nom);
        v_elel_nom.push_back(r_elel_nom);
        v_mumu_nom.push_back(r_mumu_nom);
    }






    double r_elel_var_up[n_bins];
    double r_elel_var_down[n_bins];
    double r_mumu_var_up[n_bins];
    double r_mumu_var_down[n_bins];
    double afb_var_up[n_bins];
    double afb_var_down[n_bins];
    for(int j=0; j<n_vars; j++){
        Double_t afb, r_elel, r_mumu;
        t_vars[j]->SetBranchAddress("AFB", &afb); 
        t_vars[j]->SetBranchAddress("r_elel_back", &r_elel); 
        t_vars[j]->SetBranchAddress("r_mumu_back", &r_mumu); 
        for(int i=0; i < n_bins; i++){
            t_vars[j]->GetEntry(i);
            double afb_diff = afb - v_AFB_nom[i];
            afb_var_up[i] += pow((std::max(0., afb_diff)),2.);
            afb_var_down[i] += pow(std::max(0., -afb_diff),2.);

            double r_elel_diff = r_elel - v_elel_nom[i];
            r_elel_var_up[i] += pow(std::max(0., r_elel_diff),2.);
            r_elel_var_down[i] += pow(std::max(0., -r_elel_diff),2.);

            double r_mumu_diff = r_mumu - v_mumu_nom[i];
            r_mumu_var_up[i] += pow(std::max(0., r_mumu_diff),2.);
            r_mumu_var_down[i] += pow(std::max(0., -r_mumu_diff),2.);
        }
    }

    printf("Systematic Errors  from files: \n");
    printf("%s \n", f_nominal.c_str());
    for(j=0; j<n_vars; j++){
        printf("%s \n", f_vars[j].c_str());
    }
    for(int i=0; i < n_bins; i++){
        r_elel_var_up[i] = sqrt(r_elel_var_up[i]);
        r_elel_var_down[i] = sqrt(r_elel_var_down[i]);
        r_mumu_var_up[i] = sqrt(r_mumu_var_up[i]);
        r_mumu_var_down[i] = sqrt(r_mumu_var_down[i]);

        afb_var_up[i] = sqrt(afb_var_up[i]);
        afb_var_down[i] = sqrt(afb_var_down[i]);
        printf("Bin %i: AFB = %.3f + %.3f - %.3f      R_elel = %.3f + %.3f - %.3f R_mumu = %.3f + %.3f - %.3f \n",
                i, v_AFB_nom[i], afb_var_up[i], afb_var_down[i], v_mumu_nom[i], r_mumu_var_up[i], r_mumu_var_down[i], v_mumu_nom[i], r_mumu_var_up[i], r_mumu_var_down[i]);
    }
    return;


}
        
    
    


