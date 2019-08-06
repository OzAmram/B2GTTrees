

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


const TString emu_fout_name("combine/templates/feb12_emu_templates.root");
TFile * emu_fout;




void convert_emu_qcd_to_param_hist(TH1F *h, FILE *f_log){
    //convert a hist to a parametric hist 
    RooArgList *bin_list = new RooArgList();

    char h_name[40];
    sprintf(h_name, "%s_param", h->GetName());
    RooRealVar *Rqcd;
    for(int j=1; j <= h->GetNbinsX(); j++){

        double content = h->GetBinContent(j);
        if(content<0) printf("Bin %i Content is %.0f \n", j, content);
        double error = h->GetBinError(j);
        //printf("Bin %.1f error %.1f \n", content,error);
        char bin_name[40];
        char form_name[40];
        sprintf(bin_name, "%s_bin%i",h_name, j); 
        sprintf(form_name, "%s_form%i",h_name, j); 
        //prevent underflowing by fixing super small bins
        content = max(content, 0.001);
        error = max(error, 0.0001);
        if (content < error){
            content = error/2.;
            error = 0.1*content;
        }
        else if(content < 2.5 * error){
            error = 0.3*content;
        }
        RooRealVar *bin = new RooRealVar(bin_name, bin_name, content, 0., 10000.);
        fprintf(f_log, "%s param %.4f %.4f \n", bin_name, content, error);
        //bin->Print();
        bin_list->add(*bin);
    }
    bin_list->Print();

    RooParametricHist *p= new RooParametricHist (h_name, h_name, *var, *bin_list, *h);
    char norm_name[40];
    sprintf(norm_name, "%s_norm", h_name);
    RooAddition *n = new RooAddition(norm_name, norm_name, *bin_list);

    
    w->import(*p);
    w->import(*n,RooFit::RecycleConflictNodes());
}


void make_emu_data_templates(int year){
    char title[100];
    sprintf(title, "emu%i_data_obs", year%2000);
    auto h_emu_data = new TH2F(title, "Data template of (x_f, cost_r)",
            n_xf_bins, xf_bins, n_cost_bins, cost_bins);
    h_emu_data->SetDirectory(0);
    bool ss = false;

    gen_emu_template(t_emu_data, h_emu_data,  true, year, m_low, m_high);
    auto h1_emu_data = convert2d(h_emu_data);


    printf("Integral of data templates were %.2f  \n", h_emu_data->Integral() ); 
    //h_emu_data->Write();
    write_roo_hist(h1_emu_data, var);
    printf("Made emu data templates \n");
}

void make_emu_qcd_templates(int year, FILE *f_log){
    char title[100];
    sprintf(title, "emu%i_qcd", year%2000);
    auto h_emu_qcd = new TH2F(title, "Combined background template",
            n_xf_bins, xf_bins, n_cost_bins, cost_bins);
    h_emu_qcd->SetDirectory(0);

    gen_emu_fakes_template(t_emu_WJets, t_emu_QCD, t_emu_WJets_contam,  h_emu_qcd, year, m_low, m_high);
    printf("Integral of qcd template are %.2f \n", h_emu_qcd->Integral()); 
    auto h1_emu_qcd = convert2d(h_emu_qcd);

    convert_emu_qcd_to_param_hist(h1_emu_qcd, f_log);
    printf("Made emu qcd templates \n");
}

void make_emu_mc_templates(int year){
    char title[100];
    sprintf(title, "emu%i_dy", year%2000);
    auto h_emu_dy = new TH2F(title, "dy template of (x_f, cost_r)",
            n_xf_bins, xf_bins, n_cost_bins, cost_bins);
    h_emu_dy->SetDirectory(0);
    sprintf(title, "emu%i_bk", year%2000);
    auto h_emu_bk = new TH2F(title, "bk template of (x_f, cost_r)",
            n_xf_bins, xf_bins, n_cost_bins, cost_bins);
    h_emu_bk->SetDirectory(0);

    bool ss = true;
    


    gen_emu_template(t_emu_back, h_emu_bk,  false, year, m_low, m_high);
    gen_emu_template(t_emu_dy, h_emu_dy,  false, year, m_low, m_high);

    auto h1_emu_bk = convert2d(h_emu_bk);
    auto h1_emu_dy = convert2d(h_emu_dy);




    printf("Integral of dy templates are %.2f \n", h_emu_dy->Integral()); 
    printf("Integral of bkg templates are %.2f \n", h_emu_bk->Integral()); 
    write_roo_hist(h1_emu_dy, var);
    printf("Made emu dy templates \n");

    write_roo_hist(h1_emu_bk, var);
    printf("Made emu bk templates \n");


}





void make_emu_templates(int year=2016){

    init_emu(year);


    emu_fout = TFile::Open(emu_fout_name, "RECREATE");
    FILE *f_log;
    char f_log_name[80];
    char dirname[40];


    for(int i=0; i<n_m_bins; i++){
        emu_fout->cd();
        snprintf(dirname, 10, "w%i", i);
        gDirectory->mkdir(dirname);
        gDirectory->cd(dirname);
        w = new RooWorkspace("w", "w");

        m_low = m_bins[i];
        m_high = m_bins[i+1];

        sprintf(f_log_name, "combine/EMu_fits/cards/y%i_mbin%i_bins.txt", year, i);
        f_log = fopen(f_log_name, "w");



        make_emu_data_templates(year);
        make_emu_qcd_templates(year,f_log);
        make_emu_mc_templates(year);
        emu_fout->cd();
        gDirectory->cd(dirname);
        w->Write();
        fclose(f_log);
    }
    printf("Templates written to %s \n", emu_fout_name.Data());
}





