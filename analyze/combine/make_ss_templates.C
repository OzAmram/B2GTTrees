
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
//#include"Minuit2/Minuit2Minimizer.h"
#include "Math/Functor.h"
#include "../TemplateMaker_systematics.C"
//#include "../TemplateMaker.C"
#include "../AFB_fit/FitUtils.C"

Double_t m_low;
Double_t m_high;

bool print = true;
const TString fout_name("combine/templates/feb4_param_test.root");
TFile * fout;
char dirname[10];


Double_t med_btag = 0.4432;

TH1F *h_elel_dy, *h_elel_bk, *h_elel_back,  *h_elel_data,  *h_elel_qcd;
TH1F *h_mumu_dy, *h_mumu_bk, *h_mumu_back,  *h_mumu_data, *h_mumu_qcd;

TH1F *h1_elel_dy, *h1_elel_bk, *h1_elel_back,  *h1_elel_data,  *h1_elel_qcd;
TH1F *h1_mumu_dy, *h1_mumu_bk, *h1_mumu_back,  *h1_mumu_data, *h1_mumu_qcd;
//void *x_axp, *y_axp, *z_axp;

//define axis globally for convinience 
TAxis *x_ax, *y_ax, *z_ax;



vector<double> v_elel_xF;
vector<double> v_elel_cost;
vector<double> v_mumu_xF;
vector<double> v_mumu_cost;
unsigned int nElEl_DataEvents;
unsigned int nMuMu_DataEvents;

int my_cost_bins = 20;
TH1F * dummy1 = new TH1F("dummy", "", 1, 0, 1);
bool do_emu_scale = false;
bool do_RC = true;

TH1F* convert_to_abs_cost_hist(TH1F *h){
    //conver a hist of cost to a hist of abs(cost)
    int new_bins = my_cost_bins/2;
    TH1F *h1 = new TH1F(h->GetName(), "", new_bins, 0.,1.);

    h1->SetDirectory(0);
    //printf("OLD: ");
    for(int j=1; j<=my_cost_bins; j++){
        float content = h->GetBinContent(j);
        float error = h->GetBinError(j);
        int new_bin = (j > new_bins) ? j - new_bins  : new_bins - (j-1);
        //printf("Old bin %i new bin %i \n", j, new_bin);
        //printf("%.1f ",content);
        float prev_cont = h1->GetBinContent(new_bin);
        float prev_err = h1->GetBinContent(new_bin);
        h1->SetBinContent(new_bin, prev_cont + content);
        h1->SetBinError(new_bin, sqrt(prev_err*prev_err + error*error));
    }
    /*
    printf("\nNEW: ");
    
    
    for(int i=1; i<= new_bins; i++){

        printf("%.1f ", h1->GetBinContent(i));
    }
    printf("\n");
    */
    
    return h1;
}

std::pair<RooParametricHist*, RooAddition*>  convert_to_param_hist(TH1F *h){
    //conver a hist of abs(cost) to a parametric hist 
    int nbins = my_cost_bins/2;
    RooArgList *bin_list = new RooArgList();
    RooRealVar *var = new RooRealVar("cost", "cost", 0.,1.);

    for(int j=1; j<=nbins; j++){
        float content = h->GetBinContent(j);
        if(content<0) printf("Bin %i Content is %.0f \n", j, content);
        float error = h->GetBinError(j);
        char bin_name[10];
        sprintf(bin_name, "%s_%i",h->GetName(), j); 
        RooRealVar *bin = new RooRealVar(bin_name, bin_name, 1., 0., 10000.);
        //bin->setError(100error);
        //bin->Print();
        bin_list->add(*bin);
    }
    bin_list->Print();

    char h_name[20];
    sprintf(h_name, "%s", h->GetName());
    RooParametricHist *p= new RooParametricHist (h_name, "", *var, *bin_list, *h);
    char norm_name[20];
    sprintf(norm_name, "%s_norm", h_name);
    RooAddition *n = new RooAddition(norm_name, norm_name, *bin_list);

    
    return std::make_pair(p, n);
}


void make_data_templates(){
    h_elel_data = new TH1F("ee_ss_data_obs", "",
            my_cost_bins, -1., 1.);
    h_elel_data->SetDirectory(0);
    h_mumu_data = new TH1F("mumu_ss_data_obs", "",
            my_cost_bins, -1., 1.);
    h_mumu_data->SetDirectory(0);

    make_m_cost_pt_xf_hist(t_elel_ss_data, dummy1, h_elel_data, dummy1, dummy1, true, FLAG_ELECTRONS, do_emu_scale, do_RC, m_low, m_high);
    make_m_cost_pt_xf_hist(t_mumu_ss_data, dummy1, h_mumu_data, dummy1, dummy1, true, FLAG_MUONS, do_emu_scale, do_RC, m_low, m_high);

    h1_elel_data = convert_to_abs_cost_hist(h_elel_data);
    h1_mumu_data = convert_to_abs_cost_hist(h_mumu_data);

    printf("Integral of data templates were %.2f %.2f \n", h_elel_data->Integral(), h_mumu_data->Integral()); 
    printf("Integral of new data templates (abs) are  %.2f %.2f \n", h1_elel_data->Integral(), h1_mumu_data->Integral()); 
    fout->cd();
    gDirectory->cd(dirname);
    h1_elel_data->Write();
    h1_mumu_data->Write();
    printf("Made data templates \n");
}

void make_qcd_templates(){
    h_elel_qcd = new TH1F("ee_ss_qcd", "",
            my_cost_bins, -1., 1.);
    h_elel_qcd->SetDirectory(0);
    h_mumu_qcd = new TH1F("mumu_ss_qcd", "",
            my_cost_bins, -1., 1.);
    h_mumu_qcd->SetDirectory(0);

    make_qcd_from_emu_m_cost_xf_hist(t_emu_ss_data, t_emu_ss_ttbar, t_emu_ss_diboson, t_emu_ss_dy, dummy1, h_mumu_qcd, dummy1, m_low, m_high);
    h_elel_qcd = (TH1F*) h_mumu_qcd->Clone(h_elel_qcd->GetName());

    //Fakerate_est_mu(t_mumu_ss_WJets, t_mumu_ss_QCD, t_mumu_ss_WJets_mc, t_mumu_ss_QCD_mc, dummy1, h_mumu_qcd, dummy1, dummy1, m_low, m_high);
    //Fakerate_est_el(t_elel_ss_WJets, t_elel_ss_QCD, t_elel_ss_WJets_mc, t_elel_ss_QCD_mc, dummy1, h_elel_qcd, dummy1, dummy1, m_low, m_high);
    h1_elel_qcd = convert_to_abs_cost_hist(h_elel_qcd);
    h1_mumu_qcd = convert_to_abs_cost_hist(h_mumu_qcd);
    printf("Integral of qcd templates are %.2f %.2f \n", h1_elel_qcd->Integral(), h1_mumu_qcd->Integral()); 

    /*
    RooParametricHist *p_elel_qcd, *p_mumu_qcd;
    RooAddition *n_elel_qcd, *n_mumu_qcd;
    std::tie(p_elel_qcd, n_elel_qcd) = convert_to_param_hist(h1_elel_qcd);
    std::tie(p_mumu_qcd, n_mumu_qcd) = convert_to_param_hist(h1_mumu_qcd);
    */

    fout->cd();
    gDirectory->cd(dirname);
    h1_elel_qcd->Write();
    h1_mumu_qcd->Write();
    /*
    p_elel_qcd->Write();
    p_mumu_qcd->Write();
    n_elel_qcd->Write();
    n_mumu_qcd->Write();
    */
    printf("Made qcd templates \n");
}

void make_mc_templates(){
    h_elel_dy = new TH1F("ee_ss_dy", "",
            my_cost_bins, -1., 1.);
    h_elel_dy->SetDirectory(0);
    h_mumu_dy = new TH1F("mumu_ss_dy", "",
            my_cost_bins, -1., 1.);
    h_mumu_dy->SetDirectory(0);

    h_elel_bk = new TH1F("ee_ss_bk", "",
            my_cost_bins, -1., 1.);
    h_elel_bk->SetDirectory(0);
    h_mumu_bk = new TH1F("mumu_ss_bk", "",
            my_cost_bins, -1., 1.);
    h_mumu_bk->SetDirectory(0);



    make_m_cost_pt_xf_hist(t_elel_ss_back, dummy1, h_elel_bk, dummy1, dummy1, false, FLAG_ELECTRONS, do_emu_scale, do_RC, m_low, m_high);
    make_m_cost_pt_xf_hist(t_elel_ss_dy, dummy1, h_elel_dy, dummy1, dummy1, false, FLAG_ELECTRONS, do_emu_scale, do_RC, m_low, m_high);

    //make_m_cost_pt_xf_hist(t_mumu_ss_back, dummy1, h_mumu_bk, dummy1, dummy1, false, FLAG_MUONS, do_emu_scale, do_RC, m_low, m_high);
    //make_m_cost_pt_xf_hist(t_mumu_ss_dy, dummy1, h_mumu_dy, dummy1, dummy1, false, FLAG_MUONS, do_emu_scale, do_RC, m_low, m_high);
    h1_elel_bk = convert_to_abs_cost_hist(h_elel_bk);
    h1_mumu_bk = convert_to_abs_cost_hist(h_mumu_bk);

    h1_elel_dy = convert_to_abs_cost_hist(h_elel_dy);
    h1_mumu_dy = convert_to_abs_cost_hist(h_mumu_dy);


    printf("Integral of dy templates are %.2f %.2f \n", h1_elel_dy->Integral(), h1_mumu_dy->Integral()); 
    printf("Integral of bkg templates are %.2f %.2f \n", h1_elel_bk->Integral(), h1_mumu_bk->Integral()); 
    fout->cd();
    gDirectory->cd(dirname);
    h1_elel_dy->Write();
    //h1_mumu_dy->Write();
    printf("Made dy templates \n");

    h1_elel_bk->Write();
    //h1_mumu_bk->Write();
    printf("Made bk templates \n");


}





void make_ss_templates(){

    init_ss();
    init_emu_ss();


    fout = TFile::Open(fout_name, "RECREATE");


    for(int i=0; i<n_m_bins; i++){
        fout->cd();
        snprintf(dirname, 10, "%i", i);
        gDirectory->mkdir(dirname);
        gDirectory->cd(dirname);

        m_low = m_bins[i];
        m_high = m_bins[i+1];


        make_data_templates();
        make_qcd_templates();
        make_mc_templates();
    }
    printf("Templates written to %s \n", fout_name.Data());
}





