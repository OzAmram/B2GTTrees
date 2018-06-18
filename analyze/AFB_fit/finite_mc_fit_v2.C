


//perform fits to Reconstructed MuMu data to extract Asym

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <iostream>
#include <cstring>
#include <algorithm>
#include <vector>

#include "TObject.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TH2.h"
#include "TProfile.h"
#include "TStyle.h"
#include "TPostScript.h"
#include "TLorentzVector.h"
#include "TVector3.h"
#include "TFitter.h"
#include "TSystem.h"
//#include"Minuit2/Minuit2Minimizer.h"
#include "Math/Functor.h"
//#include "../TemplateMaker_systematics.C"
#include "../TemplateMaker.C"
#include "root_files.h"

using namespace std;




const TString mumu_fout_name("AFB_fit/fit_results/m_bins/MuMu_fit_finite_mc_stat_fixed_june14.root");
const TString elel_fout_name("AFB_fit/fit_results/m_bins/ElEl_fit_finite_mc_stat_fixed_june14.root");

int FLAG = FLAG_ELECTRONS;
bool do_both = false;


float m_low;
float m_high;

bool print = true;


Double_t med_btag = 0.4432;

TH2F *h_elel_asym, *h_elel_sym, *h_elel_back,  *h_elel_data, *h_elel_dilu, *h_elel_fitted_sym, *h_elel_fitted_asym, *h_elel_cov;
TH2F *h_elel_mc_count, *h_elel_sym_count, *h_elel_mc, *h_elel_mc_errs;
TH2F *h_mumu_asym, *h_mumu_sym, *h_mumu_back,  *h_mumu_data, *h_mumu_dilu, *h_mumu_fitted_sym, *h_mumu_fitted_asym, *h_mumu_cov;
TH2F *h_mumu_mc_count, *h_mumu_sym_count, *h_mumu_mc, *h_mumu_mc_errs;

double m_elel_cov[5][10][3], m_mumu_cov[5][10][3];
//void *x_axp, *y_axp, *z_axp;

//define axis globally for convinience 
TAxis *x_ax, *y_ax, *z_ax;



vector<double> v_elel_xF;
vector<double> v_elel_cost;

vector<double> v_mumu_xF;
vector<double> v_mumu_cost;
unsigned int nElEl_DataEvents;
unsigned int nMuMu_DataEvents;


Double_t get_prob(Double_t xF, Double_t cost, TH2F *h){
    TAxis* x_ax =  h->GetXaxis();
    TAxis *y_ax =  h->GetYaxis();
    //binning the same in all templates
    int xbin = x_ax->FindBin(xF);
    int ybin = y_ax->FindBin(cost);

    return h->GetBinContent(xbin, ybin);
}

int gen_finite_mc_template_v2(TTree *t1, Double_t alpha, TH2F* h_sym, TH2F *h_asym, TH2F *h_cov, 
        Double_t var_low, Double_t var_high, int flag1 = FLAG_MUONS, int flag2 = FLAG_M_BINS,
        int pdf_sys = -1){
    Long64_t nEntries  =  t1->GetEntries();
    //printf("size is %i \n", nEntries);

    h_sym->Sumw2();
    h_asym->Sumw2();
    h_cov->Sumw2();

    //TH2F* h_reweights = (TH2F *) h_sym->Clone("h_rw");

    Double_t m, xF, cost, gen_weight, reweight, jet1_cmva, jet2_cmva, cost_st;
    Double_t bcdef_HLT_SF, bcdef_iso_SF, bcdef_id_SF;
    Double_t gh_HLT_SF, gh_iso_SF, gh_id_SF;
    Double_t gh_trk_SF, bcdef_trk_SF;
    Double_t el_id_SF, el_reco_SF, pu_SF, el_HLT_SF;
    Double_t jet1_pt, jet2_pt, jet1_b_weight, jet2_b_weight, jet1_eta, jet2_eta;
    Double_t mu1_pt, mu1_eta, mu2_pt, mu2_eta;
    Double_t el1_pt, el1_eta, el2_pt, el2_eta;
    Double_t mu_R_up, mu_R_down, mu_F_up, mu_F_down, 
             mu_RF_up, mu_RF_down, pdf_up, pdf_down;
    Float_t cost_pt, met_pt;
    Float_t pdf_weights[100];
    TLorentzVector *lep_p=0;
    TLorentzVector *lep_m=0;
    Double_t pt;
    Int_t nJets, pu_NtrueInt, jet1_flavour, jet2_flavour;
    t1->SetBranchAddress("m", &m);
    t1->SetBranchAddress("xF", &xF);
    t1->SetBranchAddress("cost", &cost);
    t1->SetBranchAddress("cost_st", &cost_st);
    t1->SetBranchAddress("jet1_CMVA", &jet1_cmva);
    t1->SetBranchAddress("jet2_CMVA", &jet2_cmva);
    t1->SetBranchAddress("met_pt", &met_pt);
    t1->SetBranchAddress("jet1_pt", &jet1_pt);
    t1->SetBranchAddress("jet2_pt", &jet2_pt);
    t1->SetBranchAddress("nJets", &nJets);
    t1->SetBranchAddress("gen_weight", &gen_weight);
    t1->SetBranchAddress("jet1_b_weight", &jet1_b_weight);
    t1->SetBranchAddress("jet2_b_weight", &jet2_b_weight);
    t1->SetBranchAddress("pu_SF", &pu_SF);
    t1->SetBranchAddress("mu_R_up", &mu_R_up);
    t1->SetBranchAddress("mu_R_down", &mu_R_down);
    t1->SetBranchAddress("mu_F_up", &mu_F_up);
    t1->SetBranchAddress("mu_F_down", &mu_F_down);
    t1->SetBranchAddress("mu_RF_up", &mu_RF_up);
    t1->SetBranchAddress("mu_RF_down", &mu_RF_down);
    t1->SetBranchAddress("pdf_up", &pdf_up);
    t1->SetBranchAddress("pdf_down", &pdf_down);
    t1->SetBranchAddress("pu_NtrueInt", &pu_NtrueInt);
    t1->SetBranchAddress("jet1_eta", &jet1_eta);
    t1->SetBranchAddress("jet2_eta", &jet2_eta);
    t1->SetBranchAddress("jet1_flavour", &jet1_flavour);
    t1->SetBranchAddress("jet2_flavour", &jet2_flavour);
    if(pdf_sys >=0) t1->SetBranchAddress("pdf_weights", &pdf_weights);
    int n = 0;


    if(flag1 == FLAG_MUONS){
        t1->SetBranchAddress("mu_p", &lep_p);
        t1->SetBranchAddress("mu_m", &lep_m);
        t1->SetBranchAddress("mu1_pt", &mu1_pt);
        t1->SetBranchAddress("mu1_eta", &mu1_eta);
        t1->SetBranchAddress("mu2_pt", &mu2_pt);
        t1->SetBranchAddress("mu2_eta", &mu2_eta);
        t1->SetBranchAddress("bcdef_HLT_SF", &bcdef_HLT_SF);
        t1->SetBranchAddress("bcdef_iso_SF", &bcdef_iso_SF);
        t1->SetBranchAddress("bcdef_id_SF", &bcdef_id_SF);
        t1->SetBranchAddress("bcdef_trk_SF", &bcdef_trk_SF);
        t1->SetBranchAddress("gh_HLT_SF", &gh_HLT_SF);
        t1->SetBranchAddress("gh_iso_SF", &gh_iso_SF);
        t1->SetBranchAddress("gh_id_SF", &gh_id_SF);
        t1->SetBranchAddress("gh_trk_SF", &gh_trk_SF);




        for (int i=0; i<nEntries; i++) {
            t1->GetEntry(i);
            bool no_bjets = has_no_bjets(nJets, jet1_pt, jet2_pt, jet1_cmva, jet2_cmva);
            if(flag2 == FLAG_PT_BINS){
                TLorentzVector cm = *lep_p + *lep_m;
                pt = cm.Pt();
            }
            bool pass = ((flag2 == FLAG_M_BINS && m >= var_low && m <= var_high) ||
                        (flag2 == FLAG_PT_BINS && m >= 150. && pt >= var_low && pt <= var_high))
                        && met_pt < 50.  && no_bjets;
            if(pass){
                reweight = (4./3.)*cost_st*(2. + alpha)/
                    (1. + cost_st*cost_st + alpha*(1.- cost_st*cost_st));
                n++;

                Double_t bcdef_weight = gen_weight * pu_SF * bcdef_HLT_SF * bcdef_iso_SF * bcdef_id_SF * bcdef_trk_SF;
                Double_t gh_weight = gen_weight * pu_SF * gh_HLT_SF * gh_iso_SF * gh_id_SF * gh_trk_SF;
                if (nJets >= 1){
                    bcdef_weight *= jet1_b_weight;
                    gh_weight *= jet1_b_weight;
                }
                if (nJets >= 2){
                    bcdef_weight *= jet2_b_weight;
                    gh_weight *= jet2_b_weight;
                }


                Double_t final_weight = 1000*(bcdef_weight*bcdef_lumi + gh_weight*gh_lumi);
                h_sym->Fill(xF, cost, final_weight); 
                h_sym->Fill(xF, -cost, final_weight); 

                h_asym->Fill(xF, cost, reweight * final_weight);
                h_asym->Fill(xF, -cost, -reweight * final_weight);
                //h_reweights->Fill(xF, fabs(cost), fabs(reweight));

                //h_count->Fill(xF, fabs(cost), 1);
                h_cov->Fill(xF, abs(cost), abs(reweight)*final_weight*final_weight);
            }
        }

    }
    else if (flag1 == FLAG_ELECTRONS) {
        t1->SetBranchAddress("el_p", &lep_p);
        t1->SetBranchAddress("el_m", &lep_m);
        t1->SetBranchAddress("el1_pt", &el1_pt);
        t1->SetBranchAddress("el1_eta", &el1_eta);
        t1->SetBranchAddress("el2_pt", &el2_pt);
        t1->SetBranchAddress("el2_eta", &el2_eta);
        t1->SetBranchAddress("el_id_SF", &el_id_SF);
        t1->SetBranchAddress("el_HLT_SF", &el_HLT_SF);
        t1->SetBranchAddress("el_reco_SF", &el_reco_SF);

        for (int i=0; i<nEntries; i++) {
            t1->GetEntry(i);
            bool no_bjets = has_no_bjets(nJets, jet1_pt, jet2_pt, jet1_cmva, jet2_cmva);
            if(flag2 == FLAG_PT_BINS){
                TLorentzVector cm = *lep_p + *lep_m;
                pt = cm.Pt();
            }
            bool pass = ((flag2 == FLAG_M_BINS && m >= var_low && m <= var_high) ||
                        (flag2 == FLAG_PT_BINS && m >= 150. && pt >= var_low && pt <= var_high))
                        && met_pt < 50.  && no_bjets;
            if(pass){
                reweight = (4./3.)*cost_st*(2. + alpha)/
                    (1. + cost_st*cost_st + alpha*(1.- cost_st*cost_st));
                n++;

                Double_t evt_weight = gen_weight * el_id_SF *el_reco_SF * pu_SF * el_HLT_SF;
                if (nJets >= 1){
                    evt_weight *= jet1_b_weight;
                }
                if (nJets >= 2){
                    evt_weight *= jet2_b_weight;
                }


                h_sym->Fill(xF, cost, evt_weight); 
                h_sym->Fill(xF, -cost, evt_weight); 

                h_asym->Fill(xF, cost, reweight * evt_weight);
                h_asym->Fill(xF, -cost, -reweight * evt_weight);

                h_cov->Fill(xF, abs(cost), abs(reweight)*evt_weight*evt_weight);

            }
        }

    }
    //h_reweights->Divide(h_count);
    //print_hist(h_reweights);

    printf("N sym is %i \n", n);
    float norm = h_sym -> Integral();
    h_sym->Scale(1./norm);
    h_cov->Scale(1./norm/norm);
    h_asym->Scale(1./norm);
    t1->ResetBranchAddresses();
    printf("MC templates generated from %i events \n \n", n);
    return 0;
}

void setHistParams(TVirtualFitter *minuit, TH2F *h_sym, TH2F *h_asym, int start_idx){
    char sym_str[20] = "sym_bin_%i_%i";
    char asym_str[20] = "asym_bin_%i_%i";
    char param_str[20];
    int offset = n_xf_bins * (n_cost_bins/2);
    for(int i=0; i< n_xf_bins; i++){
        for(int j=5; j< n_cost_bins; j++){
            int param_idx = start_idx + i*(n_cost_bins/2) + j;
            
            sprintf(param_str, sym_str, i, j);
            Float_t sym_start = h_sym->GetBinContent(i+1, j+1);
            Float_t sym_unc = h_sym->GetBinError(i+1, j+1);
            //sym_unc = 0.;
            minuit->SetParameter(param_idx,param_str, sym_start, sym_unc, 0.,0.);

            sprintf(param_str, asym_str, i, j);
            Float_t asym_start = h_asym->GetBinContent(i+1, j+1);
            Float_t asym_unc = h_asym->GetBinError(i+1, j+1);
            //asym_unc = 0.;
            minuit->SetParameter(param_idx + offset, param_str, asym_start, asym_unc, 0.,0.);
            //if(sym_start/sym_unc < 2.) printf("sym Low stats for bin %i %i, start=%.3e unc = %.3e \n", i,j, sym_start, sym_unc);
            //if(asym_start/asym_unc < 2.) printf("asym Low stats for bin %i %i, start=%.3e unc = %.3e \n", i,j, asym_start, asym_unc);
        }
    }
    return;
}
void update_templates(double *par, TH2F *h_sym, TH2F *h_asym, int start_idx){
    int offset = n_xf_bins * (n_cost_bins/2);
    int opp_cost = -5;
    for(int i=0; i< n_xf_bins; i++){
        for(int j=5; j< n_cost_bins; j++){
            int param_idx = start_idx + i*(n_cost_bins/2) + j;
            h_sym->SetBinContent(i+1, j+1, par[param_idx]);
            h_asym->SetBinContent(i+1, j+1, par[param_idx+ offset]);

            h_sym->SetBinContent(i+1, j+1+opp_cost, par[param_idx]);
            h_asym->SetBinContent(i+1, j+1+opp_cost, -par[param_idx+ offset]);
        }
    }

}

double get_hists_likelihood(TH2F *h_sym, TH2F *h_asym, TH2F *h_fitted_sym, TH2F *h_fitted_asym, double (&cov)[5][10][3]){
    //read inverse covariance matrix as ((a,c), (c,d))
    //Then the joint likelihood fcn is ax^2 + 2cxy + dy^2 where x and y are
    //differences in means
    double sum = 0;
    for(int i=0; i< n_xf_bins; i++){
        for(int j=5; j< n_cost_bins; j++){
            double a = cov[i][j][0];
            double c = cov[i][j][1];
            double d = cov[i][j][2];
            double sym_mean = h_fitted_sym->GetBinContent(i+1, j+1);
            double asym_mean = h_fitted_asym->GetBinContent(i+1, j+1);

            double sym_obs = h_sym->GetBinContent(i+1, j+1);
            double asym_obs = h_asym->GetBinContent(i+1, j+1);
            double x = sym_obs - sym_mean;
            double y = asym_obs - asym_mean;
            //printf("a=%.2e c=%.2e d=%.2e, x=%.2e y=%.2e \n", a,c,d,x,y);
            sum += a*x*x + 2*c*x*y + d*y*y;
        }
    }
    return sum;
}

void gen_cov_mat(TH2F *h_sym, TH2F *h_asym, TH2F *h_cov, double (&cov)[5][10][3]){
    //Generate inverse covariance matrix between symmetric and asym templates
    //Inverse of cov matrix of form ((a,c)(c, d)) = (1/(ad - c^2)) ((d, -c), (-c, a))
    //3rd elem of array is off diagonal term (symmetric)
    

    for(int i=0; i< n_xf_bins; i++){
        for(int j=5; j< n_cost_bins; j++){
            //errors need to be squared;
            double a = h_sym->GetBinError(i+1,j+1);
            a*=a;
            double d = h_asym->GetBinError(i+1,j+1);
            d*=d;
            double c = h_cov->GetBinContent(i+1,j+1);

            printf("a = %.3e , d = %.3e, c =%.3e \n",a,d,c);
            double norm = (a*d - c*c);
            cov[i][j][0] = d/norm;
            cov[i][j][1] = -c/norm;
            cov[i][j][2] = a/norm;
        }
    }
    return;
}




// fcn passes back f = - 2*ln(L), (where L is the likelihood of the event)
void fcn(int& npar, double* deriv, double& f, double par[], int flag){
    f=0;
    double lnL = 0.0;
    double lnL_hists = 0.0;
    int misses = 0;
    if(print) printf("\n \n \n ");

    double AFB = par[0];
    double r_back = par[1];


    if(FLAG == FLAG_ELECTRONS){
        update_templates(par, h_elel_fitted_sym, h_elel_fitted_asym, 2);
        //print_hist(h_elel_fitted_mc);
        for (int i=0; i<nElEl_DataEvents; i++){
            Double_t p_sym = get_prob(v_elel_xF[i], v_elel_cost[i], h_elel_sym);
            Double_t p_asym = get_prob(v_elel_xF[i],  v_elel_cost[i], h_elel_asym);
            Double_t p_back = get_prob(v_elel_xF[i], v_elel_cost[i], h_elel_back);



            if(p_back < 0){ 
                //Printf("P_back < 0! %.2e %.2f %.2f \n", p_back, v_xF[i], v_cost[i]);
                p_back = 1e-20;
            }

            double prob = r_back*p_back + (1 - r_back) * (p_sym + AFB*p_asym);
            if(prob > 1) printf("Warning prob is too big %.2f %.2f %.2f %.2f %.2f  \n", r_back, AFB, p_back, p_sym, p_asym);
            if(print && p_sym < 1e-20){
                misses++;
                printf(" Warning p_sym is 0 or negative! for bin xf: %0.2f cost: %1.2f \n", v_elel_xF[i], v_elel_cost[i]);
                if(p_back < 1e-20) printf("p_back Is also 0 or negative! \n");
                if(prob < 1e-20) printf("Warning prob is also 0 or negative! \n");
                p_sym = 1e-20;
                prob = r_back*p_back + (1 - r_back) * (p_sym + AFB*p_asym);
            }
            //if(prob < 1e-20) printf(" Warning prob is 0 or negative! for bin xf: %0.2f cost: %1.2f \n", v_xF[i], v_cost[i]);
            prob = max(prob, 1e-20);
            if(prob >0.0) lnL += log(prob);
        }
        lnL_hists = get_hists_likelihood(h_elel_sym, h_elel_asym, h_elel_fitted_sym, h_elel_fitted_asym, m_elel_cov);
        //printf("lnL = %.3f  hists = %.3f \n", -2.0*lnL, lnL_hists);
        f += -2.0 * (lnL) + lnL_hists;
        if(print) {
            printf("%i misses out of %i events \n\n\n", misses, nElEl_DataEvents);
            print = false;
        }
    }

    else{
        update_templates(par, h_mumu_fitted_sym, h_mumu_fitted_asym, 2);
        //print_hist(h_mumu_fitted_mc);
        for (int i=0; i<nMuMu_DataEvents; i++){
            Double_t p_sym = get_prob(v_mumu_xF[i], v_mumu_cost[i], h_mumu_sym);
            Double_t p_asym = get_prob(v_mumu_xF[i],  v_mumu_cost[i], h_mumu_asym);
            Double_t p_back = get_prob(v_mumu_xF[i], v_mumu_cost[i], h_mumu_back);



            if(p_back < 0){ 
                //Printf("P_back < 0! %.2e %.2f %.2f \n", p_back, v_xF[i], v_cost[i]);
                p_back = 1e-20;
            }

            double prob = r_back*p_back + (1 - r_back) * (p_sym + AFB*p_asym);
            if(prob > 1) printf("Warning prob is too big %.2f %.2f %.2f %.2f %.2f  \n", r_back, AFB, p_back, p_sym, p_asym);
            if(print && p_sym < 1e-20){
                misses++;
                printf(" Warning p_sym is 0 or negative! for bin xf: %0.2f cost: %1.2f \n", v_mumu_xF[i], v_mumu_cost[i]);
                if(p_back < 1e-20) printf("p_back Is also 0 or negative! \n");
                if(prob < 1e-20) printf("Warning prob is also 0 or negative! \n");
                p_sym = 1e-20;
                prob = r_back*p_back + (1 - r_back) * (p_sym + AFB*p_asym);
            }
            //if(prob < 1e-20) printf(" Warning prob is 0 or negative! for bin xf: %0.2f cost: %1.2f \n", v_xF[i], v_cost[i]);
            prob = max(prob, 1e-20);
            if(prob >0.0) lnL += log(prob);
        }
        lnL_hists = get_hists_likelihood(h_mumu_sym, h_mumu_asym, h_mumu_fitted_sym, h_mumu_fitted_asym, m_mumu_cov);
        //printf("lnL = %.3f  hists = %.3f \n", -2.0*lnL, lnL_hists);
        f += -2.0 * lnL + lnL_hists;
        if(print) {
            printf("%i misses out of %i events \n\n\n", misses, nMuMu_DataEvents);
            print = false;
        }
    }

}







void setup(){
    //setup global variables
    TH1::SetDefaultSumw2(kTRUE);

    const int n_rw_bins = 1000;
    Float_t rw_bins[n_rw_bins];
    Float_t rw_start = -1.5;
    Float_t rw_end = 1.5;
    Float_t bin_sep = (rw_end -rw_start)/n_rw_bins;
    printf("bin sep is %.4f \n", bin_sep);
    for(int i=0;i<=n_rw_bins; i++){
        rw_bins[i] = rw_start + i*bin_sep;
    }

    if(FLAG == FLAG_ELECTRONS){
        h_elel_mc_count = new TH2F("h_elel_mc_count", "Events in bins for MC templates",
                n_xf_bins, xf_bins, n_cost_bins, cost_bins);
        h_elel_mc_count->SetDirectory(0);
        h_elel_mc = new TH2F("h_elel_mc", "Events in bins for MC templates",
                n_xf_bins, xf_bins, n_cost_bins, cost_bins);
        h_elel_mc->SetDirectory(0);
        h_elel_cov = new TH2F("h_elel_cov", "Events in bins for MC templates",
                n_xf_bins, xf_bins, n_cost_bins, cost_bins);
        h_elel_cov->SetDirectory(0);
        h_elel_mc_errs = new TH2F("h_elel_mc_errs", "Events in bins for MC templates",
                n_xf_bins, xf_bins, n_cost_bins, cost_bins);
        h_elel_mc_errs->SetDirectory(0);
        h_elel_sym_count = new TH2F("h_elel_sym_count", "Events in bins for MC templates",
                n_xf_bins, xf_bins, n_cost_bins, cost_bins);
        h_elel_sym_count->SetDirectory(0);
        h_elel_sym = new TH2F("h_elel_sym", "Symmetric template of mc",
                n_xf_bins, xf_bins, n_cost_bins, cost_bins);
        h_elel_sym->SetDirectory(0);
        h_elel_asym = new TH2F("h_elel_asym", "Asymmetric template of mc",
                n_xf_bins, xf_bins, n_cost_bins, cost_bins);
        h_elel_asym->SetDirectory(0);
        h_elel_back = new TH2F("h_elel_back", "Combined background template",
                n_xf_bins, xf_bins, n_cost_bins, cost_bins);
        h_elel_back->SetDirectory(0);
        h_elel_data = new TH2F("h_elel_data", "Data template of (x_f, cost_r)",
                n_xf_bins, xf_bins, n_cost_bins, cost_bins);
        h_elel_data->SetDirectory(0);


        gen_finite_mc_template_v2(t_elel_mc, alpha, h_elel_sym, h_elel_asym, h_elel_cov,
                m_low, m_high, FLAG_ELECTRONS);
        h_elel_fitted_sym = (TH2F *) h_elel_sym->Clone("h_elel_fitted_sym");
        h_elel_fitted_asym = (TH2F *) h_elel_asym->Clone("h_elel_fitted_asym");
        gen_cov_mat(h_elel_sym, h_elel_asym, h_elel_cov, m_elel_cov);
        //gen_mc_template(t_elel_mc, alpha, h_elel_sym, h_elel_asym, h_elel_sym_count, m_low, m_high, FLAG_ELECTRONS);
        TTree *elel_ts[2] = {t_elel_back, t_elel_nosig};

        gen_fakes_template(t_elel_WJets, t_elel_QCD, t_elel_WJets_contam, t_elel_QCD_contam, h_elel_back, m_low, m_high, FLAG_ELECTRONS);
        gen_combined_background_template(2, elel_ts, h_elel_back, m_low, m_high, FLAG_ELECTRONS);

        nElEl_DataEvents = gen_data_template(t_elel_data, h_elel_data,  &v_elel_xF, &v_elel_cost, m_low, m_high, FLAG_ELECTRONS);
    }
    else{
        h_mumu_mc_count = new TH2F("h_mumu_mc_count", "Events in bins for MC templates",
                n_xf_bins, xf_bins, n_cost_bins, cost_bins);
        h_mumu_mc_count->SetDirectory(0);
        h_mumu_mc = new TH2F("h_mumu_mc", "Events in bins for MC templates",
                n_xf_bins, xf_bins, n_cost_bins, cost_bins);
        h_mumu_mc->SetDirectory(0);
        h_mumu_cov = new TH2F("h_mumu_cov", "Events in bins for MC templates",
                n_xf_bins, xf_bins, n_cost_bins, cost_bins);
        h_mumu_cov->SetDirectory(0);
        h_mumu_mc_errs = new TH2F("h_mumu_mc_errs", "Events in bins for MC templates",
                n_xf_bins, xf_bins, n_cost_bins, cost_bins);
        h_mumu_mc_errs->SetDirectory(0);
        h_mumu_sym_count = new TH2F("h_mumu_sym_count", "Events in bins for MC templates",
                n_xf_bins, xf_bins, n_cost_bins, cost_bins);
        h_mumu_sym_count->SetDirectory(0);
        h_mumu_sym = new TH2F("h_mumu_sym", "Symmetric template of mc",
                n_xf_bins, xf_bins, n_cost_bins, cost_bins);
        h_mumu_sym->SetDirectory(0);
        h_mumu_asym = new TH2F("h_mumu_asym", "Asymmetric template of mc",
                n_xf_bins, xf_bins, n_cost_bins, cost_bins);
        h_mumu_asym->SetDirectory(0);
        h_mumu_back = new TH2F("h_mumu_back", "Combined background template",
                n_xf_bins, xf_bins, n_cost_bins, cost_bins);
        h_mumu_back->SetDirectory(0);
        h_mumu_data = new TH2F("h_mumu_data", "Data template of (x_f, cost_r)",
                n_xf_bins, xf_bins, n_cost_bins, cost_bins);
        h_mumu_data->SetDirectory(0);

        //printf("size %i \n", (int) v_mumu_xF.size());
        v_mumu_xF.clear();
        v_mumu_cost.clear();

        gen_finite_mc_template_v2(t_mumu_mc, alpha, h_mumu_sym, h_mumu_asym, h_mumu_cov,
                m_low, m_high, FLAG_MUONS);
        h_mumu_fitted_sym = (TH2F *) h_mumu_sym->Clone("h_mumu_fitted_sym");
        h_mumu_fitted_asym = (TH2F *) h_mumu_asym->Clone("h_mumu_fitted_asym");
        gen_cov_mat(h_mumu_sym, h_mumu_asym, h_mumu_cov, m_mumu_cov);
        //gen_mc_template(t_mumu_mc, alpha, h_mumu_sym, h_mumu_asym, h_mumu_sym_count, m_low, m_high, FLAG_MUONS);
        TTree *mumu_ts[2] = {t_mumu_back, t_mumu_nosig};

        gen_fakes_template(t_mumu_WJets, t_mumu_QCD, t_mumu_WJets_contam, t_mumu_QCD_contam, h_mumu_back, m_low, m_high, FLAG_MUONS);
        gen_combined_background_template(2, mumu_ts, h_mumu_back, m_low, m_high, FLAG_MUONS);
        nMuMu_DataEvents = gen_data_template(t_mumu_data, h_mumu_data,  &v_mumu_xF, &v_mumu_cost, m_low, m_high, FLAG_MUONS);


    }
    printf("Finishing setup \n");
    return;
}

void cleanup(){
    //delete h_mc_count;
    //delete h_sym_count;
    //delete h_sym;
    //delete h_asym;
    //delete h_back;
    //delete h_data;
    v_elel_cost.clear();
    v_elel_xF.clear();
    v_mumu_cost.clear();
    v_mumu_xF.clear();
    printf("Finishing cleanup\n");
}

void finite_mc_fit_v2(){
    Double_t AFB_fit[n_m_bins], AFB_fit_err[n_m_bins], r_elel_back_fit[n_m_bins], r_elel_back_fit_err[n_m_bins], 
             r_mumu_back_fit[n_m_bins], r_mumu_back_fit_err[n_m_bins];

    init();
    TTree *tout= new TTree("T_fit_res", "Tree with Fit Results");
    tout->SetDirectory(0);

    Double_t AFB, AFB_err, r_elel_back, r_elel_back_err, r_mumu_back, r_mumu_back_err;

    tout->Branch("var_low", &m_low);
    tout->Branch("var_high", &m_high);
    tout->Branch("nElElEvents", &nElEl_DataEvents);
    tout->Branch("nMuMuEvents", &nMuMu_DataEvents);
    tout->Branch("AFB", &AFB);
    tout->Branch("AFB_err", &AFB_err);
    tout->Branch("r_elel_back", &r_elel_back);
    tout->Branch("r_elel_back_err", &r_elel_back_err);
    tout->Branch("r_mumu_back", &r_mumu_back);
    tout->Branch("r_mumu_back_err", &r_mumu_back_err);

    unsigned int nElElEvents[n_m_bins];
    unsigned int nMuMuEvents[n_m_bins];

    float *AFB_starts , *r_back_starts;
    float mumu_AFB_starts[] = {0.59, 0.608, 0.634, 0.642, 0.584, 0.605};
    float elel_AFB_starts[] = {0.64, 0.617, 0.589, 0.567, 0.535, 0.518};
    float r_mumu_back_starts[] = {0.110, 0.175, 0.183, 0.186, 0.166, 0.095};
    float r_elel_back_starts[] = {0.087, 0.156, 0.211, 0.185, 0.171, 0.128};
    if(FLAG == FLAG_ELECTRONS){
        AFB_starts = elel_AFB_starts;
        r_back_starts = r_elel_back_starts;
    }
    else{
        AFB_starts = mumu_AFB_starts;
        r_back_starts = r_mumu_back_starts;
    }


    for(int i=0; i<6; i++){
        printf("Starting loop \n");
        m_low = m_bins[i];
        m_high = m_bins[i+1];
        alpha = alphas[i];

        setup();



        float AFB_start = AFB_starts[i];
        float AFB_start_error = 0.04;
        float AFB_max = 0.75;
        float r_back_start = r_back_starts[i];
        float r_back_start_error = 0.04;
        float r_back_max = 0.6;

        int n_params = 2+ (n_cost_bins * n_xf_bins);

        TVirtualFitter * minuit = TVirtualFitter::Fitter(0,n_params);
        minuit->SetFCN(fcn);
        minuit->SetParameter(0,"AFB", AFB_start, AFB_start_error, -AFB_max, AFB_max);
        minuit->SetParameter(1,"r_back", r_back_start, r_back_start_error, 0., 0.5);
        if(FLAG == FLAG_ELECTRONS){
            setHistParams(minuit, h_elel_sym, h_elel_asym, 2);
        }
        else setHistParams(minuit, h_mumu_sym, h_mumu_sym, 2);
        Double_t arglist[100];
        arglist[0] = 10000.;
        minuit->ExecuteCommand("MIGRAD", arglist,0);
        Double_t up = 1.0;
        minuit->SetErrorDef(up);
        arglist[0] = 0.;
        //minuit->ExecuteCommand("MINOS", arglist, 0);



        AFB_fit[i] = minuit->GetParameter(0); 
        AFB_fit_err[i] = minuit->GetParError(0);
        AFB= AFB_fit[i];
        AFB_err = AFB_fit_err[i];
        if(FLAG == FLAG_ELECTRONS){
            r_elel_back_fit[i] = minuit->GetParameter(1); 
            r_elel_back_fit_err[i] = minuit->GetParError(1);
            r_elel_back= r_elel_back_fit[i];
            r_elel_back_err = r_elel_back_fit_err[i];

            nElElEvents[i] = nElEl_DataEvents;
        }
        else{
            nMuMuEvents[i] = nMuMu_DataEvents;

            r_mumu_back_fit[i] = minuit->GetParameter(1); 
            r_mumu_back_fit_err[i] = minuit->GetParError(1);


            r_mumu_back= r_mumu_back_fit[i];
            r_mumu_back_err = r_mumu_back_fit_err[i];
        }
        tout->Fill();

        cleanup();
    }
    TFile *fout;
    if(FLAG == FLAG_MUONS) fout = TFile::Open(mumu_fout_name, "RECREATE");
    if(FLAG == FLAG_ELECTRONS) fout = TFile::Open(elel_fout_name, "RECREATE");
    fout->cd();
    tout->Write();
    fout->Close();
    if(do_both){
        
        do_both = false;
        FLAG = 1 - FLAG;
        printf("Starting 2nd run \n");
        finite_mc_fit_v2();
        FLAG = 1 - FLAG;
    }
    if(FLAG == FLAG_MUONS){
        printf("MuMu fit results\n");
    }
    if(FLAG == FLAG_ELECTRONS){
        printf("ElEl fit results\n");
    }
    for(int i=0; i<6; i++){
        printf("\n Fit on M=[%.0f, %.0f], %i ElEl Events, %i MuMu Events: AFB = %0.5f +/- %0.5f r_elel_back = %0.3f +/- %0.3f r_mumu_back =%0.3f +/- %0.3f \n", 
                    m_bins[i], m_bins[i+1], nElElEvents[i], nMuMuEvents[i], AFB_fit[i], AFB_fit_err[i], r_elel_back_fit[i], r_elel_back_fit_err[i], r_mumu_back_fit[i], r_mumu_back_fit_err[i]);

    }

}


int main(){
    finite_mc_fit_v2();
    return 0;
}
