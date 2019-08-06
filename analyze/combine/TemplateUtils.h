#ifndef TemplateUtils
#define TemplateUtils
#include "RooWorkspace.h"
#include "RooRealVar.h"
#include "RooDataHist.h"
#include "RooFormulaVar.h"
#include "HiggsAnalysis/CombinedLimit/interface/RooParametricHist.h"
#include "../../Utils/TemplateMaker_systematics.C"
#include "../../Utils/root_files.h"

Double_t m_low;
Double_t m_high;

RooWorkspace *w;
RooRealVar *var = new RooRealVar("var", "var", 0.,n_cost_bins * n_xf_bins);
RooRealVar *var_ss = new RooRealVar("var_ss", "var_ss", 0.,(n_cost_bins/2) * n_xf_bins);
RooRealVar *Rqcd_ee_ss = new RooRealVar("Rqcd_ee_ss", "ee QCD normalization", 1, 0., 10.);
RooRealVar *Rqcd_mumu_ss = new RooRealVar("Rqcd_mumu_ss", "mumu QCD normalization", 1, 0., 10.);

TH1F *h_dummy = new TH1F("h_dummy", "", n_cost_bins*n_xf_bins, 0, n_cost_bins*n_xf_bins);



bool do_emu_scale = false;
bool do_RC = true;

void write_roo_hist(TH1F *h, RooRealVar *my_var){
    RooDataHist r(h->GetName(), h->GetName(), *my_var, h);
    delete h;
    w->import(r);
}
void symmetrize2d(TH2F *h_2d){
    int n_xf_bins = h_2d->GetNbinsX();
    int n_cost_bins = h_2d->GetNbinsY();

    for(int i=1; i<=n_xf_bins; i++){
        for(int j=1; j<= n_cost_bins/2; j++){
            float content = h_2d->GetBinContent(i,j);
            float error = h_2d->GetBinError(i,j);

            int opp_j = (n_cost_bins + 1) -j;
            float content2 = h_2d->GetBinContent(i,opp_j);
            float error2 = h_2d->GetBinError(i,opp_j);

            float new_content = (content + content2)/2.0;
            float new_error = pow((error*error + error2*error2), 0.5)/2.0;
            h_2d->SetBinContent(i,j, new_content);
            h_2d->SetBinContent(i,opp_j, new_content);

            h_2d->SetBinError(i,j, new_error);
            h_2d->SetBinError(i,opp_j, new_error);


        }
    }
}

TH1F* convert2d(TH2F *h_2d){
    int n_xf_bins = h_2d->GetNbinsX();
    int n_cost_bins = h_2d->GetNbinsY();

    TH1F *h_1d = new TH1F(h_2d->GetName(), "",  n_xf_bins * n_cost_bins, 0, n_xf_bins*n_cost_bins);
    for(int i=1; i<=n_xf_bins; i++){
        for(int j=1; j<= n_cost_bins; j++){
            float content = h_2d->GetBinContent(i,j);
            float error = h_2d->GetBinError(i,j);
            int gbin = (i-1)*n_cost_bins + j;
            //printf("gbin %i: i j %i %i \n", gbin, i, j);
            h_1d->SetBinContent(gbin, content);
            h_1d->SetBinError(gbin, error);
        }
    }
    return h_1d;
}

int TwoDToOneDIdx(int nYbins, int i, int j){
    nYbins = nYbins/2;
    j = ((j-1) % nYbins) +1;
    int g_idx = (i-1)*nYbins + j;
    return g_idx;
}
int mirrorIdx(int nYbins, int j){
    j = ((j-1) % nYbins);
    int m_idx = nYbins - j;
    return m_idx;
}

void TwoDToSymIdxs(int nYbins, int i, int j, int &sym1_idx, int &sym2_idx){
    sym1_idx = (i-1)*nYbins + j;
    sym2_idx = (i-1)*nYbins + (nYbins - (j -1));
}




char* replace(char* str, char* a, char* b)
{
    int len  = strlen(str);
    int lena = strlen(a), lenb = strlen(b);
    for (char* p = str; (p = strstr(p, a)); ++p) {
        if (lena != lenb) // shift end as needed
            memmove(p+lenb, p+lena,
                    len - (p - str) + lenb);
        memcpy(p, b, lenb);
    }
    return str;
}
#endif
