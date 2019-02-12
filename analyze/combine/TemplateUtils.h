#ifndef TemplateUtils
#define TemplateUtils
#include "RooWorkspace.h"
#include "HiggsAnalysis/CombinedLimit/interface/RooParametricHist.h"
#include "../TemplateMaker_systematics.C"
#include "../AFB_fit/FitUtils.C"

Double_t m_low;
Double_t m_high;

RooWorkspace *w;
RooRealVar *var = new RooRealVar("var", "var", 0.,n_cost_bins * n_xf_bins);
RooRealVar *Rqcd_ee_ss = new RooRealVar("Rqcd_ee_ss", "ee QCD normalization", 1, 0., 10.);
RooRealVar *Rqcd_mumu_ss = new RooRealVar("Rqcd_mumu_ss", "mumu QCD normalization", 1, 0., 10.);



bool do_emu_scale = false;
bool do_RC = true;

void write_roo_hist(TH1F *h){
    RooDataHist r(h->GetName(), h->GetName(), *var, h);
    w->import(r);
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
