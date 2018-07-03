



#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <iostream>
#include <cstring>
#include <algorithm>
#include <vector>

#include "TROOT.h"
#include "TFile.h"
#include "TObject.h"
#include "TH1F.h"
#include "TProfile.h"
#include "TStyle.h"
#include "TCanvas.h"
#include "TPostScript.h"
#include "TLorentzVector.h"
#include "TVector3.h"
#include "TFitter.h"
#include "TSystem.h"
#include "Math/Functor.h"
#include "../../analyze/TemplateMaker.C"
#include "../tdrstyle.C"
#include "../CMS_lumi.C"
#include "root_files.h"



float m_low;
float m_high;
int FLAG = FLAG_ELECTRONS;

bool print = true;


Double_t med_btag = 0.4432;

TH2F *h_asym, *h_sym, *h_back,  *h_data, *h_mc;
TH2F *h_mc_count, *h_sym_count;

int n_xf_bins = 5;
Float_t xf_bins[] = {0., 0.02, 0.04, 0.07, 0.10, 1.0};
//int n_cost_bins = 8;
//Float_t cost_bins[] = {-1.0, -.75, -.5, -.25, 0., 0.25, 0.5,  0.75, 1.0};
int n_cost_bins = 10;
Float_t cost_bins[] = {-1.0, -.8, -.6, -.4, -.2, 0., 0.2, 0.4, 0.6, 0.8, 1.0};
//int n_cost_bins = 12;
//Float_t cost_bins[] = {-1.0, -0.8333, -0.6667, -0.5, -0.3333, -0.1667, 0., 0.1667, 0.3333, 0.5, 0.6667, 0.8333, 1.0};
//int n_cost_bins = 14;
//Float_t cost_bins[] = {-1.0, -.857, -.714, -.571, -.429, -0.286, -.143,  0., 0.143, .286, 0.429, 0.571, 0.714, 0.857, 1.0};
int n_m_bins = 6;
int n_pt_bins = 6;

void setup(){
    //setup global variables
    //TH1::SetDefaultSumw2(kTRUE);
    printf("Starting setup \n");
    h_mc_count = new TH2F("h_mc_count", "Events in bins for MC templates",
            n_xf_bins, xf_bins, n_cost_bins, cost_bins);
    h_mc_count->SetDirectory(0);
    h_sym_count = new TH2F("h_sym_count", "Events in bins for MC templates",
            n_xf_bins, xf_bins, n_cost_bins, cost_bins);
    h_sym_count->SetDirectory(0);
    h_sym = new TH2F("h_sym", "Symmetric template of mc (xF, cost_r) xF>0.15",
            n_xf_bins, xf_bins, n_cost_bins, cost_bins);
    h_sym->SetDirectory(0);
    h_asym = new TH2F("h_asym", "Asymmetric template of mc (xF cost_r) xF > 0.15",
            n_xf_bins, xf_bins, n_cost_bins, cost_bins);
    h_asym->SetDirectory(0);
    h_back = new TH2F("h_back", "Combined background template",
            n_xf_bins, xf_bins, n_cost_bins, cost_bins);

    h_back->SetDirectory(0);
    h_data = new TH2F("h_data", "Data template of (x_f, cost_r) xF > 0.15",
            n_xf_bins, xf_bins, n_cost_bins, cost_bins);

    h_data->SetDirectory(0);
    printf("Generating templates \n");

    gen_mc_template(t_elel_mc, alpha, h_sym, h_asym, h_sym_count, m_low, m_high, FLAG);
    TTree *ts[2] = {t_elel_back, t_elel_nosig};

    gen_fakes_template(t_elel_WJets, t_elel_QCD, t_elel_WJets_contam, t_elel_QCD_contam, h_back, m_low, m_high, FLAG);
    gen_combined_background_template(2, ts, h_back, m_low, m_high, FLAG);

    nDataEvents = gen_data_template(t_elel_data, h_data, &v_xF, &v_cost, m_low, m_high, FLAG);
    printf("\n\n\n Printing MC counts in each bin:\n");
    for(int i=1; i<=n_xf_bins; i++){
        for(int j=1; j<=n_cost_bins; j++){
            printf("%.0f ", h_qcd->GetBinContent(i,j));
        }
        printf("\n");
    }

    printf("Finishing setup \n");
    return;
}

void draw_template(){
    for(int i=0; i<n_m_bins; i++){
        printf("Starting loop \n");
        m_low = m_bins[i];
        m_high = m_bins[i+1];
        alpha = alphas[i];

        setup();

        printf("Integrals are %f %f %f %f  \n", h_data->Integral(), h_sym->Integral(), 
                                               h_asym->Integral(), h_back->Integral() );
    }
}
