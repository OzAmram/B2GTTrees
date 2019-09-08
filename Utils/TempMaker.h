
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <cmath>
#include <cfloat>
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
#include "TH2.h"
#include "TH3.h"
#include "TTree.h"
#include "TFitter.h"
#include "TSystem.h"
#include "Math/Functor.h"
#include "ScaleFactors.C"
#include "BTagUtils.C"
#include "HistMaker.C"


class TempMaker{
    public:
        TempMaker(TTree *t, bool isdata = false, int year = 2016);
        ~TempMaker();
        void setup();
        void setup_systematic(const string &s_label);
        void getEvent(int i);
        void doCorrections();
        Double_t getEvtWeight();
        void finish();



        int year = 2016;
        bool is_data = false;
        bool is_gen_level = false;
        bool is_one_iso = false;

        bool do_electrons = false;
        bool do_muons = false;
        bool do_emu = false;
        bool do_RC = true;

        int iso_lep = -1;        
        float el_lumi = 0.;

        string sys_label = string("");
        TTree *t_in;

        Double_t m, xF, cost, gen_weight, reweight_a, reweight_s, reweight_alpha, jet1_cmva, jet2_cmva, cost_st, gen_cost;
        Double_t evt_weight;
        Double_t era1_HLT_SF, era1_iso_SF, era1_id_SF;
        Double_t era2_HLT_SF, era2_iso_SF, era2_id_SF;
        Double_t el_id_SF, el_reco_SF, pu_SF, pu_SF_up, pu_SF_down, el_HLT_SF;
        Double_t jet1_pt, jet2_pt, jet1_b_weight, jet2_b_weight, jet1_eta, jet2_eta;
        Double_t mu1_pt, mu1_eta, mu2_pt, mu2_eta;
        Double_t el1_pt, el1_eta, el2_pt, el2_eta;
        Double_t mu_R_up, mu_R_down, mu_F_up, mu_F_down, 
                 mu_RF_up, mu_RF_down, pdf_up, pdf_down;
        Double_t mu_p_SF, mu_m_SF, mu_p_SF_up, mu_m_SF_up, mu_p_SF_down, mu_m_SF_down, alphaS_up, alphaS_down;
        Float_t cost_pt, met_pt, el1_charge, el2_charge, mu1_charge, mu2_charge;

        Float_t pdf_weights[60];
        TLorentzVector *lep_p=0;
        TLorentzVector *lep_m=0;
        TLorentzVector *gen_lep_p=0;
        TLorentzVector *gen_lep_m=0;
        Double_t pt;
        Int_t nJets, pu_NtrueInt, jet1_flavour, jet2_flavour;

        Int_t has_no_bjets = 1;
        bool not_cosmic = true;

        Long64_t nEntries;


        //systematics flags

        Double_t one = 1.0;
        Double_t *systematic = &one;
        int sys_shift = 0;

        // do_sys: 0 = nominal, 1= var up, -1 = var down
        int do_pdf_sys = 0;
        int do_btag_sys = 0;
        int do_pileup_sys = 0;

        int do_muHLT_sys = 0;
        int do_muID_sys = 0;
        int do_muISO_sys = 0;
        int do_muRC_sys = 0;

        int do_elID_sys = 0;
        int do_elHLT_sys = 0;
        int do_elRECO_sys = 0;
        int do_elScale_sys = 0;
        int do_elSmear_sys = 0;
        float elp_rescale, elm_rescale;
};


//global systematics stuff so not tied to single instance of class
el_SFs el_SF;
mu_SFs era1_SFs, era2_SFs;
pileup_systematics pu_sys;

BTag_readers b_reader;
BTag_effs btag_effs;

void setup_all_SFs(int year){
    printf("Setting up SF's \n");
    setup_btag_SFs(&b_reader, &btag_effs, year);
    setup_el_SF(&el_SF, year);
    setup_mu_SFs(&era1_SFs, &era2_SFs,  year);
    setup_pileup_systematic(&pu_sys, year); 
}


