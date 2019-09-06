#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <iostream>
#include <cstring>
#include <algorithm>
#include "TFile.h"
#include "TRandom3.h"
#include "TTree.h"
#include "TLorentzVector.h"
#include "ScaleFactors.C"
#include "HistMaker.C"
#include "roccor_Run2_v3/RoccoR.cc"


#define GEN_SIZE 4000
#define MU_SIZE 100
#define EL_SIZE 100
#define JET_SIZE 60
#define MAX_SAMPLES 20

const double root2 = sqrt(2);
double Ebeam = 6500.;
double Pbeam = sqrt(Ebeam*Ebeam - 0.938*0.938);
const float mu_mass = 0.1056; // in GEV


bool is_empty_line(const char *s) {
    while (*s != '\0') {
        if (!isspace(*s))
            return false;
        s++;
    }
    return true;
}



class NTupleReader{
    public:
        NTupleReader(const char *fin_name, const char *fout_name, bool is_data_);
        ~NTupleReader();

        void setupSFs();
        void setupRC();
        void setupOutputTree(char treeName[100]);
        void getEvent(int i);
        void fillEvent();
        void fillEventSFs();
        void fillEventRC();
        bool getNextFile();
        void finish();
        bool parseGenParts(bool PRINT);




        unsigned int nFiles;
        Double_t normalization = 1.0;
        Double_t norms[MAX_SAMPLES]; // computed normalizations to apply to each event in a sample (based on xsection and total weight)
        FILE *root_files;
        TFile *fout;

        RoccoR rc;
        TRandom *rand;
        mu_SFs era1, era2;
        pileup_systematics pu_sys;
        el_SFs el_SF;

        const float mu_iso_cut = 0.15; //tight PF based iso

        int year = 2016;


        int nJobs =1;
        int iJob = 0;
        int fileCount =0;
        bool do_samesign = false;
        bool do_muons = false;
        bool do_electrons = false;
        bool do_emu = false;
        bool is_data = false;
        bool do_SFs = false;
        bool do_RC = false;
        bool RC_from_gen = false;

        unsigned int nEvents=0;
        unsigned int nSignal = 0;
        unsigned int nQQ=0;
        unsigned int nQQb=0;
        unsigned int nQGlu=0;
        unsigned int nGluGlu=0;
        unsigned int nTauTau=0;
        unsigned int nFailedID=0;
        int event_idx = 0;

        TTree *outTrees[10];
        int nOutTrees = 0;
        TFile *fin;
        TTree *tin;
        Long64_t tin_nEntries;

        float bjet_med_cut = 0.;

        Double_t cm_m, xF, cost, cost_r, cost_st, mu1_pt, mu2_pt, mu1_eta, mu2_eta, jet1_pt, jet2_pt, jet1_eta, jet2_eta, 
                 gen_weight, jet1_csv, jet1_btag, jet2_csv, jet2_btag, gen_m;
        Double_t mu_p_SF, mu_m_SF, mu_p_SF_alt, mu_m_SF_alt, mu_p_SF_up, mu_p_SF_down, mu_m_SF_up, mu_m_SF_down;
        Double_t era1_HLT_SF, era1_iso_SF, era1_id_SF, era2_HLT_SF, era2_iso_SF, era2_id_SF,
                 era1_trk_SF, era2_trk_SF,
                 jet1_b_weight, jet2_b_weight, pu_SF, pu_SF_up, pu_SF_down;

        Double_t el1_pt, el2_pt, el1_eta, el2_eta;
        Double_t el_id_SF, el_reco_SF, el_HLT_SF;
        Float_t elp_scale_stat_up, elp_scale_stat_down, elp_scale_gain_up, elp_scale_gain_down, elp_scale_syst_up, elp_scale_syst_down, 
                elm_scale_stat_up, elm_scale_stat_down, elm_scale_gain_up, elm_scale_gain_down, elm_scale_syst_up, elm_scale_syst_down, 
             elp_smear_up, elp_smear_down, elm_smear_up, elm_smear_down;


        Double_t mu_R_up, mu_R_down, mu_F_up, mu_F_down, mu_RF_up, mu_RF_down, alpha_up, alpha_down;
        Int_t nJets, jet1_flavour, jet2_flavour, pu_NtrueInt, has_nobjets;
        Bool_t is_tau_event;
        Float_t met_pt, mu1_charge, mu2_charge, el1_charge, el2_charge;
        TLorentzVector cm, gen_cm;
        TLorentzVector mu_p, mu_m, gen_mu_p_vec, gen_mu_m_vec;
        TLorentzVector el_p, el_m, gen_el_p_vec, gen_el_m_vec;
        TLorentzVector el, mu;
        Float_t scale_Weights[10], pdf_weights[60], alpha_weights[2];




        UInt_t el_size, mu_size, gen_size, jet_size, met_size, alphas_size, scale_size, pdf_size;
        Int_t gen_id[GEN_SIZE], gen_status[GEN_SIZE];
        Int_t  gen_Mom0ID[GEN_SIZE], gen_Mom0Status[GEN_SIZE], gen_Mom1ID[GEN_SIZE], gen_Mom1Status[GEN_SIZE];
        Int_t  gen_Dau0ID[GEN_SIZE], gen_Dau0Status[GEN_SIZE], gen_Dau1ID[GEN_SIZE], gen_Dau1Status[GEN_SIZE];
        Float_t gen_Pt[GEN_SIZE], gen_Eta[GEN_SIZE], gen_Phi[GEN_SIZE], gen_E[GEN_SIZE];

        Float_t mu_Pt[MU_SIZE], mu_Eta[MU_SIZE], mu_Phi[MU_SIZE], mu_E[MU_SIZE], 
                mu_Charge[MU_SIZE], mu_IsHighPtMuon[MU_SIZE], mu_IsTightMuon[MU_SIZE];
        Float_t mu_SumChargedHadronPt[MU_SIZE], mu_SumNeutralHadronPt[MU_SIZE], mu_SumPUPt[MU_SIZE], mu_SumPhotonPt[MU_SIZE],
                mu_NumberTrackerLayers[MU_SIZE];

        Float_t mu_PFIso[MU_SIZE];
        Float_t jet_Pt[JET_SIZE], jet_Eta[JET_SIZE], jet_Phi[JET_SIZE], jet_E[JET_SIZE],
                jet_btag[JET_SIZE], jet_hadronflavour[JET_SIZE], jet_genPt[JET_SIZE];

        Float_t el_Pt[EL_SIZE], el_Eta[EL_SIZE], el_Phi[EL_SIZE], el_E[EL_SIZE], 
                el_Charge[EL_SIZE], el_SCEta[EL_SIZE];

        Int_t el_IDMedium[EL_SIZE], el_IDMedium_NoIso[EL_SIZE];

        Float_t el_ScaleCorr[EL_SIZE], el_ScaleCorrStatUp[EL_SIZE], el_ScaleCorrStatDown[EL_SIZE],
                                       el_ScaleCorrGainUp[EL_SIZE], el_ScaleCorrGainDown[EL_SIZE],
                                       el_ScaleCorrSystUp[EL_SIZE], el_ScaleCorrSystDown[EL_SIZE],
            el_ScaleSmearDown[EL_SIZE], el_ScaleSmearUp[EL_SIZE];






        Float_t evt_Gen_Weight;
        Int_t HLT_IsoMu24, HLT_IsoTkMu24, HLT_IsoMu27, HLT_IsoTkMu27, HLT_El27, HLT_El35, HLT_El32, HLT_doubleEl23;


        int elp_index, elm_index;
        bool good_sign, opp_sign, good_trigger, dimuon_id, emu_ids,
             mu_iso0, mu_iso1, dielec_id, el_iso0, el_iso1;

        bool signal_event, failed_match;

        Double_t quark_dir_eta;
 private:
        void applyRC();
};
        


