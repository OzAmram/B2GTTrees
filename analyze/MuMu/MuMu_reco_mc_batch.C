#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <iostream>
#include <cstring>
#include <algorithm>
#include "TFile.h"
#include "TRandom3.h"
#include "../ScaleFactors.C"

#define GEN_SIZE 4000
#define MU_SIZE 100
#define JET_SIZE 20
#define MAX_SAMPLES 20

const double root2 = sqrt(2);
double Ebeam = 6500.;
double Pbeam = sqrt(Ebeam*Ebeam - 0.938*0.938);

char *filename("DY_files_unbinned_oct23.txt");
const TString fout_name("output_files/MuMu_DY_zpeak_unbinned_june14.root");
const bool PRINT=false;

const bool data_2016 = true;




bool is_empty_line(const char *s) {
    while (*s != '\0') {
        if (!isspace(*s))
            return false;
        s++;
    }
    return true;
}

void compute_norms(Double_t *norms, unsigned int *nFiles){
    Double_t sample_weight = 0;
    Double_t sample_xsec = 0;
    unsigned int sample_i=0;

    char lines[300];
    FILE *root_files = fopen(filename, "r");
    while(fgets(lines, 300, root_files)){
        if(lines[0] == '#' || is_empty_line(lines)) continue; // comment line
        if(lines[0] == '!'){
            //sample header line
            if(sample_xsec >0 && sample_weight >0){
                //end of old sample, record normalization
                norms[sample_i] = sample_xsec/sample_weight;
                printf("sample %i had xsec %f and weight %e and got normalization %e \n", sample_i, sample_xsec, sample_weight, norms[sample_i]);
                sample_weight = 0;
                sample_xsec = 0;
            }
            int sample_idx;
            float xsec;
            int nparams = sscanf(lines, "! idx = %i xsec = %f \n", &sample_idx, &xsec);
            if(nparams < 2 || sample_idx >= MAX_SAMPLES){
                printf("ERROR: Unable to parse sample header. Exiting");
                exit(EXIT_FAILURE);
            }
            sample_i = sample_idx;
            sample_xsec = xsec;
        }
        else{//root file

            char * end;
            //remove trailing whitespace
            end = lines + strlen(lines) - 1;
            while(end > lines && isspace((char)*end)) end--;
            // Write new null terminator
            *(end+1) = 0;

            TFile *f1=  TFile::Open(lines);
            (*nFiles)++;
            f1->cd("EventCounter");
            TDirectory *subdir = gDirectory;
            TH1D *t1 = (TH1D *)subdir->Get("totweight");
            sample_weight += t1->GetSumOfWeights();
            f1->Close();
        }
    }
    norms[sample_i] = sample_xsec/sample_weight;
    printf("sample %i had xsec %f and weight %e and got normalization %e \n", sample_i, sample_xsec, sample_weight, norms[sample_i]);
    fclose(root_files);

}




void MuMu_reco_mc_batch()
{

    Double_t norms[MAX_SAMPLES]; // computed normalizations to apply to each event in a sample (based on xsection and total weight)
    unsigned int nFiles = 0;
    printf("Computing normalizations for each sample \n");
    compute_norms(norms, &nFiles);
    printf("Done with normalizations \n\n\n");

    mu_SFs runs_bcdef, runs_gh;
    pileup_SFs pu_SFs;
    BTag_readers b_reader;
    BTag_effs btag_effs;
    //RoccoR  rc("rcdata.2016.v3"); //directory path as input for now; initialize only once, contains all variations
    TRandom *rand = new TRandom3();


    //separate SFs for runs BCDEF and GH
    setup_SFs(&runs_bcdef, &runs_gh, &b_reader, &btag_effs, &pu_SFs);
    printf("Retrieved Scale Factors \n\n");


    TFile *fout = TFile::Open(fout_name, "RECREATE");
    TTree *t_signal= new TTree("T_data", "Tree with asym events (qq bar, qg)");
    //t_signal->SetDirectory(0);
    Double_t cm_m, xF, cost_r, cost_st, mu1_pt, mu2_pt, mu1_eta, mu2_eta, jet1_pt, jet2_pt, jet1_eta, jet2_eta, deltaC, 
             gen_weight, jet1_csv, jet1_cmva, jet2_csv, jet2_cmva, gen_m;
    Double_t bcdef_HLT_SF, bcdef_iso_SF, bcdef_id_SF, gh_HLT_SF, gh_iso_SF, gh_id_SF,
             bcdef_trk_SF, gh_trk_SF,
             jet1_b_weight, jet2_b_weight, pu_SF;
    Double_t mu_R_up, mu_R_down, mu_F_up, mu_F_down, mu_RF_up, mu_RF_down, pdf_up, pdf_down;
    Int_t nJets, jet1_flavour, jet2_flavour, pu_NtrueInt;
    Bool_t is_tau_event;
    Float_t met_pt;
    TLorentzVector mu_p, mu_m, cm, gen_mu_p_vec, gen_mu_m_vec;
    Float_t scale_Weights[10], pdf_weights[100];

    t_signal->Branch("m", &cm_m, "m/D");
    t_signal->Branch("gen_m", &gen_m, "m/D");
    t_signal->Branch("xF", &xF, "xF/D");
    t_signal->Branch("cost", &cost_r, "cost/D");
    t_signal->Branch("cost_st", &cost_st, "cost_st/D");
    t_signal->Branch("mu1_pt", &mu1_pt, "mu1_pt/D");
    t_signal->Branch("mu2_pt", &mu2_pt, "mu2_pt/D");
    t_signal->Branch("mu1_eta", &mu1_eta, "mu1_eta/D");
    t_signal->Branch("mu2_eta", &mu2_eta, "mu2_eta/D");
    t_signal->Branch("mu_m", "TLorentzVector", &mu_m);
    t_signal->Branch("mu_p", "TLorentzVector", &mu_p);
    t_signal->Branch("gen_mu_m", "TLorentzVector", &gen_mu_m_vec);
    t_signal->Branch("gen_mu_p", "TLorentzVector", &gen_mu_p_vec);
    t_signal->Branch("jet1_pt", &jet1_pt, "jet1_pt/D");
    t_signal->Branch("jet1_eta", &jet1_eta, "jet1_eta/D");
    t_signal->Branch("jet1_CMVA", &jet1_cmva, "jet1_CMVA/D");
    t_signal->Branch("jet2_pt", &jet2_pt, "jet2_pt/D");
    t_signal->Branch("jet2_eta", &jet2_eta, "jet2_eta/D");
    t_signal->Branch("jet2_CMVA", &jet2_cmva, "jet2_CMVA/D");
    t_signal->Branch("met_pt", &met_pt, "met_Pt/F");
    t_signal->Branch("deltaC", &deltaC, "deltaC/D");
    t_signal->Branch("pu_SF", &pu_SF);
    t_signal->Branch("gen_weight", &gen_weight, "gen_weight/D");
    t_signal->Branch("bcdef_HLT_SF", &bcdef_HLT_SF);
    t_signal->Branch("bcdef_iso_SF", &bcdef_iso_SF);
    t_signal->Branch("bcdef_id_SF", &bcdef_id_SF);
    t_signal->Branch("bcdef_trk_SF", &bcdef_trk_SF);
    t_signal->Branch("gh_HLT_SF", &gh_HLT_SF);
    t_signal->Branch("gh_iso_SF", &gh_iso_SF);
    t_signal->Branch("gh_id_SF", &gh_id_SF);
    t_signal->Branch("gh_trk_SF", &gh_trk_SF);
    t_signal->Branch("jet1_b_weight", &jet1_b_weight);
    t_signal->Branch("jet2_b_weight", &jet2_b_weight);
    t_signal->Branch("mu_R_up", &mu_R_up);
    t_signal->Branch("mu_R_down", &mu_R_down);
    t_signal->Branch("mu_F_up", &mu_F_up);
    t_signal->Branch("mu_F_down", &mu_F_down);
    t_signal->Branch("mu_RF_up", &mu_RF_up);
    t_signal->Branch("mu_RF_down", &mu_RF_down);
    t_signal->Branch("pdf_up", &pdf_up);
    t_signal->Branch("pdf_down", &pdf_down);
    t_signal->Branch("pdf_weights", &pdf_weights, "pdf_weights[100]/F");
    t_signal->Branch("nJets", &nJets, "nJets/I");
    t_signal->Branch("jet1_flavour", &jet1_flavour, "jet1_flavour/I");
    t_signal->Branch("jet2_flavour", &jet2_flavour, "jet2_flavour/I");
    t_signal->Branch("is_tau_event", &is_tau_event);
    t_signal->Branch("pu_NtrueInt", &pu_NtrueInt);


    TTree *t_back = new TTree("T_back", "Tree for events with no asym (qq, gg)");
    //t_back->SetDirectory(0);
    t_back->Branch("m", &cm_m, "m/D");
    t_back->Branch("gen_m", &gen_m, "m/D");
    t_back->Branch("xF", &xF, "xF/D");
    t_back->Branch("cost", &cost_r, "cost/D");
    t_back->Branch("cost_st", &cost_st, "cost_st/D");
    t_back->Branch("mu1_pt", &mu1_pt, "mu1_pt/D");
    t_back->Branch("mu2_pt", &mu2_pt, "mu2_pt/D");
    t_back->Branch("mu1_eta", &mu1_pt, "mu1_eta/D");
    t_back->Branch("mu2_eta", &mu2_pt, "mu2_eta/D");
    t_back->Branch("mu_m", "TLorentzVector", &mu_m);
    t_back->Branch("mu_p", "TLorentzVector", &mu_p);
    t_back->Branch("gen_mu_m", "TLorentzVector", &gen_mu_m_vec);
    t_back->Branch("gen_mu_p", "TLorentzVector", &gen_mu_p_vec);
    t_back->Branch("jet1_pt", &jet1_pt, "jet1_pt/D");
    t_back->Branch("jet1_eta", &jet1_eta, "jet1_eta/D");
    t_back->Branch("jet1_CMVA", &jet1_cmva, "jet1_CMVA/D");
    t_back->Branch("jet2_pt", &jet2_pt, "jet2_pt/D");
    t_back->Branch("jet2_eta", &jet2_eta, "jet2_eta/D");
    t_back->Branch("jet2_CMVA", &jet2_cmva, "jet2_CMVA/D");
    t_back->Branch("met_pt", &met_pt, "met_Pt/F");
    t_back->Branch("deltaC", &deltaC, "deltaC/D");
    t_back->Branch("pu_SF", &pu_SF);
    t_back->Branch("gen_weight", &gen_weight, "gen_weight/D");
    t_back->Branch("bcdef_HLT_SF", &bcdef_HLT_SF);
    t_back->Branch("bcdef_iso_SF", &bcdef_iso_SF);
    t_back->Branch("bcdef_id_SF", &bcdef_id_SF);
    t_back->Branch("bcdef_trk_SF", &bcdef_trk_SF);
    t_back->Branch("gh_HLT_SF", &gh_HLT_SF);
    t_back->Branch("gh_iso_SF", &gh_iso_SF);
    t_back->Branch("gh_id_SF", &gh_id_SF);
    t_back->Branch("gh_trk_SF", &gh_trk_SF);
    t_back->Branch("jet1_b_weight", &jet1_b_weight);
    t_back->Branch("jet2_b_weight", &jet2_b_weight);
    t_back->Branch("nJets", &nJets, "nJets/I");
    t_back->Branch("jet1_flavour", &jet1_flavour, "jet1_flavour/I");
    t_back->Branch("jet2_flavour", &jet2_flavour, "jet2_flavour/I");
    t_back->Branch("is_tau_event", &is_tau_event);
    t_back->Branch("pu_NtrueInt", &pu_NtrueInt);



    TH1F *mu_p_dr = new TH1F("mu_p_dr", "DeltaR, reco Muon and gen Muon", 20,0,1);
    TH1F *mu_m_dr = new TH1F("mu_m_dr", "DeltaR, reco Muon and gen Muon", 20,0,1);

    unsigned int nEvents=0;
    unsigned int nSignal = 0;
    unsigned int nQQ=0;
    unsigned int nQQb=0;
    unsigned int nQGlu=0;
    unsigned int nGluGlu=0;
    unsigned int nTauTau=0;
    unsigned int nFailedID=0;
    unsigned int mismatch=0;
    unsigned int nNonIso = 0;

    Double_t normalization = 1.0;

    FILE *root_files = fopen(filename, "r");
    char lines[300];
    while(fgets(lines, 300, root_files)){
        if(lines[0] == '#' || is_empty_line(lines)) continue; //comment line
        else if(lines[0] == '!'){//sample header
            int sample_idx;
            float xsec;
            int nparams = sscanf(lines, "! idx = %i xsec = %f \n", &sample_idx, &xsec);
            if(nparams < 2 || sample_idx >= MAX_SAMPLES){
                printf("ERROR: Unable to parse sample header. Exiting");
                exit(EXIT_FAILURE);
            }
            normalization = norms[sample_idx];
            printf("Moving on to sample %i which has normalization %e \n", sample_idx, normalization);
        }
        else if(normalization > 0) {//root file



            char * end;
            //remove trailing whitespace
            end = lines + strlen(lines) - 1;
            while(end > lines && isspace((char)*end)) end--;
            // Write new null terminator
            *(end+1) = 0;

            printf("Opening file: %s \n", lines);
            TFile *f1=  TFile::Open(lines);

            f1->cd("EventCounter");
            TDirectory *subdir = gDirectory;
            TH1D *mc_pileup = (TH1D *)subdir->Get("pileup");
            mc_pileup->Scale(1./mc_pileup->Integral());
            pu_SFs.pileup_ratio->Divide(pu_SFs.data_pileup, mc_pileup);

            f1->cd("B2GTTreeMaker");
            TTree *t1 = (TTree *)gDirectory->Get("B2GTree");


            UInt_t mu_size, gen_size, jet_size, met_size;
            Int_t gen_id[GEN_SIZE], gen_status[GEN_SIZE];
            Int_t  gen_Mom0ID[GEN_SIZE], gen_Mom0Status[GEN_SIZE], gen_Mom1ID[GEN_SIZE], gen_Mom1Status[GEN_SIZE];
            Int_t  gen_Dau0ID[GEN_SIZE], gen_Dau0Status[GEN_SIZE], gen_Dau1ID[GEN_SIZE], gen_Dau1Status[GEN_SIZE];
            Float_t gen_Pt[GEN_SIZE], gen_Eta[GEN_SIZE], gen_Phi[GEN_SIZE], gen_E[GEN_SIZE];

            Float_t mu_Pt[MU_SIZE], mu_Eta[MU_SIZE], mu_Phi[MU_SIZE], mu_E[MU_SIZE], 
                    mu_Charge[MU_SIZE], mu_IsTightMuon[MU_SIZE];
            Float_t mu_SumChargedHadronPt[MU_SIZE], mu_SumNeutralHadronPt[MU_SIZE], mu_SumPUPt[MU_SIZE], mu_SumPhotonPt[MU_SIZE],
                    mu_NumberTrackerLayers[MU_SIZE];


            Float_t jet_Pt[JET_SIZE], jet_Eta[JET_SIZE], jet_Phi[JET_SIZE], jet_E[JET_SIZE],
                    jet_CSV[JET_SIZE], jet_CMVA[JET_SIZE], jet_partonflavour[JET_SIZE];


            Float_t evt_Gen_Weight;

            Int_t HLT_IsoMu, HLT_IsoTkMu;
            t1->SetBranchAddress("mu_size", &mu_size); //number of muons in the event
            t1->SetBranchAddress("mu_Pt", &mu_Pt);
            t1->SetBranchAddress("mu_Eta", &mu_Eta);
            t1->SetBranchAddress("mu_Phi", &mu_Phi);
            t1->SetBranchAddress("mu_E", &mu_E);
            t1->SetBranchAddress("mu_Charge", &mu_Charge);

            t1->SetBranchAddress("mu_IsTightMuon", &mu_IsTightMuon);
            t1->SetBranchAddress("mu_SumChargedHadronPt", &mu_SumChargedHadronPt);
            t1->SetBranchAddress("mu_SumNeutralHadronPt", &mu_SumNeutralHadronPt);
            t1->SetBranchAddress("mu_SumPUPt", &mu_SumPUPt);
            t1->SetBranchAddress("mu_SumPhotonPt", &mu_SumPhotonPt);
            //t1->SetBranchAddress("mu_NumberTrackerLayers", &mu_NumberTrackerLayers);

            t1->SetBranchAddress("jetAK4CHS_size", &jet_size);
            t1->SetBranchAddress("jetAK4CHS_Pt", &jet_Pt);
            t1->SetBranchAddress("jetAK4CHS_Eta", &jet_Eta);
            t1->SetBranchAddress("jetAK4CHS_Phi", &jet_Phi);
            t1->SetBranchAddress("jetAK4CHS_E", &jet_E);
            t1->SetBranchAddress("jetAK4CHS_CSVv2", &jet_CSV);
            t1->SetBranchAddress("jetAK4CHS_CMVAv2", &jet_CMVA);
            t1->SetBranchAddress("jetAK4CHS_PartonFlavour", &jet_partonflavour);

            t1->SetBranchAddress("HLT_IsoMu24", &HLT_IsoMu);
            t1->SetBranchAddress("HLT_IsoTkMu24", &HLT_IsoTkMu);



            t1->SetBranchAddress("gen_size", &gen_size); //number of muons in the event
            t1->SetBranchAddress("gen_Pt", &gen_Pt);
            t1->SetBranchAddress("gen_Eta", &gen_Eta);
            t1->SetBranchAddress("gen_Phi", &gen_Phi);
            t1->SetBranchAddress("gen_E", &gen_E);
            t1->SetBranchAddress("gen_ID", &gen_id);
            t1->SetBranchAddress("gen_Status", &gen_status);
            t1->SetBranchAddress("evt_Gen_Weight", &evt_Gen_Weight);
            t1->SetBranchAddress("pu_NtrueInt",&pu_NtrueInt);

            t1->SetBranchAddress("scale_Weights", &scale_Weights);
            t1->SetBranchAddress("pdf_Weights", &pdf_weights);


            t1->SetBranchAddress("gen_Mom0ID", &gen_Mom0ID);
            t1->SetBranchAddress("gen_Mom1ID", &gen_Mom1ID);
            t1->SetBranchAddress("gen_Mom0Status", &gen_Mom0Status);
            t1->SetBranchAddress("gen_Mom1Status", &gen_Mom1Status);

            t1->SetBranchAddress("gen_Dau0ID", &gen_Dau0ID);
            t1->SetBranchAddress("gen_Dau1ID", &gen_Dau1ID);
            t1->SetBranchAddress("gen_Dau0Status", &gen_Dau0Status);
            t1->SetBranchAddress("gen_Dau1Status", &gen_Dau1Status);

            t1->SetBranchAddress("met_size", &met_size);
            t1->SetBranchAddress("met_Pt", &met_pt);

            Long64_t nEntries =  t1->GetEntries();

            char out_buff[10000];
            bool print_out = true;

            for (int i=0; i<nEntries; i++) {
                t1->GetEntry(i);
                if(mu_size > MU_SIZE || gen_size >GEN_SIZE) printf("WARNING: MU_SIZE OR GEN_SIZE TOO LARGE \n");
                if(met_size != 1) printf("WARNING: Met size not equal to 1\n");
                bool good_trigger = HLT_IsoMu || HLT_IsoTkMu;
                if(good_trigger &&
                        mu_size >= 2 && ((abs(mu_Charge[0] - mu_Charge[1])) > 0.01) &&
                        mu_IsTightMuon[0] && mu_IsTightMuon[1] &&
                        mu_Pt[0] > 26. &&  mu_Pt[1] > 10. &&
                        abs(mu_Eta[0]) < 2.4 && abs(mu_Eta[1]) < 2.4){ 

                    //See https://twiki.cern.ch/twiki/bin/viewauth/CMS/SWGuideMuonIdRun2 for iso cuts
                    float iso_0 = (mu_SumChargedHadronPt[0] + max(0., mu_SumNeutralHadronPt[0] + mu_SumPhotonPt[0] - 0.5 * mu_SumPUPt[0]))/mu_Pt[0];
                    float iso_1 = (mu_SumChargedHadronPt[1] + max(0., mu_SumNeutralHadronPt[1] + mu_SumPhotonPt[1] - 0.5 * mu_SumPUPt[1]))/mu_Pt[1];
                    float tight_iso = 0.15;
                    float loose_iso = 0.25;
                    nNonIso++;
                    //only want events with 2 oppositely charged muons
                    double mu_p_mcSF, mu_m_mcSF;
                    if(mu_Charge[0] >0){
                        mu_p.SetPtEtaPhiE(mu_Pt[0], mu_Eta[0], mu_Phi[0], mu_E[0]);
                        mu_m.SetPtEtaPhiE(mu_Pt[1], mu_Eta[1], mu_Phi[1], mu_E[1]);
                        //mu_p_mcSF = rc.kScaleAndSmearMC(1, mu_Pt[0], mu_Eta[0], mu_Phi[0], 5, rand->Rndm(), rand->Rndm(), 0, 0);
                        //mu_m_mcSF = rc.kScaleAndSmearMC(-1, mu_Pt[1], mu_Eta[1], mu_Phi[1], 5, rand->Rndm(), rand->Rndm(), 0, 0);
                    }
                    else{
                        mu_m.SetPtEtaPhiE(mu_Pt[0], mu_Eta[0], mu_Phi[0], mu_E[0]);
                        mu_p.SetPtEtaPhiE(mu_Pt[1], mu_Eta[1], mu_Phi[1], mu_E[1]);
                        //mu_m_mcSF = rc.kScaleAndSmearMC(-1, mu_Pt[0], mu_Eta[0], mu_Phi[0], 5, rand->Rndm(), rand->Rndm(), 0, 0);
                        //mu_p_mcSF = rc.kScaleAndSmearMC(1, mu_Pt[1], mu_Eta[1], mu_Phi[1], 5, rand->Rndm(), rand->Rndm(), 0, 0);
                    }
                    //printf("charges are %1f %1f \n", mu_Charge[0], mu_Charge[1]);
                    //printf("(pt, eta, phi,E) = %4.2f %4.2f %4.2f %4.2f \n", mu_Pt[0], mu_Eta[0], mu_Phi[0], mu_E[0]);
                    //printf("(pt, eta, phi,E) = %4.2f %4.2f %4.2f %4.2f \n", mu_Pt[1], mu_Eta[1], mu_Phi[1], mu_E[1]);
                    //mu_p.Print();
                    //mu_m.Print();
                    //printf ("Momentum SFs are %.3f %.3f for Pts %.0f %.0f \n", mu_p_mcSF, mu_m_mcSF, mu_p.Pt(), mu_m.Pt());
                   


                    cm = mu_p + mu_m;
                    cm_m = cm.M();
                    //met and cmva cuts to reduce ttbar background
                    if (iso_0 < tight_iso && iso_1 < tight_iso && cm_m >= 50.){
                        if(PRINT) sprintf(out_buff + strlen(out_buff),"Event %i \n", i);

                        //GEN LEVEL
                        //
                        //Pythia 8 status numbering convention
                        //see:
                        //http://home.thep.lu.se/~torbjorn/pythia81html/EventRecord.html
                        int FINAL_STATE = 1;
                        int EVENT_PARTICLE = 11;
                        int BEAM_PARTICLE=4;
                        int INCIDENT_PARTICLE = 21;
                        int INTERMED_PARTICLE = 22;
                        int OUTGOING = 23;
                        //Particle ID's
                        int ELECTRON = 11; 
                        int MUON = 13;
                        int PHOTON = 22;
                        int Z=23;
                        int GLUON = 21;
                        int TAU = 15;
                        int PROTON = 2212;

                        int inc_1 =-1;
                        int inc_2 =-1;
                        int gen_mu_p=-1;
                        int gen_mu_m=-1;
                        int gen_tau_p=-1;
                        int gen_tau_m=-1;
                        int gen_e_p=-1;
                        int gen_e_m=-1;
                        int intermed=-1;

                        float quark_dir_eta;

                        bool signal_event = true;//whether it is an event with an asym or not

                        is_tau_event = false;

                        for(int k=0; k<gen_size; k++){
                            if(gen_status[k] == INCIDENT_PARTICLE && 
                                    (abs(gen_id[k]) <=6  || gen_id[k] == GLUON) && 
                                    (abs(gen_Dau0ID[k]) == MUON || gen_Dau0ID[k] == Z || 
                                     gen_Dau0ID[k] == PHOTON || abs(gen_Dau0ID[k]) == TAU)
                              ){
                                //record index of 2 initial state particles
                                if(inc_1 == -1) inc_1 = k;
                                else if(inc_2 == -1) inc_2 = k;
                                else{
                                    print_out = true;
                                    printf("WARNING: More than 2 incident particles in event\n\n");
                                }

                            }
                            //record 2 scattered muons
                            if(abs(gen_id[k]) == MUON && 
                                    (gen_Mom0ID[k] == Z || gen_Mom0ID[k] == PHOTON || abs(gen_Mom0ID[k]) == TAU 
                                      || abs(gen_Mom0ID[k]) == ELECTRON || (gen_status[k] == OUTGOING && gen_Mom0ID[k] != PROTON))) {
                                if(gen_id[k] == MUON){
                                    if(gen_mu_m == -1) gen_mu_m = k;
                                    else{
                                        if(abs(gen_Mom0ID[k]) != TAU) printf("WARNING: More than one mu_m\n\n");
                                        if(PRINT) sprintf(out_buff + strlen(out_buff), "Extra mu_m detected\n");
                                        print_out = true;
                                    }
                                }
                                if(gen_id[k] == -MUON){
                                    if(gen_mu_p == -1) gen_mu_p = k;
                                    else{
                                        if(abs(gen_Mom0ID[k]) != TAU) printf("WARNING: More than one mu_p\n\n");
                                        if(PRINT) sprintf(out_buff + strlen(out_buff), "Extra mu_p detected\n");
                                        print_out = true;
                                    }
                                }
                            }
                            //record tau's
                            if(abs(gen_id[k]) == TAU && 
                                    ( (gen_Mom0ID[k] == Z || gen_Mom0ID[k] == PHOTON || abs(gen_Mom0ID[k]) == ELECTRON ||
                                      gen_status[k] == OUTGOING)) ){
                                if(gen_id[k] == TAU){
                                    if(gen_tau_m == -1) gen_tau_m = k;
                                    else{
                                        printf("WARNING: More than one tau_m\n\n");
                                        if(PRINT) sprintf(out_buff + strlen(out_buff), "Extra mu_m detected\n");
                                        //print_out = true;
                                    }
                                }
                                if(gen_id[k] == -TAU){
                                    if(gen_tau_p == -1) gen_tau_p = k;
                                    else{
                                        printf("WARNING: More than one tau_p\n\n");
                                        if(PRINT) sprintf(out_buff + strlen(out_buff), "Extra mu_p detected\n");
                                        //print_out = true;
                                    }
                                }
                            }
                            if(PRINT){
                                if( (abs(gen_id[k]) <=6 || gen_id[k] == GLUON)){
                            //    if( (abs(gen_id[k]) <=6 || gen_id[k] == GLUON) && 
                            //            (gen_Dau0ID[k] == PHOTON &&
                            //             gen_Dau0ID[k] == Z || abs(gen_Dau0ID[k]) == MUON || abs(gen_Dau0ID[k]) == TAU)){
                                    if(PRINT) sprintf(out_buff + strlen(out_buff),"Parton (ID = %i stat = %i): \n"
                                            "    Mom1 ID: %i Mom1 Stat: %i Mom2 ID: %i Mom2 Stat %i \n"
                                            "    Dau1 ID: %i Dau1 Stat: %i Dau2 ID: %i Dau2 Stat %i \n",
                                            gen_id[k], gen_status[k], 
                                            gen_Mom0ID[k], gen_Mom0Status[k], gen_Mom1ID[k], gen_Mom1Status[k],
                                            gen_Dau0ID[k], gen_Dau0Status[k], gen_Dau1ID[k], gen_Dau1Status[k]);
                                }


                                if(gen_id[k] == -13){
                                    if(PRINT) sprintf(out_buff + strlen(out_buff),"mu_p(stat = %i): \n"
                                            "   Mom1 ID: %i Mom1 Stat: %i Mom2 ID: %i Mom2 Stat %i \n"
                                            "   Dau1 ID: %i Dau1 Stat: %i Dau2 ID: %i Dau2 Stat %i \n",
                                            gen_status[k], 
                                            gen_Mom0ID[k], gen_Mom0Status[k], gen_Mom1ID[k], gen_Mom1Status[k],
                                            gen_Dau0ID[k], gen_Dau0Status[k], gen_Dau1ID[k], gen_Dau1Status[k]);
                                    if(PRINT) sprintf(out_buff + strlen(out_buff),"(pt, eta, phi,E) = %4.2f %4.2f %4.2f %4.2f \n", 
                                            gen_Pt[k], gen_Eta[k], gen_Phi[k], gen_E[k]);
                                }
                                if(gen_id[k] == 13){

                                    if(PRINT) sprintf(out_buff + strlen(out_buff),"mu_m (stat = %i): \n"
                                            "   Mom1 ID: %i Mom1 Stat: %i Mom2 ID: %i Mom2 Stat %i \n"
                                            "   Dau1 ID: %i Dau1 Stat: %i Dau2 ID: %i Dau2 Stat %i \n",
                                            gen_status[k], 
                                            gen_Mom0ID[k], gen_Mom0Status[k], gen_Mom1ID[k], gen_Mom1Status[k],
                                            gen_Dau0ID[k], gen_Dau0Status[k], gen_Dau1ID[k], gen_Dau1Status[k]);
                                    if(PRINT) sprintf(out_buff + strlen(out_buff),"(pt, eta, phi,E) = %4.2f %4.2f %4.2f %4.2f \n", 
                                            gen_Pt[k], gen_Eta[k], gen_Phi[k], gen_E[k]);

                                }
                                if(gen_id[k] == 23){
                                    if(PRINT) sprintf(out_buff + strlen(out_buff),"Z (ID = %i, status = %i): \n"
                                            "   Mom1 ID: %i Mom1 Stat: %i Mom2 ID: %i Mom2 Stat %i \n"
                                            "   Dau1 ID: %i Dau1 Stat: %i Dau2 ID: %i Dau2 Stat %i \n",
                                            gen_id[k], gen_status[k],
                                            gen_Mom0ID[k], gen_Mom0Status[k], gen_Mom1ID[k], gen_Mom1Status[k],
                                            gen_Dau0ID[k], gen_Dau0Status[k], gen_Dau1ID[k], gen_Dau1Status[k]);
                                }
                            }
                        }


                        if(gen_mu_p != -1 && gen_mu_m != -1) {
                            if(PRINT) sprintf(out_buff + strlen(out_buff),"mu_p: \n"
                                    "   Mom1 ID: %i Mom1 Stat: %i Mom2 ID: %i Mom2 Stat %i \n"
                                    "   Dau1 ID: %i Dau1 Stat: %i Dau2 ID: %i Dau2 Stat %i \n",
                                    gen_Mom0ID[gen_mu_p], gen_Mom0Status[gen_mu_p], gen_Mom1ID[gen_mu_p], 
                                    gen_Mom1Status[gen_mu_p],
                                    gen_Dau0ID[gen_mu_p], gen_Dau0Status[gen_mu_p], gen_Dau1ID[gen_mu_p], 
                                    gen_Dau1Status[gen_mu_p]);

                            if(PRINT) sprintf(out_buff + strlen(out_buff),"mu_m: \n"
                                    "   Mom1 ID: %i Mom1 Stat: %i Mom2 ID: %i Mom2 Stat %i \n"
                                    "   Dau1 ID: %i Dau1 Stat: %i Dau2 ID: %i Dau2 Stat %i \n",
                                    gen_Mom0ID[gen_mu_m], gen_Mom0Status[gen_mu_m], gen_Mom1ID[gen_mu_m], 
                                    gen_Mom1Status[gen_mu_m],
                                    gen_Dau0ID[gen_mu_m], gen_Dau0Status[gen_mu_m], gen_Dau1ID[gen_mu_m], 
                                    gen_Dau1Status[gen_mu_m]);

                        }
                        else {
                            printf("WARNING: Unable to identify MuMu pair in event %i, skipping \n", i);
                            nFailedID ++;
                            print_out = true;
                            continue;
                        }
                        if((inc_1 == -1) || (inc_2 == -1)){
                            printf("WARNING: Unable to identify initial state particles in event %i, skipping \n", i);
                            nFailedID ++;
                            print_out = true;
                            continue;
                        }

                        else{ 
                            //printf("%i %i \n", inc_1, inc_2);
                            int inc_id1 = gen_id[inc_1];
                            int inc_id2 = gen_id[inc_2];
                            if(abs(gen_Dau0ID[inc_1]) == TAU){
                                if(gen_tau_p == -1 || gen_tau_m == -1){
                                    printf("Didn't record tau's :( \n");
                                    nFailedID ++;
                                    continue;
                                }
                                nTauTau++;
                                is_tau_event = true;
                            }
                            if((abs(inc_id1) <= 6 && abs(inc_id2) <= 6) && (inc_id1 * inc_id2 < 0)){ //a quark and anti quark
                                //qq-bar
                                signal_event = true;
                                nQQb++;
                                if(inc_id1>0) quark_dir_eta = gen_Eta[inc_1];
                                else if(inc_id2>0) quark_dir_eta = gen_Eta[inc_2];
                            }
                            else if(((abs(inc_id1) <= 6) && (inc_id2 == 21)) ||
                                    ((abs(inc_id2) <= 6) && (inc_id1 == 21))){ //qglu
                                signal_event = true;
                                int q_dir;
                                if(inc_id1 == 21){
                                    if(inc_id2 <0) quark_dir_eta = gen_Eta[inc_1];//qbar-glu, want glu dir
                                    else quark_dir_eta= gen_Eta[inc_2];//q-glu ,want q dir
                                }
                                else if(inc_id2 == 21) {
                                    if(inc_id1 <0) quark_dir_eta = gen_Eta[inc_2];//qbar-glu, want glu dir
                                    else quark_dir_eta= gen_Eta[inc_1];//q-glu ,want q dir
                                }
                                nQGlu++;
                            }
                            else if((abs(inc_id1) <= 6) && (abs(inc_id2) <= 6) && (inc_id1 * inc_id2 >0)){ //2 quarks
                                if(PRINT) sprintf(out_buff + strlen(out_buff),"QQ Event \n");
                                signal_event = false;
                                nQQ++;
                            }
                            else if((inc_id1 == 21) && (inc_id2 == 21)){ //gluglu
                                signal_event = false;
                                nGluGlu++;
                                if(PRINT) sprintf(out_buff + strlen(out_buff), "Glu Glu event \n");
                            }
                            else {
                                printf("WARNING: not qqbar, qq, qg, or gg event");
                                printf("First particle was %i second particle was %i \n \n ", inc_id1, inc_id2);
                                nFailedID ++;
                            }
                        }
                        gen_mu_p_vec.SetPtEtaPhiE(gen_Pt[gen_mu_p], gen_Eta[gen_mu_p], gen_Phi[gen_mu_p], gen_E[gen_mu_p]);
                        gen_mu_m_vec.SetPtEtaPhiE(gen_Pt[gen_mu_m], gen_Eta[gen_mu_m], gen_Phi[gen_mu_m], gen_E[gen_mu_m]);
                        TLorentzVector gen_cm = gen_mu_p_vec + gen_mu_m_vec;
                        gen_m = gen_cm.M();
                        //RECO LEVEL
                        xF = abs(2.*cm.Pz()/13000.); 

                        // compute Colins soper angle with formula
                        double mu_p_pls = (mu_p.E()+mu_p.Pz())/root2;
                        double mu_p_min = (mu_p.E()-mu_p.Pz())/root2;
                        double mu_m_pls = (mu_m.E()+mu_m.Pz())/root2;
                        double mu_m_min = (mu_m.E()-mu_m.Pz())/root2;
                        double qt2 = cm.Px()*cm.Px()+cm.Py()*cm.Py();
                        double cm_m2 = cm.M2();
                        //cost_p = cos(theta)_r (reconstructed collins soper angle, sign
                        //may be 'wrong' if lepton pair direction is not the same as inital
                        //quark direction)
                        double cost = 2*(mu_m_pls*mu_p_min - mu_m_min*mu_p_pls)/sqrt(cm_m2*(cm_m2 + qt2));

                        double cost_tau; 
                        //cos(theta) from generator level taus
                        if(is_tau_event){
                            TLorentzVector tau_p, tau_m, tau_cm;
                            tau_p.SetPtEtaPhiE(gen_Pt[gen_tau_p], gen_Eta[gen_tau_p], gen_Phi[gen_tau_p], gen_E[gen_tau_p]);
                            tau_m.SetPtEtaPhiE(gen_Pt[gen_tau_m], gen_Eta[gen_tau_m], gen_Phi[gen_tau_m], gen_E[gen_tau_m]);
                            tau_cm = tau_p + tau_m;
                            double tau_p_pls = (tau_p.E()+tau_p.Pz())/root2;
                            double tau_p_min = (tau_p.E()-tau_p.Pz())/root2;
                            double tau_m_pls = (tau_m.E()+tau_m.Pz())/root2;
                            double tau_m_min = (tau_m.E()-tau_m.Pz())/root2;
                            double tau_qt2 = tau_cm.Px()*tau_cm.Px()+tau_cm.Py()*tau_cm.Py();
                            double tau_cm_m2 = tau_cm.M2();
                            cost_tau = 2*(tau_m_pls*tau_p_min - tau_m_min*tau_p_pls)/sqrt(tau_cm_m2*(tau_cm_m2 + tau_qt2));
                        }


                        /*
                        TLorentzVector p1(0., 0., Pbeam, Ebeam);
                        TLorentzVector p2(0., 0., -Pbeam, Ebeam);

                        if(cm.Pz() < 0. ){
                            TLorentzVector p = p1;
                            p1 = p2;
                            p2 = p;
                        }

                        TVector3 beta = -cm.BoostVector();
                        mu_m.Boost(beta);
                        mu_p.Boost(beta);
                        p1.Boost(beta);
                        p2.Boost(beta);

                        // Now calculate the direction of the new z azis

                        TVector3 p1u = p1.Vect();
                        p1u.SetMag(1.0);
                        TVector3 p2u = p2.Vect();
                        p2u.SetMag(1.0);
                        TVector3 pzu = p1u - p2u;
                        pzu.SetMag(1.0);
                        mu_m.RotateUz(pzu); 
                        double cost_r_b = mu_m.CosTheta();
                        deltaC = std::abs(cost_r_b) - std::abs(cost);
                        */
                        //printf("cost_r, cost_r_b, cost_r_b2: %0.2f %0.2f %0.2f \n", cost_r, cost_r_b, cost_r_b2);

                        if(PRINT){
                            sprintf(out_buff + strlen(out_buff),  "1st Particle %i:(pt, eta, phi,E) = %4.2f %4.2f %4.2f %4.2f \n", 
                                    gen_id[inc_1], gen_Pt[inc_1], gen_Eta[inc_1], gen_Phi[inc_1], gen_E[inc_1]);
                            sprintf(out_buff + strlen(out_buff),"2nd Particle %i:(pt, eta, phi,E) = %4.2f %4.2f %4.2f %4.2f \n", 
                                    gen_id[inc_2], gen_Pt[inc_2], gen_Eta[inc_2], gen_Phi[inc_2], gen_E[inc_2]);
                        }


                        //reconstruction sign flip
                        if(cm.Pz() < 0.) cost_r = -cost;
                        else cost_r = cost;


                        gen_weight = evt_Gen_Weight * normalization;
                        mu1_pt = mu_Pt[0];
                        mu2_pt = mu_Pt[1];
                        mu1_eta = mu_Eta[0];
                        mu2_eta = mu_Eta[1];

                        //pick out 2 highest pt jets with eta < 2.4
                        nJets =0;
                        for(int j=0; j < jet_size; j++){
                            if(jet_Pt[j] > 20. && std::abs(jet_Eta[j]) < 2.4){
                                if(nJets == 1){
                                    jet2_pt = jet_Pt[j];
                                    jet2_eta = jet_Eta[j];
                                    jet2_cmva = jet_CMVA[j];
                                    jet2_flavour = jet_partonflavour[j];
                                    jet2_b_weight = get_btag_weight(jet_Pt[j], jet_Eta[j],jet_partonflavour[j],btag_effs, b_reader);
                                    nJets =2;
                                    break;
                                }
                                else if(nJets ==0){
                                    jet1_pt = jet_Pt[j];
                                    jet1_eta = jet_Eta[j];
                                    jet1_cmva = jet_CMVA[j];
                                    jet1_flavour = jet_partonflavour[j];
                                    jet1_b_weight = get_btag_weight(jet_Pt[j], jet_Eta[j],jet_partonflavour[j],btag_effs, b_reader);
                                    nJets = 1;
                                }
                            }
                        }


                        //get muon cut SFs

                        bcdef_HLT_SF = get_HLT_SF(mu1_pt, mu1_eta, mu2_pt, mu2_eta, runs_bcdef.HLT_SF, runs_bcdef.HLT_MC_EFF);
                        gh_HLT_SF = get_HLT_SF(mu1_pt, mu1_eta, mu2_pt, mu2_eta, runs_gh.HLT_SF, runs_gh.HLT_MC_EFF);

                        bcdef_iso_SF = get_SF(mu1_pt, mu1_eta, runs_bcdef.ISO_SF) * get_SF(mu2_pt, mu2_eta, runs_bcdef.ISO_SF);
                        bcdef_id_SF = get_SF(mu1_pt, mu1_eta, runs_bcdef.ID_SF) * get_SF(mu2_pt, mu2_eta, runs_bcdef.ID_SF);

                        gh_iso_SF = get_SF(mu1_pt, mu1_eta, runs_gh.ISO_SF) * get_SF(mu2_pt, mu2_eta, runs_gh.ISO_SF);
                        gh_id_SF = get_SF(mu1_pt, mu1_eta, runs_gh.ID_SF) * get_SF(mu2_pt, mu2_eta, runs_gh.ID_SF);


                        pu_SF = get_pileup_SF(pu_NtrueInt, pu_SFs.pileup_ratio);

                        bcdef_trk_SF = get_Mu_trk_SF(abs(mu1_eta), runs_bcdef.TRK_SF) * get_Mu_trk_SF(abs(mu2_eta), runs_bcdef.TRK_SF);
                        gh_trk_SF = get_Mu_trk_SF(abs(mu1_eta), runs_gh.TRK_SF) * get_Mu_trk_SF(abs(mu2_eta), runs_gh.TRK_SF);


                        mu_R_up = scale_Weights[2];
                        mu_R_down = scale_Weights[4];
                        mu_F_up = scale_Weights[0];
                        mu_F_down = scale_Weights[1];
                        mu_RF_up = scale_Weights[3];
                        mu_RF_down = scale_Weights[5];

                        Float_t pdf_avg, pdf_std_dev;

                        get_pdf_avg_std_dev(pdf_weights, &pdf_avg, &pdf_std_dev);

                        pdf_up = pdf_avg + pdf_std_dev;
                        pdf_down = pdf_avg - pdf_std_dev;

                        if(signal_event){
                            //cost_st = cos(theta)_* correct angle obtained from 'cheating' and
                            //looking at initial quark direction
                            if(!is_tau_event){
                                if(quark_dir_eta < 0){
                                    if(PRINT) sprintf(out_buff + strlen(out_buff), "Sign flip\n");
                                    cost_st = -cost;
                                }
                                else cost_st = cost;
                            }
                            else{
                                if(quark_dir_eta < 0){
                                    if(PRINT) sprintf(out_buff + strlen(out_buff), "Sign flip\n");
                                    cost_st = -cost_tau;
                                }
                                else cost_st = cost_tau;
                            }


                            //anti-symmetric template has computed weight and is 
                            //anti-symmetrized by flipping the sign of the weight and the bin
                            //in c_r
                            nSignal++;
                            t_signal->Fill();
                        }
                        else{
                            cost_st = cost_r;
                            t_back->Fill();
                        }

                        nEvents++;
                        if(PRINT && print_out){
                            sprintf(out_buff + strlen(out_buff), "\n\n");
                            fputs(out_buff, stdout);
                            print_out = false;
                        }
                        if(PRINT) memset(out_buff, 0, 10000);

                    }
                } 
            }

            f1->Close();
            printf("moving on to next file, currently %i events %i Taus %i fails \n\n", nEvents, nTauTau, nFailedID);
        }
    }
    fclose(root_files);
    printf("There were %i qqbar, %i qGlu (%i of them tautau) in %i kept events in %i files."
            "There were also %i background events (%i qq and %i gg)"
            "There were %i Failed ID's \n" , 
            nQQb, nQGlu, nTauTau, nSignal, nFiles, nQQ + nGluGlu, nQQ, nGluGlu, nFailedID);
    //printf("Ran on MC data and produced templates with %i events\n", nEvents);

    fout->cd();




    t_signal->Write();
    t_back->Write();


    printf("Writing output to file at %s \n", fout_name.Data());

    fout->Close();


    return;
}
