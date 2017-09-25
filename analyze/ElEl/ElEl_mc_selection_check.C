
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <iostream>
#include <cstring>
#include <algorithm>
#include "TFile.h"
#include "../ScaleFactors.C"

#define GEN_SIZE 4000
#define EL_SIZE 100
#define JET_SIZE 20
#define MAX_SAMPLES 20

const double root2 = sqrt(2);
double Ebeam = 6500.;
double Pbeam = sqrt(Ebeam*Ebeam - 0.938*0.938);

char *filename("DY_files_aug29.txt");
const TString fout_name("output_files/ElEl_selection_check_sep13.root");
const double alpha = 0.05;
const bool PRINT=false;





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




void ElEl_mc_selection_check()
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
    el_SFs el_SF;
    //separate SFs for runs BCDEF and GH
    setup_SFs(&runs_bcdef, &runs_gh, &b_reader, &btag_effs, &pu_SFs);
    setup_el_SF(&el_SF);
    printf("Retrieved Scale Factors \n\n");


    TFile *fout = TFile::Open(fout_name, "RECREATE");
    TTree *t_signal= new TTree("T_data", "Tree with asym events (qq bar, qg)");
    Double_t cm_m, xF, cost_r, cost_st, el1_pt, el2_pt, el1_eta, el2_eta, jet1_pt, jet2_pt, jet1_eta, jet2_eta, deltaC, 
             gen_weight, reweight, jet1_csv, jet1_cmva, jet2_csv, jet2_cmva;
    Double_t el_id_SF, jet1_b_weight, jet2_b_weight, pu_SF;
    Int_t nJets, jet1_flavour, jet2_flavour;
    Bool_t is_tau_event;
    Float_t met_pt;
    TLorentzVector el_p, el_m, cm, q1, q2;

    t_signal->Branch("m", &cm_m, "m/D");
    t_signal->Branch("xF", &xF, "xF/D");
    t_signal->Branch("cost", &cost_r, "cost/D");
    t_signal->Branch("cost_st", &cost_st, "cost_st/D");
    t_signal->Branch("el1_pt", &el1_pt, "el1_pt/D");
    t_signal->Branch("el2_pt", &el2_pt, "el2_pt/D");
    t_signal->Branch("el1_eta", &el1_eta, "el1_eta/D");
    t_signal->Branch("el2_eta", &el2_eta, "el2_eta/D");
    t_signal->Branch("el_m", "TLorentzVector", &el_m);
    t_signal->Branch("el_p", "TLorentzVector", &el_p);
    t_signal->Branch("jet1_pt", &jet1_pt, "jet1_pt/D");
    t_signal->Branch("jet1_eta", &jet1_eta, "jet1_eta/D");
    t_signal->Branch("jet1_CMVA", &jet1_cmva, "jet1_CMVA/D");
    t_signal->Branch("jet2_pt", &jet2_pt, "jet2_pt/D");
    t_signal->Branch("jet2_eta", &jet2_eta, "jet2_eta/D");
    t_signal->Branch("jet2_CMVA", &jet2_cmva, "jet2_CMVA/D");
    t_signal->Branch("met_pt", &met_pt, "met_Pt/F");
    t_signal->Branch("deltaC", &deltaC, "deltaC/D");
    t_signal->Branch("gen_weight", &gen_weight, "gen_weight/D");
    t_signal->Branch("pu_SF", &pu_SF);
    t_signal->Branch("reweight", &reweight, "reweight/D");
    t_signal->Branch("jet1_b_weight", &jet1_b_weight);
    t_signal->Branch("jet2_b_weight", &jet2_b_weight);
    t_signal->Branch("el_id_SF", &el_id_SF);
    t_signal->Branch("nJets", &nJets, "nJets/I");
    t_signal->Branch("jet1_flavour", &jet1_flavour, "jet1_flavour/I");
    t_signal->Branch("jet2_flavour", &jet2_flavour, "jet2_flavour/I");
    t_signal->Branch("is_tau_event", &is_tau_event);


    TTree *t_back = new TTree("T_back", "Tree for events with no asym (qq, gg)");
    t_back->Branch("m", &cm_m, "m/D");
    t_back->Branch("xF", &xF, "xF/D");
    t_back->Branch("cost", &cost_r, "cost/D");
    t_back->Branch("cost_st", &cost_st, "cost_st/D");
    t_back->Branch("el1_pt", &el1_pt, "el1_pt/D");
    t_back->Branch("el2_pt", &el2_pt, "el2_pt/D");
    t_back->Branch("el1_eta", &el1_pt, "el1_eta/D");
    t_back->Branch("el2_eta", &el2_pt, "el2_eta/D");
    t_back->Branch("el_m", "TLorentzVector", &el_m);
    t_back->Branch("el_p", "TLorentzVector", &el_p);
    t_back->Branch("jet1_pt", &jet1_pt, "jet1_pt/D");
    t_back->Branch("jet1_eta", &jet1_eta, "jet1_eta/D");
    t_back->Branch("jet1_CMVA", &jet1_cmva, "jet1_CMVA/D");
    t_back->Branch("jet2_pt", &jet2_pt, "jet2_pt/D");
    t_back->Branch("jet2_eta", &jet2_eta, "jet2_eta/D");
    t_back->Branch("jet2_CMVA", &jet2_cmva, "jet2_CMVA/D");
    t_back->Branch("met_pt", &met_pt, "met_Pt/F");
    t_back->Branch("deltaC", &deltaC, "deltaC/D");
    t_back->Branch("gen_weight", &gen_weight, "gen_weight/D");
    t_back->Branch("pu_SF", &pu_SF);
    t_back->Branch("reweight", &reweight, "reweight/D");
    t_back->Branch("el_id_SF", &el_id_SF);
    t_back->Branch("jet1_b_weight", &jet1_b_weight);
    t_back->Branch("jet2_b_weight", &jet2_b_weight);
    t_back->Branch("nJets", &nJets, "nJets/I");
    t_back->Branch("jet1_flavour", &jet1_flavour, "jet1_flavour/I");
    t_back->Branch("jet2_flavour", &jet2_flavour, "jet2_flavour/I");
    t_back->Branch("is_tau_event", &is_tau_event);


    TH1F *el_p_dr = new TH1F("el_p_dr", "DeltaR, reco elon and gen elon", 20,0,1);
    TH1F *el_m_dr = new TH1F("el_m_dr", "DeltaR, reco elon and gen elon", 20,0,1);
    el_p_dr->SetDirectory(0);
    el_m_dr->SetDirectory(0);


    unsigned int nEvents=0;
    unsigned int nSignal = 0;
    unsigned int nPileUp = 0;
    unsigned int nQQ=0;
    unsigned int nQQb=0;
    unsigned int nQGlu=0;
    unsigned int nGluGlu=0;
    unsigned int nTauTau=0;
    unsigned int nMuMu=0;
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

            UInt_t el_size, gen_size, jet_size, met_size;

            Int_t gen_id[GEN_SIZE], gen_status[GEN_SIZE];
            Int_t  gen_Mom0ID[GEN_SIZE], gen_Mom0Status[GEN_SIZE], gen_Mom1ID[GEN_SIZE], gen_Mom1Status[GEN_SIZE];
            Int_t  gen_Dau0ID[GEN_SIZE], gen_Dau0Status[GEN_SIZE], gen_Dau1ID[GEN_SIZE], gen_Dau1Status[GEN_SIZE];
            Float_t gen_Pt[GEN_SIZE], gen_Eta[GEN_SIZE], gen_Phi[GEN_SIZE], gen_E[GEN_SIZE];

            Float_t el_Pt[EL_SIZE], el_Eta[EL_SIZE], el_Phi[EL_SIZE], el_E[EL_SIZE], 
                    el_Charge[EL_SIZE];

            Int_t el_IDMedium[EL_SIZE];


            Float_t jet_Pt[JET_SIZE], jet_Eta[JET_SIZE], jet_Phi[JET_SIZE], jet_E[JET_SIZE],
                    jet_CSV[JET_SIZE], jet_CMVA[JET_SIZE], jet_partonflavour[JET_SIZE];

            Float_t evt_Gen_Weight;

            Int_t HLT_Ele23_WPLoose_Gsf, pu_NtrueInt;
            t1->SetBranchAddress("el_size", &el_size); //number of els in the event
            t1->SetBranchAddress("el_Pt", &el_Pt);
            t1->SetBranchAddress("el_Eta", &el_Eta);
            t1->SetBranchAddress("el_Phi", &el_Phi);
            t1->SetBranchAddress("el_E", &el_E);
            t1->SetBranchAddress("el_Charge", &el_Charge);
            t1->SetBranchAddress("el_IDMedium", &el_IDMedium);
            t1->SetBranchAddress("HLT_Ele23_WPLoose_Gsf", &HLT_Ele23_WPLoose_Gsf);


            t1->SetBranchAddress("jetAK4CHS_size", &jet_size);
            t1->SetBranchAddress("jetAK4CHS_Pt", &jet_Pt);
            t1->SetBranchAddress("jetAK4CHS_Eta", &jet_Eta);
            t1->SetBranchAddress("jetAK4CHS_Phi", &jet_Phi);
            t1->SetBranchAddress("jetAK4CHS_E", &jet_E);
            t1->SetBranchAddress("jetAK4CHS_CSVv2", &jet_CSV);
            t1->SetBranchAddress("jetAK4CHS_CMVAv2", &jet_CMVA);
            t1->SetBranchAddress("jetAK4CHS_PartonFlavour", &jet_partonflavour);



            t1->SetBranchAddress("gen_size", &gen_size); //number of muons in the event
            t1->SetBranchAddress("gen_Pt", &gen_Pt);
            t1->SetBranchAddress("gen_Eta", &gen_Eta);
            t1->SetBranchAddress("gen_Phi", &gen_Phi);
            t1->SetBranchAddress("gen_E", &gen_E);
            t1->SetBranchAddress("gen_ID", &gen_id);
            t1->SetBranchAddress("gen_Status", &gen_status);
            t1->SetBranchAddress("evt_Gen_Weight", &evt_Gen_Weight);
            t1->SetBranchAddress("pu_NtrueInt",&pu_NtrueInt);


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

            char out_buff[30000];
            bool print_out = false;

            for (int i=0; i<nEntries; i++) {
                t1->GetEntry(i);
                if(el_size > EL_SIZE || gen_size >GEN_SIZE) printf("WARNING: EL_SIZE OR GEN_SIZE TOO LARGE \n");
                if(met_size != 1) printf("WARNING: Met size not equal to 1\n");
                bool good_trigger = HLT_Ele23_WPLoose_Gsf;
                if(good_trigger &&
                        el_size >= 2 && ((abs(el_Charge[0] - el_Charge[1])) > 0.01) &&
                        el_IDMedium[0] && el_IDMedium[1] &&
                        el_Pt[0] > 26. &&  el_Pt[1] > 10. &&
                        abs(el_Eta[0]) < 2.4 && abs(el_Eta[1]) < 2.4){ 

                    //only want events with 2 oppositely charged leptons
                    if(el_Charge[0] >0){
                        el_p.SetPtEtaPhiE(el_Pt[0], el_Eta[0], el_Phi[0], el_E[0]);
                        el_m.SetPtEtaPhiE(el_Pt[1], el_Eta[1], el_Phi[1], el_E[1]);
                    }
                    else{
                        el_m.SetPtEtaPhiE(el_Pt[0], el_Eta[0], el_Phi[0], el_E[0]);
                        el_p.SetPtEtaPhiE(el_Pt[1], el_Eta[1], el_Phi[1], el_E[1]);
                    }

                    cm = el_p + el_m;
                    cm_m = cm.M();
                    if (cm_m >= 150.){
                        if(PRINT) sprintf(out_buff + strlen(out_buff),"\n \n Event %i \n", i);

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
                        int TAU = 15;
                        int PHOTON = 22;
                        int Z=23;
                        int GLUON = 21;
                        int PROTON = 2212;

                        int inc_1 =-1;
                        int inc_2 =-1;
                        int gen_el_p=-1;
                        int gen_el_m=-1;
                        int gen_pileup_el_p=-1;
                        int gen_pileup_el_m=-1;
                        int gen_tau_p=-1;
                        int gen_tau_m=-1;
                        int gen_e_p=-1;
                        int gen_e_m=-1;
                        int intermed=-1;

                        float quark_dir_eta;

                        bool signal_event = true;//whether it is an event with an asym or not
                        bool pileup_event = false; //whether it is an event with electrons coming from pileup, not real DY

                        is_tau_event = false;

                        for(int k=0; k<gen_size; k++){
                            if(gen_status[k] == INCIDENT_PARTICLE && 
                                    (abs(gen_id[k]) <=6  || gen_id[k] == GLUON) && 
                                    (abs(gen_Dau0ID[k]) == ELECTRON || gen_Dau0ID[k] == MUON || gen_Dau0ID[k] == Z || 
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

                            //record 2 scattered electrons
                            if(abs(gen_id[k]) == ELECTRON && 
                                    (gen_Mom0ID[k] == Z || gen_Mom0ID[k] == PHOTON || abs(gen_Mom0ID[k]) == TAU 
                                     || abs(gen_Mom0ID[k]) == MUON || (gen_status[k] == OUTGOING && gen_Mom0ID[k] != PROTON))) {
                                if(gen_id[k] == ELECTRON){
                                    if(gen_el_m == -1) gen_el_m = k;
                                    else{
                                        if(abs(gen_Mom0ID[k]) != TAU || abs(gen_Mom0ID[k]) != MUON) printf("WARNING: More than one el_m\n\n");
                                        if(PRINT) sprintf(out_buff + strlen(out_buff), "Extra el_m detected\n");
                                        print_out = true;
                                    }
                                }
                                if(gen_id[k] == -ELECTRON){
                                    if(gen_el_p == -1) gen_el_p = k;
                                    else{
                                        if(abs(gen_Mom0ID[k]) != TAU ||abs(gen_Mom0ID[k]) != MUON) printf("WARNING: More than one el_p\n\n");
                                        if(PRINT) sprintf(out_buff + strlen(out_buff), "Extra el_p detected\n");
                                        print_out = true;
                                    }
                                }
                            }
                            if(abs(gen_id[k]) == ELECTRON && gen_Mom0ID[k] == PROTON){
                                if(gen_id[k] == ELECTRON) gen_pileup_el_m = k;
                                if(gen_id[k] == -ELECTRON) gen_pileup_el_p = k;
                            }
                            //record tau's
                            if(abs(gen_id[k]) == TAU && 
                                    ( (gen_Mom0ID[k] == Z || gen_Mom0ID[k] == PHOTON ||
                                       gen_status[k] == OUTGOING)) ){
                                if(gen_id[k] == TAU){
                                    if(gen_tau_m == -1) gen_tau_m = k;
                                    else{
                                        printf("WARNING: More than one tau_m\n\n");
                                        if(PRINT) sprintf(out_buff + strlen(out_buff), "Extra tau_m detected\n");
                                        //print_out = true;
                                    }
                                }
                                if(gen_id[k] == -TAU){
                                    if(gen_tau_p == -1) gen_tau_p = k;
                                    else{
                                        printf("WARNING: More than one tau_p\n\n");
                                        if(PRINT) sprintf(out_buff + strlen(out_buff), "Extra tau_p detected\n");
                                        //print_out = true;
                                    }
                                }
                            }
                            if(PRINT){
                                /*
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
                                */


                                if(gen_id[k] == -ELECTRON){
                                    if(PRINT) sprintf(out_buff + strlen(out_buff),"el_p(stat = %i): \n"
                                            "   Mom1 ID: %i Mom1 Stat: %i Mom2 ID: %i Mom2 Stat %i \n"
                                            "   Dau1 ID: %i Dau1 Stat: %i Dau2 ID: %i Dau2 Stat %i \n",
                                            gen_status[k], 
                                            gen_Mom0ID[k], gen_Mom0Status[k], gen_Mom1ID[k], gen_Mom1Status[k],
                                            gen_Dau0ID[k], gen_Dau0Status[k], gen_Dau1ID[k], gen_Dau1Status[k]);
                                    if(PRINT) sprintf(out_buff + strlen(out_buff),"(pt, eta, phi,E) = %4.2f %4.2f %4.2f %4.2f \n", 
                                            gen_Pt[k], gen_Eta[k], gen_Phi[k], gen_E[k]);
                                }
                                if(gen_id[k] == ELECTRON){

                                    if(PRINT) sprintf(out_buff + strlen(out_buff),"el_m (stat = %i): \n"
                                            "   Mom1 ID: %i Mom1 Stat: %i Mom2 ID: %i Mom2 Stat %i \n"
                                            "   Dau1 ID: %i Dau1 Stat: %i Dau2 ID: %i Dau2 Stat %i \n",
                                            gen_status[k], 
                                            gen_Mom0ID[k], gen_Mom0Status[k], gen_Mom1ID[k], gen_Mom1Status[k],
                                            gen_Dau0ID[k], gen_Dau0Status[k], gen_Dau1ID[k], gen_Dau1Status[k]);
                                    if(PRINT) sprintf(out_buff + strlen(out_buff),"(pt, eta, phi,E) = %4.2f %4.2f %4.2f %4.2f \n", 
                                            gen_Pt[k], gen_Eta[k], gen_Phi[k], gen_E[k]);

                                }
                                /*
                                if(gen_id[k] == Z){
                                    if(PRINT) sprintf(out_buff + strlen(out_buff),"Z (ID = %i, status = %i): \n"
                                            "   Mom1 ID: %i Mom1 Stat: %i Mom2 ID: %i Mom2 Stat %i \n"
                                            "   Dau1 ID: %i Dau1 Stat: %i Dau2 ID: %i Dau2 Stat %i \n",
                                            gen_id[k], gen_status[k],
                                            gen_Mom0ID[k], gen_Mom0Status[k], gen_Mom1ID[k], gen_Mom1Status[k],
                                            gen_Dau0ID[k], gen_Dau0Status[k], gen_Dau1ID[k], gen_Dau1Status[k]);
                                }
                                */
                            }
                        }


                        if(gen_el_p != -1 && gen_el_m != -1) {
                            if(PRINT) sprintf(out_buff + strlen(out_buff),"el_p(stat = %i): \n"
                                    "   Mom1 ID: %i Mom1 Stat: %i Mom2 ID: %i Mom2 Stat %i \n"
                                    "   Dau1 ID: %i Dau1 Stat: %i Dau2 ID: %i Dau2 Stat %i \n",
                                    gen_status[gen_el_p],
                                    gen_Mom0ID[gen_el_p], gen_Mom0Status[gen_el_p], gen_Mom1ID[gen_el_p], 
                                    gen_Mom1Status[gen_el_p],
                                    gen_Dau0ID[gen_el_p], gen_Dau0Status[gen_el_p], gen_Dau1ID[gen_el_p], 
                                    gen_Dau1Status[gen_el_p]);

                            if(PRINT) sprintf(out_buff + strlen(out_buff),"el_m (stat = %i): \n"
                                    "   Mom1 ID: %i Mom1 Stat: %i Mom2 ID: %i Mom2 Stat %i \n"
                                    "   Dau1 ID: %i Dau1 Stat: %i Dau2 ID: %i Dau2 Stat %i \n",
                                    gen_status[gen_el_m],
                                    gen_Mom0ID[gen_el_m], gen_Mom0Status[gen_el_m], gen_Mom1ID[gen_el_m], 
                                    gen_Mom1Status[gen_el_m],
                                    gen_Dau0ID[gen_el_m], gen_Dau0Status[gen_el_m], gen_Dau1ID[gen_el_m], 
                                    gen_Dau1Status[gen_el_m]);

                        }
                        else if(gen_pileup_el_p != -1 && gen_pileup_el_m != -1){
                            pileup_event = true;
                            nPileUp++;
                            print_out = true;
                        }
                        else {
                            printf("WARNING: Unable to identify ElEl pair in event %i, skipping \n", i);
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
                            if(PRINT) sprintf(out_buff + strlen(out_buff),"inc1 (id %i stat = %i): \n"
                                    "   Mom1 ID: %i Mom1 Stat: %i Mom2 ID: %i Mom2 Stat %i \n"
                                    "   Dau1 ID: %i Dau1 Stat: %i Dau2 ID: %i Dau2 Stat %i \n",
                                    gen_id[inc_1], gen_status[inc_1],
                                    gen_Mom0ID[inc_1], gen_Mom0Status[inc_1], gen_Mom1ID[inc_1], 
                                    gen_Mom1Status[inc_1],
                                    gen_Dau0ID[inc_1], gen_Dau0Status[inc_1], gen_Dau1ID[inc_1], 
                                    gen_Dau1Status[inc_1]);
                            if(PRINT) sprintf(out_buff + strlen(out_buff),"inc2 (id %i, stat = %i): \n"
                                    "   Mom1 ID: %i Mom1 Stat: %i Mom2 ID: %i Mom2 Stat %i \n"
                                    "   Dau1 ID: %i Dau1 Stat: %i Dau2 ID: %i Dau2 Stat %i \n",
                                    gen_id[inc_2], gen_status[inc_2],
                                    gen_Mom0ID[inc_2], gen_Mom0Status[inc_2], gen_Mom1ID[inc_2], 
                                    gen_Mom1Status[inc_2],
                                    gen_Dau0ID[inc_2], gen_Dau0Status[inc_2], gen_Dau1ID[inc_2], 
                                    gen_Dau1Status[inc_2]);
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
                            if(abs(gen_Dau0ID[inc_1]) == MUON){
                                nMuMu++;
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
                        //RECO LEVEL
                        xF = abs(2.*cm.Pz()/13000.); 

                        // compute Colins soper angle with formula
                        double el_p_pls = (el_p.E()+el_p.Pz())/root2;
                        double el_p_min = (el_p.E()-el_p.Pz())/root2;
                        double el_m_pls = (el_m.E()+el_m.Pz())/root2;
                        double el_m_min = (el_m.E()-el_m.Pz())/root2;
                        double qt2 = cm.Px()*cm.Px()+cm.Py()*cm.Py();
                        double cm_m2 = cm.M2();
                        //cost_p = cos(theta)_r (reconstructed collins soper angle, sign
                        //may be 'wrong' if lepton pair direction is not the same as inital
                        //quark direction)
                        double cost = 2*(el_m_pls*el_p_min - el_m_min*el_p_pls)/sqrt(cm_m2*(cm_m2 + qt2));

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
                           el_m.Boost(beta);
                           el_p.Boost(beta);
                           p1.Boost(beta);
                           p2.Boost(beta);

                        // Now calculate the direction of the new z azis

                        TVector3 p1u = p1.Vect();
                        p1u.SetMag(1.0);
                        TVector3 p2u = p2.Vect();
                        p2u.SetMag(1.0);
                        TVector3 pzu = p1u - p2u;
                        pzu.SetMag(1.0);
                        el_m.RotateUz(pzu); 
                        double cost_r_b = el_m.CosTheta();
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
                        el1_pt = el_Pt[0];
                        el2_pt = el_Pt[1];
                        el1_eta = el_Eta[0];
                        el2_eta = el_Eta[1];

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


                        /*
                           if(jet_size >=2) nJets = 2;
                           else nJets = jet_size;
                           if(jet_size >=1 && jet_Pt[0] > 20. && std::abs(jet_Eta[0]) < 2.4){


                           } 
                           if(jet_size >=2 && jet_Pt[1] > 20. && std::abs(jet_Eta[1]) < 2.4){
                           jet2_pt = jet_Pt[1];
                           jet2_cmva = jet_CMVA[1];
                           jet2_flavour = jet_partonflavour[1];
                           jet2_b_weight = get_btag_weight(jet_Pt[1], jet_Eta[1],jet_partonflavour[1],btag_effs, b_reader);
                           }
                           */


                        //get el cut SFs

                        el_id_SF = get_el_SF(el1_pt, el1_eta, el_SF.h) * get_el_SF(el2_pt, el2_eta, el_SF.h);
                        pu_SF = get_pileup_SF(pu_NtrueInt, pu_SFs.pileup_ratio);


                        if(signal_event && !pileup_event){
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
                            TLorentzVector gen_el_p_vec, gen_el_m_vec;
                            gen_el_p_vec.SetPtEtaPhiE(gen_Pt[gen_el_p], gen_Eta[gen_el_p], gen_Phi[gen_el_p], gen_E[gen_el_p]);
                            gen_el_m_vec.SetPtEtaPhiE(gen_Pt[gen_el_m], gen_Eta[gen_el_m], gen_Phi[gen_el_m], gen_E[gen_el_m]);
                            float dr_p = gen_el_p_vec.DeltaR(el_p);
                            float dr_m = gen_el_m_vec.DeltaR(el_m);
                            el_p_dr->Fill(dr_p, gen_weight);
                            el_m_dr->Fill(dr_m, gen_weight);
                            bool evt_missmatched = false;
                            for(int j=2; j < el_size; j++){
                                TLorentzVector el_j;
                                el_j.SetPtEtaPhiE(el_Pt[j], el_Eta[j], el_Phi[j], el_E[j]);
                                float dr_j;
                                if(el_Charge[j] <0){
                                    dr_j = gen_el_m_vec.DeltaR(el_j);
                                    if (dr_j < dr_m && !evt_missmatched){ 
                                        evt_missmatched = true;
                                        mismatch++;
                                    }
                                }
                                else{
                                    dr_j = gen_el_p_vec.DeltaR(el_j);
                                    if(dr_j < dr_p && !evt_missmatched){
                                        evt_missmatched= true;
                                        mismatch++;
                                    }
                                } 
                            }
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
                        if(PRINT) memset(out_buff, 0, strlen(out_buff));

                    }
                } 
            }

            f1->Close();
            printf("moving on to next file, currently %i events %i Taus %i fails %i mismatch \n\n", nEvents, nTauTau, nFailedID, mismatch);
        }
    }
    fclose(root_files);
    printf("There were %i qqbar, %i qGlu (%i of them tautau, %i mumu) in %i kept events in %i files. "
            "There were also %i background events (%i qq and %i gg) "
            "There were also %i pileup events \n"
            "There were %i Failed ID's \n" , 
            nQQb, nQGlu, nTauTau, nMuMu, nSignal, nFiles, nQQ + nGluGlu, nQQ, nGluGlu, nPileUp, nFailedID);
    //printf("Ran on MC data and produced templates with %i events\n", nEvents);
    fout->cd();




    t_signal->Write();
    t_back->Write();
    el_p_dr->Write();
    el_m_dr->Write();


    printf("Writing output to file at %s \n", fout_name.Data());

    fout->Close();

    return;
}
