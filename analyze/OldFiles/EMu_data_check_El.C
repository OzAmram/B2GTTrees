
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <iostream>
#include <cstring>
#include <algorithm>
#include "TFile.h"
#include "../TemplateMaker.C"

#define MU_SIZE 200
#define EL_SIZE 200
#define JET_SIZE 20

const double root2 = sqrt(2);
const char* filename("SingleElectron_files_aug29.txt");
const TString fout_name("output_files/EMu_SingleElectron_data_nov3.root");

const bool data_2016 = true;

bool is_empty_line(const char *s) {
    while (*s != '\0') {
        if (!isspace(*s))
            return false;
        s++;
    }
    return true;
}

void EMu_data_check_El()
{



    TTree *tout= new TTree("T_data", "Tree with reco events");
    tout->SetDirectory(0);
    Double_t cm_m, xF, cost_r, mu1_pt, el1_pt, jet1_pt, jet2_pt, gen_weight,
             jet1_cmva, jet2_cmva, mu1_eta, el1_eta, jet1_eta, jet2_eta;
    Float_t met_pt;
    Int_t nJets;
    TLorentzVector mu, el, cm, q1, q2;
    tout->Branch("mu1_pt", &mu1_pt, "mu1_pt/D");
    tout->Branch("mu1_eta", &mu1_eta, "mu1_eta/D");
    tout->Branch("el1_pt", &el1_pt, "el1_pt/D");
    tout->Branch("el1_eta", &el1_eta, "el1_eta/D");
    tout->Branch("el", "TLorentzVector", &el);
    tout->Branch("mu", "TLorentzVector", &mu);
    tout->Branch("jet1_pt", &jet1_pt, "jet1_pt/D");
    tout->Branch("jet1_eta", &jet1_eta, "jet1_eta/D");
    tout->Branch("jet1_CMVA", &jet1_cmva, "jet1_CMVA/D");
    tout->Branch("jet2_pt", &jet2_pt, "jet2_pt/D");
    tout->Branch("jet2_eta", &jet2_eta, "jet2_eta/D");
    tout->Branch("jet2_CMVA", &jet2_cmva, "jet2_CMVA/D");
    tout->Branch("nJets", &nJets, "nJets/I");
    tout->Branch("met_pt", &met_pt, "met_Pt/F");



    printf("Opeing file list...\n");
    FILE *root_files = fopen(filename, "r");
    if(root_files == NULL){
        printf("Unable to open file list %s\n", filename);
        return;
    }
    char lines[300];
    int nFiles =0;
    unsigned int nEvents=0;
    while(fgets(lines, 300, root_files)){
        if(lines[0] == '#' || is_empty_line(lines)) continue; // comment line
        nFiles++;
        char * cur_file;

        char * end;
        //remove trailing whitespace
        end = lines + strlen(lines) - 1;
        while(end > lines && isspace((char)*end)) end--;
        // Write new null terminator
        *(end+1) = 0;

        printf("Opening file: %s \n", lines);
        TFile *f1=  TFile::Open(lines);
        f1->cd("B2GTTreeMaker");
        TDirectory *subdir = gDirectory;
        TTree *t1 = (TTree *)subdir->Get("B2GTree");

        UInt_t mu_size, met_size, jet_size, el_size;
        Float_t mu_Pt[MU_SIZE], mu_Eta[MU_SIZE], mu_Phi[MU_SIZE], mu_E[MU_SIZE], 
                mu_IsTightMuon[MU_SIZE], mu_Charge[MU_SIZE];

        Float_t el_Pt[EL_SIZE], el_Eta[EL_SIZE], el_Phi[EL_SIZE], el_E[EL_SIZE],
                el_Charge[EL_SIZE];

        Int_t el_IDMedium[EL_SIZE];


        Float_t mu_SumChargedHadronPt[MU_SIZE], mu_SumNeutralHadronPt[MU_SIZE], mu_SumPUPt[MU_SIZE], mu_SumPhotonPt[MU_SIZE];


        Float_t jet_Pt[JET_SIZE], jet_Eta[JET_SIZE], jet_Phi[JET_SIZE], jet_E[JET_SIZE],
                jet_CSV[JET_SIZE], jet_CMVA[JET_SIZE];

        Int_t HLT_IsoMu, HLT_IsoTkMu, HLT_El;
        t1->SetBranchAddress("mu_size", &mu_size); //number of muons in the event
        t1->SetBranchAddress("mu_Pt", &mu_Pt);
        t1->SetBranchAddress("mu_Eta", &mu_Eta);
        t1->SetBranchAddress("mu_Phi", &mu_Phi);
        t1->SetBranchAddress("mu_E", &mu_E);
        t1->SetBranchAddress("mu_Charge", &mu_Charge);

        t1->SetBranchAddress("el_size", &el_size); //number of els in the event
        t1->SetBranchAddress("el_Pt", &el_Pt);
        t1->SetBranchAddress("el_Eta", &el_Eta);
        t1->SetBranchAddress("el_Phi", &el_Phi);
        t1->SetBranchAddress("el_E", &el_E);
        t1->SetBranchAddress("el_Charge", &el_Charge);
        t1->SetBranchAddress("el_IDMedium", &el_IDMedium);
        t1->SetBranchAddress("HLT_Ele27_WPTight_Gsf", &HLT_El);

        t1->SetBranchAddress("mu_IsTightMuon", &mu_IsTightMuon);
        t1->SetBranchAddress("mu_SumChargedHadronPt", &mu_SumChargedHadronPt);
        t1->SetBranchAddress("mu_SumNeutralHadronPt", &mu_SumNeutralHadronPt);
        t1->SetBranchAddress("mu_SumPUPt", &mu_SumPUPt);
        t1->SetBranchAddress("mu_SumPhotonPt", &mu_SumPhotonPt);


        t1->SetBranchAddress("jetAK4CHS_size", &jet_size);
        t1->SetBranchAddress("jetAK4CHS_Pt", &jet_Pt);
        t1->SetBranchAddress("jetAK4CHS_Eta", &jet_Eta);
        t1->SetBranchAddress("jetAK4CHS_Phi", &jet_Phi);
        t1->SetBranchAddress("jetAK4CHS_E", &jet_E);
        t1->SetBranchAddress("jetAK4CHS_CSVv2", &jet_CSV);
        t1->SetBranchAddress("jetAK4CHS_CMVAv2", &jet_CMVA);

        t1->SetBranchAddress("HLT_IsoMu24", &HLT_IsoMu);
        t1->SetBranchAddress("HLT_IsoTkMu24", &HLT_IsoTkMu);

        t1->SetBranchAddress("met_size", &met_size);
        t1->SetBranchAddress("met_Pt", &met_pt);

        unsigned int nEntries =  t1->GetEntries();
        printf("there are %i entries in this tree\n", nEntries);
        for (int i=0; i<nEntries; i++) {
            t1->GetEntry(i);
            if(met_size != 1) printf("WARNING: Met size not equal to 1\n");
            if(mu_size > MU_SIZE) printf("Warning: too many muons\n");
            bool good_trigger = HLT_El;
            if(good_trigger &&
                        mu_size >= 1 && el_size >=1 && 
                        ((abs(mu_Charge[0] - el_Charge[0])) > 0.01) &&
                        mu_IsTightMuon[0] && el_IDMedium[0] && 
                        el_Pt[0] > 29. && mu_Pt[0] > 10. &&
                        abs(mu_Eta[0]) < 2.4 && abs(el_Eta[0]) < 2.5){ 

                //See https://twiki.cern.ch/twiki/bin/viewauth/CMS/SWGuideMuonIdRun2 for iso cuts
                float iso_0 = (mu_SumChargedHadronPt[0] + max(0., mu_SumNeutralHadronPt[0] + mu_SumPhotonPt[0] - 0.5 * mu_SumPUPt[0]))/mu_Pt[0];
                const float tight_iso = 0.15;
                const float loose_iso = 0.25;

                mu.SetPtEtaPhiE(mu_Pt[0], mu_Eta[0], mu_Phi[0], mu_E[0]);
                el.SetPtEtaPhiE(el_Pt[0], el_Eta[0], el_Phi[0], el_E[0]);


                cm = el + mu;
                cm_m = cm.M();

                nJets =0;
                for(int j=0; j < jet_size; j++){
                    if(jet_Pt[j] > 20. && std::abs(jet_Eta[j]) < 2.4){
                        if(nJets == 1){
                            jet2_pt = jet_Pt[j];
                            jet2_eta = jet_Eta[j];
                            jet2_cmva = jet_CMVA[j];
                            nJets =2;
                            break;
                        }
                        else if(nJets ==0){
                            jet1_pt = jet_Pt[j];
                            jet1_eta = jet_Eta[j];
                            jet1_cmva = jet_CMVA[j];
                            nJets = 1;
                        }
                    }
                }
                bool no_bjets = has_no_bjets(nJets, jet1_pt, jet2_pt, jet1_cmva, jet2_cmva);

                if (no_bjets && (met_pt < 50) && (cm_m >=150.) && iso_0 < tight_iso){


                    nJets =0;
                    for(int j=0; j < jet_size; j++){
                        if(jet_Pt[j] > 20. && std::abs(jet_Eta[j]) < 2.4){
                            if(nJets == 1){
                                jet2_pt = jet_Pt[j];
                                jet2_eta = jet_Eta[j];
                                jet2_cmva = jet_CMVA[j];
                                nJets =2;
                                break;
                            }
                            else if(nJets ==0){
                                jet1_pt = jet_Pt[j];
                                jet1_eta = jet_Eta[j];
                                jet1_cmva = jet_CMVA[j];
                                nJets = 1;
                            }
                        }
                    }
                    el1_pt = el_Pt[0];
                    el1_eta = el_Eta[0];
                    mu1_pt = mu_Pt[0];
                    mu1_eta = mu_Eta[0];
                    tout->Fill();

                    nEvents++;

                }

            }
        }
        f1->Close();
        printf("moving on to next file, currently %i events \n\n", nEvents);
    }
    printf("Ran on data from %i Files and produced template with %i Events \n", 
            nFiles, nEvents );
    printf("Writing out put to %s \n", fout_name.Data());

    TFile *fout = TFile::Open(fout_name, "RECREATE");
    fout->cd();
    tout->Write();

    fout->Close();
    return;
}
