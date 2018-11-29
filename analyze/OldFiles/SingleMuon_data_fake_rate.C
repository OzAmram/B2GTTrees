
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <iostream>
#include <cstring>
#include <algorithm>
#include "TFile.h"
#define MU_SIZE 200
#define JET_SIZE 20

const double root2 = sqrt(2);
const char* filename("SingleMuon_files_test.txt");
const TString fout_name("FakeRate/SingleMuon_data_fake_rate_sep1_test.root");
const double alpha = 0.05;

const bool data_2016 = true;

bool is_empty_line(const char *s) {
    while (*s != '\0') {
        if (!isspace(*s))
            return false;
        s++;
    }
    return true;
}

bool has_no_bjets(Int_t nJets, Double_t jet1_pt, Double_t jet2_pt, 
        Double_t jet1_cmva, Double_t jet2_cmva){
    Double_t med_btag = 0.4432;
    if(nJets ==0) return true;
    else if(nJets == 1){
        if(jet1_pt < 20.) return true;
        else return jet1_cmva < med_btag;
    }
    else{
        return (jet1_pt < 20. || jet1_cmva < med_btag) && (jet2_pt < 20. || jet2_cmva < med_btag);
    }
}

void SingleMuon_data_fake_rate()
{



    TTree *tout= new TTree("T_data", "Tree with reco events");
    tout->SetDirectory(0);
    Double_t cm_m, xF, cost_r, mu1_pt, mu2_pt, mu1_eta, mu2_eta, jet1_pt, jet2_pt,
             jet1_cmva, jet1_eta, jet2_cmva, jet2_eta;
    Int_t nJets;
    Float_t met_pt;
    TLorentzVector mu_p, mu_m, cm, q1, q2;

    Float_t pt_bins[] = {0,35, 45, 55,65,80, 120, 200, 400};
    Float_t eta_bins[] = {0, 0.9, 2.4};
    int n_eta_bins = 2;
    int n_pt_bins = 8;

    TH2D *h_pass_HLT = new TH2D("h_pass_HLT", "Rate of passing ISO cut for single muons with HLT",  n_eta_bins, eta_bins, n_pt_bins, pt_bins);
    TH2D *h_total_HLT = new TH2D("h_total_HLT", "Total number of single muons with HLT",  n_eta_bins, eta_bins, n_pt_bins, pt_bins);

    TH2D *h_pass_noHLT = new TH2D("h_pass_noHLT", "Rate of passing ISO cut for single muons with no HLT",  n_eta_bins, eta_bins, n_pt_bins, pt_bins);
    TH2D *h_total_noHLT = new TH2D("h_total_noHLT", "Total number of single muons with no HLT",  n_eta_bins, eta_bins, n_pt_bins, pt_bins);


    printf("Opeing file list...\n");
    FILE *root_files = fopen(filename, "r");
    if(root_files == NULL){
        printf("Unable to open file list %s\n", filename);
        return;
    }
    char lines[300];
    int nFiles =0;
    unsigned int nEvents=0;
    unsigned int nPass=0;
    unsigned int nTrkIso=0;
    unsigned int nLoose=0;
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

        UInt_t mu_size, met_size, jet_size;
        Float_t mu_Pt[MU_SIZE], mu_Eta[MU_SIZE], mu_Phi[MU_SIZE], mu_E[MU_SIZE], 
                mu_IsTightMuon[MU_SIZE], mu_Charge[MU_SIZE];

        Float_t mu_SumChargedHadronPt[MU_SIZE], mu_SumNeutralHadronPt[MU_SIZE], mu_SumPUPt[MU_SIZE], mu_SumPhotonPt[MU_SIZE];


        Float_t jet_Pt[JET_SIZE], jet_Eta[JET_SIZE], jet_Phi[JET_SIZE], jet_E[JET_SIZE],
                jet_CSV[JET_SIZE], jet_CMVA[JET_SIZE];

        Int_t HLT_IsoMu, HLT_IsoTkMu, evt_NIsoTrk;
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

        t1->SetBranchAddress("evt_NIsoTrk", &evt_NIsoTrk);

        if(data_2016){
            t1->SetBranchAddress("jetAK4CHS_size", &jet_size);
            t1->SetBranchAddress("jetAK4CHS_Pt", &jet_Pt);
            t1->SetBranchAddress("jetAK4CHS_Eta", &jet_Eta);
            t1->SetBranchAddress("jetAK4CHS_Phi", &jet_Phi);
            t1->SetBranchAddress("jetAK4CHS_E", &jet_E);
            t1->SetBranchAddress("jetAK4CHS_CSVv2", &jet_CSV);
            t1->SetBranchAddress("jetAK4CHS_CMVAv2", &jet_CMVA);

            t1->SetBranchAddress("HLT_IsoMu24", &HLT_IsoMu);
            t1->SetBranchAddress("HLT_IsoTkMu24", &HLT_IsoTkMu);
        }

        else{
            t1->SetBranchAddress("jetAK4Puppi_size", &jet_size);
            t1->SetBranchAddress("jetAK4Puppi_Pt", &jet_Pt);
            t1->SetBranchAddress("jetAK4Puppi_Eta", &jet_Eta);
            t1->SetBranchAddress("jetAK4Puppi_Phi", &jet_Phi);
            t1->SetBranchAddress("jetAK4Puppi_E", &jet_E);
            t1->SetBranchAddress("jetAK4Puppi_CSVv2", &jet_CSV);
            t1->SetBranchAddress("jetAK4Puppi_CMVAv2", &jet_CMVA);

            t1->SetBranchAddress("HLT_IsoMu20", &HLT_IsoMu);
            t1->SetBranchAddress("HLT_IsoTkMu20", &HLT_IsoTkMu);
        }
        t1->SetBranchAddress("met_size", &met_size);
        t1->SetBranchAddress("met_Pt", &met_pt);

        unsigned int nEntries =  t1->GetEntries();
        printf("there are %i entries in this tree\n", nEntries);
        for (int i=0; i<nEntries; i++) {
            t1->GetEntry(i);
            if(met_size != 1) printf("WARNING: Met size not equal to 1\n");
            if(mu_size > MU_SIZE) printf("Warning: too many muons\n");
            bool good_trigger = HLT_IsoMu || HLT_IsoTkMu;
            if( mu_size >= 1 && mu_Pt[0] > 26. && mu_IsTightMuon[0] && abs(mu_Eta[0]) < 2.4
                    && (mu_size == 1 || (mu_Pt[1] < 10. || !mu_IsTightMuon[1] || abs(mu_Eta[1]) > 2.4))){
                //Want events with only 1 muon
                //See https://twiki.cern.ch/twiki/bin/viewauth/CMS/SWGuideMuonIdRun2 for iso cuts

                //get jets
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
                float iso_0 = (mu_SumChargedHadronPt[0] + max(0., mu_SumNeutralHadronPt[0] + mu_SumPhotonPt[0] - 0.5 * mu_SumPUPt[0]))/mu_Pt[0];
                const float tight_iso = 0.15;
                const float loose_iso = 0.25;

                if(no_bjets && met_pt < 50){
                    if(good_trigger){
                        if(iso_0 < tight_iso){
                            nPass++;
                            h_pass_HLT->Fill(abs(mu_Eta[0]), mu_Pt[0], 1);
                        }
                        else if (iso_0 < loose_iso){
                            nLoose++;
                        }
                        h_total_HLT->Fill(abs(mu_Eta[0]), mu_Pt[0], 1);
                    }
                    else{
                        if(iso_0 < tight_iso){
                            nPass++;
                            h_pass_noHLT->Fill(abs(mu_Eta[0]), mu_Pt[0], 1);
                        }
                        h_total_noHLT->Fill(abs(mu_Eta[0]), mu_Pt[0], 1);


                    }
                    nEvents++;

                }
            }
        }
        f1->Close();
        printf("moving on to next file, currently %i events \n\n", nEvents);
    }
    printf("Ran on data from %i Files and produced template with %i pass in %i Events (%.0f%%)\n", 
            nFiles, nPass, nEvents, 100*((float) nPass)/ ((float) nEvents));
    Double_t HLT_pass = (h_pass_HLT->Integral()) / (h_total_HLT->Integral());
    Double_t noHLT_pass = (h_pass_noHLT->Integral()) / (h_total_noHLT->Integral());
    Double_t loose_frac = nLoose/ (h_total_HLT->Integral() + nPass);
    printf("Loose Frac: %.0f%% HLT rate: %.0f%%. noHLT rate: %.0f%% \n", 100*loose_frac,  100*HLT_pass, 100*noHLT_pass);

    TH2D* h_rate_HLT = (TH2D *) h_pass_HLT->Clone("h_rate_HLT");
    h_rate_HLT->Divide(h_total_HLT);
    TH2D* h_rate_noHLT = (TH2D *) h_pass_noHLT->Clone("h_rate_noHLT");
    h_rate_noHLT->Divide(h_total_noHLT);

    TFile *fout = TFile::Open(fout_name, "RECREATE");
    fout->cd();

    h_rate_HLT->Write();
    h_total_HLT->Write();

    h_rate_noHLT->Write();
    h_total_noHLT->Write();
    printf("Writing out put to %s \n", fout_name.Data());

    fout->Close();
    return;
}
