

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
const char* filename("SingleMuon_files_sep25.txt");
const TString fout_name("FakeRate/root_files/SingleMuon_data_fake_rate_v2_nov28.root");
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

bool in_Z_window(Double_t m){
    double_t Z_mass_low = 91.2 -7.;
    double_t Z_mass_high = 91.2 + 7.;
    return (m > Z_mass_low ) && (m < Z_mass_high);
}


void SingleMuon_data_fake_rate_v2()
{



    TTree *tout= new TTree("T_data", "Tree with reco events");
    tout->SetDirectory(0);
    Double_t cm_m, xF, cost_r, mu1_pt, mu2_pt, mu1_eta, mu2_eta, jet1_pt, jet2_pt,
             jet1_cmva, jet1_eta, jet2_cmva, jet2_eta;
    Int_t nJets;
    Float_t met_pt;
    TLorentzVector mu_p, mu_m, cm, q1, q2;
    Double_t mu_pt, mu_eta;
    Bool_t pass;
    
    tout->Branch("mu_pt", &mu_pt);
    tout->Branch("mu_eta", &mu_eta);
    tout->Branch("pass", &pass);

    Float_t pt_bins[] = {0, 25, 35, 50, 90, 1000};
    int n_pt_bins = 5;
    Float_t eta_bins[] = {0, 0.9, 2.4};
    int n_eta_bins = 2;

    TH2D *h_pass = new TH2D("h_pass", "Rate of passing ISO cut",  n_eta_bins, eta_bins, n_pt_bins, pt_bins);
    TH2D *h_total = new TH2D("h_total", "Total number of muons",  n_eta_bins, eta_bins, n_pt_bins, pt_bins);



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
            bool good_trigger = HLT_IsoMu || HLT_IsoTkMu;
            if( mu_size >= 3 && mu_Pt[0] > 26. && mu_IsTightMuon[0] && abs(mu_Eta[0]) < 2.4
                    && mu_Pt[1] > 10. && mu_IsTightMuon[1] && abs(mu_Eta[1]) < 2.4
                    && mu_Pt[2] > 10. && mu_IsTightMuon[2] && abs(mu_Eta[2]) < 2.4){
                //Want events with 3 muons, 2 from Z and 1 extra
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

                TLorentzVector mu0, mu1, mu2;
                mu0.SetPtEtaPhiE(mu_Pt[0], mu_Eta[0], mu_Phi[0], mu_E[0]);
                mu1.SetPtEtaPhiE(mu_Pt[1], mu_Eta[1], mu_Phi[1], mu_E[1]);
                mu2.SetPtEtaPhiE(mu_Pt[2], mu_Eta[2], mu_Phi[2], mu_E[2]);

                //mu+ and mu- from Z, extra muon
                int mu_p, mu_m, mu_extra;
                mu_p = -1;
                mu_m = -1;
                mu_extra = -1;

                Double_t m01 = (mu0 + mu1).M();
                Double_t m02 = (mu0 + mu2).M();
                Double_t m12 = (mu1 + mu2).M();


                bool m01_in_Z = in_Z_window(m01);
                bool m02_in_Z = in_Z_window(m02);
                bool m12_in_Z = in_Z_window(m12);

                float iso[3];
                iso[0] = (mu_SumChargedHadronPt[0] + max(0., mu_SumNeutralHadronPt[0] + mu_SumPhotonPt[0] - 0.5 * mu_SumPUPt[0]))/mu_Pt[0];
                iso[1] = (mu_SumChargedHadronPt[1] + max(0., mu_SumNeutralHadronPt[1] + mu_SumPhotonPt[1] - 0.5 * mu_SumPUPt[1]))/mu_Pt[1];
                iso[2] = (mu_SumChargedHadronPt[2] + max(0., mu_SumNeutralHadronPt[2] + mu_SumPhotonPt[2] - 0.5 * mu_SumPUPt[2]))/mu_Pt[2];
                const float tight_iso = 0.15;
                const float loose_iso = 0.25;

                if(m01_in_Z && !m02_in_Z && !m12_in_Z && mu_Charge[0] * mu_Charge[1] < 0 && iso[0] < tight_iso && iso[1] < tight_iso){
                    mu_extra = 2;
                    if(mu_Charge[0] >0){
                        mu_p = 0;
                        mu_m = 1;
                    }
                    else{
                        mu_p = 1;
                        mu_m =0;
                    }
                }
                else if(!m01_in_Z && m02_in_Z && !m12_in_Z && mu_Charge[0] * mu_Charge[2] < 0 && iso[0] < tight_iso && iso[2] < tight_iso){
                    mu_extra = 1;
                    if(mu_Charge[0] >0){
                        mu_p = 0;
                        mu_m = 2;
                    }
                    else{
                        mu_p = 2;
                        mu_m =0;
                    }
                }
                else if(!m01_in_Z && !m02_in_Z && m12_in_Z  && mu_Charge[1] * mu_Charge[2] < 0 && iso[1] < tight_iso && iso[2] < tight_iso){
                    mu_extra = 0;
                    if(mu_Charge[1] >0){
                        mu_p = 1;
                        mu_m = 2;
                    }
                    else{
                        mu_p = 2;
                        mu_m =2;
                    }
                }


                if( (mu_p != -1) && no_bjets && met_pt < 25.){
                    mu_pt = mu_Pt[mu_extra];
                    mu_eta = mu_Eta[mu_extra];
                    pass = iso[mu_extra] < tight_iso;
                    nEvents++;
                    if(iso[mu_extra] < tight_iso){
                        nPass++;
                        h_pass->Fill(abs(mu_Eta[mu_extra]), mu_Pt[mu_extra], 1);
                    }
                    h_total->Fill(abs(mu_Eta[mu_extra]), mu_Pt[mu_extra], 1);
                    tout->Fill();
                }

            }
        }
        f1->Close();
        printf("moving on to next file, currently %i events \n\n", nEvents);
    }
    printf("Ran on data from %i Files and produced template with %i pass in %i Events (%.0f%%)\n", 
            nFiles, nPass, nEvents, 100*((float) nPass)/ ((float) nEvents));

    TH2D* h_rate = (TH2D *) h_pass->Clone("h_rate");
    h_rate->Divide(h_total);

    TFile *fout = TFile::Open(fout_name, "RECREATE");
    fout->cd();
    tout->Write();
    h_rate->Write();
    h_total->Write();

    printf("Writing out put to %s \n", fout_name.Data());

    fout->Close();
    return;
}
