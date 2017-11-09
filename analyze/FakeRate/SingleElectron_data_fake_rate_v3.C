


#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <iostream>
#include <cstring>
#include <algorithm>
#include "TFile.h"
#define EL_SIZE 200
#define JET_SIZE 20

const double root2 = sqrt(2);
const char* filename("SinglePhoton_files_nov2.txt");
const TString fout_name("FakeRate/SingleElectron_data_fake_rate_v3_nov6.root");


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

void SingleElectron_data_fake_rate_v3()
{



    TTree *tout= new TTree("T_data", "Tree with reco events");
    tout->SetDirectory(0);
    Double_t cm_m, xF, cost_r, mu1_pt, mu2_pt, mu1_eta, mu2_eta, jet1_pt, jet2_pt,
             jet1_cmva, jet1_eta, jet2_cmva, jet2_eta;
    Int_t nJets;
    Float_t met_pt;
    TLorentzVector mu_p, mu_m, cm, q1, q2;

    Float_t pt_bins[] = {0,20, 30, 45, 70,100, 200, 1000};
    int n_pt_bins = 7;
    Float_t eta_bins[] = {0, 0.9, 2.4};
    int n_eta_bins = 2;

    TH2D *h_pass = new TH2D("h_pass", "Rate of passing ISO cut for single electrons",  n_eta_bins, eta_bins, n_pt_bins, pt_bins);
    TH2D *h_total = new TH2D("h_total", "Total number of single electrons",  n_eta_bins, eta_bins, n_pt_bins, pt_bins);



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

        UInt_t el_size, jet_size, met_size;

        Float_t el_Pt[EL_SIZE], el_Eta[EL_SIZE], el_Phi[EL_SIZE], el_E[EL_SIZE],
                el_Charge[EL_SIZE];

        Int_t el_IDMedium[EL_SIZE], el_IDMedium_NoIso[EL_SIZE];

        Float_t jet_Pt[JET_SIZE], jet_Eta[JET_SIZE], jet_Phi[JET_SIZE], jet_E[JET_SIZE],
                jet_CSV[JET_SIZE], jet_CMVA[JET_SIZE], jet_partonflavour[JET_SIZE];

        Int_t HLT_Photon22, HLT_Photon30, HLT_Photon36;
        t1->SetBranchAddress("el_size", &el_size); //number of els in the event
        t1->SetBranchAddress("el_Pt", &el_Pt);
        t1->SetBranchAddress("el_Eta", &el_Eta);
        t1->SetBranchAddress("el_Phi", &el_Phi);
        t1->SetBranchAddress("el_E", &el_E);
        t1->SetBranchAddress("el_Charge", &el_Charge);
        t1->SetBranchAddress("el_IDMedium", &el_IDMedium);
        t1->SetBranchAddress("el_IDMedium_NoIso", &el_IDMedium_NoIso);
        //t1->SetBranchAddress("HLT_Ele27_WPTight_Gsf", &HLT_El);

        t1->SetBranchAddress("HLT_Photon22", &HLT_Photon22);
        t1->SetBranchAddress("HLT_Photon30", &HLT_Photon30);
        t1->SetBranchAddress("HLT_Photon36", &HLT_Photon36);


        t1->SetBranchAddress("jetAK4CHS_size", &jet_size);
        t1->SetBranchAddress("jetAK4CHS_Pt", &jet_Pt);
        t1->SetBranchAddress("jetAK4CHS_Eta", &jet_Eta);
        t1->SetBranchAddress("jetAK4CHS_Phi", &jet_Phi);
        t1->SetBranchAddress("jetAK4CHS_E", &jet_E);
        t1->SetBranchAddress("jetAK4CHS_CSVv2", &jet_CSV);
        t1->SetBranchAddress("jetAK4CHS_CMVAv2", &jet_CMVA);

        t1->SetBranchAddress("met_size", &met_size);
        t1->SetBranchAddress("met_Pt", &met_pt);


        unsigned int nEntries =  t1->GetEntries();
        printf("there are %i entries in this tree\n", nEntries);
        for (int i=0; i<nEntries; i++) {
            t1->GetEntry(i);
            if(el_size > EL_SIZE) printf("WARNING: MU_SIZE TOO LARGE \n");
            if(met_size != 1) printf("WARNING: Met size not equal to 1\n");
            bool good_trigger = HLT_Photon22 || HLT_Photon30 || HLT_Photon36;
            if( good_trigger &&
                    el_size >= 1 && el_Pt[0] > 10. && el_IDMedium_NoIso[0] && abs(el_Eta[0]) < 2.5
                    && (el_size == 1 || (el_size >= 2  && !el_IDMedium_NoIso[1])) ){ 
                bool no_other_els = true;

                for (int j=1; j < el_size; j++){
                    if(el_IDMedium_NoIso[j]) no_other_els = false;
                    //if(el_IDMedium[j] && !el_IDMedium_NoIso[j] ) printf("HI\n");
                }

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
                float iso_0 = el_IDMedium[0];
                bool no_bjets = has_no_bjets(nJets, jet1_pt, jet2_pt, jet1_cmva, jet2_cmva);

                if( no_other_els && no_bjets && met_pt < 50){
                    nEvents++;
                    if(iso_0){
                        nPass++;
                        h_pass->Fill(abs(el_Eta[0]), el_Pt[0], 1);
                    }
                    h_total->Fill(abs(el_Eta[0]), el_Pt[0], 1);
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

    h_rate->Write();
    h_total->Write();

    printf("Writing out put to %s \n", fout_name.Data());

    fout->Close();
    return;
}
