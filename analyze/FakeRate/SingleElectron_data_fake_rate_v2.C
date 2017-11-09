


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
const char* filename("SingleElectron_files_aug29.txt");
const TString fout_name("FakeRate/SingleElectron_data_fake_rate_oct23.root");
const double alpha = 0.05;


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

void SingleElectron_data_fake_rate_v2()
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

        Int_t HLT_El;
        t1->SetBranchAddress("el_size", &el_size); //number of els in the event
        t1->SetBranchAddress("el_Pt", &el_Pt);
        t1->SetBranchAddress("el_Eta", &el_Eta);
        t1->SetBranchAddress("el_Phi", &el_Phi);
        t1->SetBranchAddress("el_E", &el_E);
        t1->SetBranchAddress("el_Charge", &el_Charge);
        t1->SetBranchAddress("el_IDMedium", &el_IDMedium);
        t1->SetBranchAddress("el_IDMedium_NoIso", &el_IDMedium_NoIso);
        t1->SetBranchAddress("HLT_Ele27_WPTight_Gsf", &HLT_El);


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
            bool good_trigger = HLT_El;
            if( good_trigger &&
                    el_size >= 3 && el_Pt[0] > 29. && el_IDMedium_NoIso[0] && abs(el_Eta[0]) < 2.4
                    && el_Pt[1] > 10. && el_IDMedium_NoIso[1] && abs(el_Eta[1]) < 2.4
                    && el_Pt[2] > 10. && el_IDMedium_NoIso[2] && abs(el_Eta[2]) < 2.4){
                //Want events with 3 elons, 2 from Z and 1 extra
                //See https://twiki.cern.ch/twiki/bin/viewauth/CMS/SWGuideelonIdRun2 for iso cuts

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
                double_t Z_mass_low = 91.2 -7;
                double_t Z_mass_high = 91.2 + 7;

                TLorentzVector el0, el1, el2;
                el0.SetPtEtaPhiE(el_Pt[0], el_Eta[0], el_Phi[0], el_E[0]);
                el1.SetPtEtaPhiE(el_Pt[1], el_Eta[1], el_Phi[1], el_E[1]);
                el2.SetPtEtaPhiE(el_Pt[2], el_Eta[2], el_Phi[2], el_E[2]);

                //el+ and el- from Z, extra elon
                int el_p, el_m, el_extra;
                el_p = -1;
                el_m = -1;
                el_extra = -1;

                Double_t m01 = (el0 + el1).M();
                Double_t m02 = (el0 + el2).M();
                Double_t m12 = (el1 + el2).M();

                Int_t iso[3];
                iso[0] = el_IDMedium[0]; 
                iso[1] = el_IDMedium[1];
                iso[2] = el_IDMedium[2];

                if(m01 > Z_mass_low && m01 < Z_mass_high && el_Charge[0] * el_Charge[1] < 0 && iso[0] && iso[1]){
                    el_extra = 2;
                    if(el_Charge[0] >0){
                        el_p = 0;
                        el_m = 1;
                    }
                    else{
                        el_p = 1;
                        el_m =0;
                    }
                }
                else if(m02 > Z_mass_low && m02 < Z_mass_high && el_Charge[0] * el_Charge[2] < 0 && iso[0]  && iso[2]){
                    el_extra = 1;
                    if(el_Charge[0] >0){
                        el_p = 0;
                        el_m = 2;
                    }
                    else{
                        el_p = 2;
                        el_m =0;
                    }
                }
                else if(m12 > Z_mass_low && m12 < Z_mass_high && el_Charge[1] * el_Charge[2] < 0 && iso[1]  && iso[2] ){
                    el_extra = 0;
                    if(el_Charge[1] >0){
                        el_p = 1;
                        el_m = 2;
                    }
                    else{
                        el_p = 2;
                        el_m =2;
                    }
                }


                if( (el_p != -1) && no_bjets && met_pt < 50){
                    nEvents++;
                    if(iso[el_extra]){
                        nPass++;
                        h_pass->Fill(abs(el_Eta[el_extra]), el_Pt[el_extra], 1);
                    }
                    h_total->Fill(abs(el_Eta[el_extra]), el_Pt[el_extra], 1);
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
