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
const TString fout_name("output_files/SingleElectron_data_sep22.root");

const bool data_2016 = true;

bool is_empty_line(const char *s) {
    while (*s != '\0') {
        if (!isspace(*s))
            return false;
        s++;
    }
    return true;
}

void ElEl_reco_data_batch()
{



    TTree *tout= new TTree("T_data", "Tree with reco events");
    tout->SetDirectory(0);
    Double_t cm_m, xF, cost_r, el1_pt, el2_pt, el1_eta, el2_eta, jet1_pt, jet2_pt,
             jet1_cmva, jet1_eta, jet2_cmva, jet2_eta;
    Int_t nJets;
    Float_t met_pt;
    TLorentzVector el_p, el_m, cm, q1, q2;
    tout->Branch("m", &cm_m, "m/D");
    tout->Branch("xF", &xF, "xF/D");
    tout->Branch("cost", &cost_r, "cost/D");
    tout->Branch("el1_pt", &el1_pt, "el1_pt/D");
    tout->Branch("el2_pt", &el2_pt, "el2_pt/D");
    tout->Branch("el1_eta", &el1_eta, "el1_eta/D");
    tout->Branch("el2_eta", &el2_eta, "el2_eta/D");
    tout->Branch("el_m", "TLorentzVector", &el_m);
    tout->Branch("el_p", "TLorentzVector", &el_p);
    tout->Branch("jet1_pt", &jet1_pt, "jet1_pt/D");
    tout->Branch("jet1_eta", &jet1_eta, "jet1_eta/D");
    tout->Branch("jet1_CMVA", &jet1_cmva, "jet1_CMVA/D");
    tout->Branch("jet2_pt", &jet2_pt, "jet2_pt/D");
    tout->Branch("jet2_eta", &jet2_eta, "jet2_eta/D");
    tout->Branch("jet2_CMVA", &jet2_cmva, "jet2_CMVA/D");
    tout->Branch("met_pt", &met_pt, "met_Pt/F");
    tout->Branch("nJets", &nJets, "nJets/I");



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
        if(lines[0] == '#' || lines[0] == '!' || is_empty_line(lines)) continue; // comment line
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

        Int_t el_IDMedium[EL_SIZE];

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
            if(good_trigger &&
                    el_size >= 2 && ((abs(el_Charge[0] - el_Charge[1])) > 0.01) &&
                    el_IDMedium[0] && el_IDMedium[1] &&
                    el_Pt[0] > 29. &&  el_Pt[1] > 10. &&
                    abs(el_Eta[0]) < 2.4 && abs(el_Eta[1]) < 2.4){ 
            /*
            if(good_trigger &&
                    el_size >=2 && 
                    (el_IDMedium[0] && el_Pt[0] > 26. && abs(el_Eta[0]) < 2.4) &&
                    ( (el_IDMedium[1] && el_Pt[1] > 10. && abs(el_Eta[1]) < 2.4)  ||
                      (el_size >=2 && el_IDMedium[2] && el_Pt[2] > 10. && abs(el_Eta[2]) < 2.4))){
                      */




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
                if (cm_m >=150.){
                    xF = abs(2.*cm.Pz()/13000.); 

                    // compute Colins soper angle with formula
                    double el_p_pls = (el_p.E()+el_p.Pz())/root2;
                    double el_p_min = (el_p.E()-el_p.Pz())/root2;
                    double el_m_pls = (el_m.E()+el_m.Pz())/root2;
                    double el_m_min = (el_m.E()-el_m.Pz())/root2;
                    double qt2 = cm.Px()*cm.Px()+cm.Py()*cm.Py();
                    //cost_p = cos(theta)_r (reconstructed collins soper angle, sign
                    //may be 'wrong' if lepton pair direction is not the same as inital
                    //quark direction
                    double cost = 2*(el_m_pls*el_p_min - el_m_min*el_p_pls)/(cm_m*sqrt(cm_m*cm_m + qt2));
                    if(cm.Pz() < 0.) cost_r = -cost;
                    else cost_r = cost;


                    el1_pt = el_Pt[0];
                    el2_pt = el_Pt[1];
                    el1_eta = el_Eta[0];
                    el2_eta = el_Eta[1];

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