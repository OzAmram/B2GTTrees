
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <iostream>
#include <cstring>
#include <algorithm>
#include "TFile.h"

#define GEN_SIZE 300
#define MU_SIZE 100
#define JET_SIZE 20
#define MAX_SAMPLES 10

const double root2 = sqrt(2);
double Ebeam = 6500.;
double Pbeam = sqrt(Ebeam*Ebeam - 0.938*0.938);

char *filename("ttbar_files_may9.txt");
const TString fout_name("SFs/BTag_efficiency_may24.root");
const double alpha = 0.05;
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





void Btag_efficiency()
{

    unsigned int nFiles = 0;

    TFile *fout = TFile::Open(fout_name, "RECREATE");
    TTree *tout= new TTree("T_data", "Tree with reco events");
    Double_t bcdef_HLT_SF, bcdef_iso_SF, bcdef_id_SF, gh_HLT_SF, gh_iso_SF, gh_id_SF;

    Double_t Eta_bins[] = {0, 0.9, 1.2, 2.1, 2.4}; 
    Int_t nEta_bins = 4;
    Double_t Pt_bins[] = {0,10,20,30,40,50,60,80,100,120,160,200,250,300,400,500};
    Int_t nPt_bins = 15;
    TH2D *b_num = new TH2D("b_num", "Efficiency for b-tagging b-jets", nPt_bins, Pt_bins, nEta_bins, Eta_bins);
    TH2D *b_denom = new TH2D("b_demon", "Total number of b-jets", nPt_bins, Pt_bins, nEta_bins, Eta_bins);
    TH2D *c_num = new TH2D("c_num", "Efficiency for b-tagging c-jets", nPt_bins, Pt_bins, nEta_bins, Eta_bins);
    TH2D *c_denom = new TH2D("c_denom", "Efficiency for b-tagging c-jets", nPt_bins, Pt_bins, nEta_bins, Eta_bins);
    TH2D *udsg_num = new TH2D("udsg_num", "Efficiency for b-tagging udsg-jets", nPt_bins, Pt_bins, nEta_bins, Eta_bins);
    TH2D *udsg_denom = new TH2D("udsg_denom", "Efficiency for b-tagging udsg-jets", nPt_bins, Pt_bins, nEta_bins, Eta_bins);

    unsigned int nJets=0;
    unsigned int nB=0;
    unsigned int nC=0;
    unsigned int nUDSG=0;


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
        }
        else {//root file



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


            UInt_t mu_size, jet_size, met_size;

            Float_t mu_Pt[MU_SIZE], mu_Eta[MU_SIZE], mu_Phi[MU_SIZE], mu_E[MU_SIZE], 
                    mu_Charge[MU_SIZE], mu_IsTightMuon[MU_SIZE];
            Float_t mu_SumChargedHadronPt[MU_SIZE], mu_SumNeutralHadronPt[MU_SIZE], mu_SumPUPt[MU_SIZE], mu_SumPhotonPt[MU_SIZE];


            Float_t jet_Pt[JET_SIZE], jet_Eta[JET_SIZE], jet_Phi[JET_SIZE], jet_E[JET_SIZE],
                    jet_CSV[JET_SIZE], jet_CMVA[JET_SIZE], jet_genpartonpt[JET_SIZE], jet_partonflavour[JET_SIZE];

            Float_t evt_Gen_Weight;
            Float_t met_pt;

            Int_t HLT_IsoMu, HLT_IsoTkMu;
            TLorentzVector mu_p, mu_m, cm, q1, q2;
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

            t1->SetBranchAddress("jetAK4CHS_size", &jet_size);
            t1->SetBranchAddress("jetAK4CHS_Pt", &jet_Pt);
            t1->SetBranchAddress("jetAK4CHS_Eta", &jet_Eta);
            t1->SetBranchAddress("jetAK4CHS_Phi", &jet_Phi);
            t1->SetBranchAddress("jetAK4CHS_E", &jet_E);
            t1->SetBranchAddress("jetAK4CHS_CSVv2", &jet_CSV);
            t1->SetBranchAddress("jetAK4CHS_CMVAv2", &jet_CMVA);
            t1->SetBranchAddress("jetAK4CHS_PartonFlavour", &jet_partonflavour);
            t1->SetBranchAddress("jetAK4CHS_GenPartonPt", &jet_genpartonpt);
            

            t1->SetBranchAddress("HLT_IsoMu24", &HLT_IsoMu);
            t1->SetBranchAddress("HLT_IsoTkMu24", &HLT_IsoTkMu);

            t1->SetBranchAddress("evt_Gen_Weight", &evt_Gen_Weight);


            t1->SetBranchAddress("met_size", &met_size);
            t1->SetBranchAddress("met_Pt", &met_pt);

            Long64_t nEntries =  t1->GetEntries();

            char out_buff[10000];
            bool print_out = false;

            int B =5;
            int C =4;

            Double_t medium_wp = 0.4432;

            for (int i=0; i<nEntries; i++) {
                t1->GetEntry(i);
                //printf("Got entry \n");
                if(mu_size > MU_SIZE) printf("WARNING: MU_SIZE TOO LARGE \n");
                
                if(met_size != 1) printf("WARNING: Met size not equal to 1\n");
                if(jet_size >= 2){
                    for(int j=0; j < jet_size; j++){
                            int flavour = std::abs(jet_partonflavour[j]);
                            if (flavour == B){

                                nB++;
                                if(jet_CMVA[j] >= medium_wp) b_num->Fill(jet_Pt[j], std::abs(jet_Eta[j]));
                                b_denom->Fill(jet_Pt[j], std::abs(jet_Eta[j]));
                            }
                            else if (flavour == C){
                                nC++;
                                if(jet_CMVA[j] >= medium_wp) c_num->Fill(jet_Pt[j], std::abs(jet_Eta[j]));
                                c_denom->Fill(jet_Pt[j], std::abs(jet_Eta[j]));
                            }
                            else if (jet_genpartonpt[j] > 8.){
                                nUDSG++;
                                if(jet_CMVA[j] >= medium_wp) udsg_num->Fill(jet_Pt[j], std::abs(jet_Eta[j]));
                                udsg_denom->Fill(jet_Pt[j], std::abs(jet_Eta[j]));
                            }
                            nJets++;
                        }

                }
                 
            }

            f1->Close();
            printf("moving on to next file, currently %i Jets \n\n", nJets);
        }
    }
    fclose(root_files);
    printf("There were %i Bs %i Cs and %i UDSGs in %i files.\n",
            nB, nC, nUDSG, nFiles);
    fout->cd();

    TH2D* b_eff = (TH2D *) b_num->Clone("b_eff");
    b_eff->Divide(b_denom);
    b_eff->Write();

    TH2D* c_eff = (TH2D *) c_num->Clone("c_eff");
    c_eff->Divide(c_denom);
    c_eff->Write();

    TH2D* udsg_eff = (TH2D *) udsg_num->Clone("udsg_eff");
    udsg_eff->Divide(udsg_denom);
    udsg_eff->Write();

    printf("Writing output to file at %s \n", fout_name.Data());

    fout->Close();

    return;
}
