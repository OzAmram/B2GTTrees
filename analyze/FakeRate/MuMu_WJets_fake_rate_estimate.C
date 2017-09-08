#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <iostream>
#include <cstring>
#include <algorithm>
#include "TFile.h"
#include "../TemplateMaker.C"
#define MU_SIZE 200
#define JET_SIZE 20

const double root2 = sqrt(2);
const char* filename("SingleMuon_files_aug17.txt");
const TString fout_name("output_files/MuMu_WJets_est_sep6.root");

const bool data_2016 = true;

bool is_empty_line(const char *s) {
    while (*s != '\0') {
        if (!isspace(*s))
            return false;
        s++;
    }
    return true;
}

typedef struct{
    TH2D *noHLT;
    TH2D *HLT;
    Double_t noHLT_avg;
    Double_t HLT_avg;
} FakeRate;



void setup_fakerate(FakeRate *FR){
    TFile *f0 = TFile::Open("FakeRate/SingleMuon_data_fake_rate_aug22.root");
    TH2D *h1 = (TH2D *) gDirectory->Get("h_rate_HLT")->Clone();
    TH2D *h2 = (TH2D *) gDirectory->Get("h_rate_noHLT")->Clone();
    h1->SetDirectory(0);
    h2->SetDirectory(0);
    FR->HLT = h1;
    FR->noHLT = h2;
    FR->noHLT_avg = 0.13; 
    FR->HLT_avg = 0.91; 
}

Double_t get_fakerate_prob(Double_t pt, Double_t eta, TH2D *h){
    if (pt >= 400) pt = 380;

    TAxis* x_ax =  h->GetXaxis();
    TAxis *y_ax =  h->GetYaxis();
    int xbin = x_ax->FindBin(std::abs(eta));
    int ybin = y_ax->FindBin(pt);

    Double_t prob = h->GetBinContent(xbin, ybin);
    if(prob < 0.001 || prob > 1) printf("Warning: %.2f Rate for pt %.0f, eta %1.1f! \n", prob, pt, eta);
    //printf("Efficiency is %f \n", eff);
    return prob;
}

void MuMu_WJets_fake_rate_estimate()
{
    FakeRate FR;
    setup_fakerate(&FR);


    TTree *tout= new TTree("T_data", "Tree with reco events");
    tout->SetDirectory(0);
    Double_t cm_m, xF, cost_r, mu1_pt, mu2_pt, mu1_eta, mu2_eta, jet1_pt, jet2_pt,
             jet1_cmva, jet1_eta, jet2_cmva, jet2_eta, evt_fakerate;
    Double_t mu_fakerate;
    Int_t nJets, iso_muon;
    Bool_t double_muon_trig;
    Float_t met_pt;
    TLorentzVector mu_p, mu_m, cm, q1, q2;
    tout->Branch("m", &cm_m, "m/D");
    tout->Branch("xF", &xF, "xF/D");
    tout->Branch("cost", &cost_r, "cost/D");
    tout->Branch("mu1_pt", &mu1_pt, "mu1_pt/D");
    tout->Branch("mu2_pt", &mu2_pt, "mu2_pt/D");
    tout->Branch("mu1_eta", &mu1_eta, "mu1_eta/D");
    tout->Branch("mu2_eta", &mu2_eta, "mu2_eta/D");
    tout->Branch("mu_fakerate", &mu_fakerate);
    tout->Branch("mu_m", "TLorentzVector", &mu_m);
    tout->Branch("mu_p", "TLorentzVector", &mu_p);
    tout->Branch("jet1_pt", &jet1_pt, "jet1_pt/D");
    tout->Branch("jet1_eta", &jet1_eta, "jet1_eta/D");
    tout->Branch("jet1_CMVA", &jet1_cmva, "jet1_CMVA/D");
    tout->Branch("jet2_pt", &jet2_pt, "jet2_pt/D");
    tout->Branch("jet2_eta", &jet2_eta, "jet2_eta/D");
    tout->Branch("jet2_CMVA", &jet2_cmva, "jet2_CMVA/D");
    tout->Branch("met_pt", &met_pt, "met_Pt/F");
    tout->Branch("evt_fakerate", &evt_fakerate);
    tout->Branch("nJets", &nJets, "nJets/I");
    tout->Branch("double_muon_trig", &double_muon_trig);
    tout->Branch("iso_muon", &iso_muon);



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

        UInt_t mu_size, met_size, jet_size;
        Float_t mu_Pt[MU_SIZE], mu_Eta[MU_SIZE], mu_Phi[MU_SIZE], mu_E[MU_SIZE], 
                mu_IsTightMuon[MU_SIZE], mu_Charge[MU_SIZE];

        Float_t mu_SumChargedHadronPt[MU_SIZE], mu_SumNeutralHadronPt[MU_SIZE], mu_SumPUPt[MU_SIZE], mu_SumPhotonPt[MU_SIZE];


        Float_t jet_Pt[JET_SIZE], jet_Eta[JET_SIZE], jet_Phi[JET_SIZE], jet_E[JET_SIZE],
                jet_CSV[JET_SIZE], jet_CMVA[JET_SIZE];

        Int_t HLT_IsoMu, HLT_IsoTkMu, HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL, HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL, 
              HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ, HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_DZ;
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

        t1->SetBranchAddress("HLT_IsoMu24", &HLT_IsoMu);
        t1->SetBranchAddress("HLT_IsoTkMu24", &HLT_IsoTkMu);
        t1->SetBranchAddress("HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_DZ", &HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_DZ);
        t1->SetBranchAddress("HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL", &HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL);
        t1->SetBranchAddress("HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ", &HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ);
        t1->SetBranchAddress("HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL", &HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL);
        t1->SetBranchAddress("met_size", &met_size);
        t1->SetBranchAddress("met_Pt", &met_pt);

        unsigned int nEntries =  t1->GetEntries();
        printf("there are %i entries in this tree\n", nEntries);
        for (int i=0; i<nEntries; i++) {
            t1->GetEntry(i);
            if(met_size != 1) printf("WARNING: Met size not equal to 1\n");
            if(mu_size > MU_SIZE) printf("Warning: too many muons\n");
            bool good_trigger = HLT_IsoMu || HLT_IsoTkMu;
            if(good_trigger &&
                    mu_size >= 2 && ((abs(mu_Charge[0] - mu_Charge[1])) > 0.01) &&
                    mu_IsTightMuon[0] && mu_IsTightMuon[1] &&
                    mu_Pt[0] > 26. &&  mu_Pt[1] > 10. &&
                    abs(mu_Eta[0]) < 2.4 && abs(mu_Eta[1]) < 2.4){ 
                //only want events with 2 oppositely charged muons
                //with pt above threshold
                //See https://twiki.cern.ch/twiki/bin/viewauth/CMS/SWGuideMuonIdRun2 for iso cuts
                float iso_0 = (mu_SumChargedHadronPt[0] + max(0., mu_SumNeutralHadronPt[0] + mu_SumPhotonPt[0] - 0.5 * mu_SumPUPt[0]))/mu_Pt[0];
                float iso_1 = (mu_SumChargedHadronPt[1] + max(0., mu_SumNeutralHadronPt[1] + mu_SumPhotonPt[1] - 0.5 * mu_SumPUPt[1]))/mu_Pt[1];
                const float tight_iso = 0.15;
                const float loose_iso = 0.25;
                if(mu_Charge[0] >0){
                    mu_p.SetPtEtaPhiE(mu_Pt[0], mu_Eta[0], mu_Phi[0], mu_E[0]);
                    mu_m.SetPtEtaPhiE(mu_Pt[1], mu_Eta[1], mu_Phi[1], mu_E[1]);
                }
                else{
                    mu_m.SetPtEtaPhiE(mu_Pt[0], mu_Eta[0], mu_Phi[0], mu_E[0]);
                    mu_p.SetPtEtaPhiE(mu_Pt[1], mu_Eta[1], mu_Phi[1], mu_E[1]);
                }

                cm = mu_p + mu_m;
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
                bool one_iso = (iso_0 < tight_iso) ^ (iso_1 < tight_iso);
                if (one_iso && cm_m >=150. && no_bjets && met_pt < 50.){
                    xF = abs(2.*cm.Pz()/13000.); 

                    // compute Colins soper angle with formula
                    double mu_p_pls = (mu_p.E()+mu_p.Pz())/root2;
                    double mu_p_min = (mu_p.E()-mu_p.Pz())/root2;
                    double mu_m_pls = (mu_m.E()+mu_m.Pz())/root2;
                    double mu_m_min = (mu_m.E()-mu_m.Pz())/root2;
                    double qt2 = cm.Px()*cm.Px()+cm.Py()*cm.Py();
                    //cost_p = cos(theta)_r (reconstructed collins soper angle, sign
                    //may be 'wrong' if lepton pair direction is not the same as inital
                    //quark direction
                    double cost = 2*(mu_m_pls*mu_p_min - mu_m_min*mu_p_pls)/(cm_m*sqrt(cm_m*cm_m + qt2));
                    if(cm.Pz() < 0.) cost_r = -cost;
                    else cost_r = cost;


                    mu1_pt = mu_Pt[0];
                    mu2_pt = mu_Pt[1];
                    mu1_eta = mu_Eta[0];
                    mu2_eta = mu_Eta[1];


                    iso_muon = 0; //which muon is isolated
                    if(iso_0 < tight_iso) iso_muon = 0;
                    else if(iso_1 < tight_iso) iso_muon = 1;
                    else printf("ERROR: Neither muon iso\n");

                    double_muon_trig = HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL || HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL || HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ || HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_DZ;
                    if(double_muon_trig){
                        for(int k =2; k < mu_size; k++){
                            if(mu_Pt[k] > 8.){
                                float iso_k = (mu_SumChargedHadronPt[k] + max(0., mu_SumNeutralHadronPt[k] + mu_SumPhotonPt[k] - 0.5 * mu_SumPUPt[k]))/mu_Pt[k];
                                if (iso_k < loose_iso) double_muon_trig = false; //other muon caused trigger
                            }
                        }
                    }

                    Double_t p1, p2;
                    if(double_muon_trig){
                        mu_fakerate = get_fakerate_prob(mu_Pt[iso_muon], mu_Eta[iso_muon], FR.HLT);
                        printf("Double trig ");
                    }
                    else{
                        printf("Single Trig ");
                        mu_fakerate = get_fakerate_prob(mu_Pt[iso_muon], mu_Eta[iso_muon], FR.noHLT);
                    }
                    //Double_t prob_WJets = 7.39e-02 / (2.031e-1 + 7.39e-02);
                    //printf("prob %.3f \n", prob_WJets);
                    evt_fakerate = mu_fakerate / (1-mu_fakerate);
                    printf("evt_fakerate = %.2f \n", evt_fakerate);

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
