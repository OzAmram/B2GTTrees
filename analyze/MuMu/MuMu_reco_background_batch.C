#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <iostream>
#include <cstring>
#include <algorithm>
#include "TFile.h"
#include "../ScaleFactors.C"

#define GEN_SIZE 300
#define MU_SIZE 100
#define JET_SIZE 20
#define MAX_SAMPLES 20

const double root2 = sqrt(2);
double Ebeam = 6500.;
double Pbeam = sqrt(Ebeam*Ebeam - 0.938*0.938);

char *filename("combined_back_files_aug29.txt");
const TString fout_name("output_files/MuMu_combined_back_sep13.root");
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




void MuMu_reco_background_batch()
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
    printf("Retrieved Scale Factors \n\n");

    TTree *tout= new TTree("T_data", "Tree with reco events");
    tout->SetDirectory(0);
    Double_t cm_m, xF, cost_r, mu1_pt, mu2_pt, mu1_eta, mu2_eta, jet1_pt, jet2_pt, deltaC, jet1_eta, jet2_eta, gen_weight,
             jet1_cmva, jet1_csv, jet2_cmva, jet2_csv;
    Double_t bcdef_HLT_SF, bcdef_iso_SF, bcdef_id_SF, gh_HLT_SF, gh_iso_SF, gh_id_SF,
             jet1_b_weight, jet2_b_weight, pu_SF;
    Int_t nJets, jet1_flavour, jet2_flavour, pu_NtrueInt;
    Float_t met_pt;
    TLorentzVector mu_p, mu_m, cm, q1, q2;
    tout->Branch("m", &cm_m, "m/D");
    tout->Branch("xF", &xF, "xF/D");
    tout->Branch("cost", &cost_r, "cost/D");
    tout->Branch("mu1_pt", &mu1_pt, "mu1_pt/D");
    tout->Branch("mu2_pt", &mu2_pt, "mu2_pt/D");
    tout->Branch("mu1_eta", &mu1_eta, "mu1_eta/D");
    tout->Branch("mu2_eta", &mu2_eta, "mu2_eta/D");
    tout->Branch("mu_m", "TLorentzVector", &mu_m);
    tout->Branch("mu_p", "TLorentzVector", &mu_p);
    tout->Branch("jet1_pt", &jet1_pt, "jet1_pt/D");
    tout->Branch("jet1_eta", &jet1_eta, "jet1_eta/D");
    tout->Branch("jet1_CMVA", &jet1_cmva, "jet1_CMVA/D");
    tout->Branch("jet2_pt", &jet2_pt, "jet2_pt/D");
    tout->Branch("jet2_eta", &jet2_eta, "jet2_eta/D");
    tout->Branch("jet2_CMVA", &jet2_cmva, "jet2_CMVA/D");
    tout->Branch("met_pt", &met_pt, "met_Pt/F");
    tout->Branch("deltaC", &deltaC, "deltaC/D");
    tout->Branch("gen_weight", &gen_weight, "gen_weight/D");
    tout->Branch("pu_SF", &pu_SF);
    tout->Branch("bcdef_HLT_SF", &bcdef_HLT_SF);
    tout->Branch("bcdef_iso_SF", &bcdef_iso_SF);
    tout->Branch("bcdef_id_SF", &bcdef_id_SF);
    tout->Branch("gh_HLT_SF", &gh_HLT_SF);
    tout->Branch("gh_iso_SF", &gh_iso_SF);
    tout->Branch("gh_id_SF", &gh_id_SF);
    tout->Branch("jet1_b_weight", &jet1_b_weight);
    tout->Branch("jet2_b_weight", &jet2_b_weight);
    tout->Branch("nJets", &nJets, "nJets/I");
    tout->Branch("jet1_flavour", &jet1_flavour, "jet1_flavour/I");
    tout->Branch("jet2_flavour", &jet2_flavour, "jet2_flavour/I");
    tout->Branch("pu_NtrueInt", &pu_NtrueInt);



    unsigned int nEvents=0;

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


            UInt_t mu_size, jet_size, met_size;

            Float_t mu_Pt[MU_SIZE], mu_Eta[MU_SIZE], mu_Phi[MU_SIZE], mu_E[MU_SIZE], 
                    mu_Charge[MU_SIZE], mu_IsTightMuon[MU_SIZE];
            Float_t mu_SumChargedHadronPt[MU_SIZE], mu_SumNeutralHadronPt[MU_SIZE], mu_SumPUPt[MU_SIZE], mu_SumPhotonPt[MU_SIZE];


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

            if(data_2016){
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
            t1->SetBranchAddress("evt_Gen_Weight", &evt_Gen_Weight);

            t1->SetBranchAddress("pu_NtrueInt",&pu_NtrueInt);

            t1->SetBranchAddress("met_size", &met_size);
            t1->SetBranchAddress("met_Pt", &met_pt);

            Long64_t nEntries =  t1->GetEntries();

            char out_buff[10000];
            bool print_out = false;

            for (int i=0; i<nEntries; i++) {
                t1->GetEntry(i);
                if(mu_size > MU_SIZE) printf("WARNING: MU_SIZE TOO LARGE \n");
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
                    //only want events with 2 oppositely charged muons
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
                    //met and cmva cuts to reduce ttbar background
                    if (iso_0 < tight_iso && iso_1 < tight_iso && cm_m >=150.){
                        if(PRINT) sprintf(out_buff + strlen(out_buff),"Event %i \n", i);

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


                        TLorentzVector p1(0., 0., Pbeam, Ebeam);
                        TLorentzVector p2(0., 0., -Pbeam, Ebeam);




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

                        /*
                           if(jet_size >=2) nJets = 2;
                           else nJets = jet_size;
                           if(jet_size >=1){
                           jet1_pt = jet_Pt[0];
                           jet1_cmva = jet_CMVA[0];
                           }
                           if(jet_size >=2){
                           jet2_pt = jet_Pt[1];
                           jet2_cmva = jet_CMVA[1];
                           }
                           */



                        //get muon cut SFs

                        bcdef_HLT_SF = get_HLT_SF(mu1_pt, mu1_eta, mu2_pt, mu2_eta, runs_bcdef.HLT_SF, runs_bcdef.HLT_MC_EFF);
                        gh_HLT_SF = get_HLT_SF(mu1_pt, mu1_eta, mu2_pt, mu2_eta, runs_gh.HLT_SF, runs_gh.HLT_MC_EFF);

                        bcdef_iso_SF = get_SF(mu1_pt, mu1_eta, runs_bcdef.ISO_SF) * get_SF(mu2_pt, mu2_eta, runs_bcdef.ISO_SF);
                        bcdef_id_SF = get_SF(mu1_pt, mu1_eta, runs_bcdef.ID_SF) * get_SF(mu2_pt, mu2_eta, runs_bcdef.ID_SF);

                        gh_iso_SF = get_SF(mu1_pt, mu1_eta, runs_gh.ISO_SF) * get_SF(mu2_pt, mu2_eta, runs_gh.ISO_SF);
                        gh_id_SF = get_SF(mu1_pt, mu1_eta, runs_gh.ID_SF) * get_SF(mu2_pt, mu2_eta, runs_gh.ID_SF);

                        pu_SF = get_pileup_SF(pu_NtrueInt, pu_SFs.pileup_ratio);


                        tout->Fill();
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
            printf("moving on to next file, currently %i events \n\n", nEvents);
        }
    }
    fclose(root_files);
    printf("There were %i ttbar events in %i files.\n",
            nEvents, nFiles);

    TFile *fout = TFile::Open(fout_name, "RECREATE");
    fout->cd();


    tout->Write();

    printf("Writing output to file at %s \n", fout_name.Data());

    fout->Close();

    return;
}
