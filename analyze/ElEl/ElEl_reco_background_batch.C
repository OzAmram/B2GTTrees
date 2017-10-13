#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <iostream>
#include <cstring>
#include <algorithm>
#include "TFile.h"
#include "../ScaleFactors.C"

#define GEN_SIZE 300
#define EL_SIZE 100
#define JET_SIZE 20
#define MAX_SAMPLES 20

const double root2 = sqrt(2);
double Ebeam = 6500.;
double Pbeam = sqrt(Ebeam*Ebeam - 0.938*0.938);

char *filename("combined_back_files_aug29.txt");
const TString fout_name("output_files/ElEl_combined_back_sep25.root");
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




void ElEl_reco_background_batch()
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

    TTree *tout= new TTree("T_data", "Tree with reco events");
    tout->SetDirectory(0);
    Double_t cm_m, xF, cost_r, el1_pt, el2_pt, el1_eta, el2_eta, jet1_pt, jet2_pt, deltaC, jet1_eta, jet2_eta, gen_weight,
             jet1_cmva, jet1_csv, jet2_cmva, jet2_csv;
    Double_t el_id_SF, el_reco_SF, jet1_b_weight, jet2_b_weight, pu_SF;
    Int_t nJets, jet1_flavour, jet2_flavour;
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
    tout->Branch("deltaC", &deltaC, "deltaC/D");
    tout->Branch("pu_SF", &pu_SF);
    tout->Branch("gen_weight", &gen_weight, "gen_weight/D");
    tout->Branch("el_id_SF", &el_id_SF);
    tout->Branch("el_reco_SF", &el_reco_SF);
    tout->Branch("jet1_b_weight", &jet1_b_weight);
    tout->Branch("jet2_b_weight", &jet2_b_weight);
    tout->Branch("nJets", &nJets, "nJets/I");
    tout->Branch("jet1_flavour", &jet1_flavour, "jet1_flavour/I");
    tout->Branch("jet2_flavour", &jet2_flavour, "jet2_flavour/I");



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

            UInt_t el_size, jet_size, met_size;

            Float_t el_Pt[EL_SIZE], el_Eta[EL_SIZE], el_Phi[EL_SIZE], el_E[EL_SIZE],
                    el_Charge[EL_SIZE];

            Int_t el_IDMedium[EL_SIZE];


            Float_t jet_Pt[JET_SIZE], jet_Eta[JET_SIZE], jet_Phi[JET_SIZE], jet_E[JET_SIZE],
                    jet_CSV[JET_SIZE], jet_CMVA[JET_SIZE], jet_partonflavour[JET_SIZE];

            Float_t evt_Gen_Weight;

            Int_t HLT_El, pu_NtrueInt;
            t1->SetBranchAddress("el_size", &el_size); //number of els in the event
            t1->SetBranchAddress("el_Pt", &el_Pt);
            t1->SetBranchAddress("el_Eta", &el_Eta);
            t1->SetBranchAddress("el_Phi", &el_Phi);
            t1->SetBranchAddress("el_E", &el_E);
            t1->SetBranchAddress("el_Charge", &el_Charge);
            t1->SetBranchAddress("el_IDMedium", &el_IDMedium);
            t1->SetBranchAddress("HLT_Ele27_WPTight_Gsf", &HLT_El);


            if(data_2016){
                t1->SetBranchAddress("jetAK4CHS_size", &jet_size);
                t1->SetBranchAddress("jetAK4CHS_Pt", &jet_Pt);
                t1->SetBranchAddress("jetAK4CHS_Eta", &jet_Eta);
                t1->SetBranchAddress("jetAK4CHS_Phi", &jet_Phi);
                t1->SetBranchAddress("jetAK4CHS_E", &jet_E);
                t1->SetBranchAddress("jetAK4CHS_CSVv2", &jet_CSV);
                t1->SetBranchAddress("jetAK4CHS_CMVAv2", &jet_CMVA);
                t1->SetBranchAddress("jetAK4CHS_PartonFlavour", &jet_partonflavour);

            }

            else{
                t1->SetBranchAddress("jetAK4Puppi_size", &jet_size);
                t1->SetBranchAddress("jetAK4Puppi_Pt", &jet_Pt);
                t1->SetBranchAddress("jetAK4Puppi_Eta", &jet_Eta);
                t1->SetBranchAddress("jetAK4Puppi_Phi", &jet_Phi);
                t1->SetBranchAddress("jetAK4Puppi_E", &jet_E);
                t1->SetBranchAddress("jetAK4Puppi_CSVv2", &jet_CSV);
                t1->SetBranchAddress("jetAK4Puppi_CMVAv2", &jet_CMVA);

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
                if(el_size > EL_SIZE) printf("WARNING: MU_SIZE TOO LARGE \n");
                if(met_size != 1) printf("WARNING: Met size not equal to 1\n");
                bool good_trigger = HLT_El;
                if(good_trigger &&
                        el_size >= 2 && ((abs(el_Charge[0] - el_Charge[1])) > 0.01) &&
                        el_IDMedium[0] && el_IDMedium[1] &&
                        el_Pt[0] > 29. &&  el_Pt[1] > 10. &&
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
                    //met and cmva cuts to reduce ttbar background
                    if (cm_m >=150.){
                        if(PRINT) sprintf(out_buff + strlen(out_buff),"Event %i \n", i);

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


                        TLorentzVector p1(0., 0., Pbeam, Ebeam);
                        TLorentzVector p2(0., 0., -Pbeam, Ebeam);




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
                           if(jet_size >=1){
                           jet1_pt = jet_Pt[0];
                           jet1_cmva = jet_CMVA[0];
                           }
                           if(jet_size >=2){
                           jet2_pt = jet_Pt[1];
                           jet2_cmva = jet_CMVA[1];
                           }
                           */



                        //get el cut SFs


                        el_id_SF = get_el_SF(el1_pt, el1_eta, el_SF.ID_SF) * get_el_SF(el2_pt, el2_eta, el_SF.ID_SF);
                        el_reco_SF = get_el_SF(el1_pt, el1_eta, el_SF.RECO_SF) * get_el_SF(el2_pt, el2_eta, el_SF.RECO_SF);


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