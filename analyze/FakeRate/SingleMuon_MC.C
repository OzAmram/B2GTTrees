

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <iostream>
#include <cstring>
#include <algorithm>
#include "TFile.h"
#include "TFile.h"
#include "../ScaleFactors.C"
#include "../TemplateMaker.C"

#define GEN_SIZE 300
#define MU_SIZE 100
#define JET_SIZE 20
#define MAX_SAMPLES 20

const double root2 = sqrt(2);
const char* filename("QCD_files_aug17.txt");
const TString fout_name("output_files/SingleMuon_QCD_aug22.root");

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

void SingleMuon_MC()
{

    Double_t norms[MAX_SAMPLES]; // computed normalizations to apply to each event in a sample (based on xsection and total weight)
    unsigned int nFiles = 0;
    printf("Computing normalizations for each sample \n");
    compute_norms(norms, &nFiles);
    printf("Done with normalizations \n\n\n");

    SFs runs_bcdef, runs_gh;
    BTag_readers b_reader;
    BTag_effs btag_effs;
    //separate SFs for runs BCDEF and GH
    setup_SFs(&runs_bcdef, &runs_gh, &b_reader, &btag_effs);
    printf("Retrieved Scale Factors \n\n");

    //TFile *fout = TFile::Open(fout_name, "RECREATE");
    Double_t cm_m, xF, cost_r, mu1_pt, mu2_pt, mu1_eta, mu2_eta, jet1_pt, jet2_pt, deltaC, jet1_eta, jet2_eta, gen_weight,
             jet1_cmva, jet1_csv, jet2_cmva, jet2_csv;
    Double_t bcdef_HLT_SF, bcdef_iso_SF, bcdef_id_SF, gh_HLT_SF, gh_iso_SF, gh_id_SF,
             jet1_b_weight, jet2_b_weight;
    Int_t nJets, jet1_flavour, jet2_flavour;
    Float_t met_pt;
    TLorentzVector mu_p, mu_m, cm, q1, q2;



    unsigned int nEvents=0;
    Double_t tot_HLT_weight=0;
    Double_t tot_noHLT_weight=0;

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
            f1->cd("B2GTTreeMaker");
            TDirectory *subdir = gDirectory;
            TTree *t1 = (TTree *)subdir->Get("B2GTree");


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


            t1->SetBranchAddress("met_size", &met_size);
            t1->SetBranchAddress("met_Pt", &met_pt);

            Long64_t nEntries =  t1->GetEntries();

            char out_buff[10000];
            bool print_out = false;
            printf("there are %i entries in this tree\n", nEntries);
            for (int i=0; i<nEntries; i++) {
                t1->GetEntry(i);
                if(met_size != 1) printf("WARNING: Met size not equal to 1\n");
                if(mu_size > MU_SIZE) printf("Warning: too many muons\n");
                bool good_trigger = HLT_IsoMu || HLT_IsoTkMu;
                //bool good_trigger = true;
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
                                jet2_flavour = jet_partonflavour[j];
                                jet2_b_weight = get_btag_weight(jet_Pt[j], jet_Eta[j],jet_partonflavour[j],btag_effs, b_reader);
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
                    bool no_bjets = has_no_bjets(nJets, jet1_pt, jet2_pt, jet1_cmva, jet2_cmva);
                    float iso_0 = (mu_SumChargedHadronPt[0] + max(0., mu_SumNeutralHadronPt[0] + mu_SumPhotonPt[0] - 0.5 * mu_SumPUPt[0]))/mu_Pt[0];
                    const float tight_iso = 0.15;

                    if(no_bjets && met_pt < 50){
                        mu1_pt = mu_Pt[0];
                        mu1_eta = abs(mu_Eta[0]);
                        bcdef_HLT_SF = get_HLT_SF_1mu(mu1_pt, mu1_eta, runs_bcdef.HLT_SF, runs_bcdef.HLT_MC_EFF);
                        gh_HLT_SF = get_HLT_SF_1mu(mu1_pt, mu1_eta, runs_gh.HLT_SF, runs_gh.HLT_MC_EFF);

                        //bcdef_iso_SF = get_SF(mu1_pt, mu1_eta, runs_bcdef.ISO_SF) * get_SF(mu2_pt, mu2_eta, runs_bcdef.ISO_SF);
                        bcdef_id_SF = get_SF(mu1_pt, mu1_eta, runs_bcdef.ID_SF);

                        //gh_iso_SF = get_SF(mu1_pt, mu1_eta, runs_gh.ISO_SF) * get_SF(mu2_pt, mu2_eta, runs_gh.ISO_SF);
                        gh_id_SF = get_SF(mu1_pt, mu1_eta, runs_gh.ID_SF);

                        gen_weight = evt_Gen_Weight * normalization;
                        if(good_trigger){
                            Double_t bcdef_weight = gen_weight * bcdef_HLT_SF *  bcdef_id_SF;
                            Double_t gh_weight = gen_weight * gh_HLT_SF * gh_id_SF;
                            Double_t final_weight = (bcdef_weight *bcdef_lumi + gh_weight * gh_lumi)/(bcdef_lumi + gh_lumi);
                            //printf("Gen_weight, final_weight: %.2e, %.2e \n", gen_weight, final_weight);


                            tot_HLT_weight += final_weight;
                        }
                        else{
                            Double_t bcdef_weight = gen_weight * bcdef_id_SF;
                            Double_t gh_weight = gen_weight * gh_id_SF;
                            Double_t final_weight = (bcdef_weight *bcdef_lumi + gh_weight * gh_lumi)/(bcdef_lumi + gh_lumi);
                            //printf("Gen_weight, final_weight: %.2e, %.2e \n", gen_weight, final_weight);


                            tot_noHLT_weight += final_weight;
                        }



                        nEvents++;

                    }

                }
            }
        f1->Close();
        }
        printf("moving on to next file, currently %i events \nHLT xsection is %.3e noHLT xsection is %.3e \n\n", nEvents, tot_HLT_weight, tot_noHLT_weight);
    }
    printf("Final output for %s file \n", filename);
    printf("Total xsec of selected HLT events %.3e \n", tot_HLT_weight);
    printf("Total xsec of selected noHLT events %.3e \n", tot_noHLT_weight);

    //fout->cd();
    //printf("Writing out put to %s \n", fout_name.Data());

    //fout->Close();
    return;
}
