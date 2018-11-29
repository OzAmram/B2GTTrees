



#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <iostream>
#include <cstring>
#include <algorithm>
#include "TFile.h"
#include "TFile.h"
#include "../TemplateMaker.C"

#define GEN_SIZE 300
#define EL_SIZE 100
#define JET_SIZE 20
#define MAX_SAMPLES 20

const double root2 = sqrt(2);
char *filename("diboson_files_aug7.txt");
const TString fout_name("output_files/SingleElectron_mc_fakerate_contam_v2_nov27.root");


bool is_empty_line(const char *s) {
    while (*s != '\0') {
        if (!isspace(*s))
            return false;
        s++;
    }
    return true;
}

bool in_Z_window(Double_t m){
    double_t Z_mass_low = 91.2 -7.;
    double_t Z_mass_high = 91.2 + 7.;
    return (m > Z_mass_low ) && (m < Z_mass_high);
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

void SingleElectron_mc_fake_rate_v2(int nJobs=1, int iJob =0)
{

    Double_t norms[MAX_SAMPLES]; // computed normalizations to apply to each event in a sample (based on xsection and total weight)
    unsigned int nFiles = 0;
    printf("Computing normalizations for each sample \n");
    compute_norms(norms, &nFiles);
    printf("Done with normalizations \n\n\n");


    mu_SFs runs_bcdef, runs_gh;
    BTag_readers b_reader;
    BTag_effs btag_effs;
    el_SFs el_SF;
    setup_el_SF(&el_SF);
    //separate SFs for runs BCDEF and GH
    printf("Retrieved Scale Factors \n\n");

    TFile *fout = TFile::Open(fout_name, "RECREATE");
    Float_t pt_bins[] = {0, 25, 35, 50, 90, 1000};
    int n_pt_bins = 5;
    Float_t eta_bins[] = {0, 0.9, 2.4};
    int n_eta_bins = 2;

    TH2D *h_pass = new TH2D("h_pass", "Rate of passing ISO cut for single electrons",  n_eta_bins, eta_bins, n_pt_bins, pt_bins);
    TH2D *h_total = new TH2D("h_total", "Total number of single electrons",  n_eta_bins, eta_bins, n_pt_bins, pt_bins);

    TTree *tout= new TTree("T_data", "Tree with reco events");
    Double_t cm_m, xF, cost_r, mu_pt, mu_eta, jet1_pt, jet2_pt, deltaC, jet1_eta, jet2_eta, gen_weight,
             jet1_cmva, jet1_csv, jet2_cmva, jet2_csv;
    Double_t bcdef_HLT_SF, bcdef_iso_SF, bcdef_id_SF, gh_HLT_SF, gh_iso_SF, gh_id_SF,
             jet1_b_weight, jet2_b_weight, pu_SF;
    Int_t nJets, jet1_flavour, jet2_flavour;
    Float_t met_pt;
    Double_t el_pt, el_eta;
    Bool_t pass;

    TLorentzVector mu_p, mu_m, cm, q1, q2;
    tout->Branch("el_pt", &el_pt);
    tout->Branch("el_eta", &el_eta);
    tout->Branch("gen_weight", &gen_weight);
    tout->Branch("pass", &pass);



    unsigned int nEvents=0;
    Double_t tot_HLT_weight=0;
    Double_t tot_noHLT_weight=0;

    Double_t normalization = 1.0;

    FILE *root_files = fopen(filename, "r");
    char lines[300];
    int count = 0;
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

        count++;
        if(count % nJobs != iJob) continue; 

            char * end;
            //remove trailing whitespace
            end = lines + strlen(lines) - 1;
            while(end > lines && isspace((char)*end)) end--;
            // Write new null terminator
            *(end+1) = 0;

            printf("Opening file: %s \n", lines);
            TFile *f1=  TFile::Open(lines);
            printf("opened \n");


            f1->cd("B2GTTreeMaker");
            TTree *t1 = (TTree *)gDirectory->Get("B2GTree");


            UInt_t el_size, jet_size, met_size;

            Float_t el_Pt[EL_SIZE], el_Eta[EL_SIZE], el_Phi[EL_SIZE], el_E[EL_SIZE],
                    el_Charge[EL_SIZE];

            Float_t el_ScaleCorr[EL_SIZE], el_ScaleCorrUp[EL_SIZE], el_ScaleCorrDown[EL_SIZE],
                el_ScaleSmearDown[EL_SIZE], el_ScaleSmearUp[EL_SIZE], el_SCEta[EL_SIZE];

            Int_t el_IDMedium[EL_SIZE], el_IDMedium_NoIso[EL_SIZE];

            Float_t jet_Pt[JET_SIZE], jet_Eta[JET_SIZE], jet_Phi[JET_SIZE], jet_E[JET_SIZE],
                    jet_CSV[JET_SIZE], jet_CMVA[JET_SIZE], jet_partonflavour[JET_SIZE];

            Float_t evt_Gen_Weight;
            Int_t pu_NtrueInt;

            Int_t HLT_El;
            t1->SetBranchAddress("el_size", &el_size); //number of els in the event
            t1->SetBranchAddress("el_Pt", &el_Pt);
            t1->SetBranchAddress("el_Eta", &el_Eta);
            t1->SetBranchAddress("el_Phi", &el_Phi);
            t1->SetBranchAddress("el_E", &el_E);
            t1->SetBranchAddress("el_Charge", &el_Charge);
            t1->SetBranchAddress("el_SCEta", &el_SCEta);
            t1->SetBranchAddress("el_ScaleCorr", &el_ScaleCorr);
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

            t1->SetBranchAddress("evt_Gen_Weight", &evt_Gen_Weight);
            t1->SetBranchAddress("pu_NtrueInt",&pu_NtrueInt);


            t1->SetBranchAddress("met_MuCleanOnly_size", &met_size);
            t1->SetBranchAddress("met_MuCleanOnly_Pt", &met_pt);

            Long64_t nEntries =  t1->GetEntries();

            char out_buff[10000];
            bool print_out = false;
            printf("there are %i entries in this tree\n", nEntries);
            for (int i=0; i<nEntries; i++) {
                t1->GetEntry(i);
                if(el_size > EL_SIZE) printf("WARNING: MU_SIZE TOO LARGE \n");
                if(met_size != 1) printf("WARNING: Met size not equal to 1\n");
                bool good_trigger = HLT_El;
                if( good_trigger && el_size >= 3 && 
                        el_Pt[0] * el_ScaleCorr[0] > 29. && el_IDMedium_NoIso[0] && goodElEta(el_SCEta[0])
                        && el_Pt[1] * el_ScaleCorr[1] > 15. && el_IDMedium_NoIso[1] && goodElEta(el_SCEta[1])
                        && el_Pt[2] * el_ScaleCorr[2] > 15. && el_IDMedium_NoIso[2] && goodElEta(el_SCEta[1])){
                    //Want events with 3 elons, 2 from Z and 1 extra
                    //See https://twiki.cern.ch/twiki/bin/viewauth/CMS/SWGuideelonIdRun2 for iso cuts

                    //get jets
                    nJets =0;
                    jet1_b_weight = 1;
                    jet2_b_weight = 1;
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

                    TLorentzVector el0, el1, el2;
                    el0.SetPtEtaPhiE(el_ScaleCorr[0] * el_Pt[0], el_Eta[0], el_Phi[0], el_ScaleCorr[0] *el_E[0]);
                    el1.SetPtEtaPhiE(el_ScaleCorr[1] * el_Pt[1], el_Eta[1], el_Phi[1], el_ScaleCorr[1] *el_E[1]);
                    el2.SetPtEtaPhiE(el_ScaleCorr[2] * el_Pt[2], el_Eta[2], el_Phi[2], el_ScaleCorr[2] *el_E[2]);

                    //el+ and el- from Z, extra elon
                    int el_p, el_m, el_extra;
                    el_p = -1;
                    el_m = -1;
                    el_extra = -1;

                    Double_t m01 = (el0 + el1).M();
                    Double_t m02 = (el0 + el2).M();
                    Double_t m12 = (el1 + el2).M();

                    bool m01_in_Z = in_Z_window(m01);
                    bool m02_in_Z = in_Z_window(m02);
                    bool m12_in_Z = in_Z_window(m12);

                    Int_t iso[3];
                    iso[0] = el_IDMedium[0]; 
                    iso[1] = el_IDMedium[1];
                    iso[2] = el_IDMedium[2];

                    if(m01_in_Z && !m02_in_Z && !m12_in_Z && el_Charge[0] * el_Charge[1] < 0 && iso[0] && iso[1]){
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
                    else if(!m01_in_Z && m02_in_Z && !m12_in_Z && el_Charge[0] * el_Charge[2] < 0 && iso[0]  && iso[2]){
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
                    else if(!m01_in_Z && !m02_in_Z && m12_in_Z && el_Charge[1] * el_Charge[2] < 0 && iso[1]  && iso[2] ){
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
                    



                    if( (el_p != -1) && no_bjets && met_pt < 25.){

                        //get el cut SFs


                        Double_t el_id_SF = get_el_SF(el_Pt[el_p], el_Eta[el_p], el_SF.ID_SF) * get_el_SF(el_Pt[el_m], el_Eta[el_m], el_SF.ID_SF);
                        Double_t el_reco_SF = get_el_SF(el_Pt[el_p], el_Eta[el_p], el_SF.RECO_SF) * get_el_SF(el_Pt[el_m], el_Eta[el_m], el_SF.RECO_SF);;


                        gen_weight = evt_Gen_Weight * normalization * el_id_SF * el_reco_SF;

                        nEvents++;
                        el_pt = el_Pt[el_extra];
                        el_eta = abs(el_Eta[el_extra]);
                        pass = iso[el_extra];
                        if(iso[el_extra]){
                            h_pass->Fill(abs(el_Eta[el_extra]), el_Pt[el_extra], gen_weight);
                        }
                        h_total->Fill(abs(el_Eta[el_extra]), el_Pt[el_extra], gen_weight);
                        tout->Fill();






                    }

                }
            }
            f1->Close();
        }
        printf("moving on to next file, currently %i events \n", nEvents);
    }
    printf("Final output for %s file \n. %i events", filename, nEvents);
    fout->cd();
    printf("Writing out put to %s \n", fout_name.Data());
    tout->Write();
    h_pass->Write();
    h_total->Write();


    fout->Close();
    //
    return;
}
