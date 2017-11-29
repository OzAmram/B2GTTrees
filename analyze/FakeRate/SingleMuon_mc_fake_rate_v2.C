
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <iostream>
#include <cstring>
#include <algorithm>
#include "TFile.h"
#include "../ScaleFactors.C"
#include "../TemplateMaker.C"

#define GEN_SIZE 300
#define MU_SIZE 100
#define JET_SIZE 20
#define MAX_SAMPLES 20

const double root2 = sqrt(2);
char *filename("non_QCD_files_aug29.txt");
const TString fout_name("FakeRate/root_files/SingleMu_mc_fakerate_contam_v2_nov28.root");


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

void SingleMuon_mc_fake_rate_v2()
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
    setup_el_SF(&el_SF);
    //separate SFs for runs BCDEF and GH
    setup_SFs(&runs_bcdef, &runs_gh, &b_reader, &btag_effs, &pu_SFs);
    printf("Retrieved Scale Factors \n\n");

    TFile *fout = TFile::Open(fout_name, "RECREATE");
    Float_t pt_bins[] = {0, 25, 35, 50, 90, 1000};
    int n_pt_bins = 5;
    Float_t eta_bins[] = {0, 0.9, 2.4};
    int n_eta_bins = 2;

    TH2D *h_pass = new TH2D("h_pass", "Rate of passing ISO cut for single muons",  n_eta_bins, eta_bins, n_pt_bins, pt_bins);
    TH2D *h_total = new TH2D("h_total", "Total number of single muons",  n_eta_bins, eta_bins, n_pt_bins, pt_bins);
    TTree *tout= new TTree("T_data", "Tree with reco events");
    Double_t cm_m, xF, cost_r, mu1_pt, mu2_pt, mu1_eta, mu2_eta, jet1_pt, jet2_pt, deltaC, jet1_eta, jet2_eta, gen_weight,
             jet1_cmva, jet1_csv, jet2_cmva, jet2_csv;
    Double_t bcdef_HLT_SF, bcdef_iso_SF, bcdef_id_SF, gh_HLT_SF, gh_iso_SF, gh_id_SF,
             jet1_b_weight, jet2_b_weight, pu_SF;
    Int_t nJets, jet1_flavour, jet2_flavour;
    Float_t met_pt;
    TLorentzVector mu_p, mu_m, cm, q1, q2;

    Double_t mu_pt, mu_eta;
    Bool_t pass;
    tout->Branch("mu_pt", &mu_pt);
    tout->Branch("mu_eta", &mu_eta);
    tout->Branch("gen_weight", &gen_weight);
    tout->Branch("pass", &pass);



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

            f1->cd("EventCounter");
            TDirectory *subdir = gDirectory;
            TH1D *mc_pileup = (TH1D *)subdir->Get("pileup");
            mc_pileup->Scale(1./mc_pileup->Integral());
            pu_SFs.pileup_ratio->Divide(pu_SFs.data_pileup, mc_pileup);

            f1->cd("B2GTTreeMaker");
            TTree *t1 = (TTree *)gDirectory->Get("B2GTree");


            UInt_t mu_size, met_size, jet_size;
            Float_t mu_Pt[MU_SIZE], mu_Eta[MU_SIZE], mu_Phi[MU_SIZE], mu_E[MU_SIZE], 
                    mu_IsTightMuon[MU_SIZE], mu_Charge[MU_SIZE];

            Float_t mu_SumChargedHadronPt[MU_SIZE], mu_SumNeutralHadronPt[MU_SIZE], mu_SumPUPt[MU_SIZE], mu_SumPhotonPt[MU_SIZE];


            Float_t jet_Pt[JET_SIZE], jet_Eta[JET_SIZE], jet_Phi[JET_SIZE], jet_E[JET_SIZE],
                    jet_CSV[JET_SIZE], jet_CMVA[JET_SIZE], jet_partonflavour[JET_SIZE];

            Float_t evt_Gen_Weight;

            Int_t HLT_IsoMu, HLT_IsoTkMu, pu_NtrueInt;
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

            t1->SetBranchAddress("met_size", &met_size);
            t1->SetBranchAddress("met_Pt", &met_pt);

            t1->SetBranchAddress("evt_Gen_Weight", &evt_Gen_Weight);
            t1->SetBranchAddress("pu_NtrueInt",&pu_NtrueInt);

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
                    jet1_b_weight = 1;
                    jet2_b_weight = 1;
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
                    double_t Z_mass_low = 91.2 -7;
                    double_t Z_mass_high = 91.2 + 7;

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
                        //get el cut SFs
                        mu1_pt = mu_Pt[mu_p];
                        mu2_pt = mu_Pt[mu_m];
                        mu1_eta = mu_Eta[mu_p];
                        mu2_eta = mu_Eta[mu_m];

                        /*
                        if(mu1_pt > 26 || mu2_pt > 26){
                            bcdef_HLT_SF = get_HLT_SF(mu1_pt, mu1_eta, mu2_pt, mu2_eta, runs_bcdef.HLT_SF, runs_bcdef.HLT_MC_EFF);
                            gh_HLT_SF = get_HLT_SF(mu1_pt, mu1_eta, mu2_pt, mu2_eta, runs_gh.HLT_SF, runs_gh.HLT_MC_EFF);
                        }
                        else{
                            bcdef_HLT_SF = get_HLT_SF_1mu(mu_Pt[mu_extra], mu_Eta[mu_extra], runs_bcdef.HLT_SF);
                            gh_HLT_SF = get_HLT_SF_1mu(mu_Pt[mu_extra], mu_Eta[mu_extra], runs_gh.HLT_SF);
                        }
                        */
                        bcdef_HLT_SF = 1;
                        gh_HLT_SF = 1;
                        bcdef_iso_SF = get_SF(mu1_pt, mu1_eta, runs_bcdef.ISO_SF) * get_SF(mu2_pt, mu2_eta, runs_bcdef.ISO_SF);
                        bcdef_id_SF = get_SF(mu1_pt, mu1_eta, runs_bcdef.ID_SF) * get_SF(mu2_pt, mu2_eta, runs_bcdef.ID_SF);

                        gh_iso_SF = get_SF(mu1_pt, mu1_eta, runs_gh.ISO_SF) * get_SF(mu2_pt, mu2_eta, runs_gh.ISO_SF);
                        gh_id_SF = get_SF(mu1_pt, mu1_eta, runs_gh.ID_SF) * get_SF(mu2_pt, mu2_eta, runs_gh.ID_SF);

                        Double_t bcdef_weight = bcdef_HLT_SF * bcdef_iso_SF * bcdef_id_SF;
                        Double_t gh_weight = gh_HLT_SF * gh_iso_SF * gh_id_SF;

                        pu_SF = get_pileup_SF(pu_NtrueInt, pu_SFs.pileup_ratio);
                        gen_weight = pu_SF * evt_Gen_Weight * normalization * jet1_b_weight * jet2_b_weight* (bcdef_weight *bcdef_lumi + gh_weight * gh_lumi)/(bcdef_lumi + gh_lumi);
                        if(TMath::IsNaN(gen_weight)){
                            printf("NAN ENTRY!!!! \n");
                            printf(" mu1: %.0f %.2f \n mu2: %.0f %.2f \n HLT ISO ID %.0e %.0e %.0e \n", mu1_pt, mu1_eta, mu2_pt, mu2_eta, bcdef_HLT_SF, bcdef_iso_SF, bcdef_id_SF);
                            continue;
                        }
                        //printf("%.e %.e %.e %.e %.e %.e\n", gen_weight, pu_SF, evt_Gen_Weight, normalization, bcdef_weight, gh_weight);

                        nEvents++;
                        pass = iso[mu_extra] < tight_iso;
                        mu_pt = mu_Pt[mu_extra];
                        mu_eta = mu_Eta[mu_extra];
                        if(iso[mu_extra] < tight_iso){
                            h_pass->Fill(abs(mu_Eta[mu_extra]), mu_Pt[mu_extra], gen_weight);
                        }
                        h_total->Fill(abs(mu_Eta[mu_extra]), mu_Pt[mu_extra], gen_weight);

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
