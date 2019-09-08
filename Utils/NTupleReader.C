#include "NTupleReader.h"

void compute_norms(FILE *root_files, Double_t *norms, unsigned int *nFiles){
    Double_t sample_weight = 0;
    Double_t sample_xsec = 0;
    unsigned int sample_i=0;

    char lines[300];
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
    rewind(root_files);

}

NTupleReader::NTupleReader(const char *fin_name, const char *fout_name, bool is_data_){

    is_data = is_data_;
    root_files = fopen(fin_name , "r");
    if(root_files == NULL){
        printf("Problem opening file: %s. \nExiting. \n", fin_name);
        exit(1);
    }
    if(!is_data){
        nFiles = 0;
        printf("Computing normalizations for each sample \n");
        compute_norms(root_files, norms, &nFiles);
        printf("Done with normalizations \n\n\n");
    }
    fout = TFile::Open(fout_name, "RECREATE");
}

void NTupleReader::setupSFs(){
    do_SFs = true;

    printf("getting pu SFs \n");

    setup_pileup_systematic(&pu_sys, year);

    if(do_muons || do_emu){
        printf("getting muon SFs \n");
        setup_mu_SFs(&era1, &era2, year);
    }



    if(do_electrons || do_emu){
        printf("getting electron SFs \n");
        setup_el_SF(&el_SF, year);
    }
    printf("Retrieved Scale Factors \n\n");
}

bool NTupleReader::getNextFile(){
    if(fileCount != 0) fin->Close();
    char lines[300];
    while(fgets(lines, 300, root_files)){
        if(lines[0] == '#' || is_empty_line(lines)) continue; //comment line
        else if(lines[0] == '!'){//sample header
            int sample_idx;
            float xsec;
            if(!is_data){
                int nparams = sscanf(lines, "! idx = %i xsec = %f \n", &sample_idx, &xsec);
                if(nparams < 2 || sample_idx >= MAX_SAMPLES){
                    printf("ERROR: Unable to parse sample header. Exiting");
                    exit(EXIT_FAILURE);
                }
                normalization = norms[sample_idx];
                printf("Moving on to sample %i which has normalization %e \n", sample_idx, normalization);
            }
        }
        else {//root file
            if(normalization <= 0. && !is_data) printf("WARNING NORM IS ZERO FOR THIS SAMPLE!! This is bad! \n");
            fileCount++;
            if(fileCount % nJobs != iJob) continue;



            char * end;
            //remove trailing whitespace
            end = lines + strlen(lines) - 1;
            while(end > lines && isspace((char)*end)) end--;
            // Write new null terminator
            *(end+1) = 0;

            printf("Opening file: %s   ", lines);
            fin=  TFile::Open(lines);

            if(!is_data && do_SFs){
                fin->cd("EventCounter");
                TDirectory *subdir = gDirectory;
                TH1D *mc_pileup = (TH1D *)subdir->Get("pileup");
                mc_pileup->Scale(1./mc_pileup->Integral());


                pu_sys.ratio_pileup_nom->Divide(pu_sys.data_pileup_nom, mc_pileup);
                pu_sys.ratio_pileup_up->Divide(pu_sys.data_pileup_up, mc_pileup);
                pu_sys.ratio_pileup_down->Divide(pu_sys.data_pileup_down, mc_pileup);

            }


            fin->cd("B2GTTreeMaker");
            tin = (TTree *)gDirectory->Get("B2GTree");

            tin->SetBranchAddress("jetAK4CHS_size", &jet_size);
            tin->SetBranchAddress("jetAK4CHS_Pt", &jet_Pt);
            tin->SetBranchAddress("jetAK4CHS_Eta", &jet_Eta);
            tin->SetBranchAddress("jetAK4CHS_Phi", &jet_Phi);
            tin->SetBranchAddress("jetAK4CHS_E", &jet_E);
            tin->SetBranchAddress("jetAK4CHS_HadronFlavour", &jet_hadronflavour);
            tin->SetBranchAddress("jetAK4CHS_GenJetPt", &jet_genPt);

            if(year == 2016) {
                //tin->SetBranchAddress("jetAK4CHS_CSVv2", &jet_CSV);
                tin->SetBranchAddress("jetAK4CHS_CMVAv2", &jet_btag);
            }
            else{
                tin->SetBranchAddress("jetAK4CHS_DeepCSV", &jet_btag);
            }


            tin->SetBranchAddress("met_MuCleanOnly_size", &met_size);
            tin->SetBranchAddress("met_MuCleanOnly_Pt", &met_pt);

            if(do_muons || do_emu){

                tin->SetBranchAddress("mu_size", &mu_size); //number of muons in the event
                tin->SetBranchAddress("mu_Pt", &mu_Pt);
                tin->SetBranchAddress("mu_Eta", &mu_Eta);
                tin->SetBranchAddress("mu_Phi", &mu_Phi);
                tin->SetBranchAddress("mu_E", &mu_E);
                tin->SetBranchAddress("mu_Charge", &mu_Charge);

                tin->SetBranchAddress("mu_IsTightMuon", &mu_IsTightMuon);
                tin->SetBranchAddress("mu_NumberTrackerLayers", &mu_NumberTrackerLayers);
                tin->SetBranchAddress("mu_Iso04", &mu_PFIso);
                //tin->SetBranchAddress("mu_IsHighPtMuon", &mu_IsHighPtMuon);
                //tin->SetBranchAddress("mu_TrackerIso", &mu_TrackerIso);
                //tin->SetBranchAddress("mu_SumChargedHadronPt", &mu_SumChargedHadronPt);
                //tin->SetBranchAddress("mu_SumNeutralHadronPt", &mu_SumNeutralHadronPt);
                //tin->SetBranchAddress("mu_SumPUPt", &mu_SumPUPt);
                //tin->SetBranchAddress("mu_SumPhotonPt", &mu_SumPhotonPt);
                //
                if(year == 2016){
                    tin->SetBranchAddress("HLT_IsoMu24", &HLT_IsoMu24);
                    tin->SetBranchAddress("HLT_IsoTkMu24", &HLT_IsoTkMu24);
                }
                else if(year == 2017){
                    tin->SetBranchAddress("HLT_IsoMu27", &HLT_IsoMu27);
                    tin->SetBranchAddress("HLT_IsoTkMu27", &HLT_IsoTkMu27);
                }
                else if(year == 2018){
                    tin->SetBranchAddress("HLT_IsoMu24", &HLT_IsoMu24);
                }

            }
            if(do_electrons || do_emu){

                tin->SetBranchAddress("el_size", &el_size); //number of els in the event
                tin->SetBranchAddress("el_Pt", &el_Pt);
                tin->SetBranchAddress("el_Eta", &el_Eta);
                tin->SetBranchAddress("el_Phi", &el_Phi);
                tin->SetBranchAddress("el_E", &el_E);
                tin->SetBranchAddress("el_Charge", &el_Charge);
                tin->SetBranchAddress("el_IDMedium", &el_IDMedium);
                tin->SetBranchAddress("el_SCEta", &el_SCEta);
                tin->SetBranchAddress("el_ScaleCorr", &el_ScaleCorr);
                tin->SetBranchAddress("el_ScaleCorrStatUp", &el_ScaleCorrStatUp);
                tin->SetBranchAddress("el_ScaleCorrStatDown", &el_ScaleCorrStatDown);
                tin->SetBranchAddress("el_ScaleCorrSystUp", &el_ScaleCorrSystUp);
                tin->SetBranchAddress("el_ScaleCorrSystDown", &el_ScaleCorrSystDown);
                tin->SetBranchAddress("el_ScaleCorrGainUp", &el_ScaleCorrGainUp);
                tin->SetBranchAddress("el_ScaleCorrGainDown", &el_ScaleCorrGainDown);
                tin->SetBranchAddress("el_ScaleSmearUp", &el_ScaleSmearUp);
                tin->SetBranchAddress("el_ScaleSmearDown", &el_ScaleSmearDown);

                if(year == 2016){
                    tin->SetBranchAddress("HLT_Ele27_WPTight_Gsf", &HLT_El27);
                    tin->SetBranchAddress("el_IDMedium_NoIso", &el_IDMedium_NoIso);
                }
                else{
                    tin->SetBranchAddress("el_IDMediumNoIso", &el_IDMedium_NoIso);
                    tin->SetBranchAddress("HLT_Ele32_WPTight_Gsf", &HLT_El32);
                    tin->SetBranchAddress("HLT_Ele35_WPTight_Gsf", &HLT_El35);
                    tin->SetBranchAddress("HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL", &HLT_doubleEl23);
                }

            }


            if(!is_data){
                tin->SetBranchAddress("gen_size", &gen_size); //number of muons in the event
                tin->SetBranchAddress("gen_Pt", &gen_Pt);
                tin->SetBranchAddress("gen_Eta", &gen_Eta);
                tin->SetBranchAddress("gen_Phi", &gen_Phi);
                tin->SetBranchAddress("gen_E", &gen_E);
                tin->SetBranchAddress("gen_ID", &gen_id);
                tin->SetBranchAddress("gen_Status", &gen_status);
                tin->SetBranchAddress("evt_Gen_Weight", &evt_Gen_Weight);
                tin->SetBranchAddress("pu_NtrueInt",&pu_NtrueInt);

                tin->SetBranchAddress("scale_Weights", &scale_Weights);
                tin->SetBranchAddress("pdf_Weights", &pdf_weights);
                tin->SetBranchAddress("alphas_Weights", &alpha_weights);

                tin->SetBranchAddress("alphas_size", &alphas_size);
                tin->SetBranchAddress("scale_size", &scale_size);
                tin->SetBranchAddress("pdf_size", &pdf_size);

                tin->SetBranchAddress("gen_Mom0ID", &gen_Mom0ID);
                tin->SetBranchAddress("gen_Mom1ID", &gen_Mom1ID);

                tin->SetBranchAddress("gen_Dau0ID", &gen_Dau0ID);
                tin->SetBranchAddress("gen_Dau1ID", &gen_Dau1ID);
            }



            tin_nEntries =  tin->GetEntries();
            printf("done \n");
            return true;
        }
    }
    return false;
}

void NTupleReader::setupOutputTree(char treeName[100]){
    fout->cd();
    int idx = nOutTrees;
    nOutTrees++;
    outTrees[idx] = new TTree(treeName, "");

    outTrees[idx]->Branch("year", &year, "year/I");
    outTrees[idx]->Branch("m", &cm_m, "m/D");
    outTrees[idx]->Branch("xF", &xF, "xF/D");
    outTrees[idx]->Branch("cost", &cost_r, "cost/D");
    outTrees[idx]->Branch("cost_st", &cost_st, "cost_st/D");
    outTrees[idx]->Branch("jet1_pt", &jet1_pt, "jet1_pt/D");
    outTrees[idx]->Branch("jet1_eta", &jet1_eta, "jet1_eta/D");
    outTrees[idx]->Branch("jet2_pt", &jet2_pt, "jet2_pt/D");
    outTrees[idx]->Branch("jet2_eta", &jet2_eta, "jet2_eta/D");
    outTrees[idx]->Branch("nJets", &nJets, "nJets/I");
    outTrees[idx]->Branch("met_pt", &met_pt, "met_Pt/F");


    outTrees[idx]->Branch("has_nobjets", &has_nobjets, "has_nobjets/I");
    outTrees[idx]->Branch("jet1_btag", &jet1_btag, "jet1_btag/D");
    outTrees[idx]->Branch("jet2_btag", &jet2_btag, "jet2_btag/D");
    /*
    if(year != 2016){
        outTrees[idx]->Branch("jet1_btag", &jet1_btag, "jet1_btag/D");
        outTrees[idx]->Branch("jet2_btag", &jet2_btag, "jet2_btag/D");
    }
    else{
        outTrees[idx]->Branch("jet1_CMVA", &jet1_btag, "jet1_btag/D");
        outTrees[idx]->Branch("jet2_CMVA", &jet2_btag, "jet2_btag/D");
    }
    */


    if(do_muons){
        outTrees[idx]->Branch("mu1_eta", &mu1_eta, "mu1_eta/D");
        outTrees[idx]->Branch("mu2_eta", &mu2_eta, "mu2_eta/D");
        outTrees[idx]->Branch("mu_m", "TLorentzVector", &mu_m);
        outTrees[idx]->Branch("mu_p", "TLorentzVector", &mu_p);
        outTrees[idx]->Branch("mu1_pt", &mu1_pt, "mu1_pt/D");
        outTrees[idx]->Branch("mu2_pt", &mu2_pt, "mu2_pt/D");
        outTrees[idx]->Branch("mu1_charge", &mu1_charge, "mu1_charge/F");
        outTrees[idx]->Branch("mu2_charge", &mu2_charge, "mu2_charge/F");


        if(do_RC){
            outTrees[idx]->Branch("mu_p_SF", &mu_p_SF, "mu_p_SF/D");
            outTrees[idx]->Branch("mu_m_SF", &mu_m_SF, "mu_m_SF/D");
            outTrees[idx]->Branch("mu_p_SF_up", &mu_p_SF_up, "mu_p_SF_up/D");
            outTrees[idx]->Branch("mu_m_SF_up", &mu_m_SF_up, "mu_m_SF_up/D");
            outTrees[idx]->Branch("mu_p_SF_down", &mu_p_SF_down, "mu_p_SF_down/D");
            outTrees[idx]->Branch("mu_m_SF_down", &mu_m_SF_down, "mu_m_SF_down/D");
        }
    }
    if(do_electrons){
        outTrees[idx]->Branch("el1_pt", &el1_pt, "el1_pt/D");
        outTrees[idx]->Branch("el2_pt", &el2_pt, "el2_pt/D");
        outTrees[idx]->Branch("el1_eta", &el1_eta, "el1_eta/D");
        outTrees[idx]->Branch("el2_eta", &el2_eta, "el2_eta/D");
        outTrees[idx]->Branch("el_m", "TLorentzVector", &el_m);
        outTrees[idx]->Branch("el_p", "TLorentzVector", &el_p);
        outTrees[idx]->Branch("el1_charge", &el1_charge, "el1_charge/F");
        outTrees[idx]->Branch("el2_charge", &el2_charge, "el2_charge/F");
    }
    if(do_emu){
        outTrees[idx]->Branch("el", "TLorentzVector", &el);
        outTrees[idx]->Branch("mu", "TLorentzVector", &mu);
        outTrees[idx]->Branch("mu1_eta", &mu1_eta, "mu1_eta/D");
        outTrees[idx]->Branch("mu1_pt", &mu1_pt, "mu1_pt/D");
        outTrees[idx]->Branch("mu1_charge", &mu1_charge, "mu1_charge/F");
        outTrees[idx]->Branch("el1_pt", &el1_pt, "el1_pt/D");
        outTrees[idx]->Branch("el1_eta", &el1_eta, "el1_eta/D");
        outTrees[idx]->Branch("el1_charge", &el1_charge, "el1_charge/F");
    }





    if(!is_data){

        outTrees[idx]->Branch("pu_SF", &pu_SF);
        outTrees[idx]->Branch("pu_SF_up", &pu_SF_up);
        outTrees[idx]->Branch("pu_SF_down", &pu_SF_down);
        outTrees[idx]->Branch("gen_weight", &gen_weight, "gen_weight/D");
        outTrees[idx]->Branch("gen_m", &gen_m, "gen_m/D");
        outTrees[idx]->Branch("mu_R_up", &mu_R_up);
        outTrees[idx]->Branch("mu_R_down", &mu_R_down);
        outTrees[idx]->Branch("mu_F_up", &mu_F_up);
        outTrees[idx]->Branch("mu_F_down", &mu_F_down);
        outTrees[idx]->Branch("mu_RF_up", &mu_RF_up);
        outTrees[idx]->Branch("mu_RF_down", &mu_RF_down);
        outTrees[idx]->Branch("alpha_down", &alpha_down);
        outTrees[idx]->Branch("alpha_up", &alpha_up);
        outTrees[idx]->Branch("pdf_weights", &pdf_weights, "pdf_weights[60]/F");
        outTrees[idx]->Branch("jet1_flavour", &jet1_flavour, "jet1_flavour/I");
        outTrees[idx]->Branch("jet2_flavour", &jet2_flavour, "jet2_flavour/I");
        outTrees[idx]->Branch("is_tau_event", &is_tau_event);
        outTrees[idx]->Branch("pu_NtrueInt", &pu_NtrueInt);

        if(do_muons || do_emu){
            outTrees[idx]->Branch("gen_mu_m", "TLorentzVector", &gen_mu_m_vec);
            outTrees[idx]->Branch("gen_mu_p", "TLorentzVector", &gen_mu_p_vec);
            outTrees[idx]->Branch("era1_HLT_SF", &era1_HLT_SF);
            outTrees[idx]->Branch("era1_iso_SF", &era1_iso_SF);
            outTrees[idx]->Branch("era1_id_SF", &era1_id_SF);
            outTrees[idx]->Branch("era2_HLT_SF", &era2_HLT_SF);
            outTrees[idx]->Branch("era2_iso_SF", &era2_iso_SF);
            outTrees[idx]->Branch("era2_id_SF", &era2_id_SF);

        }

        if(do_electrons || do_emu){

            outTrees[idx]->Branch("elp_scale_stat_up", &elp_scale_stat_up);
            outTrees[idx]->Branch("elp_scale_stat_down", &elp_scale_stat_down);
            outTrees[idx]->Branch("elp_scale_gain_up", &elp_scale_gain_up);
            outTrees[idx]->Branch("elp_scale_gain_down", &elp_scale_gain_down);
            outTrees[idx]->Branch("elp_scale_syst_up", &elp_scale_syst_up);
            outTrees[idx]->Branch("elp_scale_syst_down", &elp_scale_syst_down);
            outTrees[idx]->Branch("elp_smear_up", &elp_smear_up);
            outTrees[idx]->Branch("elp_smear_down", &elp_smear_down);

            outTrees[idx]->Branch("elm_scale_stat_up", &elm_scale_stat_up);
            outTrees[idx]->Branch("elm_scale_stat_down", &elm_scale_stat_down);
            outTrees[idx]->Branch("elm_scale_gain_up", &elm_scale_gain_up);
            outTrees[idx]->Branch("elm_scale_gain_down", &elm_scale_gain_down);
            outTrees[idx]->Branch("elm_scale_syst_up", &elm_scale_syst_up);
            outTrees[idx]->Branch("elm_scale_syst_down", &elm_scale_syst_down);
            outTrees[idx]->Branch("elm_smear_up", &elm_smear_up);
            outTrees[idx]->Branch("elm_smear_down", &elm_smear_down);

            outTrees[idx]->Branch("gen_el_m", "TLorentzVector", &gen_el_m_vec);
            outTrees[idx]->Branch("gen_el_p", "TLorentzVector", &gen_el_p_vec);
            outTrees[idx]->Branch("el_id_SF", &el_id_SF);
            outTrees[idx]->Branch("el_reco_SF", &el_reco_SF);
            outTrees[idx]->Branch("el_HLT_SF", &el_HLT_SF);
        }
    }
}

void NTupleReader::setupRC(){
    do_RC = true;

    printf("Getting RC  \n");
    char path[400], env_var[400];

    if (year == 2016){
        sprintf(path, "%s/src/Analysis/B2GTTrees/Utils/roccor_Run2_v3/RoccoR2016.txt", getenv("CMSSW_BASE"));
    }
    else if(year ==2017){
        sprintf(path, "%s/src/Analysis/B2GTTrees/Utils/roccor_Run2_v3/RoccoR2017.txt", getenv("CMSSW_BASE"));
    }
    else if(year == 2018){
        sprintf(path, "%s/src/Analysis/B2GTTrees/Utils/roccor_Run2_v3/RoccoR2018.txt", getenv("CMSSW_BASE"));
        //sprintf(path, "%s/src/Analysis/B2GTTrees/Utils/roccor_Run2_v3/RoccoR2017.txt", getenv("CMSSW_BASE"));
    }

   
    printf("path %s \n", path);
    rc =  RoccoR(path);
    rand = new TRandom3();
}


void NTupleReader::getEvent(int i){
    tin->GetEntry(i);
    if(year == 2016) bjet_med_cut = 0.6321;
    if(year == 2017) bjet_med_cut = 0.4941;
    if(year == 2018) bjet_med_cut = 0.4184;
    event_idx = i;
    if(mu_size > MU_SIZE || el_size > EL_SIZE ||  gen_size >GEN_SIZE) printf("WARNING: MU_SIZE EL_SIZE OR GEN_SIZE TOO LARGE \n");
    if(met_size != 1) printf("WARNING: Met size not equal to 1\n");
    if(do_muons){
        opp_sign = good_trigger = dimuon_id = mu_iso0 = mu_iso1 = false;
        if(mu_size >= 2){

            opp_sign = ((abs(mu_Charge[0] - mu_Charge[1])) > 0.01);
            good_sign = opp_sign ^ do_samesign;

            double min_pt = 26.;

            if(year == 2016){ 
                good_trigger = HLT_IsoMu24 || HLT_IsoTkMu24;
                min_pt = 26.;
            }
            else if(year == 2017) {
                good_trigger = HLT_IsoMu27;
                min_pt = 29.;
            }
            else if(year == 2018) {
                good_trigger = HLT_IsoMu24;
                min_pt = 26.;
            }

            dimuon_id = mu_IsTightMuon[0] && mu_IsTightMuon[1] &&
                mu_Pt[0] > min_pt &&  mu_Pt[1] > 15. &&
                abs(mu_Eta[0]) < 2.4 && abs(mu_Eta[1]) < 2.4;
            //See https://twiki.cern.ch/twiki/bin/viewauth/CMS/SWGuideMuonIdRun2 for iso cuts
            mu_iso0 = mu_PFIso[0] < mu_iso_cut;
            mu_iso1 = mu_PFIso[1] < mu_iso_cut;

            if(mu_Charge[0] >0){
                mu_p.SetPtEtaPhiM(mu_Pt[0], mu_Eta[0], mu_Phi[0], mu_mass);
                mu_m.SetPtEtaPhiM(mu_Pt[1], mu_Eta[1], mu_Phi[1], mu_mass);
            }
            else{
                mu_m.SetPtEtaPhiM(mu_Pt[0], mu_Eta[0], mu_Phi[0], mu_mass);
                mu_p.SetPtEtaPhiM(mu_Pt[1], mu_Eta[1], mu_Phi[1], mu_mass);
            }



            cm = mu_p + mu_m;
            cm_m = cm.M();
        }
    }

    if(do_electrons){
        opp_sign = good_trigger = dielec_id = el_iso0 = el_iso1 = false;
        if(el_size >= 2){

            opp_sign = ((abs(el_Charge[0] - el_Charge[1])) > 0.01);
            good_sign = opp_sign ^ do_samesign;

            double min_pt = 29.;

            
            if(year == 2016){
                good_trigger = HLT_El27;
                min_pt = 29.;
            }
            else if(year == 2017){
                good_trigger = HLT_El35;
                min_pt = 38.;
            }
            else if(year == 2018){
                good_trigger = HLT_El32;
                min_pt = 35.;
            }

            dielec_id = el_IDMedium_NoIso[0] && el_IDMedium_NoIso[1] &&
                el_ScaleCorr[0] * el_Pt[0] > min_pt &&  el_ScaleCorr[1] * el_Pt[1] > 15. &&
                goodElEta(el_SCEta[0]) && goodElEta(el_SCEta[1]);

            el_iso0 = el_IDMedium[0];
            el_iso1 = el_IDMedium[1];


            if(el_Charge[0] >0){
                el_p.SetPtEtaPhiE(el_ScaleCorr[0] * el_Pt[0], el_Eta[0], el_Phi[0], el_ScaleCorr[0] * el_E[0]);
                el_m.SetPtEtaPhiE(el_ScaleCorr[1] * el_Pt[1], el_Eta[1], el_Phi[1], el_ScaleCorr[1] * el_E[1]);
                elp_index = 0;
                elm_index = 1;
            }
            else{
                el_m.SetPtEtaPhiE(el_ScaleCorr[0] * el_Pt[0], el_Eta[0], el_Phi[0], el_ScaleCorr[0] * el_E[0]);
                el_p.SetPtEtaPhiE(el_ScaleCorr[1] * el_Pt[1], el_Eta[1], el_Phi[1], el_ScaleCorr[1] * el_E[1]);
                elm_index = 0;
                elp_index = 1;
            }

            cm = el_p + el_m;
            cm_m = cm.M();
        }
    }
    if(do_emu){
        opp_sign = good_trigger = emu_ids = el_iso0 = mu_iso0 = false;
        if(el_size >=1 && mu_size >= 1){
            opp_sign = ((abs(el_Charge[0] - mu_Charge[0])) > 0.01);
            good_sign = opp_sign ^ do_samesign;

            double min_pt = 26.;

            if(year == 2016){ 
                good_trigger = HLT_IsoMu24 || HLT_IsoTkMu24;
                min_pt = 26.;
            }
            else if(year == 2017) {
                good_trigger = HLT_IsoMu27;
                min_pt = 29.;
            }
            else if(year == 2018) {
                good_trigger = HLT_IsoMu24;
                min_pt = 26.;
            }


            emu_ids = el_IDMedium_NoIso[0] && mu_IsTightMuon[0] &&
                mu_Pt[0] > min_pt &&  el_ScaleCorr[1] * el_Pt[1] > 15. &&
                abs(mu_Eta[0])  && goodElEta(el_SCEta[1]);

            el_iso0 = el_IDMedium[0];
            mu_iso0 = mu_PFIso[0] < mu_iso_cut;

            mu.SetPtEtaPhiE(mu_Pt[0], mu_Eta[0], mu_Phi[0], mu_E[0]);
            el.SetPtEtaPhiE(el_ScaleCorr[0] * el_Pt[0], el_Eta[0], el_Phi[0], el_ScaleCorr[0] * el_E[0]);

            cm = el + mu;
            cm_m = cm.M();
        }
    }

}
void NTupleReader::fillEvent(){
    nEvents++;
    xF = compute_xF(cm); 
    //pick out 2 highest pt jets with eta < 2.4
    nJets =0;
    has_nobjets = 1;

    for(unsigned int j=0; j < jet_size; j++){
        if(jet_Pt[j] > 20. && std::abs(jet_Eta[j]) < 2.4){
            if(nJets == 1){
                jet2_pt = jet_Pt[j];
                jet2_eta = jet_Eta[j];
                jet2_btag = jet_btag[j];
                has_nobjets = has_nobjets && (jet2_btag < bjet_med_cut);
                if(!is_data) jet2_flavour = jet_hadronflavour[j];
                nJets =2;
                break;
            }
            else if(nJets ==0){
                jet1_pt = jet_Pt[j];
                jet1_eta = jet_Eta[j];
                jet1_btag = jet_btag[j];
                has_nobjets = has_nobjets && (jet1_btag < bjet_med_cut);
                if(!is_data) jet1_flavour = jet_hadronflavour[j];
                nJets = 1;
            }
        }
    }
    if(do_muons){
        mu1_pt = mu_Pt[0];
        mu2_pt = mu_Pt[1];
        mu1_eta = mu_Eta[0];
        mu2_eta = mu_Eta[1];
        mu1_charge = mu_Charge[0];
        mu2_charge = mu_Charge[1];

        cost = get_cost(mu_p, mu_m, false);
        if(cm.Pz() < 0.) cost_r = -cost;
        else cost_r = cost;
    }
    if(do_electrons){
        el1_pt = el_Pt[0];
        el2_pt = el_Pt[1];
        el1_eta = el_Eta[0];
        el2_eta = el_Eta[1];
        el1_charge = el_Charge[0];
        el2_charge = el_Charge[1];

        cost = get_cost(el_p, el_m, false);
        if(cm.Pz() < 0.) cost_r = -cost;
        else cost_r = cost;
    }
    if(do_emu){
            mu1_pt = mu_Pt[0];
            mu1_eta = mu_Eta[0];
            mu1_charge = mu_Charge[0];

            el1_pt = el_Pt[0];
            el1_eta = el_Eta[0];
            el1_charge = el_Charge[0];

            if(el1_charge > 0) cost = get_cost(el, mu, false);
            else cost = get_cost(mu,el, false);
            if(cm.Pz() < 0.) cost_r = -cost;
            else cost_r = cost;
    }

    if(!is_data){
        gen_weight = evt_Gen_Weight * normalization;
    }

}

void NTupleReader::fillEventSFs(){
    pu_SF = get_pileup_SF(pu_NtrueInt, pu_sys.ratio_pileup_nom);
    pu_SF_up = get_pileup_SF(pu_NtrueInt, pu_sys.ratio_pileup_up);
    pu_SF_down = get_pileup_SF(pu_NtrueInt, pu_sys.ratio_pileup_down);
    //printf("pu, pu_up, pu_down: %.2f %.2f %.2f \n", pu_SF, pu_SF_up, pu_SF_down);

    if(pdf_size <60){
        for(int i=0;i<60; i++){
            pdf_weights[i] = 1.;
        }
    }
    if(scale_size > 0){
        mu_R_up = scale_Weights[2];
        mu_R_down = scale_Weights[4];
        mu_F_up = scale_Weights[0];
        mu_F_down = scale_Weights[1];
        mu_RF_up = scale_Weights[3];
        mu_RF_down = scale_Weights[5];
    }
    else{
        mu_R_up = mu_R_down = mu_F_up = mu_F_down = mu_RF_up = mu_RF_down = 1.0;
    }

    if(alphas_size >= 2){
        alpha_up = alpha_weights[0];
        alpha_down = alpha_weights[1];
    }
    else{
        alpha_up = alpha_down = 1.;
    }

    if(do_muons){

        era1_HLT_SF = get_HLT_SF(mu1_pt, mu1_eta, mu2_pt, mu2_eta, era1.HLT_SF, era1.HLT_MC_EFF);
        era2_HLT_SF = get_HLT_SF(mu1_pt, mu1_eta, mu2_pt, mu2_eta, era2.HLT_SF, era2.HLT_MC_EFF);

        era1_id_SF = get_mu_SF(mu1_pt, mu1_eta,  year, era1.ID_SF) * get_mu_SF(mu2_pt, mu2_eta,  year, era1.ID_SF);
        era2_id_SF = get_mu_SF(mu1_pt, mu1_eta,  year, era2.ID_SF) * get_mu_SF(mu2_pt, mu2_eta,  year, era2.ID_SF);

        era2_iso_SF = get_mu_SF(mu1_pt, mu1_eta, year, era2.ISO_SF) * get_mu_SF(mu2_pt, mu2_eta, year, era2.ISO_SF);
        era1_iso_SF = get_mu_SF(mu1_pt, mu1_eta, year, era1.ISO_SF) * get_mu_SF(mu2_pt, mu2_eta, year, era1.ISO_SF);

    }
    if(do_electrons){

        elp_scale_stat_up = el_ScaleCorrStatUp[elp_index] / el_ScaleCorr[elp_index];
        elp_scale_stat_down = el_ScaleCorrStatDown[elp_index]/ el_ScaleCorr[elp_index];
        elp_scale_syst_up = el_ScaleCorrSystUp[elp_index] / el_ScaleCorr[elp_index];
        elp_scale_syst_down = el_ScaleCorrSystDown[elp_index]/ el_ScaleCorr[elp_index];
        elp_scale_gain_up = el_ScaleCorrGainUp[elp_index] / el_ScaleCorr[elp_index];
        elp_scale_gain_down = el_ScaleCorrGainDown[elp_index]/ el_ScaleCorr[elp_index];

        elp_smear_up = el_ScaleSmearUp[elp_index]/ el_ScaleCorr[elp_index];
        elp_smear_down = el_ScaleSmearDown[elp_index]/ el_ScaleCorr[elp_index];

        elm_scale_stat_up = el_ScaleCorrStatUp[elm_index] / el_ScaleCorr[elm_index];
        elm_scale_stat_down = el_ScaleCorrStatDown[elm_index]/ el_ScaleCorr[elm_index];
        elm_scale_syst_up = el_ScaleCorrSystUp[elm_index] / el_ScaleCorr[elm_index];
        elm_scale_syst_down = el_ScaleCorrSystDown[elm_index]/ el_ScaleCorr[elm_index];
        elm_scale_gain_up = el_ScaleCorrGainUp[elm_index] / el_ScaleCorr[elm_index];
        elm_scale_gain_down = el_ScaleCorrGainDown[elm_index]/ el_ScaleCorr[elm_index];

        elm_smear_up = el_ScaleSmearUp[elm_index]/ el_ScaleCorr[elm_index];
        elm_smear_down = el_ScaleSmearDown[elm_index]/ el_ScaleCorr[elm_index];

        el_id_SF = get_el_SF(el1_pt, el1_eta, el_SF.ID_SF) * get_el_SF(el2_pt, el2_eta, el_SF.ID_SF);
        el_reco_SF = get_el_SF(el1_pt, el1_eta, el_SF.RECO_SF) * get_el_SF(el2_pt, el2_eta, el_SF.RECO_SF);
        el_HLT_SF = get_el_HLT_SF(el1_pt, el1_eta, el2_pt, el2_eta, el_SF.HLT_SF, el_SF.HLT_MC_EFF);
    }
    if(do_emu){

        era1_HLT_SF = get_HLT_SF_1mu(mu1_pt, mu1_eta, era1.HLT_SF);
        era2_HLT_SF = get_HLT_SF_1mu(mu1_pt, mu1_eta, era2.HLT_SF);

        era1_iso_SF = get_mu_SF(mu1_pt, mu1_eta,year, era1.ISO_SF);
        era1_id_SF = get_mu_SF(mu1_pt, mu1_eta, year, era1.ID_SF);

        era2_iso_SF = get_mu_SF(mu1_pt, mu1_eta, year, era2.ISO_SF);
        era2_id_SF = get_mu_SF(mu1_pt, mu1_eta,  year, era2.ID_SF);


        el_id_SF = get_el_SF(el1_pt, el1_eta, el_SF.ID_SF);
        el_reco_SF = get_el_SF(el1_pt, el1_eta, el_SF.RECO_SF);
    }
}

void NTupleReader::fillEventRC(){

    int mu_p_n_TL, mu_m_n_TL;
    double mu_p_SF_alt1, mu_p_SF_alt2, mu_p_SF_alt3, mu_p_SF_alt4; 
    double mu_m_SF_alt1, mu_m_SF_alt2, mu_m_SF_alt3, mu_m_SF_alt4; 



    if(mu_Charge[0] < 0){
        mu_m_n_TL = (int) mu_NumberTrackerLayers[0];
        mu_p_n_TL = (int) mu_NumberTrackerLayers[1];
    }
    else{
        mu_m_n_TL = (int) mu_NumberTrackerLayers[1];
        mu_p_n_TL = (int) mu_NumberTrackerLayers[0];
    }

    Double_t mu_p_SF_vars[100], mu_m_SF_vars[100];
    if(is_data){
        mu_p_SF = rc.kScaleDT(1, mu_p.Pt(), mu_p.Eta(), mu_p.Phi(), 0, 0);
        mu_m_SF = rc.kScaleDT(-1, mu_m.Pt(), mu_m.Eta(), mu_m.Phi(), 0, 0);

        applyRC();
        return;

    }
    else if(!is_data && RC_from_gen){

        mu_p_SF = rc.kSpreadMC(1, mu_p.Pt(), mu_p.Eta(), mu_p.Phi(),  gen_mu_p_vec.Pt(),  0, 0);
        mu_m_SF = rc.kSpreadMC(-1, mu_m.Pt(), mu_m.Eta(), mu_m.Phi(),  gen_mu_m_vec.Pt(),  0, 0);

        mu_p_SF_alt1 = rc.kSpreadMC(1, mu_p.Pt(), mu_p.Eta(), mu_p.Phi(),  gen_mu_p_vec.Pt(),  2, 0) - mu_p_SF;
        mu_m_SF_alt1 = rc.kSpreadMC(-1, mu_m.Pt(), mu_m.Eta(), mu_m.Phi(),  gen_mu_m_vec.Pt(),  2, 0) - mu_m_SF;

        mu_p_SF_alt2 = rc.kSpreadMC(1, mu_p.Pt(), mu_p.Eta(), mu_p.Phi(),  gen_mu_p_vec.Pt(),  3, 0) - mu_p_SF;
        mu_m_SF_alt2 = rc.kSpreadMC(-1, mu_m.Pt(), mu_m.Eta(), mu_m.Phi(),  gen_mu_m_vec.Pt(),  3, 0) - mu_m_SF;;

        mu_p_SF_alt3 = rc.kSpreadMC(1, mu_p.Pt(), mu_p.Eta(), mu_p.Phi(),  gen_mu_p_vec.Pt(),  4, 0) - mu_p_SF;
        mu_m_SF_alt3 = rc.kSpreadMC(-1, mu_m.Pt(), mu_m.Eta(), mu_m.Phi(),  gen_mu_m_vec.Pt(),  4, 0) - mu_m_SF;

        mu_p_SF_alt4 = rc.kSpreadMC(1, mu_p.Pt(), mu_p.Eta(), mu_p.Phi(),  gen_mu_p_vec.Pt(),  5, 0) - mu_p_SF;
        mu_m_SF_alt4 = rc.kSpreadMC(-1, mu_m.Pt(), mu_m.Eta(), mu_m.Phi(),  gen_mu_m_vec.Pt(),  5, 0) - mu_m_SF;


        for(int k=0; k<100; k++){
            //printf("P %.1f %.1f %.1f %i %.1f %.3f %i %i \n" ,mu_p.Pt(), mu_p.Eta(), mu_p.Phi(),  gen_mu_p_vec.Pt(), rand1, 1, k);
            mu_p_SF_vars[k] = rc.kSpreadMC(1, mu_p.Pt(), mu_p.Eta(), mu_p.Phi(),  gen_mu_p_vec.Pt(), 1, k);
            //printf("M %.1f %.1f %.1f %i %.1f %.3f %i %i \n" ,mu_m.Pt(), mu_m.Eta(), mu_m.Phi(),  gen_mu_m_vec.Pt(), rand2, 1, k);
            mu_m_SF_vars[k] = rc.kSpreadMC(-1, mu_m.Pt(), mu_m.Eta(), mu_m.Phi(),  gen_mu_m_vec.Pt(), 1, k);
        }
        //printf("7 \n");
    }

    else{
        double rand1a = rand->Rndm(); 
        double rand2a = rand->Rndm(); 


        mu_p_SF = rc.kSmearMC(1, mu_p.Pt(), mu_p.Eta(), mu_p.Phi(), mu_p_n_TL, rand1a,  0, 0);
        mu_m_SF = rc.kSmearMC(-1, mu_m.Pt(), mu_m.Eta(), mu_m.Phi(), mu_m_n_TL, rand2a,  0, 0);

        mu_p_SF_alt1 = rc.kSmearMC(1, mu_p.Pt(), mu_p.Eta(), mu_p.Phi(), mu_p_n_TL, rand1a,  2, 0) - mu_p_SF;
        mu_m_SF_alt1 = rc.kSmearMC(-1, mu_m.Pt(), mu_m.Eta(), mu_m.Phi(), mu_m_n_TL, rand2a,  2, 0)- mu_m_SF;

        mu_p_SF_alt2 = rc.kSmearMC(1, mu_p.Pt(), mu_p.Eta(), mu_p.Phi(), mu_p_n_TL, rand1a,  3, 0) - mu_p_SF;
        mu_m_SF_alt2 = rc.kSmearMC(-1, mu_m.Pt(), mu_m.Eta(), mu_m.Phi(), mu_m_n_TL, rand2a,  3, 0)- mu_m_SF;

        mu_p_SF_alt3 = rc.kSmearMC(1, mu_p.Pt(), mu_p.Eta(), mu_p.Phi(), mu_p_n_TL, rand1a,  4, 0) - mu_p_SF;
        mu_m_SF_alt3 = rc.kSmearMC(-1, mu_m.Pt(), mu_m.Eta(), mu_m.Phi(), mu_m_n_TL, rand2a,  4, 0)- mu_m_SF;

        mu_p_SF_alt4 = rc.kSmearMC(1, mu_p.Pt(), mu_p.Eta(), mu_p.Phi(), mu_p_n_TL, rand1a,  5, 0) - mu_p_SF;
        mu_m_SF_alt4 = rc.kSmearMC(-1, mu_m.Pt(), mu_m.Eta(), mu_m.Phi(), mu_m_n_TL, rand2a,  5, 0)- mu_m_SF;
        

        for(int k=0; k<100; k++){
            mu_p_SF_vars[k] = rc.kSmearMC(1, mu_p.Pt(), mu_p.Eta(), mu_p.Phi(), mu_p_n_TL, rand1a,  1, k);
            mu_m_SF_vars[k] = rc.kSmearMC(-1, mu_m.Pt(), mu_m.Eta(), mu_m.Phi(), mu_m_n_TL, rand2a,  1, k);
        }
    }


    double mu_p_SF_var = get_var(mu_p_SF_vars);
    double mu_m_SF_var = get_var(mu_m_SF_vars);


    double mu_p_unc = sqrt(mu_p_SF_var + mu_p_SF_alt1*mu_p_SF_alt1 + mu_p_SF_alt2*mu_p_SF_alt2  + mu_p_SF_alt3*mu_p_SF_alt3  + mu_p_SF_alt4*mu_p_SF_alt4);
    double mu_m_unc = sqrt(mu_m_SF_var + mu_m_SF_alt1*mu_m_SF_alt1 + mu_m_SF_alt2*mu_m_SF_alt2  + mu_m_SF_alt3*mu_m_SF_alt3  + mu_m_SF_alt4*mu_m_SF_alt4);



    mu_p_SF_up = mu_p_SF + mu_p_unc;
    mu_p_SF_down = mu_p_SF - mu_p_unc;
    mu_m_SF_up = mu_m_SF + mu_m_unc;
    mu_m_SF_down = mu_m_SF - mu_m_unc;

    applyRC();

}

void NTupleReader::applyRC(){
    mu_p.SetPtEtaPhiE(mu_p.Pt() * mu_p_SF, mu_p.Eta(), mu_p.Phi(), mu_p_SF * mu_p.E());
    mu_m.SetPtEtaPhiE(mu_m.Pt() * mu_m_SF, mu_m.Eta(), mu_m.Phi(), mu_m_SF * mu_m.E());
    cost = get_cost(mu_p, mu_m);
    cm = mu_p + mu_m;
    cm_m = cm.M();
    xF = compute_xF(cm); 
}


bool NTupleReader::parseGenParts(bool PRINT = false){
    //returns false if unable to match all gen parts

    char out_buff[50000];
    bool print_out = false;


    if(PRINT) memset(out_buff, 0, 10000);
    if(PRINT) sprintf(out_buff + strlen(out_buff),"Event %i \n", event_idx);

    //GEN LEVEL
    //
    //Pythia 8 status numbering convention
    //see:
    //http://home.thep.lu.se/~torbjorn/pythia81html/EventRecord.html
    int FINAL_STATE = 1;
    int EVENT_PARTICLE = 11;
    int BEAM_PARTICLE=4;
    int INCIDENT_PARTICLE = 21;
    int INTERMED_PARTICLE = 22;
    int OUTGOING = 23;
    //Particle ID's
    int ELECTRON = 11; 
    int MUON = 13;
    int TAU = 15;
    int PHOTON = 22;
    int Z=23;
    int GLUON = 21;
    int PROTON = 2212;


    int MY_LEP;
    if(do_muons) MY_LEP = MUON;
    else MY_LEP = ELECTRON;

    int inc_1 =-1;
    int inc_2 =-1;
    int gen_lep_p=-1;
    int gen_lep_m=-1;
    int gen_tau_p=-1;
    int gen_tau_m=-1;
    int gen_e_p=-1;
    int gen_e_m=-1;
    int intermed=-1;


    signal_event = false;//whether it is an event with an asym or not
    failed_match = false;

    is_tau_event = false;

    for(unsigned int k=0; k<gen_size; k++){
        if(gen_status[k] == INCIDENT_PARTICLE && 
                (abs(gen_id[k]) <=6  || gen_id[k] == GLUON) && 
                (abs(gen_Dau0ID[k]) == MY_LEP || gen_Dau0ID[k] == Z || 
                 gen_Dau0ID[k] == PHOTON || abs(gen_Dau0ID[k]) == TAU)
          ){
            //record index of 2 initial state particles
            if(inc_1 == -1) inc_1 = k;
            else if(inc_2 == -1) inc_2 = k;
            else{
                print_out = true;
                printf("WARNING: More than 2 incident particles in event\n\n");
            }

        }
        //record 2 final leptons
        if(abs(gen_id[k]) == MY_LEP && 
                (gen_Mom0ID[k] == Z || gen_Mom0ID[k] == PHOTON || (abs(gen_Mom0ID[k]) == TAU && gen_Pt[k] > 10.)
                 || (gen_status[k] == OUTGOING && gen_Mom0ID[k] != PROTON))) {
            if(gen_id[k] == MY_LEP){
                if(gen_lep_m == -1) gen_lep_m = k;
                else{
                    if(abs(gen_Mom0ID[k]) != TAU) printf("WARNING: More than one lep_m\n\n");
                    if(PRINT) sprintf(out_buff + strlen(out_buff), "Extra lep_m detected\n");
                    print_out = true;
                }
            }
            if(gen_id[k] == -MY_LEP){
                if(gen_lep_p == -1) gen_lep_p = k;
                else{
                    if(abs(gen_Mom0ID[k]) != TAU) printf("WARNING: More than one lep_p\n\n");
                    if(PRINT) sprintf(out_buff + strlen(out_buff), "Extra lep_p detected\n");
                    print_out = true;
                }
            }
        }
        //record tau's
        if(abs(gen_id[k]) == TAU && 
                ( (gen_Mom0ID[k] == Z || gen_Mom0ID[k] == PHOTON || abs(gen_Mom0ID[k]) == TAU ||
                   gen_status[k] == OUTGOING)) ){
            if(gen_id[k] == TAU){
                if(gen_tau_m == -1) gen_tau_m = k;
                else{
                    //printf("WARNING: More than one tau_m\n\n");
                    if(PRINT) sprintf(out_buff + strlen(out_buff), "Extra tau_m detected\n");
                    //print_out = true;
                }
            }
            if(gen_id[k] == -TAU){
                if(gen_tau_p == -1) gen_tau_p = k;
                else{
                    //printf("WARNING: More than one tau_p\n\n");
                    if(PRINT) sprintf(out_buff + strlen(out_buff), "Extra tau_p detected\n");
                    //print_out = true;
                }
            }
        }
        if(PRINT){
            if( (abs(gen_id[k]) <=6 || gen_id[k] == GLUON)){
                if(PRINT) sprintf(out_buff + strlen(out_buff),"Parton (ID = %i stat = %i): \n"
                        "    Mom1 ID: %i Mom1 Stat:  Mom2 ID: %i Mom2 Stat  \n"
                        "    Dau1 ID: %i Dau1 Stat: Dau2 ID: %i Dau2 Stat  \n",
                        gen_id[k], gen_status[k], 
                        gen_Mom0ID[k],  gen_Mom1ID[k], 
                        gen_Dau0ID[k],  gen_Dau1ID[k] );
            }


            if(abs(gen_id[k]) == MUON || abs(gen_id[k]) == ELECTRON || abs(gen_id[k]) == TAU){
                if(PRINT) sprintf(out_buff + strlen(out_buff),"lep (id %i stat = %i): \n"
                        "   Mom1 ID: %i Mom1 Stat:  Mom2 ID: %i Mom2 Stat  \n"
                        "   Dau1 ID: %i Dau1 Stat:  Dau2 ID: %i Dau2 Stat  \n",
                        gen_id[k], gen_status[k], 
                        gen_Mom0ID[k],  gen_Mom1ID[k], 
                        gen_Dau0ID[k],  gen_Dau1ID[k] );

                if(PRINT) sprintf(out_buff + strlen(out_buff),"(pt, eta, phi,E) = %4.2f %4.2f %4.2f %4.2f \n", 
                        gen_Pt[k], gen_Eta[k], gen_Phi[k], gen_E[k]);
            }
            if(gen_id[k] == Z){
                if(PRINT) sprintf(out_buff + strlen(out_buff),"Z (ID = %i, status = %i): \n"
                        "   Mom1 ID: %i Mom1 Stat:  Mom2 ID: %i Mom2 Stat  \n"
                        "   Dau1 ID: %i Dau1 Stat:  Dau2 ID: %i Dau2 Stat  \n",
                        gen_id[k], gen_status[k],
                        gen_Mom0ID[k],  gen_Mom1ID[k], 
                        gen_Dau0ID[k],  gen_Dau1ID[k] );
            }
        }
    }


    if(gen_lep_p != -1 && gen_lep_m != -1) {
        if(PRINT) sprintf(out_buff + strlen(out_buff),"lep_p: \n"
                "   Mom1 ID: %i Mom1 Stat: Mom2 ID: %i Mom2 Stat \n"
                "   Dau1 ID: %i Dau1 Stat:  Dau2 ID: %i Dau2 Stat \n",
                gen_Mom0ID[gen_lep_p],  gen_Mom1ID[gen_lep_p], 
                gen_Dau0ID[gen_lep_p], gen_Dau1ID[gen_lep_p]);

        if(PRINT) sprintf(out_buff + strlen(out_buff),"lep_m: \n"
                "   Mom1 ID: %i Mom1 Stat:  Mom2 ID: %i Mom2 Stat  \n"
                "   Dau1 ID: %i Dau1 Stat:  Dau2 ID: %i Dau2 Stat \n",
                gen_Mom0ID[gen_lep_m],  gen_Mom1ID[gen_lep_m], 
                gen_Dau0ID[gen_lep_m],  gen_Dau1ID[gen_lep_m]
                );

    }
    else {
        printf("WARNING: Unable to identify lepton pair in event %i, skipping \n", event_idx);
        nFailedID ++;
        print_out = true;
        if(PRINT && print_out){
            sprintf(out_buff + strlen(out_buff), "\n\n");
            fputs(out_buff, stdout);
            print_out = false;
        }
        failed_match = true;
        return false;
    }
    if((inc_1 == -1) || (inc_2 == -1)){
        printf("WARNING: Unable to identify initial state particles in event %i, skipping \n", event_idx);
        nFailedID ++;
        print_out = true;
        if(PRINT && print_out){
            sprintf(out_buff + strlen(out_buff), "\n\n");
            fputs(out_buff, stdout);
            print_out = false;
        }
        failed_match = true;
        return false;
    }

    else{ 
        //printf("%i %i \n", inc_1, inc_2);
        int inc_id1 = gen_id[inc_1];
        int inc_id2 = gen_id[inc_2];
        if(abs(gen_Dau0ID[inc_1]) == TAU){
            if(gen_tau_p == -1 || gen_tau_m == -1){
                printf("Didn't record tau's :( \n");
                nFailedID ++;
                failed_match = true;
                return false;
            }
            nTauTau++;
            is_tau_event = true;
            signal_event = false;
        }
        if((abs(inc_id1) <= 6 && abs(inc_id2) <= 6) && (inc_id1 * inc_id2 < 0)){ //a quark and anti quark
            //qq-bar
            signal_event = true;
            nQQb++;
            if(inc_id1>0) quark_dir_eta = gen_Eta[inc_1];
            else if(inc_id2>0) quark_dir_eta = gen_Eta[inc_2];
        }
        else if(((abs(inc_id1) <= 6) && (inc_id2 == 21)) ||
                ((abs(inc_id2) <= 6) && (inc_id1 == 21))){ //qglu
            signal_event = true;
            int q_dir;
            if(inc_id1 == 21){
                if(inc_id2 <0) quark_dir_eta = gen_Eta[inc_1];//qbar-glu, want glu dir
                else quark_dir_eta= gen_Eta[inc_2];//q-glu ,want q dir
            }
            else if(inc_id2 == 21) {
                if(inc_id1 <0) quark_dir_eta = gen_Eta[inc_2];//qbar-glu, want glu dir
                else quark_dir_eta= gen_Eta[inc_1];//q-glu ,want q dir
            }
            nQGlu++;
        }
        else if((abs(inc_id1) <= 6) && (abs(inc_id2) <= 6) && (inc_id1 * inc_id2 >0)){ //2 quarks
            if(PRINT) sprintf(out_buff + strlen(out_buff),"QQ Event \n");
            signal_event = false;
            nQQ++;
        }
        else if((inc_id1 == 21) && (inc_id2 == 21)){ //gluglu
            signal_event = false;
            nGluGlu++;
            if(PRINT) sprintf(out_buff + strlen(out_buff), "Glu Glu event \n");
        }
        else {
            printf("WARNING: not qqbar, qq, qg, or gg event");
            printf("First particle was %i second particle was %i \n \n ", inc_id1, inc_id2);
            failed_match = true;
            nFailedID ++;
        }
    }

    if(PRINT){
        sprintf(out_buff + strlen(out_buff),  "1st Particle %i:(pt, eta, phi,E) = %4.2f %4.2f %4.2f %4.2f \n", 
                gen_id[inc_1], gen_Pt[inc_1], gen_Eta[inc_1], gen_Phi[inc_1], gen_E[inc_1]);
        sprintf(out_buff + strlen(out_buff),"2nd Particle %i:(pt, eta, phi,E) = %4.2f %4.2f %4.2f %4.2f \n", 
                gen_id[inc_2], gen_Pt[inc_2], gen_Eta[inc_2], gen_Phi[inc_2], gen_E[inc_2]);
    }
    float gen_cost;
    if(do_muons){
        gen_mu_p_vec.SetPtEtaPhiE(gen_Pt[gen_lep_p], gen_Eta[gen_lep_p], gen_Phi[gen_lep_p], gen_E[gen_lep_p]);
        gen_mu_m_vec.SetPtEtaPhiE(gen_Pt[gen_lep_m], gen_Eta[gen_lep_m], gen_Phi[gen_lep_m], gen_E[gen_lep_m]);
        gen_cm = gen_mu_p_vec + gen_mu_m_vec;
        gen_cost = get_cost(gen_mu_p_vec, gen_mu_m_vec, false);
        if(std::isnan(gen_cost)) gen_cost = get_cost_v2(gen_mu_p_vec, gen_mu_m_vec);
        if(std::isnan(gen_cost)) gen_cost = cost;
    }
    else{
        gen_el_p_vec.SetPtEtaPhiE(gen_Pt[gen_lep_p], gen_Eta[gen_lep_p], gen_Phi[gen_lep_p], gen_E[gen_lep_p]);
        gen_el_m_vec.SetPtEtaPhiE(gen_Pt[gen_lep_m], gen_Eta[gen_lep_m], gen_Phi[gen_lep_m], gen_E[gen_lep_m]);
        gen_cm = gen_el_p_vec + gen_el_m_vec;
        gen_cost = get_cost(gen_el_p_vec, gen_el_m_vec, false);
        if(std::isnan(gen_cost)) gen_cost = get_cost_v2(gen_el_p_vec, gen_el_m_vec);
        if(std::isnan(gen_cost)) gen_cost = cost;
    }
    if(quark_dir_eta < 0){
        cost_st = -gen_cost;
    }
    else cost_st = gen_cost;

    gen_m = gen_cm.M();

    if(PRINT && print_out){
        sprintf(out_buff + strlen(out_buff), "\n\n");
        fputs(out_buff, stdout);
        print_out = false;
    }

    if(PRINT) memset(out_buff, 0, 10000);
    return !failed_match;
}

void NTupleReader::finish(){

    fout->cd();

    printf("Finished. There were %i events from %i files \n\n", nEvents, fileCount);
    printf("Writing output to file at %s \n", fout->GetName());
    for(int i=0; i<nOutTrees; i++){
        outTrees[i]->Write();
    }

    fout->Close();
    fclose(root_files);
}




NTupleReader::~NTupleReader(){

}
