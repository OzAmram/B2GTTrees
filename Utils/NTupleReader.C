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

NTupleReader::NTupleReader(char file_list[100], char fout_name[100], bool is_data_){

    is_data = is_data_;
    root_files = fopen(file_list, "r");
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

    //separate SFs for runs BCDEF and GH
    printf("getting SFs \n");
    setup_SFs(&runs_bcdef, &runs_gh, &pu_SFs);



    if(do_electrons){
        setup_el_SF(&el_SF);
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
            int nparams = sscanf(lines, "! idx = %i xsec = %f \n", &sample_idx, &xsec);
            if(nparams < 2 || sample_idx >= MAX_SAMPLES){
                printf("ERROR: Unable to parse sample header. Exiting");
                exit(EXIT_FAILURE);
            }
            normalization = norms[sample_idx];
            printf("Moving on to sample %i which has normalization %e \n", sample_idx, normalization);
        }
        else if(normalization > 0) {//root file
            fileCount++;
            if(fileCount % nJobs != iJob) continue;



            char * end;
            //remove trailing whitespace
            end = lines + strlen(lines) - 1;
            while(end > lines && isspace((char)*end)) end--;
            // Write new null terminator
            *(end+1) = 0;

            printf("Opening file: %s \n", lines);
            fin=  TFile::Open(lines);

            fin->cd("EventCounter");
            TDirectory *subdir = gDirectory;
            TH1D *mc_pileup = (TH1D *)subdir->Get("pileup");
            mc_pileup->Scale(1./mc_pileup->Integral());
            pu_SFs.pileup_ratio->Divide(pu_SFs.data_pileup, mc_pileup);

            fin->cd("B2GTTreeMaker");
            tin = (TTree *)gDirectory->Get("B2GTree");

            tin->SetBranchAddress("jetAK4CHS_size", &jet_size);
            tin->SetBranchAddress("jetAK4CHS_Pt", &jet_Pt);
            tin->SetBranchAddress("jetAK4CHS_Eta", &jet_Eta);
            tin->SetBranchAddress("jetAK4CHS_Phi", &jet_Phi);
            tin->SetBranchAddress("jetAK4CHS_E", &jet_E);
            tin->SetBranchAddress("jetAK4CHS_CSVv2", &jet_CSV);
            tin->SetBranchAddress("jetAK4CHS_CMVAv2", &jet_CMVA);
            tin->SetBranchAddress("jetAK4CHS_PartonFlavour", &jet_partonflavour);

            tin->SetBranchAddress("met_MuCleanOnly_size", &met_size);
            tin->SetBranchAddress("met_MuCleanOnly_Pt", &met_pt);

            if(do_muons){
                tin->SetBranchAddress("HLT_IsoMu24", &HLT_IsoMu);
                tin->SetBranchAddress("HLT_IsoTkMu24", &HLT_IsoTkMu);

                tin->SetBranchAddress("mu_size", &mu_size); //number of muons in the event
                tin->SetBranchAddress("mu_TunePMuonBestTrackPt", &mu_Pt);
                tin->SetBranchAddress("mu_Eta", &mu_Eta);
                tin->SetBranchAddress("mu_Phi", &mu_Phi);
                tin->SetBranchAddress("mu_E", &mu_E);
                tin->SetBranchAddress("mu_Charge", &mu_Charge);

                tin->SetBranchAddress("mu_IsHighPtMuon", &mu_IsHighPtMuon);
                tin->SetBranchAddress("mu_TrackerIso", &mu_TrackerIso);
                tin->SetBranchAddress("mu_NumberTrackerLayers", &mu_NumberTrackerLayers);
                //tin->SetBranchAddress("mu_SumChargedHadronPt", &mu_SumChargedHadronPt);
                //tin->SetBranchAddress("mu_SumNeutralHadronPt", &mu_SumNeutralHadronPt);
                //tin->SetBranchAddress("mu_SumPUPt", &mu_SumPUPt);
                //tin->SetBranchAddress("mu_SumPhotonPt", &mu_SumPhotonPt);
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

                tin->SetBranchAddress("gen_Mom0ID", &gen_Mom0ID);
                tin->SetBranchAddress("gen_Mom1ID", &gen_Mom1ID);
                tin->SetBranchAddress("gen_Mom0Status", &gen_Mom0Status);
                tin->SetBranchAddress("gen_Mom1Status", &gen_Mom1Status);

                tin->SetBranchAddress("gen_Dau0ID", &gen_Dau0ID);
                tin->SetBranchAddress("gen_Dau1ID", &gen_Dau1ID);
                tin->SetBranchAddress("gen_Dau0Status", &gen_Dau0Status);
                tin->SetBranchAddress("gen_Dau1Status", &gen_Dau1Status);
            }



            tin_nEntries =  tin->GetEntries();
            return true;
        }
    }
    return false;
}

void NTupleReader::setupOutputTree(char treeName[100]){
    int idx = nOutTrees;
    nOutTrees++;
    outTrees[idx] = new TTree(treeName, "");

    outTrees[idx]->Branch("m", &cm_m, "m/D");
    outTrees[idx]->Branch("xF", &xF, "xF/D");
    outTrees[idx]->Branch("cost", &cost_r, "cost/D");
    outTrees[idx]->Branch("cost_st", &cost_st, "cost_st/D");
    outTrees[idx]->Branch("jet1_pt", &jet1_pt, "jet1_pt/D");
    outTrees[idx]->Branch("jet1_eta", &jet1_eta, "jet1_eta/D");
    outTrees[idx]->Branch("jet1_CMVA", &jet1_cmva, "jet1_CMVA/D");
    outTrees[idx]->Branch("jet2_pt", &jet2_pt, "jet2_pt/D");
    outTrees[idx]->Branch("jet2_eta", &jet2_eta, "jet2_eta/D");
    outTrees[idx]->Branch("jet2_CMVA", &jet2_cmva, "jet2_CMVA/D");
    outTrees[idx]->Branch("nJets", &nJets, "nJets/I");

    if(do_muons){
        outTrees[idx]->Branch("met_pt", &met_pt, "met_Pt/F");
        outTrees[idx]->Branch("mu1_eta", &mu1_eta, "mu1_eta/D");
        outTrees[idx]->Branch("mu2_eta", &mu2_eta, "mu2_eta/D");
        outTrees[idx]->Branch("mu_m", "TLorentzVector", &mu_m);
        outTrees[idx]->Branch("mu_p", "TLorentzVector", &mu_p);
        outTrees[idx]->Branch("mu1_pt", &mu1_pt, "mu1_pt/D");
        outTrees[idx]->Branch("mu2_pt", &mu2_pt, "mu2_pt/D");
    }



    if(!is_data){

        outTrees[idx]->Branch("pu_SF", &pu_SF);
        outTrees[idx]->Branch("gen_weight", &gen_weight, "gen_weight/D");
        outTrees[idx]->Branch("gen_m", &gen_m, "m/D");
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

        if(do_muons){
            outTrees[idx]->Branch("mu_p_SF", &mu_p_SF, "mu_p_SF/D");
            outTrees[idx]->Branch("mu_m_SF", &mu_m_SF, "mu_m_SF/D");
            outTrees[idx]->Branch("mu_p_SF_up", &mu_p_SF_up, "mu_p_SF_up/D");
            outTrees[idx]->Branch("mu_m_SF_up", &mu_m_SF, "mu_m_SF_up/D");
            outTrees[idx]->Branch("mu_p_SF_down", &mu_p_SF_down, "mu_p_SF_down/D");
            outTrees[idx]->Branch("mu_m_SF_down", &mu_m_SF, "mu_m_SF_down/D");
            outTrees[idx]->Branch("mu_p_SF_alt", &mu_p_SF_alt, "mu_p_SF_alt/D");
            outTrees[idx]->Branch("mu_m_SF_alt", &mu_m_SF_alt, "mu_m_SF_alt/D");
            outTrees[idx]->Branch("gen_mu_m", "TLorentzVector", &gen_mu_m_vec);
            outTrees[idx]->Branch("gen_mu_p", "TLorentzVector", &gen_mu_p_vec);
            outTrees[idx]->Branch("bcdef_HLT_SF", &bcdef_HLT_SF);
            outTrees[idx]->Branch("bcdef_iso_SF", &bcdef_iso_SF);
            outTrees[idx]->Branch("bcdef_id_SF", &bcdef_id_SF);
            outTrees[idx]->Branch("bcdef_trk_SF", &bcdef_trk_SF);
            outTrees[idx]->Branch("gh_HLT_SF", &gh_HLT_SF);
            outTrees[idx]->Branch("gh_iso_SF", &gh_iso_SF);
            outTrees[idx]->Branch("gh_id_SF", &gh_id_SF);
            outTrees[idx]->Branch("gh_trk_SF", &gh_trk_SF);
        }
    }
}

void NTupleReader::setupRC(){
    use_RC = true;
    printf("Getting RC \n");
    //RoccoR  rc("rcdata.2016.v3"); //directory path as input for now; initialize only once, contains all variations
    rc =  RoccoR("/uscms_data/d3/oamram/CMSSW_8_0_24_patch1/src/Analysis/B2GTTrees/Utils/rcdata.2016.v3");
    rand = new TRandom3();
}


void NTupleReader::getEvent(int i){
    tin->GetEntry(i);
    event_idx = i;
    if(mu_size > MU_SIZE || gen_size >GEN_SIZE) printf("WARNING: MU_SIZE OR GEN_SIZE TOO LARGE \n");
    if(met_size != 1) printf("WARNING: Met size not equal to 1\n");
    if(do_muons){
        opp_sign = good_trigger = dimuon_id = mu_iso0 = mu_iso1 = false;
        if(mu_size >= 2){

            opp_sign = ((abs(mu_Charge[0] - mu_Charge[1])) > 0.01);
            good_sign = opp_sign ^ do_samesign;
            good_trigger = HLT_IsoMu || HLT_IsoTkMu;

            dimuon_id = mu_IsHighPtMuon[0] && mu_IsHighPtMuon[1] &&
                mu_Pt[0] > 26. &&  mu_Pt[1] > 15. &&
                abs(mu_Eta[0]) < 2.4 && abs(mu_Eta[1]) < 2.4;
            //See https://twiki.cern.ch/twiki/bin/viewauth/CMS/SWGuideMuonIdRun2 for iso cuts
            mu_iso0 = mu_TrackerIso[0] < mu_iso;
            mu_iso1 = mu_TrackerIso[1] < mu_iso;

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
}
void NTupleReader::fillEvent(){
    nEvents++;
    xF = abs(2.*cm.Pz()/13000.); 
    //pick out 2 highest pt jets with eta < 2.4
    nJets =0;
    for(int j=0; j < jet_size; j++){
        if(jet_Pt[j] > 20. && std::abs(jet_Eta[j]) < 2.4){
            if(nJets == 1){
                jet2_pt = jet_Pt[j];
                jet2_eta = jet_Eta[j];
                jet2_cmva = jet_CMVA[j];
                if(!is_data) jet2_flavour = jet_partonflavour[j];
                nJets =2;
                break;
            }
            else if(nJets ==0){
                jet1_pt = jet_Pt[j];
                jet1_eta = jet_Eta[j];
                jet1_cmva = jet_CMVA[j];
                if(!is_data) jet1_flavour = jet_partonflavour[j];
                nJets = 1;
            }
        }
    }
    if(do_muons){
        mu1_pt = mu_Pt[0];
        mu2_pt = mu_Pt[1];
        mu1_eta = mu_Eta[0];
        mu2_eta = mu_Eta[1];
        cost = get_cost(mu_p, mu_m, false);
        if(cm.Pz() < 0.) cost_r = -cost;
        else cost_r = cost;
    }

    if(!is_data){
        gen_weight = evt_Gen_Weight * normalization;
    }

}

void NTupleReader::fillEventSFs(){
    pu_SF = get_pileup_SF(pu_NtrueInt, pu_SFs.pileup_ratio);

    mu_R_up = scale_Weights[2];
    mu_R_down = scale_Weights[4];
    mu_F_up = scale_Weights[0];
    mu_F_down = scale_Weights[1];
    mu_RF_up = scale_Weights[3];
    mu_RF_down = scale_Weights[5];

    alpha_up = alpha_weights[0];
    alpha_down = alpha_weights[1];

    if(do_muons){

        bcdef_HLT_SF = get_HLT_SF(mu1_pt, mu1_eta, mu2_pt, mu2_eta, runs_bcdef.HLT_SF, runs_bcdef.HLT_MC_EFF);
        gh_HLT_SF = get_HLT_SF(mu1_pt, mu1_eta, mu2_pt, mu2_eta, runs_gh.HLT_SF, runs_gh.HLT_MC_EFF);

        bcdef_iso_SF = get_SF(mu1_pt, mu1_eta, runs_bcdef.ISO_SF) * get_SF(mu2_pt, mu2_eta, runs_bcdef.ISO_SF);
        bcdef_id_SF = get_SF(mu1_pt, mu1_eta, runs_bcdef.ID_SF) * get_SF(mu2_pt, mu2_eta, runs_bcdef.ID_SF);

        gh_iso_SF = get_SF(mu1_pt, mu1_eta, runs_gh.ISO_SF) * get_SF(mu2_pt, mu2_eta, runs_gh.ISO_SF);
        gh_id_SF = get_SF(mu1_pt, mu1_eta, runs_gh.ID_SF) * get_SF(mu2_pt, mu2_eta, runs_gh.ID_SF);

        bcdef_trk_SF = get_Mu_trk_SF(abs(mu1_eta), runs_bcdef.TRK_SF) * get_Mu_trk_SF(abs(mu2_eta), runs_bcdef.TRK_SF);
        gh_trk_SF = get_Mu_trk_SF(abs(mu1_eta), runs_gh.TRK_SF) * get_Mu_trk_SF(abs(mu2_eta), runs_gh.TRK_SF);
    }
}

void NTupleReader::fillEventRC(){

    int mu_p_n_TL, mu_m_n_TL;
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
        mu_p_SF_alt = rc.kScaleDT(1, mu_p.Pt(), mu_p.Eta(), mu_p.Phi(), 2, 0);
        mu_m_SF_alt = rc.kScaleDT(-1, mu_m.Pt(), mu_m.Eta(), mu_m.Phi(), 2, 0);


        Double_t mu_p_SF_vars[100], mu_m_SF_vars[100];
        for(int k=0; k<100; k++){
            mu_p_SF_vars[k] = rc.kScaleDT(1, mu_p.Pt(), mu_p.Eta(), mu_p.Phi(), 1, k);
            mu_m_SF_vars[k] = rc.kScaleDT(-1, mu_m.Pt(), mu_m.Eta(), mu_m.Phi(), 1, k);
        }
    }
    else if(!is_data && RC_from_gen){
        double rand1 = rand->Rndm(); 
        double rand2 = rand->Rndm(); 

        mu_p_SF = rc.kScaleFromGenMC(1, mu_p.Pt(), mu_p.Eta(), mu_p.Phi(), mu_p_n_TL, gen_mu_p_vec.Pt(), rand1, 0, 0);
        mu_m_SF = rc.kScaleFromGenMC(-1, mu_m.Pt(), mu_m.Eta(), mu_m.Phi(), mu_m_n_TL, gen_mu_m_vec.Pt(), rand2, 0, 0);
        mu_p_SF_alt = rc.kScaleFromGenMC(1, mu_p.Pt(), mu_p.Eta(), mu_p.Phi(), mu_p_n_TL, gen_mu_p_vec.Pt(), rand1, 2, 0);
        mu_m_SF_alt = rc.kScaleFromGenMC(-1, mu_m.Pt(), mu_m.Eta(), mu_m.Phi(), mu_m_n_TL, gen_mu_m_vec.Pt(), rand2, 2, 0);

        for(int k=0; k<100; k++){
            mu_p_SF_vars[k] = rc.kScaleFromGenMC(1, mu_p.Pt(), mu_p.Eta(), mu_p.Phi(), mu_p_n_TL, gen_mu_p_vec.Pt(), rand1, 1, k);
            mu_m_SF_vars[k] = rc.kScaleFromGenMC(-1, mu_m.Pt(), mu_m.Eta(), mu_m.Phi(), mu_m_n_TL, gen_mu_m_vec.Pt(), rand2, 1, k);
        }
    }

    else{
        double rand1a = rand->Rndm(); 
        double rand1b = rand->Rndm(); 
        double rand2a = rand->Rndm(); 
        double rand2b = rand->Rndm(); 


        mu_p_SF = rc.kScaleAndSmearMC(1, mu_p.Pt(), mu_p.Eta(), mu_p.Phi(), mu_p_n_TL, rand1a, rand1b, 0, 0);
        mu_m_SF = rc.kScaleAndSmearMC(-1, mu_m.Pt(), mu_m.Eta(), mu_m.Phi(), mu_m_n_TL, rand2a, rand2b, 0, 0);
        mu_p_SF_alt = rc.kScaleAndSmearMC(1, mu_p.Pt(), mu_p.Eta(), mu_p.Phi(), mu_p_n_TL, rand1a, rand1b, 2, 0);
        mu_m_SF_alt = rc.kScaleAndSmearMC(-1, mu_m.Pt(), mu_m.Eta(), mu_m.Phi(), mu_m_n_TL, rand2a, rand2b, 2, 0);

        Double_t mu_p_SF_vars[100], mu_m_SF_vars[100];
        for(int k=0; k<100; k++){
            mu_p_SF_vars[k] = rc.kScaleAndSmearMC(1, mu_p.Pt(), mu_p.Eta(), mu_p.Phi(), mu_p_n_TL, rand1a, rand1b, 1, k);
            mu_m_SF_vars[k] = rc.kScaleAndSmearMC(-1, mu_m.Pt(), mu_m.Eta(), mu_m.Phi(), mu_m_n_TL, rand2a, rand2b, 1, k);
        }
    }


    double mu_p_SF_std = sqrt(get_var(mu_p_SF_vars));
    double mu_m_SF_std = sqrt(get_var(mu_m_SF_vars));
    mu_p_SF_up = mu_p_SF + mu_p_SF_std;
    mu_p_SF_down = mu_p_SF - mu_p_SF_std;
    mu_m_SF_up = mu_m_SF + mu_m_SF_std;
    mu_m_SF_down = mu_m_SF - mu_m_SF_std;

}

void NTupleReader::parseGenParts(bool PRINT = false){

    char out_buff[10000];
    bool print_out = false;

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
    int PHOTON = 22;
    int Z=23;
    int GLUON = 21;
    int TAU = 15;
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
        //record 2 scattered muons
        if(abs(gen_id[k]) == MY_LEP && 
                (gen_Mom0ID[k] == Z || gen_Mom0ID[k] == PHOTON || (abs(gen_Mom0ID[k]) == TAU && gen_Pt[k] > 10.)
                 || abs(gen_Mom0ID[k]) == ELECTRON || (gen_status[k] == OUTGOING && gen_Mom0ID[k] != PROTON))) {
            if(gen_id[k] == MY_LEP){
                if(gen_lep_m == -1) gen_lep_m = k;
                else{
                    if(abs(gen_Mom0ID[k]) != TAU) printf("WARNING: More than one mu_m\n\n");
                    if(PRINT) sprintf(out_buff + strlen(out_buff), "Extra mu_m detected\n");
                    print_out = true;
                }
            }
            if(gen_id[k] == -MY_LEP){
                if(gen_lep_p == -1) gen_lep_p = k;
                else{
                    if(abs(gen_Mom0ID[k]) != TAU) printf("WARNING: More than one mu_p\n\n");
                    if(PRINT) sprintf(out_buff + strlen(out_buff), "Extra mu_p detected\n");
                    print_out = true;
                }
            }
        }
        //record tau's
        if(abs(gen_id[k]) == TAU && 
                ( (gen_Mom0ID[k] == Z || gen_Mom0ID[k] == PHOTON || abs(gen_Mom0ID[k]) == ELECTRON ||
                   gen_status[k] == OUTGOING)) ){
            if(gen_id[k] == TAU){
                if(gen_tau_m == -1) gen_tau_m = k;
                else{
                    printf("WARNING: More than one tau_m\n\n");
                    if(PRINT) sprintf(out_buff + strlen(out_buff), "Extra mu_m detected\n");
                    //print_out = true;
                }
            }
            if(gen_id[k] == -TAU){
                if(gen_tau_p == -1) gen_tau_p = k;
                else{
                    printf("WARNING: More than one tau_p\n\n");
                    if(PRINT) sprintf(out_buff + strlen(out_buff), "Extra mu_p detected\n");
                    //print_out = true;
                }
            }
        }
        if(PRINT){
            if( (abs(gen_id[k]) <=6 || gen_id[k] == GLUON)){
                if(PRINT) sprintf(out_buff + strlen(out_buff),"Parton (ID = %i stat = %i): \n"
                        "    Mom1 ID: %i Mom1 Stat: %i Mom2 ID: %i Mom2 Stat %i \n"
                        "    Dau1 ID: %i Dau1 Stat: %i Dau2 ID: %i Dau2 Stat %i \n",
                        gen_id[k], gen_status[k], 
                        gen_Mom0ID[k], gen_Mom0Status[k], gen_Mom1ID[k], gen_Mom1Status[k],
                        gen_Dau0ID[k], gen_Dau0Status[k], gen_Dau1ID[k], gen_Dau1Status[k]);
            }


            if(gen_id[k] == -MY_LEP){
                if(PRINT) sprintf(out_buff + strlen(out_buff),"mu_p(stat = %i): \n"
                        "   Mom1 ID: %i Mom1 Stat: %i Mom2 ID: %i Mom2 Stat %i \n"
                        "   Dau1 ID: %i Dau1 Stat: %i Dau2 ID: %i Dau2 Stat %i \n",
                        gen_status[k], 
                        gen_Mom0ID[k], gen_Mom0Status[k], gen_Mom1ID[k], gen_Mom1Status[k],
                        gen_Dau0ID[k], gen_Dau0Status[k], gen_Dau1ID[k], gen_Dau1Status[k]);
                if(PRINT) sprintf(out_buff + strlen(out_buff),"(pt, eta, phi,E) = %4.2f %4.2f %4.2f %4.2f \n", 
                        gen_Pt[k], gen_Eta[k], gen_Phi[k], gen_E[k]);
            }
            if(gen_id[k] == MY_LEP){

                if(PRINT) sprintf(out_buff + strlen(out_buff),"mu_m (stat = %i): \n"
                        "   Mom1 ID: %i Mom1 Stat: %i Mom2 ID: %i Mom2 Stat %i \n"
                        "   Dau1 ID: %i Dau1 Stat: %i Dau2 ID: %i Dau2 Stat %i \n",
                        gen_status[k], 
                        gen_Mom0ID[k], gen_Mom0Status[k], gen_Mom1ID[k], gen_Mom1Status[k],
                        gen_Dau0ID[k], gen_Dau0Status[k], gen_Dau1ID[k], gen_Dau1Status[k]);
                if(PRINT) sprintf(out_buff + strlen(out_buff),"(pt, eta, phi,E) = %4.2f %4.2f %4.2f %4.2f \n", 
                        gen_Pt[k], gen_Eta[k], gen_Phi[k], gen_E[k]);

            }
            if(gen_id[k] == Z){
                if(PRINT) sprintf(out_buff + strlen(out_buff),"Z (ID = %i, status = %i): \n"
                        "   Mom1 ID: %i Mom1 Stat: %i Mom2 ID: %i Mom2 Stat %i \n"
                        "   Dau1 ID: %i Dau1 Stat: %i Dau2 ID: %i Dau2 Stat %i \n",
                        gen_id[k], gen_status[k],
                        gen_Mom0ID[k], gen_Mom0Status[k], gen_Mom1ID[k], gen_Mom1Status[k],
                        gen_Dau0ID[k], gen_Dau0Status[k], gen_Dau1ID[k], gen_Dau1Status[k]);
            }
        }
    }


    if(gen_lep_p != -1 && gen_lep_m != -1) {
        if(PRINT) sprintf(out_buff + strlen(out_buff),"mu_p: \n"
                "   Mom1 ID: %i Mom1 Stat: %i Mom2 ID: %i Mom2 Stat %i \n"
                "   Dau1 ID: %i Dau1 Stat: %i Dau2 ID: %i Dau2 Stat %i \n",
                gen_Mom0ID[gen_lep_p], gen_Mom0Status[gen_lep_p], gen_Mom1ID[gen_lep_p], 
                gen_Mom1Status[gen_lep_p],
                gen_Dau0ID[gen_lep_p], gen_Dau0Status[gen_lep_p], gen_Dau1ID[gen_lep_p], 
                gen_Dau1Status[gen_lep_p]);

        if(PRINT) sprintf(out_buff + strlen(out_buff),"mu_m: \n"
                "   Mom1 ID: %i Mom1 Stat: %i Mom2 ID: %i Mom2 Stat %i \n"
                "   Dau1 ID: %i Dau1 Stat: %i Dau2 ID: %i Dau2 Stat %i \n",
                gen_Mom0ID[gen_lep_m], gen_Mom0Status[gen_lep_m], gen_Mom1ID[gen_lep_m], 
                gen_Mom1Status[gen_lep_m],
                gen_Dau0ID[gen_lep_m], gen_Dau0Status[gen_lep_m], gen_Dau1ID[gen_lep_m], 
                gen_Dau1Status[gen_lep_m]);

    }
    else {
        printf("WARNING: Unable to identify MuMu pair in event %i, skipping \n", event_idx);
        nFailedID ++;
        print_out = true;
        if(PRINT && print_out){
            sprintf(out_buff + strlen(out_buff), "\n\n");
            fputs(out_buff, stdout);
            print_out = false;
        }
        return;
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
        return;
    }

    else{ 
        //printf("%i %i \n", inc_1, inc_2);
        int inc_id1 = gen_id[inc_1];
        int inc_id2 = gen_id[inc_2];
        if(abs(gen_Dau0ID[inc_1]) == TAU){
            if(gen_tau_p == -1 || gen_tau_m == -1){
                printf("Didn't record tau's :( \n");
                nFailedID ++;
                return;
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
    }
    else{
        gen_el_p_vec.SetPtEtaPhiE(gen_Pt[gen_lep_p], gen_Eta[gen_lep_p], gen_Phi[gen_lep_p], gen_E[gen_lep_p]);
        gen_el_m_vec.SetPtEtaPhiE(gen_Pt[gen_lep_m], gen_Eta[gen_lep_m], gen_Phi[gen_lep_m], gen_E[gen_lep_m]);
        gen_cm = gen_el_p_vec + gen_el_m_vec;
        gen_cost = get_cost(gen_el_p_vec, gen_el_m_vec, false);
    }
    if(quark_dir_eta < 0){
        cost_st = -gen_cost;
    }
    else cost_st = gen_cost;

    gen_m = gen_cm.M();

    if(PRINT) memset(out_buff, 0, 10000);
}

void NTupleReader::finish(){

    fout->cd();

    printf("Writing output to file at %s \n", fout->GetName());
    for(int i=0; i<nOutTrees; i++){
        outTrees[i]->Write();
    }

    fout->Close();
    fclose(root_files);
}




NTupleReader::~NTupleReader(){

}
