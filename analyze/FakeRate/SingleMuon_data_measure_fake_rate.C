#include "../../Utils/NTupleReader.C"


bool in_Z_window(Double_t m){
    double_t Z_mass_low = 91.2 -7.;
    double_t Z_mass_high = 91.2 + 7.;
    return (m > Z_mass_low ) && (m < Z_mass_high);
}


void SingleMuon_data_measure_fake_rate(int nJobs =1, int iJob=0, string fin="")
{


    if(fin == "") fin = string("EOS_files/2017/SingleMuon_files.txt");
    NTupleReader nt(fin.c_str(),"output_files/MuMu_fake_rate_test.root", true);
    nt.year = 2017;
    nt.nJobs = nJobs;
    nt.iJob = iJob;
    nt.do_muons = true;
    nt.do_SFs = false;
    nt.do_RC = true;

    TTree *tout= new TTree("T_data", "Tree with reco events");
    nt.nOutTrees++;
    nt.outTrees[0] = tout;
    tout->SetDirectory(0);
    Double_t mu_pt, mu_eta;
    Bool_t pass;
    tout->Branch("mu_pt", &mu_pt);
    tout->Branch("mu_eta", &mu_eta);
    tout->Branch("pass", &pass);
    int nEvents =0;

    while(nt.getNextFile()){
        for (int i=0; i<nt.tin_nEntries; i++) {
            nt.getEvent(i);

            if(nt.good_trigger &&  nt.dimuon_id && nt.mu_size >= 3
                    && nt.mu_Pt[2] > 15. && abs(nt.mu_Eta[2]) < 2.4 &&  nt.mu_IsTightMuon[2]){

                //Want events with 3 muons, 2 from Z and 1 extra

                nt.fillEvent();
                const float mu_mass = 0.1056; // in GEV

                TLorentzVector mu0, mu1, mu2;
                mu0.SetPtEtaPhiM(nt.mu_Pt[0], nt.mu_Eta[0], nt.mu_Phi[0], mu_mass);
                mu1.SetPtEtaPhiM(nt.mu_Pt[1], nt.mu_Eta[1], nt.mu_Phi[1], mu_mass);
                mu2.SetPtEtaPhiM(nt.mu_Pt[2], nt.mu_Eta[2], nt.mu_Phi[2], mu_mass);

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

                bool iso[3];
                iso[0] = nt.mu_PFIso[0] < nt.mu_iso_cut;
                iso[1] = nt.mu_PFIso[1] < nt.mu_iso_cut;
                iso[2] = nt.mu_PFIso[2] < nt.mu_iso_cut;

                if(m01_in_Z && !m02_in_Z && !m12_in_Z && nt.mu_Charge[0] * nt.mu_Charge[1] < 0 && iso[0] && iso[1]){
                    mu_extra = 2;
                    if(nt.mu_Charge[0] >0){
                        mu_p = 0;
                        mu_m = 1;
                    }
                    else{
                        mu_p = 1;
                        mu_m =0;
                    }
                }
                else if(!m01_in_Z && m02_in_Z && !m12_in_Z && nt.mu_Charge[0] * nt.mu_Charge[2] < 0 && iso[0] && iso[2] ){
                    mu_extra = 1;
                    if(nt.mu_Charge[0] >0){
                        mu_p = 0;
                        mu_m = 2;
                    }
                    else{
                        mu_p = 2;
                        mu_m =0;
                    }
                }
                else if(!m01_in_Z && !m02_in_Z && m12_in_Z  && nt.mu_Charge[1] * nt.mu_Charge[2] < 0 && iso[1] && iso[2]){
                    mu_extra = 0;
                    if(nt.mu_Charge[1] >0){
                        mu_p = 1;
                        mu_m = 2;
                    }
                    else{
                        mu_p = 2;
                        mu_m =2;
                    }
                }


                if( (mu_p != -1) && nt.has_nobjets && nt.met_pt < 25.){
                    mu_pt = nt.mu_Pt[mu_extra];
                    mu_eta = nt.mu_Eta[mu_extra];
                    pass = iso[mu_extra];
                    nEvents++;
                    tout->Fill();
                }

            }
        }
        printf("moving on to next file, currently %i events \n\n", nEvents);
    }

    nt.finish();
    return;
}
