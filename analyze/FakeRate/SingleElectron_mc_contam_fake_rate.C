#include "../../Utils/NTupleReader.C"


bool in_Z_window(Double_t m){
    double_t Z_mass_low = 91.2 -7.;
    double_t Z_mass_high = 91.2 + 7.;
    return (m > Z_mass_low ) && (m < Z_mass_high);
}


void SingleElectron_mc_contam_fake_rate(int nJobs =1, int iJob=0, string fin="")
{


    if(fin == "") fin = string("EOS_files/2017/diboson_files_may29.txt");
    NTupleReader nt(fin.c_str(),"output_files/ElEl_fake_rate_contam_test.root", false);
    nt.year = 2017;
    nt.nJobs = nJobs;
    nt.iJob = iJob;
    nt.do_electrons =  true;
    nt.do_SFs = false;

    TTree *tout= new TTree("T_data", "Tree with reco events");
    nt.nOutTrees++;
    nt.outTrees[0] = tout;
    tout->SetDirectory(0);
    Double_t el_pt, el_eta, gen_weight;
    Bool_t pass;
    tout->Branch("el_pt", &el_pt);
    tout->Branch("el_eta", &el_eta);
    tout->Branch("pass", &pass);
    tout->Branch("gen_weight", &gen_weight);
    int nEvents =0;

    while(nt.getNextFile()){
        for (int i=0; i<nt.tin_nEntries; i++) {
            nt.getEvent(i);

            if(nt.good_trigger &&  nt.dielec_id && nt.el_size >= 3
                    && nt.el_ScaleCorr[2] * nt.el_Pt[2] > 15. && goodElEta(nt.el_SCEta[2]) &&  nt.el_IDMedium_NoIso[2]){

                //Want events with 3 electrons, 2 from Z and 1 extra

                nt.fillEvent();
                gen_weight = nt.gen_weight;
                const float el_mass = 0.; // in GEV

                TLorentzVector el0, el1, el2;
                el0.SetPtEtaPhiM(nt.el_Pt[0], nt.el_Eta[0], nt.el_Phi[0], el_mass);
                el1.SetPtEtaPhiM(nt.el_Pt[1], nt.el_Eta[1], nt.el_Phi[1], el_mass);
                el2.SetPtEtaPhiM(nt.el_Pt[2], nt.el_Eta[2], nt.el_Phi[2], el_mass);

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

                bool iso[3];
                iso[0] = nt.el_IDMedium[0];
                iso[1] = nt.el_IDMedium[1];
                iso[2] = nt.el_IDMedium[2];

                if(m01_in_Z && !m02_in_Z && !m12_in_Z && nt.el_Charge[0] * nt.el_Charge[1] < 0 && iso[0] && iso[1]){
                    el_extra = 2;
                    if(nt.el_Charge[0] >0){
                        el_p = 0;
                        el_m = 1;
                    }
                    else{
                        el_p = 1;
                        el_m =0;
                    }
                }
                else if(!m01_in_Z && m02_in_Z && !m12_in_Z && nt.el_Charge[0] * nt.el_Charge[2] < 0 && iso[0] && iso[2] ){
                    el_extra = 1;
                    if(nt.el_Charge[0] >0){
                        el_p = 0;
                        el_m = 2;
                    }
                    else{
                        el_p = 2;
                        el_m =0;
                    }
                }
                else if(!m01_in_Z && !m02_in_Z && m12_in_Z  && nt.el_Charge[1] * nt.el_Charge[2] < 0 && iso[1] && iso[2]){
                    el_extra = 0;
                    if(nt.el_Charge[1] >0){
                        el_p = 1;
                        el_m = 2;
                    }
                    else{
                        el_p = 2;
                        el_m =2;
                    }
                }


                if( (el_p != -1) && nt.has_nobjets && nt.met_pt < 25.){
                    el_pt = nt.el_Pt[el_extra];
                    el_eta = nt.el_Eta[el_extra];
                    pass = iso[el_extra];
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
