#include "../../Utils/NTupleReader.C"




void MuMu_WJets_fake_rate_estimate(int nJobs =1, int iJob = 0)
{


    NTupleReader nt("EOS_files/SingleMuon_files_test.txt","output_files/test.root", true);
    nt.nJobs = nJobs;
    nt.iJob = iJob;
    nt.do_muons = true;
    nt.setupOutputTree("T_data");

    int iso_mu;
    nt.outTrees[0]->Branch("iso_mu", &iso_mu); 


    while(nt.getNextFile()){


        for (int i=0; i<nt.tin_nEntries; i++) {
            nt.getEvent(i);
            bool one_iso = nt.mu_iso0 ^ nt.mu_iso1;
            if(nt.good_trigger && nt.dimuon_id && one_iso && nt.cm_m > 15. ){
                nt.fillEvent();
                if(nt.mu_iso0) iso_mu = 0;
                else           iso_mu = 1;
                nt.outTrees[0]->Fill();


            }
        } 


        printf("moving on to next file, currently %i events \n\n", nt.nEvents);

    }
    nt.finish();

    return;
}

