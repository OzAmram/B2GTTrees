#include "../../Utils/NTupleReader.C"




void MuMu_WJets_MC(int nJobs =1, int iJob = 0)
{


    NTupleReader nt("EOS_files/DY_files_test.txt","output_files/test.root", false);
    nt.nJobs = nJobs;
    nt.iJob = iJob;
    nt.do_muons = true;
    nt.do_SFs = true;
    nt.setupSFs();
    nt.setupOutputTree("T_data");

    int iso_mu;
    nt.outTrees[0]->Branch("iso_mu", &iso_mu); 


    while(nt.getNextFile()){


        for (int i=0; i<nt.tin_nEntries; i++) {
            nt.getEvent(i);
            bool one_iso = nt.mu_iso0 ^ nt.mu_iso1;
            if(nt.good_trigger && nt.dimuon_id && one_iso && nt.cm_m > 15. ){
                nt.fillEvent();
                nt.fillEventSFs();
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

