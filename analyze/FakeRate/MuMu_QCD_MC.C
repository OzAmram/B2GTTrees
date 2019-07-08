#include "../../Utils/NTupleReader.C"




void MuMu_QCD_MC(int nJobs =1, int iJob = 0)
{


    NTupleReader nt("EOS_files/2016/non_QCD_files_may29.txt","output_files/test.root", false);
    nt.year= 2016;

    nt.nJobs = nJobs;
    nt.iJob = iJob;
    nt.do_muons = true;
    nt.do_SFs = true;
    nt.setupSFs();
    nt.setupOutputTree("T_data");



    while(nt.getNextFile()){


        for (int i=0; i<nt.tin_nEntries; i++) {
            nt.getEvent(i);
            if(nt.good_trigger && nt.dimuon_id && !nt.mu_iso0 && !nt.mu_iso1 && nt.cm_m > 15. ){
                nt.fillEvent();
                nt.fillEventSFs();
                nt.outTrees[0]->Fill();


            }
        } 


        printf("moving on to next file, currently %i events \n\n", nt.nEvents);

    }
    nt.finish();

    return;
}

