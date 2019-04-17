#include "../../Utils/NTupleReader.C"




void ElEl_reco_background_batch(int nJobs =1, int iJob = 0)
{


    NTupleReader nt("EOS_files/DY_files_test.txt","output_files/ElEl_test.root", false);
    nt.nJobs = nJobs;
    nt.iJob = iJob;
    nt.do_electrons = true;
    nt.do_SFs = true;
    nt.setupSFs();
    nt.setupOutputTree("T_data");
    nt.setupOutputTree("T_back");


    while(nt.getNextFile()){


        for (int i=0; i<nt.tin_nEntries; i++) {
            nt.getEvent(i);
            if(nt.good_trigger && nt.good_sign && nt.dielec_id &&
                    nt.el_iso0 && nt.el_iso1 && nt.cm_m > 50. ){
                nt.fillEvent();
                nt.fillEventSFs();
                nt.outTrees[0]->Fill();


            }
        } 

        printf("moving on to next file, currently %i events \n\n", nt.nEvents);


    }
    printf("Finished. Selected %i events from %i files \n", nt.nEvents, nt.fileCount);
    nt.finish();

    return;
}

