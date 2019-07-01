#include "../../Utils/NTupleReader.C"




void ElEl_reco_background_batch(int nJobs =1, int iJob = 0, string fin ="")
{


    if(fin == "") fin = string("EOS_files/2016/combined_back_files_may29.txt");
    NTupleReader nt(fin.c_str(),"output_files/ElEl_back_june25.root", false);
    nt.year = 2016;
    nt.nJobs = nJobs;
    nt.iJob = iJob;
    nt.do_electrons = true;
    nt.do_SFs = true;
    nt.setupSFs();
    nt.setupOutputTree("T_data");


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

