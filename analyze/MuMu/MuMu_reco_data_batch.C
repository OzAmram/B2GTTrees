#include "../../Utils/NTupleReader.C"




void MuMu_reco_data_batch(int nJobs =1, int iJob = 0)
{


    NTupleReader nt("EOS_files/SingleMuon_files_test.txt","output_files/MuMu_data_test.root", true);
    nt.nJobs = nJobs;
    nt.iJob = iJob;
    nt.do_muons = true;
    nt.do_SFs = false;
    nt.use_RC = true;
    nt.setupRC();
    nt.setupOutputTree("T_data");


    while(nt.getNextFile()){


        for (int i=0; i<nt.tin_nEntries; i++) {
            nt.getEvent(i);
            if(nt.good_trigger && nt.good_sign && nt.dimuon_id &&
                    nt.mu_iso0 && nt.mu_iso1 && nt.cm_m > 50. ){
                nt.fillEvent();
                nt.fillEventRC();
                nt.outTrees[0]->Fill();

            }
        } 

        printf("moving on to next file, currently %i events \n\n", nt.nEvents);


    }
    printf("Finished. There were %i events in %i files \n", nt.nEvents, nt.fileCount);
    nt.finish();

    return;
}

int main(){
    MuMu_reco_data_batch();
    return 0;
}
