#include "../../Utils/NTupleReader.C"




void MuMu_reco_mc_batch(int nJobs =1, int iJob = 0)
{


    NTupleReader nt("EOS_files/test_file.txt","output_files/MuMu_test.root", false);
    nt.do_muons = true;
    nt.do_SFs = true;
    nt.use_RC = true;
    nt.RC_from_gen = true;
    nt.setupSFs();
    nt.setupRC();
    nt.setupOutputTree("T_data");
    nt.setupOutputTree("T_back");


    while(nt.getNextFile()){


        for (int i=0; i<nt.tin_nEntries; i++) {
            nt.getEvent(i);
            if(nt.good_trigger && nt.good_sign && nt.dimuon_id &&
                    nt.mu_iso0 && nt.mu_iso1 && nt.cm_m > 50. ){
                //printf("HI THERE %i \n", i);
                nt.fillEvent();
                nt.fillEventSFs();
                nt.parseGenParts();
                nt.fillEventRC();

                if(nt.signal_event){
                    nt.nSignal++;
                    nt.outTrees[0]->Fill();
                }
                else{
                    nt.outTrees[1]->Fill();
                }


            }
        } 

        printf("moving on to next file, currently %i events %i Taus %i fails \n\n", nt.nEvents, nt.nTauTau, nt.nFailedID);


    }
    printf("There were %i qqbar, %i qGlu (%i of them tautau) in %i kept events in %i files."
            "There were also %i background events (%i qq and %i gg)"
            "There were %i Failed ID's \n" , 
            nt.nQQb, nt.nQGlu, nt.nTauTau, nt.nSignal, nt.fileCount, nt.nQQ + nt.nGluGlu, nt.nQQ, nt.nGluGlu, nt.nFailedID);
    //printf("Ran on MC data and produced templates with %i events\n", nEvents);
    nt.finish();

    return;
}

int main(){
    MuMu_reco_mc_batch();
    return 0;
}
