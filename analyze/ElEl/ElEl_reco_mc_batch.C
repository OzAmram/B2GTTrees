#include "../../Utils/NTupleReader.C"




void ElEl_reco_mc_batch(int nJobs =1, int iJob = 0, string fin = "", bool do_ss = false)
{


    if(fin == "") fin = string("EOS_files/2016/DY_files_test.txt");
    NTupleReader nt(fin.c_str(),"output_files/ElEl_dy_test.root", false);
    nt.year = 2016;
    nt.do_samesign = do_ss;


    nt.nJobs = nJobs;
    nt.iJob = iJob;
    nt.do_electrons = true;
    nt.do_SFs = true;
    nt.setupSFs();
    nt.setupOutputTree("T_sig");
    nt.setupOutputTree("T_back");


    while(nt.getNextFile()){


        for (int i=0; i<nt.tin_nEntries; i++) {
            nt.getEvent(i);
            if(nt.good_trigger && nt.good_sign && nt.dielec_id &&
                    nt.el_iso0 && nt.el_iso1 && nt.cm_m > 130. ){
                nt.fillEvent();
                nt.fillEventSFs();
                nt.parseGenParts(false);

                if(nt.signal_event && !nt.failed_match){
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
    printf("There were %i qqbar, %i qGlu (%i of them tautau) in %i kept events in %i files. "
            "There were also %i background events (%i qq and %i gg). "
            "There were %i Failed ID's \n" , 
            nt.nQQb, nt.nQGlu, nt.nTauTau, nt.nSignal, nt.fileCount, nt.nQQ + nt.nGluGlu, nt.nQQ, nt.nGluGlu, nt.nFailedID);
    //printf("Ran on MC data and produced templates with %i events\n", nEvents);
    nt.finish();

    return;
}

