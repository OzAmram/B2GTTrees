#include "../../Utils/NTupleReader.C"




void EMu_background_check_Mu(int nJobs =1, int iJob = 0, string fin = "")
{

    if(fin == "") fin = string("EOS_files/2016/non_QCD_files_may29.txt");
    NTupleReader nt(fin.c_str() ,"output_files/test.root", false);
    nt.nJobs = nJobs;
    nt.iJob = iJob;
    nt.do_emu = true;
    nt.do_SFs = true;
    nt.setupSFs();
    nt.setupOutputTree("T_data");


    while(nt.getNextFile()){


        for (int i=0; i<nt.tin_nEntries; i++) {
            nt.getEvent(i);
            if(nt.good_trigger && nt.emu_ids && nt.mu_iso0 && nt.el_iso0 && nt.cm_m > 50. ){
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

