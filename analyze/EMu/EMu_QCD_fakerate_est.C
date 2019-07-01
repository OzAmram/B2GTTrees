#include "../../Utils/NTupleReader.C"




void EMu_QCD_fakerate_est(int nJobs =1, int iJob = 0, string fin ="")
{


    if (fin == "") fin = string("EOS_files/2016/SingleMuon_files_may31.txt");
    NTupleReader nt(fin.c_str(),"output_files/test.root", true);
    nt.year = 2016;
    nt.nJobs = nJobs;
    nt.iJob = iJob;
    nt.do_emu = true;
    nt.setupOutputTree("T_data");


    while(nt.getNextFile()){


        for (int i=0; i<nt.tin_nEntries; i++) {
            nt.getEvent(i);
            if(nt.good_trigger && nt.emu_ids && !nt.mu_iso0 && !nt.el_iso0 && nt.cm_m > 50.){
                nt.fillEvent();
                nt.outTrees[0]->Fill();


            }
        } 


        printf("moving on to next file, currently %i events \n\n", nt.nEvents);

    }
    nt.finish();

    return;
}

