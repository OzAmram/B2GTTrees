#include "../../Utils/NTupleReader.C"




void EMu_WJets_fakerate_est(int nJobs =1, int iJob = 0, string fin = "")
{


    if (fin == "") fin = string("EOS_files/2016/SingleMuon_files_may31.txt");
    NTupleReader nt(fin.c_str(),"output_files/test.root", true);
    nt.year = 2016;
    nt.nJobs = nJobs;
    nt.iJob = iJob;
    nt.do_emu = true;
    nt.setupOutputTree("T_data");
    int iso_lep;
    nt.outTrees[0]->Branch("iso_lep", &iso_lep);


    while(nt.getNextFile()){


        for (int i=0; i<nt.tin_nEntries; i++) {
            nt.getEvent(i);
            bool one_iso = nt.mu_iso0 ^ nt.el_iso0;
            if(nt.good_trigger && nt.emu_ids && one_iso && nt.cm_m > 50.){
                nt.fillEvent();
                if(nt.mu_iso0) iso_lep = 0;
                else          iso_lep = 1;
                nt.outTrees[0]->Fill();


            }
        } 


        printf("moving on to next file, currently %i events \n\n", nt.nEvents);

    }
    nt.finish();

    return;
}

