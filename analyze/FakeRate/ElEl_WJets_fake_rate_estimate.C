
#include "../../Utils/NTupleReader.C"




void ElEl_WJets_fake_rate_estimate(int nJobs =1, int iJob = 0, string fin="")
{


    if (fin == "") fin = string("EOS_files/2016/SingleElectron_files_may31.txt");
    NTupleReader nt(fin.c_str(),"output_files/SingleElectron_files_test.root", true);
    nt.year = 2016;
    nt.nJobs = nJobs;
    nt.iJob = iJob;
    nt.do_electrons = true;
    nt.setupOutputTree("T_data");

    int iso_el;
    nt.outTrees[0]->Branch("iso_el", &iso_el); 


    while(nt.getNextFile()){


        for (int i=0; i<nt.tin_nEntries; i++) {
            nt.getEvent(i);
            bool one_iso = nt.el_iso0 ^ nt.el_iso1;
            if(nt.good_trigger && nt.dielec_id && one_iso && nt.cm_m > 130. ){
                nt.fillEvent();
                if(nt.el_iso0) iso_el = 0;
                else           iso_el = 1;
                nt.outTrees[0]->Fill();


            }
        } 


        printf("moving on to next file, currently %i events \n\n", nt.nEvents);

    }
    nt.finish();

    return;
}

