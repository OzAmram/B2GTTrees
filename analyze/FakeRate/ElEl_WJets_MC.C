
#include "../../Utils/NTupleReader.C"




void ElEl_WJets_MC(int nJobs =1, int iJob = 0)
{


    NTupleReader nt("EOS_files/2016/non_QCD_files_may29.txt","output_files/test.root", false);
    nt.year= 2016;

    nt.nJobs = nJobs;
    nt.iJob = iJob;
    nt.do_electrons = true;
    nt.do_SFs = true;
    nt.setupSFs();
    nt.setupOutputTree("T_data");

    int iso_el;
    nt.outTrees[0]->Branch("iso_el", &iso_el); 


    while(nt.getNextFile()){


        for (int i=0; i<nt.tin_nEntries; i++) {
            nt.getEvent(i);
            bool one_iso = nt.el_iso0 ^ nt.el_iso1;
            if(nt.good_trigger && nt.dielec_id && one_iso && nt.cm_m > 15. ){
                nt.fillEvent();
                nt.fillEventSFs();
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

