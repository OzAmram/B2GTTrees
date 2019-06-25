#include "../../Utils/NTupleReader.C"




void MuMu_WJets_fake_rate_estimate(int nJobs =1, int iJob = 0, string fin="")
{


    if(fin == "") fin = string("EOS_files/2016/SingleMuon_files_may31.txt");
    NTupleReader nt(fin.c_str(),"output_files/MuMu_data_test.root", true);
    nt.year = 2016;
    nt.nJobs = nJobs;
    nt.iJob = iJob;
    nt.do_muons = true;
    nt.do_RC = true;
    nt.setupOutputTree("T_data");

    int iso_mu;
    nt.outTrees[0]->Branch("iso_mu", &iso_mu); 


    while(nt.getNextFile()){


        for (int i=0; i<nt.tin_nEntries; i++) {
            nt.getEvent(i);
            bool one_iso = nt.mu_iso0 ^ nt.mu_iso1;
            if(nt.good_trigger && nt.dimuon_id && one_iso && nt.cm_m > 130. ){
                nt.fillEvent();
                if(nt.mu_iso0) iso_mu = 0;
                else           iso_mu = 1;
                nt.outTrees[0]->Fill();


            }
        } 


        printf("moving on to next file, currently %i events \n\n", nt.nEvents);

    }
    nt.finish();

    return;
}

