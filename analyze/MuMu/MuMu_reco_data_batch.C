#include "../../Utils/NTupleReader.C"




void MuMu_reco_data_batch(int nJobs =1, int iJob = 0, string fin="", bool do_ss=false)
{

    if(fin == "") fin = string("EOS_files/2017/SingleMuon_files.txt");
    NTupleReader nt(fin.c_str(),"output_files/MuMu_data_test.root", true);
    nt.year = 2017;
    nt.do_samesign = false;

    nt.nJobs = nJobs;
    nt.iJob = iJob;
    nt.do_muons = true;
    nt.do_RC = true;
    nt.setupRC();

    nt.setupOutputTree("T_sig");
    nt.setupOutputTree("T_WJets");
    nt.setupOutputTree("T_QCD");
    nt.setupOutputTree("T_ss");

    int iso_mu;
    nt.outTrees[1]->Branch("iso_mu", &iso_mu); 

    while(nt.getNextFile()){
        for (int i=0; i<nt.tin_nEntries; i++) {
            nt.getEvent(i);

            if(nt.good_trigger && nt.dimuon_id && nt.cm_m > 130.){
                nt.fillEvent();
                nt.fillEventRC();
                bool one_iso = nt.mu_iso0 ^ nt.mu_iso1;

                //pick the category
                if(nt.opp_sign && nt.mu_iso0 && nt.mu_iso1){ //signal region
                    nt.outTrees[0]->Fill();
                }
                else if(!nt.opp_sign && nt.mu_iso0 && nt.mu_iso1){ //samesign region
                    nt.outTrees[3]->Fill();
                }
                else if(one_iso){ //wjets control region
                    if(nt.mu_iso0) iso_mu = 0;
                    else           iso_mu = 1;
                    nt.outTrees[1]->Fill();
                }
                else if(!nt.mu_iso0 && !nt.mu_iso1){ //qcd control region
                    nt.outTrees[2]->Fill();
                }
            }

            
        } 

        printf("moving on to next file, currently %i events \n\n", nt.nEvents);


    }
    nt.finish();

    return;
}

int main(){
    MuMu_reco_data_batch();
    return 0;
}
