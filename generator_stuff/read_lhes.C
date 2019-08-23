#include "madgraph_lhe_reader.C"
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <iostream>
#include <cstring>
#include <algorithm>
#include "TFile.h"
#include "TLorentzVector.h"


void read_lhes(){
    TTree *t1 = new TTree("T_lhe", "Lhe event info for mass binned DY");
    t1->SetDirectory(0);
    //string f1("/home/ozamram/Documents/Research/Generators/POWHEG-BOX-V2/Z_ew-BMNNPV/DY_m150_may6/pwgevents.lhe");
    //string f1("/home/ozamram/Documents/Research/B2GTTrees/generator_stuff/DY_M700_events.lhe");
    //string f1("/uscms_data/d3/oamram/CMSSW_8_0_24_patch1/src/Analysis/B2GTTrees/generator_stuff/MG5/m200/cmsgrid_final.lhe");
    string f1("/home/ozamram/Documents/Research/B2GTTrees/generator_stuff/phot_induced_aug21.lhe");
    //string f2("/uscms_data/d3/oamram/CMSSW_8_0_24_patch1/src/Analysis/B2GTTrees/generator_stuff/m200/cmsgrid_final.lhe");
    //string f3("/uscms_data/d3/oamram/CMSSW_8_0_24_patch1/src/Analysis/B2GTTrees/generator_stuff/mass_binned/cmsgrid_final.lhe");
    fill_tree(f1, t1, true);
    //fill_tree(f2, t1, true);
    //fill_tree(f3, t1, true);

    TFile *fout1 = TFile::Open("root_files/phot_induced_aug21.root", "RECREATE");
    fout1->cd();
    t1->Write();
    fout1->Close();
    delete t1;




    /*
    TTree *t2 = new TTree("T_lhe", "Lhe event info for mass unbinned DY");
    t2->SetDirectory(0);
    string f2("mass_binned_100k.lhe");
    fill_tree(f2, t2, true);

    TFile *fout2 = TFile::Open("mass_binned_100k.root", "RECREATE");
    fout2->cd();
    t2->Write();
    fout2->Close();
    delete t2;
    */

    /*
    char base_str[120] = "/uscms_data/d3/oamram/condor_jobs/DY_jobs/condor_jobs_binned_SEED_%i/cmsgrid_final.lhe";
    char root_base[120] = "condor_files/mass_binned_100k_%i.root";
    char root_file[120];
    char file_str[120];
    for(int i=1; i<= 2; i++){
        sprintf(file_str, base_str, i);
        sprintf(root_file, root_base, i);
        string f(file_str);
        TTree *t1 = new TTree("T_lhe", "Lhe event info for mass unbinned DY");
        t1->SetDirectory(0);
        fill_tree(f,t1, true);
        TFile *fout = TFile::Open(root_file, "RECREATE");
        fout->cd();
        t1->Write();
        fout->Close();
        delete t1;
    }
    */


    return;
}

