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
    /*
    TTree *t1 = new TTree("T_lhe", "Lhe event info for mass binned DY");
    t1->SetDirectory(0);
    string f1("mass_binned.lhe");
    fill_tree(f1, t1, true);

    TFile *fout1 = TFile::Open("mass_binned_evts.root", "RECREATE");
    fout1->cd();
    t1->Write();
    fout1->Close();
    delete t1;
    */




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

    char base_str[80] = "/uscms_data/d3/oamram/condor_jobs/DY_jobs/condor_jobs_SEED_%i/cmsgrid_final.lhe";
    char root_base[80] = "condor_files/mass_unbinned_100k_%i.root";
    char root_file[80];
    char file_str[80];
    for(int i=1; i<= 20; i++){
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


    return;
}

