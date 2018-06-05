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
    string f1("mass_binned.lhe");
    fill_tree(f1, t1, true);

    TFile *fout1 = TFile::Open("mass_binned_evts.root", "RECREATE");
    fout1->cd();
    t1->Write();
    fout1->Close();
    delete t1;




    TTree *t2 = new TTree("T_lhe", "Lhe event info for mass unbinned DY");
    t2->SetDirectory(0);
    string f2("mass_unbinned.lhe");
    fill_tree(f2, t2, true);

    TFile *fout2 = TFile::Open("mass_unbinned_evts.root", "RECREATE");
    fout2->cd();
    t2->Write();
    fout2->Close();
    delete t2;

    return;
}

