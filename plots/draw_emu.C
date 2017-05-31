#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <iostream>
#include <cstring>
#include <algorithm>
#include <vector>

#include "TROOT.h"
#include "TFile.h"
#include "TObject.h"
#include "TH1F.h"
#include "TProfile.h"
#include "TStyle.h"
#include "TCanvas.h"
#include "TPostScript.h"
#include "TLorentzVector.h"
#include "TVector3.h"
#include "TFitter.h"
#include "TSystem.h"
#include "Math/Functor.h"



Double_t count_tree(TTree *t1,  bool is_data=false){
    //count events in the tree
    Long64_t size  =  t1->GetEntries();
    Double_t m, xF, cost, mu1_pt, mu2_pt, jet1_cmva, jet2_cmva, gen_weight;
    Float_t cost_pt, met_pt;
    t1->SetBranchAddress("met_pt", &met_pt);
    t1->SetBranchAddress("jet2_CMVA", &jet2_cmva);
    t1->SetBranchAddress("jet1_CMVA", &jet1_cmva);
    if(!is_data) t1->SetBranchAddress("gen_weight", &gen_weight);
    jet1_cmva = -1.;
    jet2_cmva = -1.;
    Double_t count = 0;

    for (int i=0; i<size; i++) {
        t1->GetEntry(i);
        //if(met_pt < 50. && jet1_cmva < 0. && jet2_cmva < 0){
            if(is_data){
                count += 1;
            }
            else{
                count += gen_weight;
            }

            
        //}
    }
    printf("%f \n", count);
    return count;
}


void draw_emu(){
    TFile *f_data = TFile::Open("../analyze/output_files/EMu_data_May18.root");
    TTree *t_data = (TTree *)f_data->Get("T_data");

                                
    TFile *f_back = TFile::Open("../analyze/output_files/ttbar_EMu_may18.root");
    TTree *t_back = (TTree *)f_back->Get("T_data");

    Double_t data_count = count_tree(t_data, true);
    Double_t ttbar_count = count_tree(t_back);
    Double_t lumi = 35.867;
    ttbar_count *= lumi * 1000;

    printf("Data count %e \n", data_count);
    printf("TTbar count %e \n", ttbar_count);
    return;
}









