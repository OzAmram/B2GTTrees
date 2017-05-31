//perform fits to Reconstructed MuMu data to extract Asym

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


void draw1(){
    TFile *f_data = TFile::Open("output_files/DYToLL_data_feb8.root");
    TH1F *h_cost = new TH1F("cost", "Cos(#theta) for M in[150,200] (data) (Gev);cos(#theta)", 50, -1,1);
    TH1F *h_xf = new TH1F("xF", "Feynman X for M in [150,200] (data); xF", 20,-1,1);

    //read event data
    TTree *t1 = (TTree *)f_data->Get("T_data"); //3d histogram of data
    Long64_t size  =  t1->GetEntries();
    Double_t m, xF, cost;
    t1->SetBranchAddress("m", &m);
    t1->SetBranchAddress("xF", &xF);
    t1->SetBranchAddress("cost", &cost);
    int nBack = 0;
    int nEvents = 0;
    int nForward = 0;
    int nF = 0;
    int nB=0;
    for (int i=0; i<size; i++) {
        t1->GetEntry(i);
        if(m >= 150 && m<=200){
            nEvents++;
            h_cost->Fill(cost);
            h_xf->Fill(xF);
            if(cost <0) nBack++;
            if(cost >0) nForward++;
            if(cost > 0 && xF > 0.1) nF++;
            if(cost < 0 && xF>0.1) nB++;
        }
    }
    float AFB = float(nForward - nBack)/float(nForward + nBack);
    float AFB_highx = float(nF - nB)/float(nF + nB);
    printf("high x count %i, AFB %0.3f \n", nF+nB, AFB_highx);

    printf("Running through %i events. Counting AFB is %0.3f \n", 
            nEvents, AFB);
    printf("Running fit on %i events. There are %i back and %i forward \n", nEvents, nBack, nForward);

    TCanvas *c_cost = new TCanvas("c_cost", "Histograms", 200, 10, 900, 700);
    h_cost->Draw();
    h_cost->SetFillColor(35);
    c_cost->Update();

    TCanvas *c_xF = new TCanvas("c_xF", "Histograms", 200, 10, 900, 700);
    h_xf->Draw();
    h_xf->SetFillColor(35);
    c_xF->Update();
    

 
}

    
    
