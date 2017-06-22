

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

#include "tdrstyle.C"
#include "CMS_lumi.C"


void draw_AFB(){
    setTDRStyle();
    writeExtraText = true;       // if extra text
    extraText  = "Preliminary";  // default extra text is "Preliminary"
    lumi_sqrtS = "13 TeV";

    const int n = 6;
    Double_t x[n] = {175,225,300,425,600,900};
    Double_t y[n] = {0.653, 0.560,0.621,0.633,0.582,0.594};
    Double_t x_errs[n] = {25,25,50,75,100,200};
    Double_t y_errs[n] = {0.018, 0.027,0.028,0.037,0.059,0.091};
        
    TCanvas *c1 = new TCanvas("c1", "AFB as function of mass", 200,10, 800,600);
    TGraphErrors *g1 = new TGraphErrors(n, x,y,x_errs,y_errs);
    //g1->SetTitle("Drell-Yan Forward Backward Asymmetry");
    g1->SetMinimum(0);

    g1->SetMarkerStyle(20);
    //g1->GetYaxis()->SetTitle("A_{FB}");
    g1->GetYaxis()->SetTitle("Forward Backward Asymmetry");
    g1->GetXaxis()->SetTitle("DiMuon Mass (GeV)");
    g1->Draw("A  P");
    //g1->SetTitle("A_{FB} as a function of mass");
    c1->Update();

    writeExtraText = true;
    extraText = "Preliminary";
    //lumi_sqrtS = "";       // used with iPeriod = 0, e.g. for simulation-only plots (default is an empty string)
    int iPeriod = 4; 
    CMS_lumi( c1, iPeriod, 0 );
}

