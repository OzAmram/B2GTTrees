

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
    Double_t x[n] = {175,225,300,425,800,600};
    Double_t y[n] = {0.613, 0.557,0.637,0.626,0.614,0.583};
    Double_t x_errs[n] = {25,25,50,75,100,100};
    Double_t y_errs[n] = {0.021, 0.03,0.031,0.045,0.059,0.112};
        
    TCanvas *c1 = new TCanvas("c1", "AFB as function of mass", 200,10, 600,800);
    TGraphErrors *g1 = new TGraphErrors(n,x,y,x_errs,y_errs);
    g1->SetMinimum(0);
    g1->GetYaxis()->SetTitle("AFB");
    g1->GetXaxis()->SetTitle("DiMuon Mass (GeV)");

    g1->SetTitle("Forward Backward Asymmetry");
    g1->Draw("AP");
    return;
}

