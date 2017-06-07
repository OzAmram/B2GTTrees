
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
#include "../analyze/TemplateMaker.C"


void draw_m_cost(){
    TFile *f = TFile::Open("../analyze/output_files/ttbar_background_jun05.root");
    TTree *t = (TTree *)f->Get("T_data");

    TH1F *h_m = new TH1F("h_m", "ttbar Background", 30, 150, 1000);

    TH1F *h_cost = new TH1F("back_cost", "ttbar Background", 40, -1.,1.);


    make_m_cost_hist(t, h_m, h_cost, false);

    TCanvas *c1 = new TCanvas("c1", "ZZ back M", 100,200, 900, 700);
    c1->cd();
    h_m->Print();
    h_m->Draw("hist");
    h_m->SetFillColor(kBlue);
    h_m->SetMarkerStyle(21);
    h_m->SetMarkerColor(kBlue);
    c1->Update();

    TCanvas *c2 = new TCanvas("c2", "ZZ back cost", 900, 700);
    c2->cd();
    h_cost->Print();
    h_cost->Draw("hist");
    h_cost->SetFillColor(kBlue);
    h_cost->SetMarkerStyle(21);
    h_cost->SetMarkerColor(kBlue);
    c2->Update();

}
