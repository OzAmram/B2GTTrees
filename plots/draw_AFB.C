

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

    Double_t ysm[n] = {0.618, 0.618,0.617,0.611,0.617,0.617};
    Double_t ysm_errs[n] = {0.002, 0.002,0.003,0.003,0.003,0.002};

    Double_t yratio[n], yratio_errs[n];
    for(i=0; i<n; i++){
        yratio[i] = y[i]/ysm[i];
        yratio_errs[i] = sqrt(pow(y_errs[i]/ysm[i],2) + 
                              pow((y[i]*ysm_errs[i]/ysm[i]/ysm[i]), 2));
    }


    TGraphErrors *gdata = new TGraphErrors(n, x,y,x_errs,y_errs);
    TGraph *gsm = new TGraph(n, x,ysm);
    TGraphErrors *gratio = new TGraphErrors(n, x, yratio, x_errs, yratio_errs);
    //g1->SetTitle("Drell-Yan Forward Backward Asymmetry");

    gsm->SetLineWidth(4);
    gsm->SetMarkerStyle(kFullSquare);
    gsm->SetMarkerColor(kBlue);
    gsm->SetLineColor(kBlue);

    gdata->SetMarkerStyle(kFullSquare);
    gdata->SetMarkerColor(kBlack);
    gdata->SetLineWidth(2);
    gdata->SetLineColor(kBlack);

    //TMultiGraph *mg = new TMultiGraph();
    //mg->Add(gdata);
    //mg->Add(gsm);

    gsm->SetMinimum(0);
    gsm->SetMaximum(0.75);
    gdata->SetMinimum(0);
    gdata->SetMaximum(0.75);

    //gsm->SetMarkerStyle(20);
    //g1->GetYaxis()->SetTitle("A_{FB}");
    TCanvas *c1 = new TCanvas("c1", "AFB as function of mass", 200,10, 800,600);
    TPad *pad1 = new TPad("pad1", "pad1", 0.,0.3,0.98,1.);
    pad1->SetTopMargin(0.08);
    pad1->SetBottomMargin(0.02);
    pad1->Draw();
    pad1->cd();
    gsm->Draw("A P C");
    gdata->Draw("P same");
    gsm->GetYaxis()->SetTitle("Forward Backward Asymmetry");
    gsm->GetYaxis()->SetTitleSize(20);
    gsm->GetYaxis()->SetTitleFont(43);
    gsm->GetYaxis()->SetTitleOffset(1.1);
    //gdata->GetXaxis()->SetTitle("DiMuon Mass (GeV)");
    gsm->GetXaxis()->SetLabelSize(0);
    gsm->GetXaxis()->SetRangeUser(100,1000);

    gStyle->SetLegendBorderSize(0);
    TLegend *leg2 = new TLegend(0.5, 0.65, 0.75, 0.8);
    leg2->AddEntry(gsm, "Standard Model A_{FB} from POWHEG", "l");
    leg2->AddEntry(gdata, "Data", "p");
    leg2->Draw();
    pad1->Update();

    c1->cd();
    TPad *pad2 = new TPad("pad2", "pad2", 0.,0,.98,0.3);
    pad2->SetTopMargin(0.02);
    pad2->SetBottomMargin(0.4);
    pad2->Draw();
    pad2->cd();
    gratio->Draw("A P");
    TLine *l1 = new TLine(100,1,1000,1);
    l1->SetLineStyle(2);
    l1->Draw();

    gratio->GetYaxis()->SetTitle("Data/SM");
    gratio->GetYaxis()->SetRangeUser(0.85,1.15);
    gratio->GetYaxis()->SetNdivisions(505);
    gratio->GetYaxis()->SetTitleSize(20);
    gratio->GetYaxis()->SetTitleFont(43);
    gratio->GetYaxis()->SetTitleOffset(1.2);
    gratio->GetYaxis()->SetLabelFont(43); // Absolute font size in pixel (precision 3)
    gratio->GetYaxis()->SetLabelSize(15);
    // X axis ratio plot settings
    gratio->GetXaxis()->SetRangeUser(100,1000);
    gratio->GetXaxis()->SetTitle("M_{#mu#mu} (GeV)");
    gratio->GetXaxis()->SetTitleSize(30);
    gratio->GetXaxis()->SetTitleFont(43);
    gratio->GetXaxis()->SetTitleOffset(3.);
    gratio->GetXaxis()->SetLabelFont(43); // Absolute font size in pixel (precision 3)
    gratio->GetXaxis()->SetLabelSize(20);

    writeExtraText = true;
    extraText = "Preliminary";
    //lumi_sqrtS = "";       // used with iPeriod = 0, e.g. for simulation-only plots (default is an empty string)
    int iPeriod = 4; 
    CMS_lumi( pad1, iPeriod, 0 );
}

