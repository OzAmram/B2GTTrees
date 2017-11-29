
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
#include "../../analyze/TemplateMaker.C"

/*
void SetErrors(TH1D *h_pass, TH1D *h_total){
    int nBins = h_pass->GetXaxis()->GetNbins();
    for (int i=0; i < nBins; i++){
        //binomial distribution divided by n
        Double_t n = h_total->GetBinContent(i);
        Double_t p = h_pass->GetBinContent(i);
        Double_t err = sqrt(p*(1-p)/n);
        h_pass->SetBinError(i, err);
    }
    return;
}
*/
        
    


void draw_mc_fakerate(){
    TFile *f = TFile::Open("../analyze/FakeRate/root_files/SingleEl_mc_fakerate_contam_v2_nov28.root");

    TH2D* h_pass = (TH2D *)f->Get("h_pass"); 
    TH2D* h_total = (TH2D *)f->Get("h_total"); 

    TH2D * h_rate = (TH2D *) h_pass->Clone("h_rate");
    h_rate->Divide(h_total);
    h_pass->Scale(1000*tot_lumi);
    h_total->Scale(1000*tot_lumi);

    TH1D *pass_barrel = h_pass->ProjectionY("pass_barrel", 1,1);
    TH1D *pass_endcap = h_pass->ProjectionY("pass_endcap", 2,2);

    TH1D *rate_barrel = h_rate->ProjectionY("rate_barrel", 1,1);
    TH1D *rate_endcap = h_rate->ProjectionY("rate_endcap", 2,2);

    TH1D *total_barrel = h_total->ProjectionY("total_bar", 1,1);
    TH1D *total_endcap = h_total->ProjectionY("total_endcap", 2,2);
    //SetErrors(pass_barrel, total_barrel);
    //SetErrors(pass_endcap, total_endcap);

    

    TCanvas *c1 = new TCanvas("c1", "Histograms", 200, 10, 900, 700);
    pass_barrel->SetTitle("Muons fakepass");
    pass_barrel->SetStats(0);
    pass_barrel->SetLineWidth(3);
    pass_barrel ->Draw("E1");
    pass_barrel->GetXaxis()->SetTitle("p_t (Gev)");

    //TCanvas *c2 = new TCanvas("c2", "Histograms", 200, 10, 900, 700);
    //pass_endcap->SetTitle("Fake-pass Muons with Trigger: End Cap");
    pass_endcap->SetStats(0);
    pass_endcap->SetLineWidth(3);
    pass_endcap->Draw("E1 same");
    pass_endcap->SetLineColor(kRed);
    TLegend *leg1 = new TLegend(0.5, 0.65, 0.75, 0.8);
    leg1->AddEntry(pass_barrel, "Barrel", "f");
    leg1->AddEntry(pass_endcap, "Endcap",  "f");
    leg1->Draw();
    
    TCanvas *c2 = new TCanvas("c2", "Histograms", 200, 10, 900, 700);
    rate_barrel->SetTitle("Muons fakerate");
    rate_barrel->SetStats(0);
    rate_barrel->SetLineWidth(3);
    rate_barrel ->Draw("E1");
    rate_barrel->GetXaxis()->SetTitle("p_t (Gev)");

    //TCanvas *c2 = new TCanvas("c2", "Histograms", 200, 10, 900, 700);
    //rate_endcap->SetTitle("Fake-rate Muons with Trigger: End Cap");
    rate_endcap->SetStats(0);
    rate_endcap->SetLineWidth(3);
    rate_endcap->Draw("E1 same");
    rate_endcap->SetLineColor(kRed);
    TLegend *leg3 = new TLegend(0.5, 0.65, 0.75, 0.8);
    leg3->AddEntry(rate_barrel, "Barrel", "f");
    leg3->AddEntry(rate_endcap, "Endcap",  "f");
    leg3->Draw();



    TCanvas *c3 = new TCanvas("c3", "Histograms", 200, 10, 900, 700);
    total_barrel->SetTitle("Number of barrel events");
    total_barrel->SetStats(0);
    total_barrel->SetLineWidth(3);
    total_barrel ->Draw("hist");
    total_barrel->GetXaxis()->SetTitle("p_t (Gev)");
    c3->SetLogy();

    //TCanvas *c4 = new TCanvas("c4", "Histograms", 200, 10, 900, 700);
    total_endcap->SetTitle("Number of events");
    total_endcap->SetStats(0);
    total_endcap->SetLineWidth(3);
    total_endcap->Draw("hist");
    total_endcap->GetXaxis()->SetTitle("p_t (Gev)");
    TLegend *leg2 = new TLegend(0.5, 0.65, 0.75, 0.8);
    leg2->AddEntry(total_barrel, "Barrel",  "f");
    leg2->AddEntry(total_endcap, "Endcap" , "f");
    leg2->Draw();
    //c4->SetLogy();
}

    
    
