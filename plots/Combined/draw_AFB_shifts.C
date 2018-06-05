

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
#include "../tdrstyle.C"
#include "../CMS_lumi.C"

void draw_AFB_shifts(){
    setTDRStyle();
    Double_t x[11] = {0., 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0};

    Double_t shifts[6][11]= 
    {-0.141,  -0.124,  -0.107,  -0.089,  -0.071,  -0.053,  -0.035,  -0.017,  0.002 , 0.021 , 0.040,
	-0.143,  -0.126,  -0.108,  -0.090,  -0.073,  -0.054,  -0.036,  -0.017,  0.001 , 0.021 , 0.040 , 
	-0.166,  -0.147,  -0.127,  -0.107,  -0.086,  -0.066,  -0.045,  -0.023,  -0.002,  0.021,  0.044,
	-0.115,  -0.102,  -0.088,  -0.075,  -0.060,  -0.046,  -0.032,  -0.017,  -0.002,  0.014,  0.030,
	-0.124,  -0.109,  -0.094,  -0.079,  -0.064,  -0.049,  -0.033,  -0.017,  -0.001,  0.015,  0.032,
	-0.104,  -0.092,  -0.080,  -0.067,  -0.054,  -0.042,  -0.028,  -0.015,  -0.001,  0.013,  0.028}; 

    Double_t AFB_meas[6] = {0.611, 0.612, 0.614, 0.608, 0.559, 0.532};
    //correct factor for AFB
	Double_t corr[6][11];
	for(int i=0; i<6; i++){
		for(int j=0;j<11; j++){
            corr[i][j] = (AFB_meas[i] + shifts[i][j])/AFB_meas[i];
        }
    }

	TGraph *g[6];
	g[0] = new TGraph(11, x, corr[0]);
	g[1] = new TGraph(11, x, corr[1]);
	g[2] = new TGraph(11, x, corr[2]);
	g[3] = new TGraph(11, x, corr[3]);
	g[4] = new TGraph(11, x, corr[4]);
	g[5] = new TGraph(11, x, corr[5]);
		

    g[0]->SetMarkerColor(kBlue);
    g[0]->SetLineColor(kBlue);
    g[0]->SetLineWidth(2);
    g[0]->GetYaxis()->SetRangeUser(0.5,1.3);


    g[1]->SetMarkerColor(kGreen+3);
    g[1]->SetLineColor(kGreen+3);
    g[1]->SetLineWidth(2);

    g[2]->SetMarkerColor(kMagenta +3);
    g[2]->SetLineColor(kMagenta +3);
    g[2]->SetLineWidth(2);

    g[3]->SetMarkerColor(kRed-6);
    g[3]->SetLineColor(kRed-6);
    g[3]->SetLineWidth(2);

    g[4]->SetMarkerColor(kOrange+7);
    g[4]->SetLineColor(kOrange+7);
    g[4]->SetLineWidth(2);

    g[5]->SetLineWidth(2);


    TCanvas *c_m = new TCanvas("c_m", "Histograms", 200, 10, 900, 700);
    g[0]->Draw("ALP");
    g[1]->Draw("LP same");
    g[2]->Draw("LP same");
    g[3]->Draw("LP same");
    g[4]->Draw("LP same");
    g[5]->Draw("LP same");

    gStyle->SetLegendBorderSize(0);
    TLegend *leg1 = new TLegend(0.5, 0.65, 0.75, 0.8);
    leg1->AddEntry(g[0], "M150", "p");
    leg1->AddEntry(g[1], "M200", "p");
    leg1->AddEntry(g[2], "M250", "p");
    leg1->AddEntry(g[3], "M350", "p");
    leg1->AddEntry(g[4], "M500", "p");
    leg1->AddEntry(g[5], "M700", "p");
    leg1->Draw();



    g[0]->GetXaxis()->SetTitle("u-type quark fraction");
    g[0]->GetYaxis()->SetTitle("AFB Correction Factor");
    int iPeriod = 4; 
    CMS_lumi(c_m, iPeriod, 11 );
    c_m->Update();
    
    printf("Performing fits: \n");
    TF1 *f[6];
    for(int i=0; i<6; i++){
        f[i] = new TF1("f", "[0]*x + [1]", 0,1);
        g[i]->Fit(f[i], "Q");
        Double_t m = f[i]->GetParameter(0);
        Double_t b = f[i]->GetParameter(1);
        printf("Bin %i: m=%.3f b=%.3f \n", i, m, b);
    }
}


