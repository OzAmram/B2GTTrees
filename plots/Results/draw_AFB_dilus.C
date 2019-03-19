

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

void draw_AFB_dilus(){
    setTDRStyle();
    Double_t x[11] = {0., 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0};

    Double_t meas0[11] ={0.750, 0.750, 0.750, 0.750, 0.734, 0.700, 0.668, 0.639, 0.612, 0.587, 0.564};
    Double_t meas1[11] ={0.706, 0.688, 0.670, 0.653, 0.637, 0.622, 0.608, 0.594, 0.581, 0.568, 0.556};
    Double_t meas2[11] ={0.728, 0.713, 0.697, 0.681, 0.666, 0.651, 0.636, 0.621, 0.607, 0.593, 0.579};
    Double_t meas3[11] ={0.750, 0.747, 0.721, 0.696, 0.671, 0.648, 0.625, 0.603, 0.581, 0.561, 0.542};
    Double_t meas4[11] ={0.714, 0.695, 0.677, 0.660, 0.644, 0.628, 0.613, 0.599, 0.586, 0.573, 0.560};
    Double_t meas5[11] ={0.641, 0.628, 0.614, 0.601, 0.588, 0.575, 0.563, 0.552, 0.541, 0.529, 0.519};

	TGraph *g[6];
	g[0] = new TGraph(7, &x[4], &meas0[4]);
	g[1] = new TGraph(11, x, meas1);
	g[2] = new TGraph(11, x, meas2);
	g[3] = new TGraph(10, &x[1], &meas3[1]);
	g[4] = new TGraph(11, x, meas4);
	g[5] = new TGraph(11, x, meas5);
		

    g[0]->SetMarkerColor(kBlue);
    g[0]->SetLineColor(kBlue);
    g[0]->SetLineWidth(2);


    g[1]->SetMarkerColor(kGreen+3);
    g[1]->SetLineColor(kGreen+3);
    g[1]->SetLineWidth(2);
    g[1]->GetYaxis()->SetRangeUser(0.4,0.8);

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
    g[1]->Draw("ALP");
    g[0]->Draw("LP same");
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
    g[0]->GetYaxis()->SetTitle("AFB");
    int iPeriod = 4; 
    CMS_lumi(c_m, iPeriod, 11 );
    c_m->Update();
    
    printf("Performing fits: \n");
    TF1 *f[6];
    for(int i=0; i<6; i++){
        f[i] = new TF1("f", "[0]*x + [1]", 0,1);
        g[i]->Fit(f[i], "QN");
        Double_t m = f[i]->GetParameter(0);
        Double_t b = f[i]->GetParameter(1);
        printf("Bin %i: m=%.3f b=%.3f \n", i, m, b);
    }
}


