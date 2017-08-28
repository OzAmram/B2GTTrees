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


void draw_project(){
    TFile *f_mc = TFile::Open("output_files/DYToLL_MC_mar2.root");

    TH3F* f_sym = (TH3F *)f_mc->Get("f_sym"); //3d histogram of data
    TH3F* f_asym = (TH3F *)f_mc->Get("f_asym"); //3d histogram of data

    TH1D *sym_x = f_sym->ProjectionX("sym_x", 3,4, 0, -1);
    TH1D *sym_cost = f_sym->ProjectionZ("sym_cost", 0,-1, 3,4);

    TH1D *asym_x = f_asym->ProjectionX("asym_x", 3,4, 0, -1);
    TH1D *asym_cost = f_asym->ProjectionZ("asym_cost", 0,-1, 3,4);

    TH2D *sym_2d = (TH2D *)f_sym->Project3D("xz");
    TH2D *asym_2d = (TH2D *)f_asym->Project3D("xz");

    TCanvas *c1 = new TCanvas("c1", "Histograms", 200, 10, 900, 700);
    sym_x->Draw();

    TCanvas *c2 = new TCanvas("c2", "Histograms", 200, 10, 900, 700);
    sym_cost->Draw();
    
    TCanvas *c3 = new TCanvas("c3", "Histograms", 200, 10, 900, 700);
    asym_x->Draw();

    TCanvas *c4 = new TCanvas("c4", "Histograms", 200, 10, 900, 700);
    asym_cost->Draw();

    TCanvas *c5 = new TCanvas("c5", "Histograms", 200, 10, 900, 700);
    sym_2d->Draw("colz");
 
    TCanvas *c6 = new TCanvas("c6", "Histograms", 200, 10, 900, 700);
    asym_2d->Draw("colz");
}

    
    
