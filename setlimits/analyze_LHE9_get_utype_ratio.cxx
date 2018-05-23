
//! \file analyze_LHE
//!
//! Code to read and analyze LHE files generated by Powheg/MadGraph
// get dilution factor for u and d quarks


#include <stdlib.h> 
#include <stdio.h>
#include <string.h>
#include <stddef.h>
#include <math.h>
#include <algorithm>
#include <vector>
#include <iostream>
#include <iomanip>
#include <sstream>
#include <fstream>
#include <sys/time.h>

using namespace std;

#include "TROOT.h"
#include "TFile.h"
#include "TString.h"
#include "TObject.h"
#include "TH1F.h"
#include "TF1.h"
#include "TProfile.h"
#include "TStyle.h"
#include "TCanvas.h"
#include "TPostScript.h"
#include "TLorentzVector.h"
#include "TVector3.h"
#include "TFitter.h"
#include "TSystem.h"

// Global definitions 

struct particle {         //!< individual particle structure
    int id;                 //! pdg particle id
    int status;             //! decay status
    int mother[2];          //! particle mother indices
    int color[2];           //!< color info
    double p[5];            //! 4-vector + mass
    float vtime;            //! lifetime
    float spin;             //! spin (helicity) info
} ;




// fcn passes back f = - 2*ln(L), the function to be minimized.


// Main program  

int main(int argc, char *argv[])
{
    // Local variables 
    int i, j, k, iqk=0, iqb=0, itq=0, itb=0, ig1=0, ig2=0, ijt, njet=0, npflavor=0, nevent, nfit;
    int nqqb=0, ngg=0, itbq, itbb, itqk, itqb, itlp, itnu, iqk_extra=0, iqb_extra=0;
    int nup, idrup; 
    float xwgtup, scalup, aqedup, aqcdup;
    double alpha;
    static char infile[80], header[80], inp1[80], outfile0[80], outfile1[80], outfile2[80];
    FILE *ifp;
    particle in;

    struct timeval now0, now1;
    struct timezone timz;
    long deltas, deltaus;
    double deltat;

    // Construct file name to analyze

    strcpy(infile, "/uscms_data/d3/oamram/condor_jobs/condor_output/Zp_M2000/Zp_events_M2000_kL.80_bin6.lhe");
    TString fout_name("dilus/SM/M700_dilus.root");
    bool write_out = false;
    const double m_low = 700;
    const double m_high = 10000.;
    const double pt_low = 0.;
    const double pt_high = 100000.;
    char hist_name[200];
    sprintf (hist_name, "Dilepton Angular Distribution. M in [%.0f, %.0f]. Pt in [%.0f, %.0f]; cos(#theta^{*})", 
            m_low, m_high, pt_low, pt_high);


    // Make sure file is available 

    ifp = fopen(infile, "r");
    if (ifp==NULL) {
        printf("can't find %s\n",infile);
        return 1;
    }

    //  Book Histograms
    int n_xf_bins = 5;
    Float_t xf_bins[] = {0., 0.02, 0.04, 0.07, 0.10, 1.0};
    TH1F * utype_Ni= new TH1F("utype_Ni", "", n_xf_bins, xf_bins);
    TH1F * dtype_Ni= new TH1F("utype_Nc", "", n_xf_bins, xf_bins);
    TH1F * utype_Nc= new TH1F("dtype_Ni", "", n_xf_bins, xf_bins);
    TH1F * dtype_Nc= new TH1F("dtype_Nc", "", n_xf_bins, xf_bins);


    double root2 = sqrt(2.);
    double Ebeam = 6500.;
    double Pbeam = sqrt(Ebeam*Ebeam - 0.938*0.938);

    nevent=0; nfit=0;
    int event_num=0;
    int nFail =0;
    double nF = 0, nB = 0;

    // Scan throught file first for the beginning of an event record	

    while (fscanf(ifp,"%s", inp1) !=EOF) {
        if(strcmp(inp1,"<event>") == 0) {
            event_num++;

            // event found


            // read event info

            fscanf(ifp,"%d %d %f %f %f %f", &nup, &idrup, &xwgtup, &scalup, &aqedup, &aqcdup);

            // Make a vector to hold all particles for each event

            std::vector<particle> up;

            iqk = -1; iqb = -1; itq = -1; itb = -1; ig1 = -1; ig2 = -1; ijt = -1;
            iqk_extra = -1;
            iqb_extra = -1; //in case event with 2 quarks or 2 anti quarks
            njet = 0; npflavor = 0; nqqb = 0; ngg = 0;

            // Now loop over all particles

            for(i=0; i<nup; ++i) {
                fscanf(ifp,"%d %d %d %d %d %d %lf %lf %lf %lf %lf %f %f", &in.id, &in.status, &in.mother[0], &in.mother[1],
                        &in.color[0], &in.color[1], &in.p[0], &in.p[1], &in.p[2], &in.p[3], &in.p[4], &in.vtime, &in.spin);	
                up.push_back(in);
            } 

            //  Do event by event analysis here	

            for(i=0; i<nup; ++i) {
                if(up[i].mother[0] == 0 && up[i].mother[1] == 0) {
                    if(up[i].id == 1 || (up[i].id == 2 || (up[i].id == 3 || (up[i].id == 4 || up[i].id == 5)))) {
                        if(iqk == -1) iqk = i; 
                        else{ 
                            iqk_extra = i;
                            //printf("Found double quark event \n");
                        }
                        continue;
                    }
                    if(up[i].id == -1 || (up[i].id == -2 || (up[i].id == -3|| (up[i].id == -4 || up[i].id == -5)))) {
                        if (iqb == -1) iqb = i; 
                        else{ 
                            iqb_extra = i;
                            //printf("Found double quark event \n");
                        }
                        continue;
                    }
                    if(ig1 == -1) {ig1 = i;} else {ig2 = i;}
                    continue;
                }
                if(up[i].id == 13) {itq = i; continue;}
                if(up[i].id == -13) {itb = i; continue;}

            }



            if(((iqk == -1) && (ig1 == -1 && iqb_extra == -1)) ||
               ((iqb == -1) && (ig1 == -1 && iqk_extra ==-1)) ||
                ((iqk == -1 && iqb == -1) && (ig2 == -1)) )   {
                printf("unable to determine initial state for event %i \n", event_num);
                printf("indices are %i %i %i %i \n", iqk, iqb, ig1, ig2);
                printf("printing particles \n");
                for(i=0; i<nup; i++){
                    printf("%d %d %d %d \n", up[i].id, up[i].status, up[i].mother[0], up[i].mother[1]);
                }
                nFail++;
                continue;
            }
            if(itq == -1 || itb == -1){
                printf("no mus \n");
                nFail++;
                continue;
            }


            TLorentzVector tq(up[itq].p[0], up[itq].p[1], up[itq].p[2], up[itq].p[3]);
            TLorentzVector tb(up[itb].p[0], up[itb].p[1], up[itb].p[2], up[itb].p[3]);
            TLorentzVector p1(0., 0., Pbeam, Ebeam);
            TLorentzVector p2(0., 0., -Pbeam, Ebeam);
            if(iqk > -1) {
                if(up[iqk].p[2] < 0.) {
                    TLorentzVector p = p1;
                    p1 = p2;
                    p2 = p;
                }
            } else {
                if(up[iqb].p[2] > 0.) {
                    TLorentzVector p = p1;
                    p1 = p2;
                    p2 = p;
                }
            }

            // Create vectors for the parents
            TLorentzVector cm = tq + tb;
            Double_t pt = cm.Pt();

            if(cm.M() >= m_low && cm.M() <= m_high && pt >= pt_low && pt <= pt_high){
                ++nevent;
                // Invariant mass of pair
                double mtt = cm.M();

                double xF = abs(2.*cm.Pz()/13000.); 
                // Boost everything into the t-tbar rest frame

                TVector3 beta = -cm.BoostVector();
                tq.Boost(beta);
                tb.Boost(beta);
                p1.Boost(beta);
                p2.Boost(beta);

                // Now calculate the direction of the new z azis

                TVector3 p1u = p1.Vect();
                p1u.SetMag(1.0);
                TVector3 p2u = p2.Vect();
                p2u.SetMag(1.0);
                TVector3 pzu = p1u - p2u;
                pzu.SetMag(1.0);
                tq.RotateUz(pzu);
                double costcs = tq.CosTheta();
                if(costcs > 0.) {++nF;} else {++nB;}
                double cost2 = costcs*costcs;
                alpha = 0.;

                if(p1.Pz() * cm.Pz() <0 ){ //wrong sign in reco (cm and quark directions opposite)
                    if(up[iqk].id % 2 == 0 || up[iqb].id % 2 == 0){ //even id means u-type (unfilled will be -1)
                        utype_Ni ->Fill(xF);
                    }
                    else{
                        dtype_Ni ->Fill(xF);
                    }
                }
                else{
                    if(up[iqk].id % 2 == 0 || up[iqb].id % 2 == 0){ //even id means u-type (unfilled will be -1)
                        utype_Nc ->Fill(xF);
                    }
                    else{
                        dtype_Nc ->Fill(xF);
                    }
                }


                // Save the kinematic variables for later use


            }
        }
    }

    /* close input file */
    fclose(ifp);
    printf("nf, nb = (%.0f, %.0f) \n", nF, nB);
    double AFB = ((nF - nB))/((nF+nB));
    double dAFB = (1.-AFB*AFB)/sqrt((nF+nB));

    printf( "A_FB count = %.3f +- %.3f\n" 
            "%i events failed to id \n",
            AFB, dAFB, nFail);
    printf("total events %i %i fails \n \n", nevent, nFail++);
    Double_t u_tot = utype_Ni->Integral() + utype_Nc->Integral();
    Double_t d_tot = dtype_Ni->Integral() + dtype_Nc->Integral();
    printf("Total number of events: u-type %.0f d-type %.0f \n", u_tot, d_tot);
    printf("Fraction of u: %.3f \n", u_tot/(d_tot + u_tot));

    TH1F *utype_num = (TH1F *)utype_Nc->Clone("utype_num");
    TH1F *utype_denom = (TH1F *)utype_Nc->Clone("utype_num");
    utype_num->Add(utype_Ni, -1);
    utype_denom->Add(utype_Ni, 1);
    TH1F *utype_dilu = (TH1F *)utype_num->Clone("utype_dilu");
    utype_dilu->Divide(utype_denom);

    TH1F *dtype_num = (TH1F *)dtype_Nc->Clone("dtype_num");
    TH1F *dtype_denom = (TH1F *)dtype_Nc->Clone("dtype_num");
    dtype_num->Add(dtype_Ni, -1);
    dtype_denom->Add(dtype_Ni, 1);
    TH1F *dtype_dilu = (TH1F *)dtype_num->Clone("dtype_dilu");
    dtype_dilu->Divide(dtype_denom);

    TH1F *tot_num = (TH1F *) utype_num->Clone("tot_num");
    TH1F *tot_denom = (TH1F *) utype_denom->Clone("tot_denom");
    tot_num->Add(dtype_num);
    tot_denom->Add(dtype_denom);
    TH1F *mix_dilu = (TH1F *)tot_num->Clone("mix_dilu");
    mix_dilu->Divide(tot_denom);


    if(write_out){
        TFile *fout = TFile::Open(fout_name, "RECREATE");
        fout->cd();
        utype_dilu->Write();
        dtype_dilu->Write();
        mix_dilu->Write();
        fout->Close();
        printf("Wrote out dilus to %s \n", fout_name.Data());
    }

    /*
    printf("dtype correct vs. incorrect \n");
    for(int i=0; i<= n_xf_bins; i++){
        printf("%.0f ", dtype_Nc->GetBinContent(i));
    }
    printf("\n");
    for(int i=0; i<= n_xf_bins; i++){
        printf("%.0f ", dtype_Ni->GetBinContent(i));
    }
    printf("\n");
    */

    return 0;
} // MAIN__ 

