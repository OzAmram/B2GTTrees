#include "find_limit.C"

void test2(){
    Double_t AFB_measured[6] = {0.615, 0.600, 0.612, 0.608, 0.552, 0.558};
    Double_t AFB_SM[6] =  {0.620,  0.615, 0.603, 0.590, 0.586, 0.588};
    Double_t AFB_unc[6] = {0.014, 0.020,0.021,0.028,0.045, 0.063};


    Double_t AFB_Zp[6];
    Double_t Zp_mass = 2000.;
    for(int i = 0; i<6; i++){
        AFB_Zp[i] =  get_AFB(Zp_mass, i);
    }

    /*
    int n_trials = 50000;
    TH1D *h_gaus = new TH1D("h_gaus", "gauss", 100.,-5.,5.);
    TRandom3 *r1 = new TRandom3();
    for(int i=0; i<n_trials; i++){
        h_gaus->Fill(r1->Gaus(0.,1.));
    }
    h_gaus->Scale(1./n_trials);

    Double_t int1 = get_pval(h_gaus, 0.);
    Double_t int2 = get_pval(h_gaus, 1.);
    Double_t int3 = get_pval(h_gaus, 2.);

    printf("Integrals are %.4f %.4f %.4f \n", int1, int2,int3);

    TH1D * h_dist = new TH1D("h_dist", "Distribution of test statistic", 200,-100,100);
    TRandom3 *r3 = new TRandom3();
    Double_t x[6]; //randomly generated x vals
    for(int i=0; i<n_trials; i++){
        for(int j=0; j<6; j++){
            x[j] = r3->Gaus(AFB_Zp[j], AFB_unc[j]);
        }
        Double_t test_val = test_stat(x, AFB_Zp);
        h_dist->Fill(test_val);
    }
    
    TCanvas *can2 = new TCanvas("c2", "Histograms", 200, 10, 900, 700);
    h_dist->SetFillColor(kBlue);
    h_dist->Scale(1./n_trials);
    h_dist->Draw("hist");

    Double_t t_obs = test_stat(AFB_measured, AFB_Zp);

    Double_t pval = get_pval(h_dist, t_obs);
    */

    Double_t pval1 = test_Zp(1000.);
    Double_t pval2 = test_Zp(2000.);
    Double_t pval3 = test_Zp(3000.);
    Double_t pval4 = test_Zp(10000.);

    printf("P-values are %.3e %.3e %.3e %.3e \n", pval1, pval2, pval3, pval4);
    //printf("P-values are %.3e %.3e %.3e %.3e \n", pval, pval2, pval3, pval4);

    
    return;

}
