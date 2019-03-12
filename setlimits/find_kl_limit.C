
//#include "madgraph_lhe_reader.C"
#include "read_AFBs_from_file.C"

//Double_t AFB_SM[6] =  {0.620,  0.615, 0.603, 0.590, 0.586, 0.588};
Double_t AFB_SM[6] = {0.60, 0.612, 0.603, 0.602, 0.60, 0.60};
Double_t AFB_unc[6] = {0.014, 0.020,0.020,0.027,0.044, 0.059};
Double_t AFB_measured[6] = {0.609, 0.578, 0.600, 0.574, 0.579, 0.532};

Double_t test_stat(Double_t *x, Double_t *params){
    //6 x values are asymmetries in the 6 different mass bins
    //6 params are asymmetries of Zprime model in different mass bins
    //6 hard coded SM AFB values
    //6 hard coded AFB uncertainties (become widths of Gaussians)



    //
    Double_t log_L = 0; 
    for(int i=0; i<6; i++){
        //Null hypothesis is Zprime, Alternate is SM
        Double_t chi_sm = pow((x[i] - AFB_SM[i])/AFB_unc[i],2);
        Double_t chi_zp = pow((x[i] - params[i])/AFB_unc[i],2);
        log_L += (chi_zp - chi_sm);
    }
    Double_t ret_val = log_L;
    //printf("%.2e \n", ret_val);
    return ret_val;
}


Double_t get_pval(TH1D *h, Double_t x){
    //integrate the histogram from x to infinity (To find p-value)
    //
    TAxis *ax = h->GetXaxis();
    Int_t xbin = ax->FindBin(x);
    Int_t maxbin = ax->FindBin(ax->GetXmax());
    return h->Integral(xbin, maxbin)/ h->Integral();
}


Double_t test_Zp(FILE *f1, int M_Zp, Double_t cpl, Double_t *AFB_test){
    //hard code measured AFB's + uncertainties and SM predictions for AFB
    //return pvalue for a given Zp mass


    Double_t AFB_Zp[6];

    read_AFBs(f1, AFB_Zp, M_Zp, cpl);
    /*
    printf("AFBs are: ");
    for(int i = 0; i<6; i++){
        //AFB_Zp[i] =  get_AFB(M_Zp, cpl, i);
        printf("%.2f ", AFB_Zp[i]);
    }
    printf("\n");
    */

    int n_trials = 50000;

    TH1D * h_dist = new TH1D("h_dist1", "Distribution of test statistic", 400,-50,50);
    TRandom3 *r3 = new TRandom3();
    Double_t x[6]; //randomly generated x vals
    //Find distribution of test statistic under Null hypothesis (ZPrime is
    //there)
    for(int i=0; i<n_trials; i++){
        for(int j=0; j<6; j++){
            x[j] = r3->Gaus(AFB_Zp[j], AFB_unc[j]);
        }
        Double_t test_val = test_stat(x, AFB_Zp);
        h_dist->Fill(test_val);
    }
    //float bin_width = 100./400.;
    //h_dist->Scale(1./n_trials/bin_width);
    //printf("h_dist mean and std dev are %.2f %.2f \n", h_dist->GetMean(), h_dist->GetStdDev()); 
    //TF1 *f = new TF1("chi2", "ROOT::Math::chisquared_pdf(x,6)", 0., 50);
    //h_dist->Draw("hist");
    //f->Draw("same");

    Double_t t_obs = test_stat(AFB_test, AFB_Zp);
    //prin2.2"T_obs is %.2f \n", t_obs);

    Double_t pval = get_pval(h_dist, t_obs);
    //if pval < 0.05 we will reject Null hypothesis (of Zprime existiing)

    delete h_dist; 
    
    return pval;
}


void find_kl_limit(){
    Double_t *AFB_test = AFB_measured;
    int m_start = 1000;
    int m_max = 3300;
    int m_step = 20;
    int m;
    Double_t kl_start = 0.4;
    Double_t kl_min = 0.05;
    Double_t kl_max = 2.2;
    Double_t kl_step = 0.05;
    Double_t alpha = 0.05;
    Double_t pval = 0.;
    Double_t kl;
    vector<Double_t> limits;
    vector<int> ms;
    FILE *f1;

    f1 = fopen("AFBs.txt", "r");

    m=2300.;
    kl = 1.0;

    /*
    pval = test_Zp(f1, m, kl, AFB_test);
    Double_t AFB_Zp[6];
    read_AFBs(f1, AFB_Zp, m, kl);
    Double_t t_obs = test_stat(AFB_test, AFB_Zp);
    Double_t pval2 = ROOT::Math::chisquared_cdf(t_obs, 6);
    printf("test stat %.2f numeric %.3f analytic %.3f \n", t_obs, pval, pval2);
    */

    for(m = m_start; m <=m_max; m+=m_step){
        for(kl = kl_start; kl >= kl_min && kl <= kl_max; kl-=kl_step){
            if(m > 2490. || m<1990.) m_step = 50;
            else m_step = 20;
            //printf("%.2f kl \n", kl);
            pval = test_Zp(f1, m, kl, AFB_test);
            if (pval > alpha && kl < kl_start) break;
            else if(pval  > alpha){
                printf("Started too low, changing kl start from %.2f to %.2f\n", kl_start, kl_start +2*kl_step);
                kl_start = kl_start + 2*kl_step;
                kl = kl_start;
                printf("%i %.2f \n", m, kl);
            }
        }
        printf("Limit found!\n"
               "For M = %i, kl=%.2f, pval is %0.3f \n", m, kl, pval);
        limits.push_back(kl);
        ms.push_back(m);
    }
    printf("\n\nLimits summary: \n");
    for(int i=0; i<limits.size(); i++){
        //limit is 1 step less than where we cross the boundary (last point we
        //can exlude)
        
        printf("M=%i, kl < %.2f \n", ms[i], limits[i] + kl_step);

    }

    return;
}
   




