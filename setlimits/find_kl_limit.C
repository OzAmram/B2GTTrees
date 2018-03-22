
#include "madgraph_lhe_reader.C"

//Double_t AFB_SM[6] =  {0.620,  0.615, 0.603, 0.590, 0.586, 0.588};
Double_t AFB_SM[6] = {0.60, 0.612, 0.603, 0.602, 0.60, 0.60};
Double_t AFB_unc[6] = {0.014, 0.020,0.021,0.028,0.045, 0.063};
//Double_t AFB_measured[6] = {0.615, 0.600, 0.612, 0.608, 0.552, 0.558};
Double_t *AFB_measured = AFB_SM;

Double_t test_stat(Double_t *x, Double_t *params){
    //6 x values are asymmetries in the 6 different mass bins
    //6 params are asymmetries of Zprime model in different mass bins
    //6 hard coded SM AFB values
    //6 hard coded AFB uncertainties (become widths of Gaussians)



    //
    Double_t log_L = 0; 
    for(int i=0; i<6; i++){
        //Null hypothesis is Zprime, Alternate is SM
        Double_t term1 = pow((x[i] - AFB_SM[i]),2);
        Double_t term2 = pow((x[i] - params[i]),2);
        log_L += (-term1 + term2)/(2*pow(AFB_unc[i],2));
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
    return h->Integral(xbin, maxbin);
}


Double_t test_Zp(int M_Zp, Double_t cpl){
    //hard code measured AFB's + uncertainties and SM predictions for AFB
    //return pvalue for a given Zp mass


    Double_t AFB_Zp[6];
    //printf("AFBs are: ");
    for(int i = 0; i<6; i++){
        AFB_Zp[i] =  get_AFB(M_Zp, cpl, i);
        //printf("%.2f ", AFB_Zp[i]);
    }
    //printf("\n");

    int n_trials = 50000;

    TH1D * h_dist = new TH1D("h_dist1", "Distribution of test statistic", 400,-200,200);
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
    h_dist->Scale(1./n_trials);
    //printf("h_dist mean and std dev are %.2f %.2f \n", h_dist->GetMean(), h_dist->GetStdDev()); 

    Double_t t_obs = test_stat(AFB_measured, AFB_Zp);
    //printf("T_obs is %.2f \n", t_obs);

    Double_t pval = get_pval(h_dist, t_obs);
    //if pval < 0.05 we will reject Null hypothesis (of Zprime existiing)

    //del h_dist; 
    
    return pval;
}


void find_kl_limit(){
    int m_start = 1000;
    int m_max = 2000;
    int m_step = 50;
    Double_t kl_start = 0.3;
    Double_t kl_min = 0.05;
    Double_t kl_step = 0.05;
    Double_t alpha = 0.05;
    Double_t pval = 0.;
    Double_t m;
    Double_t kl;
    vector<double> limits;
    for(m = m_start; m <=m_max; m+=m_step){
        for(kl = kl_start; kl >= kl_min; kl-=kl_step){
            //printf("%.2f kl \n", kl);
            pval = test_Zp(m, kl);
            if (pval > alpha && kl < kl_start) break;
            else if(pval  > alpha){
                printf("Started too low, changing kl start from %.2f to %.2f\n", kl_start, kl_start +2*kl_step);
                kl_start = kl_start + 2*kl_step;
                kl = kl_start;
            }
        }
        printf("Limit found!\n"
               "For M = %.0f, kl=%.2f, pval is %0.3f \n", m, kl, pval);
        limits.push_back(kl);
    }
    printf("\n\nLimits summary: \n");
    for(int i=0; i<limits.size(); i++){
        //limit is 1 step less than where we cross the boundary (last point we
        //can exlude)
        printf("M=%i, kl < %.2f \n", m_start + i*m_step, limits[i] + kl_step);
    }

    return;
}
   



