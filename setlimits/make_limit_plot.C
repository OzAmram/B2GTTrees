#include "../plots/tdrstyle.C"
#include "../plots/CMS_lumi.C"
#include "find_kl_limit.C"

Double_t get_kl_limit(FILE *f1, int M_Zp, Double_t kl_start, Double_t *AFB_test){
    Double_t kl_min = 0.05;
    Double_t kl_max = 2.2;
    Double_t kl_step = 0.05;
    Double_t alpha = 0.05;
    Double_t pval = 0.;
    Double_t kl;


    for(kl = kl_start; kl >= kl_min && kl <= kl_max; kl-=kl_step){
        //printf("%.2f kl \n", kl);
        pval = test_Zp(f1, M_Zp, kl, AFB_test);
        if (pval > alpha && kl < kl_start) return kl;
        else if(pval  > alpha){
            //printf("Started too low, changing kl start from %.2f to %.2f\n", kl_start, kl_start +2*kl_step);
            kl_start = kl_start + 2*kl_step;
            kl = kl_start;
            //printf("%i %.2f \n", m, kl);
        }
    }

    //cant find limit, return minimum
    printf("Couldn't find limit for M=%i, returning kl_min \n", M_Zp);
    return kl_min;
}
Double_t get_var(int n_vals, Double_t *vals){
    Double_t mean, var;
    int n_entries = n_vals;
    for(int i =0; i< n_vals; i++){
        //printf("%.2f \n", vals[i]);
        if(vals[i] < 0. || vals[i] > 1. || std::isnan((float)vals[i])) n_entries--;
        else{
            //printf("val %.2f \n", vals[i]);
            mean += vals[i];
        }
        //printf("%.3e \n", vals[i]);
    }
    mean = mean / n_entries;
    //printf("mean %.3f n_entries %i\n", mean, n_entries);

    for(int  i=0; i< n_vals; i++){
        if(vals[i] < 0. || vals[i] > 1. || std::isnan((float)vals[i])) continue;
        else var += pow(vals[i] - mean, 2);
    }
    var = var/(n_entries -1);
    //printf("std %.3f \n\n\n", sqrt(var));
    return sqrt(var);
}
TGraph *makeAFillGraph(int len, Double_t *x, Double_t *y1, Double_t *y2, int linecolor=1, int fillcolor=0, int fillstyle=0){
    vector<Double_t> g_x;
    vector<Double_t> g_y;
    for(int i =0; i<len; i++){
        printf("filling in graph for m=%.2f \n", x[i]);
        g_x.push_back(x[i]);
        g_y.push_back(y1[i]);
    }
    for(int i =len -1; i>=0; i--){
        g_x.push_back(x[i]);
        g_y.push_back(y2[i]);
    }
    TGraph * gr = new TGraph(2*len, g_x.data(), g_y.data());
    gr->SetLineColor(linecolor);
    gr->SetFillColor(fillcolor);
    gr->SetFillStyle(fillstyle);
    return gr;
}

void make_limit_plot(){
    int n_m_bins = 100;
    int n_trials = 200;
    int m_start = 1000;
    //int m_max = 1500;
    int m_max = 3500;
    int m_step = 40;
    int m;
    int i=0;
    Double_t kl_start = 2.2;
    Double_t kl_min = 0.05;
    Double_t kl_step = 0.05;
    FILE *f1 = fopen("AFBs.txt", "r");
    TRandom3 *r3 = new TRandom3();

    Double_t masses[n_m_bins], kl_limit[n_m_bins], kl_expected_lim[n_m_bins], kl_expected_lim_up[n_m_bins], 
             kl_expected_lim_upup[n_m_bins], kl_expected_lim_down[n_m_bins],    
             kl_expected_lim_downdown[n_m_bins], kl_lim_trials[n_m_bins][n_trials],
             kl_limit_stds[n_m_bins];

    for(m = m_start; m <=m_max; m+=m_step){
        printf("trying m=%i \n", m);
        if(m > 2490 || m<1990) m_step = 50;
        if(m > 2990) m_step = 100;
        else m_step = 40;
        Double_t meas_lim = get_kl_limit(f1, m, kl_start, AFB_measured);
        Double_t exp_lim  = get_kl_limit(f1, m, kl_start, AFB_SM);
        if(meas_lim <= 0.09 || exp_lim <= 0.09 || meas_lim >=2.2 || exp_lim >= 2.2){
            printf("bad lim for m=%i exp = %.2f meas = %.2f \n", m, exp_lim, meas_lim); 
            continue;
        }
        if(i>0) meas_lim = max(meas_lim, kl_limit[i-1]);
        if(i>0) exp_lim = max(exp_lim, kl_expected_lim[i-1]);
        masses[i] = m;
        kl_limit[i] = meas_lim;
        kl_expected_lim[i] = exp_lim;
        //generate random expected limits
        Double_t AFB_rand[6];
        for(int j=0; j<n_trials; j++){
            for(int k=0; k<6; k++){
                AFB_rand[k] = r3->Gaus(AFB_SM[k], AFB_unc[k]);
            }
            kl_lim_trials[i][j] = get_kl_limit(f1, m, kl_expected_lim[i]+2*kl_step, AFB_rand);
        }
        kl_limit_stds[i] = get_var(n_trials, kl_lim_trials[i]);
        kl_expected_lim_up[i] = kl_expected_lim[i] + kl_limit_stds[i];
        kl_expected_lim_upup[i] = kl_expected_lim[i] + 2*kl_limit_stds[i];
        kl_expected_lim_down[i] = max(kl_expected_lim[i] - kl_limit_stds[i], 0.);
        kl_expected_lim_downdown[i] = max(kl_expected_lim[i] - 2*kl_limit_stds[i], 0.);
        i++;
    }
    setTDRStyle();
    TGraph *meas = new TGraph(i, masses, kl_limit);    
    TGraph *exp = new TGraph(i, masses, kl_expected_lim);    
    TGraph *one_sig = makeAFillGraph(i, masses, kl_expected_lim_down, kl_expected_lim_up, 0, kGreen+1, 1001);
    TGraph *two_sig = makeAFillGraph(i, masses, kl_expected_lim_downdown, kl_expected_lim_upup, 0, kOrange, 1001);

    meas->SetLineColor(1);
    meas->SetLineStyle(1);
    meas->SetLineWidth(4);

    exp->SetLineColor(1);
    exp->SetLineStyle(2);
    exp->SetLineWidth(4);
    auto legbb = new TLegend(0.45,0.55,0.90,0.85);
	legbb->SetFillStyle(1001);
	legbb->SetFillColor(0);    
	legbb->SetBorderSize(0);  
	legbb->SetTextSize(0.040);
	legbb->SetTextFont(42);
	legbb->AddEntry(meas,"Observed","l");
	legbb->AddEntry(exp,"Expected","l");
	legbb->AddEntry(one_sig,"#pm 1 std. deviation","f");
	legbb->AddEntry(two_sig,"#pm 2 std. deviation","f");

    TCanvas *c1 = new TCanvas("c1", "", 1200, 800);
    two_sig->SetTitle("Z' Limits");
    two_sig->GetYaxis()->SetTitle("#kappa_{L}(g_{Z'}/g_{Z})");
    two_sig->GetYaxis()->CenterTitle();
    two_sig->GetYaxis()->SetTitleOffset(1.0);
    two_sig->GetXaxis()->SetTitle("Z' mass (GeV)");

    two_sig->Draw("Af");
    one_sig->Draw("fsames");
    exp->Draw("Csames");
    meas->Draw("Csames");
    legbb->Draw();
    
    writeExtraText = true;
    extraText = "Preliminary";
    //lumi_sqrtS = "";       // used with iPeriod = 0, e.g. for simulation-only plots (default is an empty string)
    int iPeriod = 4; 
    CMS_lumi( c1, iPeriod, 0 );

}
