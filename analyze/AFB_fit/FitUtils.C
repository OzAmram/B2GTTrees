#include "root_files.h"

double get_chi_sq(TH2F *h_data, TH2F *h_sym, TH2F* h_asym, TH2F* h_back, double AFB, double r_back){
    h_data->Scale(1./h_data->Integral());
    float sum=0;

    for(int i=1; i<n_xf_bins; i++){
        for(int j=1; j<n_cost_bins; j++){
            float sym = h_sym->GetBinContent(i,j);
            float asym = h_asym->GetBinContent(i,j);
            float back = h_back->GetBinContent(i,j);
    
            float fit = r_back*back + (1 - r_back) * (sym + AFB*asym);
            float data = h_data->GetBinContent(i,j);
            float unc = h_data->GetBinError(i,j);

            //printf("fit %.2e data %.2e unc %.2e \n", fit,data,unc);
            sum += pow((fit-data)/unc, 2);
           
        }
    }
    return sum;
}

double get_chi_sq_v2(TH2F *h_data, TH2F *h_sym, TH2F* h_asym, TH2F* h_back, TH2F *h_qcd, double AFB, double r_back, double r_qcd){
    h_data->Scale(1./h_data->Integral());
    float sum=0;

    for(int i=1; i<n_xf_bins; i++){
        for(int j=1; j<n_cost_bins; j++){
            float sym = h_sym->GetBinContent(i,j);
            float asym = h_asym->GetBinContent(i,j);
            float back = h_back->GetBinContent(i,j);
            float qcd = h_qcd->GetBinContent(i,j);
    
            float fit = r_back*back + r_qcd*qcd + (1 - r_back - r_qcd) * (sym + AFB*asym);
            float data = h_data->GetBinContent(i,j);
            float unc = h_data->GetBinError(i,j);

            //printf("fit %.2e data %.2e unc %.2e \n", fit,data,unc);
            sum += pow((fit-data)/unc, 2);
           
        }
    }
    return sum;
}
