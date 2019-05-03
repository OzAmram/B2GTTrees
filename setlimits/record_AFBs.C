
#include "madgraph_lhe_reader.C"



void record_AFB(FILE *fout, int M_Zp, Double_t cpl){


    Double_t AFB_Zp[6];
    //printf("AFBs are: ");
    for(int i = 0; i<6; i++){
        AFB_Zp[i] =  get_AFB(M_Zp, cpl, i, true, true);
        //if can't get a good AFB, don't record
        if(AFB_Zp[i] > 1.0 || AFB_Zp[i] <= 0.01) return; 

    }
        printf("Writing out AFBs for M %i kL %.2f \n", M_Zp, cpl);
        fprintf(fout, "M_Zp=%i kL=%.2f: %.3f %.3f %.3f %.3f %.3f %.3f \n", M_Zp, cpl, 
                AFB_Zp[0],AFB_Zp[1],AFB_Zp[2], AFB_Zp[3], AFB_Zp[4], AFB_Zp[5]);
}


void record_AFBs(){
    int m_start = 2000;
    int m_max = 2400;
    int m_step = 20;

    FILE *fout = fopen("test_AFBs.txt", "a");
    int m;
    Double_t kl_start = 1.0;
    Double_t kl_min = 0.05;
    Double_t kl_max = 2.3;
    Double_t kl_step = 0.05;
    Double_t alpha = 0.05;
    Double_t pval = 0.;
    Double_t kl;
    vector<Double_t> limits;
    vector<int> ms;

    //for(m = m_start; m <=m_max; m+=m_step){
        //for(kl = kl_start; kl <= kl_max; kl+=kl_step){
            //if(m%100 == 0) continue;
            //if(m > 2490) m_step = 50;
            //printf("%.2f kl \n", kl);
            m = 2500;
            kl = 1.0;
            record_AFB(fout, m, kl);
        //}
    //}
    fclose(fout);

    return;
}
   




