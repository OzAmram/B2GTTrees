
#include "madgraph_lhe_reader.C"



void record_ZChi_AFBs(){

    FILE *fout = fopen("ZChi_AFBs.txt", "a");
    int m = 2500;
    const int n_bins = 6;
    Double_t AFB_ZChi[6];
    //files for events for different mass bins (150,200,250,350,500,700, inf)
    double mbins[] = {150., 200., 250., 350., 500., 700., 100000.};
    char files[n_bins][140];
    sprintf(files[0], "/home/ozamram/Documents/Research/Generators/MG5_aMC_v2_6_2/ZChi_to_mumuj_M2500/Events/run_04/events.lhe");
    sprintf(files[1], "/home/ozamram/Documents/Research/Generators/MG5_aMC_v2_6_2/ZChi_to_mumuj_M2500/Events/run_05/events.lhe");
    sprintf(files[2], "/home/ozamram/Documents/Research/Generators/MG5_aMC_v2_6_2/ZChi_to_mumuj_M2500/Events/run_06/events.lhe");
    sprintf(files[3], "/home/ozamram/Documents/Research/Generators/MG5_aMC_v2_6_2/ZChi_to_mumuj_M2500/Events/run_07/events.lhe");
    sprintf(files[4], "/home/ozamram/Documents/Research/Generators/MG5_aMC_v2_6_2/ZChi_to_mumuj_M2500/Events/run_08/events.lhe");
    sprintf(files[5], "/home/ozamram/Documents/Research/Generators/MG5_aMC_v2_6_2/ZChi_to_mumuj_M2500/Events/run_09/events.lhe");

    char file1[] = "/home/ozamram/Documents/Research/Generators/MG5_aMC_v2_6_2/ZChi_test/Events/run_02/unweighted_events.lhe";
    cout << get_AFB(file1, 700, 10000,-1, true) << endl;

    /*
    for(int i=0; i<n_bins; i++){
        AFB_ZChi[i] = get_AFB(files[i], mbins[i], mbins[i+1],-1, true);
    }

    fprintf(fout, "M_Zchi=%i: %.3f %.3f %.3f %.3f %.3f %.3f \n", m, 
            AFB_ZChi[0],AFB_ZChi[1],AFB_ZChi[2], AFB_ZChi[3], AFB_ZChi[4], AFB_ZChi[5]);
            */
    fclose(fout);

    return;
}
   




