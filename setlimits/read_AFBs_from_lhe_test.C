
#include "madgraph_lhe_reader.C"




void read_AFBs_from_lhe_test(){
    //char fname[] = "/home/ozamram/Documents/Research/Generators/MG5_aMC_v2_6_2/Z_Chi_test/Events/run_02/events.lhe";
    char fname[] = "/home/ozamram/Documents/Research/Generators/MG5_aMC_v2_6_2/ZChi_test2/Events/run_01/events.lhe";
    Double_t AFB_uncorr = get_AFB(fname, 2950,3050 , -1,true);
    Double_t AFB_corr = get_AFB(fname, 2950, 3050, 5,true);

    printf("AFB corrected for changing dilution is %.3f, not corrected is %.3f \n", AFB_corr, AFB_uncorr);

}
   




