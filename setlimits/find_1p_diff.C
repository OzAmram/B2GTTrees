
#include "compute_AFB.C"


void find_1p_diff(){
    Double_t ZP_M = 3000.;
    Double_t SM_afb = 0.6;
    Double_t ZP_afb = 0.5;
    Double_t M_step = 10.;
    SM_afb = SM_AFB(1000., 0);

    while(abs(SM_afb-ZP_afb) > 0.01){
        ZP_M += M_step;

        ZP_afb = Zp_AFB(1000.,ZP_M, 0);
    }
    printf("Limit found. ZP of M = %.0f has AFB diff of %.4f at E = 1000 GeV \n", ZP_M, abs(SM_afb - ZP_afb));
    printf("Has AFB diff of %.4f at E = 700 GeV, diff of %.4f at E = 500 GeV \n", abs(Zp_AFB(700, ZP_M, 0) - SM_AFB(700,0)), abs(Zp_AFB(500, ZP_M, 0) - SM_AFB(500,0)));
    return ;
}


