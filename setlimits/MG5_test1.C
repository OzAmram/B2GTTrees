#include "madgraph_lhe_reader.C"


void MG5_test1(){
    for(int i=0; i<6; i++){
        Double_t AFB = get_AFB(1950, 0.2, i, true );
        //printf("AFB for %i bin is %.3f \n", i, AFB);
    }
    return;
}
