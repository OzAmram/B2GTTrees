#include "read_AFBs_from_file.C"

void read_test(){
    FILE *f1;
    //setup_AFB_file(f1);
    f1 = fopen("AFBs.txt", "r");
    if (f1 == NULL) printf("ERROR opening file! \n");
    Double_t AFB[6];
    read_AFBs(f1, AFB, 2200, 0.7);

    printf("Printing AFBs: \n");
    for(int i=0; i<6; i++){
        printf(" %.3lf", AFB[i]);
    }
    printf("\n");

    return;
}
