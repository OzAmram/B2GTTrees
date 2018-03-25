#include <stdio.h>
#include <cstdio>
#include <stdlib.h>
#include <math.h>
#include <iostream>
#include <cstring>
#include <algorithm>
#include <vector>


void read_AFBs(FILE *f1, Double_t *AFBs, int M_Zp, Double_t kL){
    rewind(f1);
    char match_str[80];
    sprintf(match_str, "M_Zp=%i kL=%.2f", M_Zp, kL);
    char line[120];
    while(fgets(line, 120, f1)){
        //printf("%s \n", line);
        if(strstr(line,match_str) != NULL) {
            //printf("Line found! %s \n", line);
            sscanf(line, "M_Zp=%i kL=%lf: %lf %lf %lf %lf %lf %lf \n", &M_Zp, &kL, 
                    &AFBs[0], &AFBs[1], &AFBs[2], &AFBs[3], &AFBs[4], &AFBs[5]);
            return;
        }
    }
    printf("Couldn't find AFBs for M_Zp=%i kL=%.2f \n", M_Zp, kL);
    AFBs[0] = AFBs[1] = AFBs[2] = AFBs[3] = AFBs[4] = AFBs[5] = -1.;
    return;
}

