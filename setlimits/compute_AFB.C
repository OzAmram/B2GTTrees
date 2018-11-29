
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <iostream>
#include <cstring>
#include <algorithm>
#include <vector>

#include "TROOT.h"
#include "TFile.h"
#include "TObject.h"
#include "TH1F.h"
#include "TProfile.h"
#include "TStyle.h"
#include "TCanvas.h"
#include "TPostScript.h"
#include "TLorentzVector.h"
#include "TVector3.h"
#include "TFitter.h"
#include "TSystem.h"
#include "TComplex.h"
#include "TMath.h"
#include "TRandom3.h"


Double_t M_Z0 = 91.1876;
Double_t width_Z0 = 2.4952;
Double_t alpha = 1./127.9; //fine structure constant at Z mass
Double_t e2_nom = alpha * 4 * TMath::Pi();
Double_t x = 0.231; //sin^2(theta_w)
Double_t gz2 = e2_nom/(x*(1-x));
Double_t g_theta = sqrt(5 * e2_nom / (3 *(1 - x)));
Double_t gz = sqrt(gz2);

Double_t e2_renorm(Double_t s){
    Double_t scale  = TMath::Log(s/M_Z0/M_Z0);
    Double_t denom = 1. - (e2_nom * scale / (12. * TMath::Pi() * TMath::Pi()));
    return e2_nom/denom;
}

TComplex SM_Amplitude(Double_t s, int hand1, int hand2, int qtype){
    //s is inv. mass, hand1 = 0 = LEFT, =1=RIGHT, qtype = 0= u,c,t, =1 = d,s,b
    Double_t Q, C_V, C_A;
    
    if(qtype == 0){ 
        //u c t quark type
        Q = 2./3.;
        C_V = gz*(-0.25+ 2.*x /3.);
        C_A = gz*(-0.25);
    }
    else if(qtype == 1){ 
        Q = -1./3.;
        C_V = gz*(0.25 - x/3.);
        C_A = gz*(0.25);
    }
    Double_t C_L = C_V + C_A;
    Double_t C_R = C_V - C_A;

    Double_t C_V_lep = gz*(0.25 - x);
    Double_t C_A_lep = gz*(0.25);
    Double_t C_L_lep = C_V_lep + C_A_lep;
    Double_t C_R_lep = C_V_lep - C_A_lep;


    Double_t C_q = C_R*hand1 + C_L*(1-hand1);

    Double_t C_lep = C_R_lep*hand2 + C_L_lep*(1-hand2);

    TComplex term1 = TComplex(-Q*e2_nom, 0);
    TComplex term2_num = TComplex(s*C_q*C_lep,0);
    TComplex term2_denom = TComplex(s - M_Z0*M_Z0, M_Z0*width_Z0);

    TComplex term2 = term2_num/term2_denom;

    return term1 + term2;
}

Double_t SM_AFB(Double_t E, int qtype){
    // E is energy to compute AFB, qtype = 0 = u,c,t, =1= d,s,b
    TComplex A_LL = SM_Amplitude(E*E, 0,0,qtype); 
    TComplex A_LR = SM_Amplitude(E*E, 0,1,qtype); 
    TComplex A_RL = SM_Amplitude(E*E, 1,0,qtype); 
    TComplex A_RR = SM_Amplitude(E*E, 1,1,qtype); 

    Double_t num = 3. * (A_LL.Rho2() + A_RR.Rho2() - A_LR.Rho2() - A_RL.Rho2());
    Double_t denom = 4. * (A_LL.Rho2() + A_RR.Rho2() + A_LR.Rho2() + A_RL.Rho2());
    
    Double_t AFB = num/denom;

    if(fabs(AFB) >= 0.75) printf("Warning AFB too large for E=%.0f! \n", E);
    return AFB;
}

TComplex Zp_Amplitude(Double_t s, Double_t M_Zp, int hand1, int hand2, int qtype, int type=0){
    //type 0 = SSM, 1 = E6 Z_x (chi)
    TComplex SM_Amp = SM_Amplitude(s,hand1, hand2, qtype);

    Double_t Q, C_V, C_A;
    if(qtype == 0){ 
        //u c t quark type
        Q = 2./3.;
        if(type == 0){
            C_V = gz*(-0.25+ 2.*x /3.);
            C_A = gz*(-0.25);
        }
        else{
            C_V = 0.;
            C_A = -g_theta/2./sqrt(10);
        }

    }
    else{ 
        Q = -1./3.;
        if(type == 0){
            C_V = gz*(0.25 - x/3.);
            C_A = gz*(0.25);
        }
        else{
            C_V = -g_theta/sqrt(10);
            C_A = g_theta/2./sqrt(10);
        }
    }


    Double_t C_L = C_V + C_A;
    Double_t C_R = C_V - C_A;

    Double_t C_V_lep, C_A_lep;
    if(type == 0){
        C_V_lep = gz*(0.25 - x);
        C_A_lep = gz*(0.25);
    }
    else{
        C_V_lep = g_theta/sqrt(10);
        C_A_lep = g_theta/2./sqrt(10);
    }

    Double_t C_L_lep = C_V_lep + C_A_lep;
    Double_t C_R_lep = C_V_lep - C_A_lep;


    Double_t C_q = C_R*hand1 + C_L*(1-hand1);
    Double_t C_lep = C_R_lep*hand2 + C_L_lep*(1-hand2);

    Double_t alpha_new = e2_renorm(M_Zp*M_Zp) / 4. / TMath::Pi();
    //printf("New old alpha %.1f %.1f \n", 1./alpha_new, 1./alpha);
    Double_t width_Zp = 5./2. * alpha_new * M_Zp / (1-x);

    TComplex term3_num = TComplex(s*C_q*C_lep,0);
    TComplex term3_denom = TComplex(s - M_Zp*M_Zp, M_Zp*width_Zp);

    TComplex term3 = term3_num/term3_denom;
    return SM_Amp + term3;
}

Double_t Zp_AFB(Double_t E, Double_t M_Zp, int qtype, int type =0){
    // E is energy to compute AFB, qtype = 0 = u,c,t, =1= d,s,b
    // type is type of Z' model (0 = SSM, 1 = E6 Z_x
    TComplex A_LL = Zp_Amplitude(E*E, M_Zp, 0,0,qtype, type); 
    TComplex A_LR = Zp_Amplitude(E*E, M_Zp, 0,1,qtype, type); 
    TComplex A_RL = Zp_Amplitude(E*E, M_Zp, 1,0,qtype, type); 
    TComplex A_RR = Zp_Amplitude(E*E, M_Zp, 1,1,qtype, type); 

    Double_t num = 3. * (A_LL.Rho2() + A_RR.Rho2() - A_LR.Rho2() - A_RL.Rho2());
    Double_t denom = 4. * (A_LL.Rho2() + A_RR.Rho2() + A_LR.Rho2() + A_RL.Rho2());
    
    Double_t AFB = num/denom;

    if(fabs(AFB) >= 0.75) printf("Warning AFB too large for E=%.0f! \n", E);
    return AFB;
}



Double_t SM_xsec(Double_t E, Double_t c_s, int qtype){
    TComplex A_LL = SM_Amplitude(E*E, 0,0,qtype); 
    TComplex A_LR = SM_Amplitude(E*E, 0,1,qtype); 
    TComplex A_RL = SM_Amplitude(E*E, 1,0,qtype); 
    TComplex A_RR = SM_Amplitude(E*E, 1,1,qtype); 

    Double_t norm = 1./128./TMath::Pi()/E/E;

    Double_t term1 = (A_LL.Rho2() + A_RR.Rho2()) * pow(1. + c_s, 2);
    Double_t term2 = (A_LR.Rho2() + A_LR.Rho2()) * pow(1. - c_s, 2);

    return norm * (term1 + term2);
}


Double_t Zp_xsec(Double_t E, Double_t c_s, Double_t M_Zp, int qtype, int type =0){
    TComplex A_LL = Zp_Amplitude(E*E, M_Zp, 0,0,qtype, type); 
    TComplex A_LR = Zp_Amplitude(E*E, M_Zp, 0,1,qtype, type); 
    TComplex A_RL = Zp_Amplitude(E*E, M_Zp, 1,0,qtype, type); 
    TComplex A_RR = Zp_Amplitude(E*E, M_Zp, 1,1,qtype, type); 

    Double_t norm = 1./128./TMath::Pi()/E/E;

    Double_t term1 = (A_LL.Rho2() + A_RR.Rho2()) * pow(1. + c_s, 2);
    Double_t term2 = (A_LR.Rho2() + A_LR.Rho2()) * pow(1. - c_s, 2);

    return norm * (term1 + term2);
}



