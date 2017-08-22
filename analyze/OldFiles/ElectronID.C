#include <stdio.h>
#include <stdlib.h>
#include <math.h>

bool get_el_id(Float_t eta, Float_t full5x5, Float_t dEtaInSeed, Float_t dPhiIn, 
        Float_t HoE, Float_t iso, Float_t ooEmooP, Float_t missHits, Float_t matchConVeto){
    //see
    //https://twiki.cern.ch/twiki/bin/view/CMS/CutBasedElectronIdentificationRun2
    //using medium ID criteria
    Float_t full5x5_med, dEtaInSeed_med, dPhiIn_med, HoE_med, iso_med, ooEmooP_med, missHits_med;
    eta = std::abs(eta);

    Float_t eta_barrel = 1.479;
    if(eta < eta_barrel){
        full5x5_med = 0.00998;
        dEtaInSeed_med = 0.00311;
        dPhiIn_med = 0.103;
        HoE_med = 0.253;
        iso_med = 0.0695;
        ooEmooP_med = 0.134;
        missHits_med =1;
    }
    else{
        full5x5_med = 0.0298;
        dEtaInSeed_med = 0.00609;
        dPhiln_med = 0.045;
        HoE_med = 0.0878;
        iso_med = 0.0821;
        ooEmooP_med = 0.13;
        missHits_med =1;
    }
    return full5x5 < full5x5_med && 
        std::abs(dEtaInSeed) < dEtaInSeed_med &&
        std::abs(dPhiIn) < dPhiIn_med &&
        HoE < HoE_med &&
        iso < iso_med &&
        ooEmooP < ooEmooP_med &&
        missHits <= missHits_med &&
        matchConVeto > 0.1;
}
               

