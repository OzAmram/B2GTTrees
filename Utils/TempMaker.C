#include "TempMaker.h"

TempMaker::TempMaker(TTree *t, bool isdata=false, int y  = 2016){
    is_data = isdata;
    year = y;
    t_in = t;
}

void TempMaker::setup(){

    nEntries  =  t_in->GetEntries();
    if(year == 2016) el_lumi=el_lumi16;
    if(year == 2017) el_lumi=el_lumi17;
    if(year == 2018) el_lumi=el_lumi18;


    t_in->SetBranchAddress("m", &m);
    t_in->SetBranchAddress("xF", &xF);
    t_in->SetBranchAddress("cost", &cost);
    t_in->SetBranchAddress("jet1_CMVA", &jet1_cmva);
    t_in->SetBranchAddress("jet2_CMVA", &jet2_cmva);
    t_in->SetBranchAddress("met_pt", &met_pt);
    t_in->SetBranchAddress("jet1_pt", &jet1_pt);
    t_in->SetBranchAddress("jet2_pt", &jet2_pt);
    t_in->SetBranchAddress("jet1_eta", &jet1_eta);
    t_in->SetBranchAddress("jet2_eta", &jet2_eta);
    t_in->SetBranchAddress("has_no_bjets", &has_no_bjets);

    if(do_muons){
        t_in->SetBranchAddress("mu_p", &lep_p);
        t_in->SetBranchAddress("mu_m", &lep_m);
        t_in->SetBranchAddress("mu1_pt", &mu1_pt);
        t_in->SetBranchAddress("mu1_eta", &mu1_eta);
        t_in->SetBranchAddress("mu2_pt", &mu2_pt);
        t_in->SetBranchAddress("mu2_eta", &mu2_eta);
        t_in->SetBranchAddress("mu_p_SF", &mu_p_SF);
        t_in->SetBranchAddress("mu_m_SF", &mu_m_SF);
        t_in->SetBranchAddress("mu1_charge", &mu1_charge);
        t_in->SetBranchAddress("mu2_charge", &mu2_charge);
    }

    if(do_electrons){

        t_in->SetBranchAddress("el_p", &lep_p);
        t_in->SetBranchAddress("el_m", &lep_m);
        t_in->SetBranchAddress("el1_pt", &el1_pt);
        t_in->SetBranchAddress("el1_eta", &el1_eta);
        t_in->SetBranchAddress("el2_pt", &el2_pt);
        t_in->SetBranchAddress("el2_eta", &el2_eta);
        t_in->SetBranchAddress("el1_charge", &el1_charge);
        t_in->SetBranchAddress("el2_charge", &el2_charge);
    }
    if(is_one_iso){
        if(do_muons) t_in->SetBranchAddress("iso_mu", &iso_lep);
        if(do_electrons) t_in->SetBranchAddress("iso_el", &iso_lep);
    }


    if(!is_data){
        if(is_gen_level){ 
            t_in->SetBranchAddress("cost_st", &cost_st);
            t_in->SetBranchAddress("pdf_weights", &pdf_weights);
        }
        t_in->SetBranchAddress("nJets", &nJets);
        t_in->SetBranchAddress("gen_weight", &gen_weight);
        t_in->SetBranchAddress("pu_SF", &pu_SF);
        t_in->SetBranchAddress("pu_SF_up", &pu_SF_up);
        t_in->SetBranchAddress("pu_SF_down", &pu_SF_down);
        t_in->SetBranchAddress("mu_R_up", &mu_R_up);
        t_in->SetBranchAddress("mu_R_down", &mu_R_down);
        t_in->SetBranchAddress("mu_F_up", &mu_F_up);
        t_in->SetBranchAddress("mu_F_down", &mu_F_down);
        t_in->SetBranchAddress("mu_RF_up", &mu_RF_up);
        t_in->SetBranchAddress("mu_RF_down", &mu_RF_down);
        t_in->SetBranchAddress("pu_NtrueInt", &pu_NtrueInt);
        t_in->SetBranchAddress("jet1_flavour", &jet1_flavour);
        t_in->SetBranchAddress("jet2_flavour", &jet2_flavour);
        t_in->SetBranchAddress("alpha_up", &alphaS_up);
        t_in->SetBranchAddress("alpha_down", &alphaS_down);

        if(do_muons){
            if(is_gen_level){
                t_in->SetBranchAddress("gen_mu_p", &gen_lep_m);
                t_in->SetBranchAddress("gen_mu_m", &gen_lep_p);
            }
            t_in->SetBranchAddress("era1_HLT_SF", &era1_HLT_SF);
            t_in->SetBranchAddress("era1_iso_SF", &era1_iso_SF);
            t_in->SetBranchAddress("era1_id_SF", &era1_id_SF);
            t_in->SetBranchAddress("era2_HLT_SF", &era2_HLT_SF);
            t_in->SetBranchAddress("era2_iso_SF", &era2_iso_SF);
            t_in->SetBranchAddress("era2_id_SF", &era2_id_SF);
            t_in->SetBranchAddress("mu_p_SF_up", &mu_p_SF_up);
            t_in->SetBranchAddress("mu_m_SF_up", &mu_m_SF_up);
            t_in->SetBranchAddress("mu_p_SF_down", &mu_p_SF_down);
            t_in->SetBranchAddress("mu_m_SF_down", &mu_m_SF_down);
        }

        if(do_electrons){

            if(is_gen_level){
                t_in->SetBranchAddress("gen_el_p", &gen_lep_p);
                t_in->SetBranchAddress("gen_el_m", &gen_lep_m);
            }
            t_in->SetBranchAddress("el_id_SF", &el_id_SF);
            t_in->SetBranchAddress("el_HLT_SF", &el_HLT_SF);
            t_in->SetBranchAddress("el_reco_SF", &el_reco_SF);
            if(do_elScale_sys || do_elSmear_sys){
                if(do_elScale_sys >0){
                    if(sys_label.find("elScaleStat") != string::npos){
                        t_in->SetBranchAddress("elp_scale_stat_up", &elp_rescale);
                        t_in->SetBranchAddress("elm_scale_stat_up", &elm_rescale);
                    }
                    if(sys_label.find("elScaleSyst") != string::npos){
                        t_in->SetBranchAddress("elp_scale_syst_up", &elp_rescale);
                        t_in->SetBranchAddress("elm_scale_syst_up", &elm_rescale);
                    }
                    if(sys_label.find("elScaleGain") != string::npos){
                        t_in->SetBranchAddress("elp_scale_gain_up", &elp_rescale);
                        t_in->SetBranchAddress("elm_scale_gain_up", &elm_rescale);
                    }
                }
                else if(do_elScale_sys < 0){
                    if(sys_label.find("elScaleStat") != string::npos){
                        t_in->SetBranchAddress("elp_scale_stat_down", &elp_rescale);
                        t_in->SetBranchAddress("elm_scale_stat_down", &elm_rescale);
                    }
                    if(sys_label.find("elScaleSyst") != string::npos){
                        t_in->SetBranchAddress("elp_scale_syst_down", &elp_rescale);
                        t_in->SetBranchAddress("elm_scale_syst_down", &elm_rescale);
                    }
                    if(sys_label.find("elScaleGain") != string::npos){
                        t_in->SetBranchAddress("elp_scale_gain_down", &elp_rescale);
                        t_in->SetBranchAddress("elm_scale_gain_down", &elm_rescale);
                    }
                }
                else if(do_elSmear_sys < 0){
                    t_in->SetBranchAddress("elp_smear_down", &elp_rescale);
                    t_in->SetBranchAddress("elm_smear_down", &elm_rescale);
                }
                else if(do_elSmear_sys > 0){
                    t_in->SetBranchAddress("elp_smear_up", &elp_rescale);
                    t_in->SetBranchAddress("elm_smear_up", &elm_rescale);
                }

            }

        }
    }
}


void TempMaker::setup_systematic(const string &s_label){
    sys_label = s_label;

    if(!sys_label.empty()){
        //doing some systematic
        if(sys_label.find("Up") != string::npos){
            sys_shift = 1;
        }
        else if(sys_label.find("Down") != string::npos){
            sys_shift = -1;
        }

        else{
            printf("systematic label not empty, but doesn't have Up or Down \n");
        }
    }

    if(sys_shift !=0){
        if(sys_label.find("BTAG") != string::npos) do_btag_sys = sys_shift;
        else if(sys_label.find("Pu") != string::npos) do_pileup_sys = sys_shift;
        else if(sys_label.find("muHLT") != string::npos) do_muHLT_sys = sys_shift;
        else if(sys_label.find("muID") != string::npos) do_muID_sys = sys_shift;
        else if(sys_label.find("muISO") != string::npos) do_muISO_sys = sys_shift;
        else if(sys_label.find("muRC") != string::npos) do_muRC_sys = sys_shift;

        else if(sys_label.find("elID") != string::npos) do_elID_sys = sys_shift;
        else if(sys_label.find("elHLT") != string::npos) do_elHLT_sys = sys_shift;
        else if(sys_label.find("elRECO") != string::npos) do_elRECO_sys = sys_shift;
        else if(sys_label.find("elScale") != string::npos) do_elScale_sys = sys_shift;
        else if(sys_label.find("elSmear") != string::npos) do_elSmear_sys = sys_shift;


        else if(sys_label.find("RENORM") != string::npos && sys_shift > 0) systematic = &mu_R_up;
        else if(sys_label.find("RENORM") != string::npos && sys_shift < 0) systematic = &mu_R_down;
        else if(sys_label.find("FAC") != string::npos && sys_shift > 0) systematic = &mu_F_up;
        else if(sys_label.find("FAC") != string::npos && sys_shift < 0) systematic = &mu_F_down;
        else if(sys_label.find("alphaS") != string::npos && sys_shift < 0) systematic = &alphaS_down;
        else if(sys_label.find("alphaS") != string::npos && sys_shift > 0) systematic = &alphaS_up;
        else if(sys_label.find("alpha") != string::npos) systematic = &one;
        else if(sys_label.find("pdf") != string::npos){
            if(sys_shift > 0) sscanf(sys_label.c_str(), "_pdf%iUp", &do_pdf_sys);
            else sscanf(sys_label.c_str(), "_pdf%iDown", &do_pdf_sys);
            printf("Doing pdf sys %i \n", do_pdf_sys);
        }

        else printf("COULDN'T PARSE SYSTEMATIC %s !!! \n \n", sys_label.c_str());

    }
}


void TempMaker::getEvent(int i){
    t_in->GetEntry(i);
    not_cosmic = notCosmic(*lep_p, *lep_m);




}

void TempMaker::doCorrections(){
    if(do_muons){
        if(!do_RC || do_muRC_sys != 0){

            Double_t mu_p_corr, mu_m_corr;
            if(!do_RC){
                //undo RC 
                mu_p_corr = 1./mu_p_SF;
                mu_m_corr = 1./mu_m_SF;
            }
            else if(do_muRC_sys != 0){
                if(do_muRC_sys > 0){
                    mu_p_corr = mu_p_SF_up/mu_p_SF;
                    mu_m_corr = mu_m_SF_up/mu_m_SF;
                }
                else if(do_muRC_sys < 0){
                    mu_p_corr = mu_p_SF_down/mu_p_SF;
                    mu_m_corr = mu_m_SF_down/mu_p_SF;
                }
            }
            lep_p->SetPtEtaPhiE(lep_p->Pt() * mu_p_corr, lep_p->Eta(), lep_p->Phi(), mu_p_corr * lep_p->E());
            lep_m->SetPtEtaPhiE(lep_m->Pt() * mu_m_corr, lep_m->Eta(), lep_m->Phi(), mu_m_corr * lep_m->E());


        }
    }

    if(do_electrons){

        if(do_elScale_sys || do_elSmear_sys){
            lep_p->SetPtEtaPhiE(lep_p->Pt() * (elp_rescale), lep_p->Eta(), lep_p->Phi(), lep_p->E() *(elp_rescale));
            lep_m->SetPtEtaPhiE(lep_m->Pt() * (elm_rescale), lep_m->Eta(), lep_m->Phi(), lep_m->E() *(elm_rescale));
        }
    }
    TLorentzVector cm = *lep_p + *lep_m;
    m = cm.M();
    xF = compute_xF(cm);
    cost = get_cost(*lep_p, *lep_m);

    if(std::isnan(cost_st)){
        printf("Gen level cost is Nan, using reco level \n");
        cost_st = cost;
    }

}

Double_t TempMaker::getEvtWeight(){
    if(is_data) return 1.;


    if(do_pileup_sys == -1) pu_SF = pu_SF_down;
    if(do_pileup_sys == 1) pu_SF = pu_SF_up;
    if(do_pdf_sys != 0){
        if(sys_shift > 0) *systematic = pdf_weights[do_pdf_sys-1];
        if(sys_shift < 0) *systematic = 2. - pdf_weights[do_pdf_sys-1];
    }
    if(abs(*systematic) > 10.0 || abs(*systematic) < 0.01 || std::isnan(*systematic)){
        //printf("sys is %.4f  \n", *systematic);
        *systematic = 1.;
    }
    double base_weight = gen_weight * (*systematic) * pu_SF;
    jet1_b_weight = get_btag_weight(jet1_pt, jet1_eta, (Float_t) jet1_flavour , btag_effs, b_reader, do_btag_sys);
    jet2_b_weight = get_btag_weight(jet2_pt, jet2_eta, (Float_t) jet2_flavour , btag_effs, b_reader, do_btag_sys);
    if (nJets >= 1){
        base_weight *= jet1_b_weight;
    }
    if (nJets >= 2){
        base_weight *= jet2_b_weight;
    }


    if(do_muons){
        if(do_muHLT_sys)   era1_HLT_SF = get_HLT_SF(mu1_pt, mu1_eta, mu2_pt, mu2_eta, era1_SFs.HLT_SF, era1_SFs.HLT_MC_EFF, do_muHLT_sys);
        if(do_muHLT_sys)   era2_HLT_SF = get_HLT_SF(mu1_pt, mu1_eta, mu2_pt, mu2_eta, era2_SFs.HLT_SF, era2_SFs.HLT_MC_EFF, do_muHLT_sys);

        if(do_muISO_sys)   era1_iso_SF = get_mu_SF(mu1_pt, mu1_eta, year, era1_SFs.ISO_SF,  do_muISO_sys) * get_mu_SF(mu2_pt, mu2_eta, year, era1_SFs.ISO_SF,  do_muISO_sys);
        if(do_muISO_sys)   era2_iso_SF = get_mu_SF(mu1_pt, mu1_eta, year, era2_SFs.ISO_SF,  do_muISO_sys) * get_mu_SF(mu2_pt, mu2_eta, year, era2_SFs.ISO_SF,  do_muISO_sys);

        if(do_muID_sys)    era1_id_SF = get_mu_SF(mu1_pt, mu1_eta, year, era1_SFs.ID_SF,  do_muID_sys) * get_mu_SF(mu2_pt, mu2_eta, year, era1_SFs.ID_SF,  do_muID_sys);
        if(do_muID_sys)    era2_id_SF = get_mu_SF(mu1_pt, mu1_eta, year, era2_SFs.ID_SF,  do_muID_sys) * get_mu_SF(mu2_pt, mu2_eta, year, era2_SFs.ID_SF,  do_muID_sys);


        Double_t era1_weight = base_weight * era1_HLT_SF * era1_iso_SF * era1_id_SF ;
        Double_t era2_weight = base_weight * era2_HLT_SF * era2_iso_SF * era2_id_SF ;


        if(year==2016) evt_weight = 1000*(era1_weight*bcdef_lumi16 + era2_weight*gh_lumi16);
        if(year==2017) evt_weight = 1000*(era1_weight*mu_lumi17);
        if(year==2018) evt_weight = 1000*(era1_weight*mu_lumi18);
    }
    if(do_electrons){
        if(do_elID_sys) el_id_SF = get_el_SF(el1_pt, el1_eta, el_SF.ID_SF, do_elID_sys) * get_el_SF(el2_pt, el2_eta, el_SF.ID_SF, do_elID_sys);
        if(do_elRECO_sys) el_reco_SF = get_el_SF(el1_pt, el1_eta, el_SF.RECO_SF, do_elRECO_sys) * get_el_SF(el2_pt, el2_eta, el_SF.RECO_SF, do_elRECO_sys);
        if(do_elHLT_sys) el_HLT_SF = get_el_HLT_SF(el1_pt, el1_eta, el2_pt, el2_eta, el_SF.HLT_SF, el_SF.HLT_MC_EFF, do_elHLT_sys);




        evt_weight = base_weight * el_id_SF *el_reco_SF * el_HLT_SF * 1000. * el_lumi;

    }
    return evt_weight;

}

void TempMaker::finish(){
    t_in->ResetBranchAddresses();
}

TempMaker::~TempMaker(){}
