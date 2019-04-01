#include "pythia_lhe_reader.C"

Double_t AFB_utype(Double_t *x, Double_t *params){
    Double_t E = x[0];
    return SM_AFB(E, 0);
}

Double_t AFB_dtype(Double_t *x, Double_t *params){
    Double_t E = x[0];
    return SM_AFB(E, 1);
}

Double_t Zp_AFB_utype(Double_t *x, Double_t *params){
    Double_t E = x[0];
    return Zp_AFB(E, 2000., 0);
}

Double_t Zp_AFB_dtype(Double_t *x, Double_t *params){
    Double_t E = x[0];
    return Zp_AFB(E, 2000., 1);
}

Double_t Zp_E6_AFB_utype(Double_t *x, Double_t *params){
    Double_t E = x[0];
    return Zp_AFB(E, 2000., 1, 1);
}

Double_t SM_xsec_slice(Double_t *x, Double_t *params){
    Double_t E = x[0];
    return SM_xsec(E,0., 1);
}

Double_t Zp_xsec_slice(Double_t *x, Double_t *params){
    Double_t E = x[0];
    return Zp_xsec(E,0., 2000., 1);
}

Double_t Zp_xsec_cost_1(Double_t *x, Double_t *params){
    Double_t cost = x[0];
    return Zp_xsec(500.,cost, 1500., 1);
}

Double_t Zp_xsec_cost_2(Double_t *x, Double_t *params){
    Double_t cost = x[0];
    return Zp_xsec(500.,cost, 600., 1);
}
void test1(){
    TComplex c1 = SM_Amplitude(100.*100., 0,0,0);
    TComplex c2 = SM_Amplitude(200.*200., 0,0,0);

    printf("Real %.2e Im %.2e \n", c1.Re(), c1.Im());
    printf("Real %.2e Im %.2e \n", c2.Re(), c2.Im());

    Double_t AFB1 = SM_AFB(600., 1);
    Double_t AFB2 = SM_AFB(600., 0);

    Double_t AFB3 = Zp_AFB(600.,1500., 1);
    Double_t AFB4 = Zp_AFB(600.,1500., 0);

    Double_t AFB5 = SM_AFB(91.18, 1);
    Double_t AFB6 = SM_AFB(91.18, 0);

    Double_t AFB7 = SM_AFB(60., 1);
    Double_t AFB8 = SM_AFB(60., 0);

    //Double_t alpha_new = e2_renorm(1000.*1000.) / 4. / TMath::Pi();
    //printf("New old alpha %.1f %.1f \n", 1./alpha_new, 1./alpha);

    printf("AFBs are SM 225 %.3f %.3f \n"
            "Zp 225 %.3f %.3f \n" 
            "91.18 %.2f %.2f \n"
            "60 %.2f %.2f \n", 
            AFB1, AFB2, AFB3, AFB4, AFB5, AFB6, AFB7, AFB8);


    TF1 *AFB_u = new TF1("AFB", AFB_utype, 1.,2000.);
    TF1 *AFB_d = new TF1("AFB_d", AFB_dtype, 1.,2000.);

    TF1 *Zp_AFB_u = new TF1("Zp_AFB_u", Zp_AFB_utype, 1.,2000.);
    TF1 *Zp_AFB_d = new TF1("Zp_AFB_d", Zp_AFB_dtype, 1.,2000.);

    TF1 *Zp_E6_AFB_u = new TF1("Zp_AFB_E6_u", Zp_E6_AFB_utype, 1.,2000.);
    AFB_u->SetLineColor(kBlue);
    Zp_AFB_u->SetLineColor(kGreen);
    Zp_E6_AFB_u->SetLineColor(kRed);
    TCanvas *can = new TCanvas("c1", "Histograms", 200, 10, 900, 700);
    AFB_u->Draw();
    AFB_u->GetXaxis()->SetTitle("Dilepton Mass (GeV)");
    Zp_AFB_u->Draw("same");
    Zp_E6_AFB_u->Draw("same");

    gStyle->SetLegendBorderSize(0);
    TLegend *leg1 = new TLegend(0.5, 0.65, 0.75, 0.8);
    leg1->AddEntry(AFB_u, "SM u-type AFB", "l");
    leg1->AddEntry(Zp_AFB_u, "SM + SSM Z' (2 Tev) u-type AFB", "l");
    leg1->AddEntry(Zp_E6_AFB_u, "SM + Z_{#chi} (2 Tev) u-type AFB", "l");
    leg1->Draw();



    TCanvas *can2 = new TCanvas("c2", "Histograms", 200, 10, 900, 700);
    TF1 *sm_xsec = new TF1("xsec1", SM_xsec_slice, 1., 2000.);
    TF1 *zp_xsec = new TF1("xsec2", Zp_xsec_slice, 1., 2000.);
    sm_xsec->SetLineColor(kBlue);
    zp_xsec->SetLineColor(kGreen);

    sm_xsec->Draw();
    zp_xsec->Draw("same");
    can2->SetLogy();

    /*
    TCanvas *can3 = new TCanvas("c3", "Histograms", 200, 10, 900, 700);
    TF1 *zp1_cost = new TF1("zp1_cost", Zp_xsec_cost_1, -1., 1.);
    TF1 *zp2_cost = new TF1("zp2_cost", Zp_xsec_cost_2, -1., 1.);
    zp1_cost->SetLineColor(kBlue);
    zp2_cost->SetLineColor(kGreen);

    zp1_cost->Draw();
    zp2_cost->Draw("same");
    */

    Double_t Zp_AFB1 = get_AFB(2500., 0);
    Double_t Zp_AFB2 = get_AFB(3000., 0);
    Double_t Zp_AFB3 = get_AFB(2500., 5);
    Double_t Zp_AFB4 = get_AFB(3000., 5);
    printf("AFBs are: 150-200 %.2f %.2f \n" "700+ %.2f %.2f \n", Zp_AFB1, Zp_AFB2, Zp_AFB3, Zp_AFB4);
    
}


