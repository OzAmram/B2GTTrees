
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

Double_t Zp_E6_AFB_utype(Double_t *x, Double_t *params){
    Double_t E = x[0];
    return Zp_AFB(E, 2000., 0, 1);
}

void draw_AFBs(){
    Double_t m_min = 10.;
    Double_t m_max = 2000.;
    TF1 *AFB_u = new TF1("AFB (u-type)", AFB_utype, m_min,m_max);

    TF1 *Zp_AFB_u = new TF1("Zp_AFB_u", Zp_AFB_utype, m_min,m_max);

    TF1 *Zp_E6_AFB_u = new TF1("Zp_AFB_E6_u", Zp_E6_AFB_utype, m_min,m_max);
    AFB_u->SetLineColor(kBlue);
    Zp_AFB_u->SetLineColor(kGreen);
    Zp_E6_AFB_u->SetLineColor(kRed);
    TCanvas *can = new TCanvas("c1", "Histograms", 200, 10, 900, 700);
    AFB_u->Draw();
    Zp_AFB_u->Draw("same");
    //Zp_E6_AFB_u->Draw("same");
    AFB_u->GetYaxis()->SetRangeUser(-0.75, 0.75);
    AFB_u->GetXaxis()->SetTitle("Dilepton Mass (GeV)");
    AFB_u->GetYaxis()->SetTitle("AFB");
    //can->SetLogx();

    gStyle->SetLegendBorderSize(0);
    TLegend *leg1 = new TLegend(0.5, 0.65, 0.75, 0.8);
    leg1->AddEntry(AFB_u, "SM", "l");
    leg1->AddEntry(Zp_AFB_u, "SM + Z'_{SSM} (1.5 Tev)", "l");
    //leg1->AddEntry(Zp_E6_AFB_u, "SM + Z'_{#chi} (1.5 Tev)", "l");
    leg1->Draw();
}
