#include "draw_generator_cmp.C"




void draw_phot_pdf_cmp(){
    gStyle->SetOptStat(0);
    TFile *f_mad= TFile::Open("/home/ozamram/Documents/Research/B2GTTrees/generator_stuff/root_files/photOnly_lux_may2.root");
    //TFile * f_mad = TFile::Open("/home/ozamram/Documents/Research/B2GTTrees/generator_stuff/root_files/phot_induced_aug21.root");
    TTree *t_mad = (TTree *)f_mad->Get("T_lhe");

    TFile *f_pwg = TFile::Open("/home/ozamram/Documents/Research/B2GTTrees/generator_stuff/root_files/photInd_lux_may2.root");
    TTree *t_pwg = (TTree *) f_pwg->Get("T_lhe");

    char title[80] = "LO vs. NLO QED, Lux (M > 150)";
    TH1F *h_pwg_cost_st = new TH1F("h_pwg_cost_st", title, 20, -1., 1.);
    TH1F *h_pwg_cost_r = new TH1F("h_pwg_cost_r", title, 20, -1., 1.);
    TH1F *h_pwg_pt = new TH1F("h_pwg_pt", title, 20, 0., 200.);
    TH1F *h_pwg_xf = new TH1F("h_pwg_xf", title, 20, 0., 0.5);

    TH1F *h_mad_cost_st = new TH1F("h_mad_cost_st", title, 20, -1., 1.);
    TH1F *h_mad_cost_r = new TH1F("h_mad_cost_r", title, 20, -1., 1.);
    TH1F *h_mad_pt = new TH1F("h_mad_pt", title, 20, 0., 200.);
    TH1F *h_mad_xf = new TH1F("h_mad_xf", title, 20, 0., 0.5);

    float m_low = 150.;
    float m_high = 10000.;
    bool phot_ind = true;

    make_gen_cost(t_mad,  h_mad_cost_st, h_mad_cost_r, h_mad_pt, h_mad_xf, m_low, m_high, phot_ind);
    make_gen_cost(t_pwg,  h_pwg_cost_st, h_pwg_cost_r, h_pwg_pt, h_pwg_xf, m_low, m_high, phot_ind);

    h_mad_cost_st->Scale(1./h_mad_cost_st->Integral());
    h_pwg_cost_st->Scale(1./h_pwg_cost_st->Integral());

    h_mad_cost_st->SetLineColor(kRed);
    h_pwg_cost_st->SetLineColor(kBlue);

    h_mad_cost_st->SetLineWidth(3);
    h_pwg_cost_st->SetLineWidth(3);

    h_mad_cost_r->Scale(1./h_mad_cost_r->Integral());
    h_pwg_cost_r->Scale(1./h_pwg_cost_r->Integral());

    h_mad_cost_r->SetLineColor(kRed);
    h_pwg_cost_r->SetLineColor(kBlue);

    h_mad_cost_r->SetLineWidth(3);
    h_pwg_cost_r->SetLineWidth(3);

    h_mad_pt->Scale(1./h_mad_pt->Integral());
    h_pwg_pt->Scale(1./h_pwg_pt->Integral());

    h_mad_pt->SetLineColor(kRed);
    h_pwg_pt->SetLineColor(kBlue);

    h_mad_pt->SetLineWidth(3);
    h_pwg_pt->SetLineWidth(3);

    h_mad_xf->Scale(1./h_mad_xf->Integral());
    h_pwg_xf->Scale(1./h_pwg_xf->Integral());

    h_mad_xf->SetLineColor(kRed);
    h_pwg_xf->SetLineColor(kBlue);

    h_mad_xf->SetLineWidth(3);
    h_pwg_xf->SetLineWidth(3);

    
    make_ratio_plot("lo_vs_nlo_lux_cost_cmp.pdf", h_mad_cost_st, "lo",h_pwg_cost_st, "nlo", "lo/nlo", "cos(#theta_{*})", false, true);
    make_ratio_plot("lo_vs_nlo_lux_xf_cmp.pdf", h_mad_xf, "lo",h_pwg_xf, "nlo", "lo/nlo", "x_{F}", false, true);
    make_ratio_plot("lo_vs_nlo_lux_pt_cmp.pdf", h_mad_pt, "lo",h_pwg_pt, "nlo", "lo/nlo", "p_{T}", false, true);
    make_ratio_plot("lo_vs_nlo_cost_r_cmp.pdf", h_mad_cost_r, "lo",h_pwg_cost_r, "nlo", "lo/nlo", "cos(#theta_{r})", false, true);

    return;
}

