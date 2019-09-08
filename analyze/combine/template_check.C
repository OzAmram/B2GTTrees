#include "make_templates.C"


void template_check(){
    
        int year = 2017;
        init(year);
        setup_all_SFs(year);
        string sys_label = string("");

        char title[100];
        h_mumu_sym = new TH2F(title, "Symmetric template of mc",
                n_xf_bins, xf_bins, n_cost_bins, cost_bins);
        h_mumu_sym->SetDirectory(0);
        sprintf(title, "mumu%i_alpha%s", year %2000, sys_label.c_str());
        h_mumu_alpha = new TH2F(title, "Gauge boson polarization template of mc",
                n_xf_bins, xf_bins, n_cost_bins, cost_bins);
        h_mumu_alpha->SetDirectory(0);
        sprintf(title, "mumu%i_asym%s", year %2000, sys_label.c_str());
        h_mumu_asym = new TH2F(title, "Asymmetric template of mc",
                n_xf_bins, xf_bins, n_cost_bins, cost_bins);
        h_mumu_asym->SetDirectory(0);

        TH2F * h_mumu_plain = new TH2F("mumu_plain", "", n_xf_bins, xf_bins, n_cost_bins, cost_bins);

        h_elel_sym = new TH2F(title, "Symmetric template of mc",
                n_xf_bins, xf_bins, n_cost_bins, cost_bins);
        h_elel_sym->SetDirectory(0);
        sprintf(title, "elel%i_alpha%s", year %2000, sys_label.c_str());
        h_elel_alpha = new TH2F(title, "Gauge boson polarization template of mc",
                n_xf_bins, xf_bins, n_cost_bins, cost_bins);
        h_elel_alpha->SetDirectory(0);
        sprintf(title, "elel%i_asym%s", year %2000, sys_label.c_str());
        h_elel_asym = new TH2F(title, "Asymmetric template of mc",
                n_xf_bins, xf_bins, n_cost_bins, cost_bins);
        h_elel_asym->SetDirectory(0);

        TH2F * h_elel_plain = new TH2F("elel_plain", "", n_xf_bins, xf_bins, n_cost_bins, cost_bins);

        double m_low = 200.;
       
        double m_high = 250.;
        Double_t alpha_denom = alphas_denom[1];

        bool ss = false;

        gen_mc_template(t_elel_mc, alpha_denom, h_elel_sym, h_elel_asym, h_elel_alpha, year, m_low, m_high, FLAG_ELECTRONS, do_RC, "");
        //h_elel_sym->Print("all");
        //h_elel_asym->Print("all");
        //h_elel_alpha->Print("all");


        gen_mc_template(t_mumu_mc, alpha_denom, h_mumu_sym, h_mumu_asym, h_mumu_alpha, year, m_low, m_high, FLAG_MUONS, do_RC, "");
        TTree *mumu_ts[1] = {t_mumu_mc};
        printf("Making mumu back \n");
        gen_combined_background_template(1, mumu_ts, h_mumu_plain, year, m_low, m_high, FLAG_MUONS,  do_RC, ss, "");
        auto h1_mumu_back = convert2d(h_mumu_plain);

        TTree *elel_ts[1] = {t_elel_mc};
        printf("Making elel back \n");
        gen_combined_background_template(1, elel_ts, h_elel_plain, year, m_low, m_high, FLAG_ELECTRONS,  do_RC, ss, "");
        auto h1_elel_back = convert2d(h_elel_plain);



        TH2F *h_elel_dy, *h_mumu_dy;
        float afb = 0.59;

        printf("trying to add");

        double norm = 3./4./(2.+alpha_denom);

        h_mumu_alpha->Scale(norm * alpha_denom);

        auto h_mumu_pl = *h_mumu_sym + *h_mumu_asym;
        auto h_mumu_mn = *h_mumu_sym - *h_mumu_asym;
        h_mumu_pl.Scale(0.5);
        h_mumu_mn.Scale(0.5);
        h_mumu_dy = (TH2F *) h_mumu_pl.Clone("h_mumu_dy");
        h_mumu_dy->Add(&h_mumu_pl, &h_mumu_mn, (norm + afb), (norm - afb));


        auto h1_mumu_lo = convert2d(h_mumu_dy);
        auto h1_mumu_alpha = convert2d(h_mumu_alpha);
        TH1F *h1_mumu_templ = (TH1F *) h1_mumu_lo->Clone("h1_mumu_templ");
        h1_mumu_templ->Add(h1_mumu_alpha);
        
        printf("MuMu Back is %.2f.Comb Temp is %.2f  LO is %.2f h_alpha is %.2f  \n", 
                h1_mumu_back->Integral(),  h1_mumu_templ->Integral(), h1_mumu_lo->Integral(), h1_mumu_alpha->Integral());


        h1_mumu_back->SetLineColor(kGreen +3);
        h1_mumu_templ->SetLineColor(kBlack);
        h1_mumu_lo->SetLineColor(kRed+1);


        TCanvas *c_mumu = new TCanvas("c_mumu", "Histograms", 200, 10, 900, 700);
        h1_mumu_lo->Draw("hist");
        h1_mumu_templ->Draw("hist same ");
        h1_mumu_back->Draw("hist same");

        TLegend *leg = new TLegend(0.15, 0.3, 0.3, 0.15);
        leg->AddEntry(h1_mumu_templ, "Template", "f");
        leg->AddEntry(h1_mumu_lo, "LO Template", "f");
        leg->AddEntry(h1_mumu_back, "Straight Hist", "f");
        leg->Draw();

        c_mumu->Print("mumu17_template_check.png");


        h_elel_alpha->Scale(norm * alpha_denom);

        auto h_elel_pl = *h_elel_sym + *h_elel_asym;
        auto h_elel_mn = *h_elel_sym - *h_elel_asym;
        h_elel_pl.Scale(0.5);
        h_elel_mn.Scale(0.5);
        h_elel_dy = (TH2F *) h_elel_pl.Clone("h_elel_dy");
        h_elel_dy->Add(&h_elel_pl, &h_elel_mn, (norm + afb), (norm - afb));


        auto h1_elel_lo = convert2d(h_elel_dy);
        auto h1_elel_alpha = convert2d(h_elel_alpha);
        TH1F *h1_elel_templ = (TH1F *) h1_elel_lo->Clone("h1_elel_templ");
        h1_elel_templ->Add(h1_elel_alpha);
        
        printf("elel Back is %.2f.Comb Temp is %.2f  LO is %.2f h_alpha is %.2f  \n", 
                h1_elel_back->Integral(),  h1_elel_templ->Integral(), h1_elel_lo->Integral(), h1_elel_alpha->Integral());

        h1_elel_back->SetLineColor(kGreen +3);
        h1_elel_templ->SetLineColor(kBlack);
        h1_elel_lo->SetLineColor(kRed+1);



        TCanvas *c_elel = new TCanvas("c_elel", "Histograms", 200, 10, 900, 700);
        h1_elel_lo->Draw("hist");
        h1_elel_templ->Draw("hist same ");
        h1_elel_back->Draw("hist same");

        leg->Draw();

        c_elel->Print("elel17_template_check.png");


        /*
        */

}




