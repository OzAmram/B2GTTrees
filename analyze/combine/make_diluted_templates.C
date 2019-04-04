#include "make_templates.C"

void do_dilu_rw(TH2F *h_asym, TH1F* h_rw, TH2F *h_asym_rw){
    for(int i=1; i<= h_asym->GetNbinsX(); i++){
        float rw = h_rw->GetBinContent(i);
        for(int j=1; j<= h_asym->GetNbinsY(); j++){
            float content = h_asym->GetBinContent(i,j);
            float new_content = content*rw;
            h_asym_rw->SetBinContent(i,j, new_content);
        }
    }
    return;
}



void make_diluted_templates(){
    init();
    setup_all_SFs();

    const TString fout_name("combine/templates/mar18_dilu_templates.root");
    char dirname[40];

    int i_start=0;
    int i_max = n_m_bins;

    fout = TFile::Open(fout_name, "RECREATE");
    FILE *f_log;
    char f_log_name[80];


    for(int i=i_start; i<i_max; i++){
        fout->cd();
        snprintf(dirname, 10, "w%i", i);
        gDirectory->mkdir(dirname);
        gDirectory->cd(dirname);
        w = new RooWorkspace("w", "w");

        m_low = m_bins[i];
        m_high = m_bins[i+1];

        h_mumu_sym = new TH2F((string("mumu_sym") ).c_str(), "Symmetric template of mc",
                n_xf_bins, xf_bins, n_cost_bins, cost_bins);
        h_mumu_sym->SetDirectory(0);
        h_mumu_asym = new TH2F((string("mumu_asym") ).c_str(), "Asymmetric template of mc",
                n_xf_bins, xf_bins, n_cost_bins, cost_bins);
        h_mumu_asym->SetDirectory(0);

        gen_mc_template(t_mumu_mc, alpha, h_mumu_sym, h_mumu_asym, m_low, m_high, FLAG_MUONS, FLAG_M_BINS, do_RC);

        h_elel_sym = new TH2F((string("ee_sym") ).c_str(), "Symmetric template of mc",
                n_xf_bins, xf_bins, n_cost_bins, cost_bins);
        h_elel_sym->SetDirectory(0);
        h_elel_asym = new TH2F((string("ee_asym") ).c_str(), "Asymmetric template of mc",
                n_xf_bins, xf_bins, n_cost_bins, cost_bins);
        h_elel_asym->SetDirectory(0);

        gen_mc_template(t_elel_mc, alpha, h_elel_sym, h_elel_asym,  m_low, m_high, FLAG_ELECTRONS, FLAG_M_BINS, do_RC);



        char f_in[100];
        sprintf(f_in, "../setlimits/dilus/M%i_dilus.root", (int) m_low);
        printf("Getting dilus from file %s \n", f_in);
        auto f_dilu = (TFile*) TFile::Open(f_in);
        TH1F *utype_dilu = (TH1F *) f_dilu ->Get("utype_dilu");
        TH1F *dtype_dilu = (TH1F *) f_dilu ->Get("dtype_dilu");
        TH1F *mix_dilu = (TH1F *) f_dilu ->Get("mix_dilu");
        TH1F *h_dilu_ratio = new TH1F("dilu_ratio", "", n_xf_bins, xf_bins);

        for(int l=0; l<=10; l++){
            Float_t x = 0.1*((float) l);
            for(int i=1; i<=n_xf_bins; i++){
                Float_t u_dilu = utype_dilu->GetBinContent(i);
                Float_t d_dilu = dtype_dilu->GetBinContent(i);
                Float_t m_dilu = mix_dilu->GetBinContent(i);
                Float_t ratio = (u_dilu*x + d_dilu*(1-x))/m_dilu;
                h_dilu_ratio->SetBinContent(i,ratio);
            }
            TH2F *h_mumu_asym_rw = (TH2F *) h_mumu_asym->Clone("mumu_asym_rw");
            TH2F *h_elel_asym_rw = (TH2F *) h_elel_asym->Clone("elel_asym_rw");

            do_dilu_rw(h_mumu_asym, h_dilu_ratio, h_mumu_asym_rw);
            do_dilu_rw(h_elel_asym, h_dilu_ratio, h_elel_asym_rw);


            auto h_mumu_pl = *h_mumu_sym + *h_mumu_asym_rw;
            auto h_mumu_mn = *h_mumu_sym - *h_mumu_asym_rw;
            h_mumu_pl.Scale(0.5);
            h_mumu_mn.Scale(0.5);
            h1_mumu_pl = convert2d(&h_mumu_pl);
            h1_mumu_mn = convert2d(&h_mumu_mn);
            char mumu_pl_name_up[100], mumu_mn_name_up[100];
            char mumu_pl_name_down[100], mumu_mn_name_down[100];
            sprintf(mumu_pl_name_up, "mumu_fpl_Dilu%iUp", l);
            sprintf(mumu_mn_name_up, "mumu_fmn_Dilu%iUp", l);
            h1_mumu_pl->SetName(mumu_pl_name_up);
            h1_mumu_mn->SetName(mumu_mn_name_up);

            sprintf(mumu_pl_name_down, "mumu_fpl_Dilu%iDown", l);
            sprintf(mumu_mn_name_down, "mumu_fmn_Dilu%iDown", l);
            TH1F* h1_mumu_pl_dwn = (TH1F*)  h1_mumu_pl->Clone();
            TH1F* h1_mumu_mn_dwn = (TH1F*) h1_mumu_mn->Clone();

            h1_mumu_pl_dwn->SetName(mumu_pl_name_down);
            h1_mumu_mn_dwn->SetName(mumu_mn_name_down);
            write_roo_hist(h1_mumu_pl);
            write_roo_hist(h1_mumu_mn);
            write_roo_hist(h1_mumu_pl_dwn);
            write_roo_hist(h1_mumu_mn_dwn);


            auto h_elel_pl = *h_elel_sym + *h_elel_asym_rw;
            auto h_elel_mn = *h_elel_sym - *h_elel_asym_rw;
            h_elel_pl.Scale(0.5);
            h_elel_mn.Scale(0.5);
            h1_elel_pl = convert2d(&h_elel_pl);
            h1_elel_mn = convert2d(&h_elel_mn);
            char elel_pl_name_up[100], elel_mn_name_up[100];
            char elel_pl_name_down[100], elel_mn_name_down[100];
            sprintf(elel_pl_name_up, "ee_fpl_Dilu%iUp", l);
            sprintf(elel_mn_name_up, "ee_fmn_Dilu%iUp", l);
            h1_elel_pl->SetName(elel_pl_name_up);
            h1_elel_mn->SetName(elel_mn_name_up);

            sprintf(elel_pl_name_down, "ee_fpl_Dilu%iDown", l);
            sprintf(elel_mn_name_down, "ee_fmn_Dilu%iDown", l);
            TH1F* h1_elel_pl_dwn = (TH1F*) h1_elel_pl->Clone();
            TH1F* h1_elel_mn_dwn = (TH1F*) h1_elel_mn->Clone();

            h1_elel_pl_dwn->SetName(elel_pl_name_down);
            h1_elel_mn_dwn->SetName(elel_mn_name_down);
            write_roo_hist(h1_elel_pl);
            write_roo_hist(h1_elel_mn);
            write_roo_hist(h1_elel_mn_dwn);
            write_roo_hist(h1_elel_pl_dwn);

        }
        fout->cd();
        gDirectory->cd(dirname);
        w->Write();
    }

    fout->Close();
    printf("Templates written to %s \n", fout_name.Data());
}
