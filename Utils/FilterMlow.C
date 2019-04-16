
void FilterMlow(){
    TFile *f_old = (TFile*) TFile::Open("output_files/MuMu_WJets_mc_mlow_feb25.root");
    TTree *t_old = (TTree *) f_old ->Get("T_data");


    TFile *fout = (TFile*) TFile::Open("output_files/MuMu_WJets_dy_mc_filt_mlow_feb25.root", "RECREATE");
    TTree *t_new = t_old->CopyTree("m<100.");


    fout->cd();
    t_new->Write();
}

