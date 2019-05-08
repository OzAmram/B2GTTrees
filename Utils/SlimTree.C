



void SlimTree(){
    TFile *f_old = (TFile*) TFile::Open("analyze/output_files/MuMu_combined_back_may7.root");
    TTree *t_old = (TTree *) f_old ->Get("T_data");
    //TTree *t_old2 = (TTree *) f_old ->Get("T_back");


    string fout_name("analyze/output_files/MuMu_comb_back_slim_may7.root");
    TFile *fout = TFile::Open(fout_name.c_str(), "RECREATE");
    TTree *t_new = t_old->CopyTree("m>140.");
    //TTree *t_new2 = t_old2->CopyTree("m>140.");


    fout->cd();
    t_new->Write();
    //t_new2->Write();
}

