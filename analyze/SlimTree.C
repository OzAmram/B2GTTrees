



void SlimTree(){
    TFile *f_old = (TFile*) TFile::Open("output_files/MuMu_DY_unbinned_sep11.root");
    TTree *t_old = (TTree *) f_old ->Get("T_data");
    TTree *t_old2 = (TTree *) f_old ->Get("T_back");


    string fout_name("output_files/MuMu_DY_unbinned_slim_sep11.root");
    TFile *fout = TFile::Open(fout_name.c_str(), "RECREATE");
    TTree *t_new = t_old->CopyTree("m>130.");
    TTree *t_new2 = t_old2->CopyTree("m>130.");


    fout->cd();
    t_new->Write();
    t_new2->Write();
}

