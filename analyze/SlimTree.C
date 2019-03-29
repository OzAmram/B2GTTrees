



void SlimTree(){
    TFile *f_old = (TFile*) TFile::Open("output_files/MuMu_dy_mar18.root");
    TTree *t_old = (TTree *) f_old ->Get("T_data");
    TTree *t_old2 = (TTree *) f_old ->Get("T_back");


    string fout_name("output_files/MuMu_dy_slim_mar18.root");
    TFile *fout = TFile::Open(fout_name.c_str(), "RECREATE");
    TTree *t_new = t_old->CopyTree("m>140.");
    TTree *t_new2 = t_old2->CopyTree("m>140.");


    fout->cd();
    t_new->Write();
    t_new2->Write();
}

