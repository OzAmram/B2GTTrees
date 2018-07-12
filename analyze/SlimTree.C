



void SlimTree(){
    TFile *f_old = (TFile*) TFile::Open("output_files/MuMu_DY_july10.root");
    TTree *t_old = (TTree *) f_old ->Get("T_data");
    TTree *t_old2 = (TTree *) f_old ->Get("T_back");


    string fout_name("output_files/MuMu_DY_slim_july10.root");
    TFile *fout = TFile::Open(fout_name.c_str(), "RECREATE");
    TTree *t_new = t_old->CopyTree("m>120.");
    TTree *t_new2 = t_old2->CopyTree("m>120.");


    fout->cd();
    t_new->Write();
    t_new2->Write();
}

