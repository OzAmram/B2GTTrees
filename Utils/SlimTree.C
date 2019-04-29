



void SlimTree(){
    TFile *f_old = (TFile*) TFile::Open("analyze/output_files/ElEl_comb_back_april18.root");
    TTree *t_old = (TTree *) f_old ->Get("T_data");
    //TTree *t_old2 = (TTree *) f_old ->Get("T_back");


    string fout_name("analyze/output_files/ElEl_comb_back_slim_april18.root");
    TFile *fout = TFile::Open(fout_name.c_str(), "RECREATE");
    TTree *t_new = t_old->CopyTree("m>140.");
    //TTree *t_new2 = t_old2->CopyTree("m>140.");


    fout->cd();
    t_new->Write();
    //t_new2->Write();
}

