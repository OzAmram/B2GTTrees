



void SlimTree(){
    TFile *f_old = (TFile*) TFile::Open("output_files/SingleMuon_data_july10.root");
    TTree *t_old = (TTree *) f_old ->Get("T_data");
    //TTree *t_old2 = (TTree *) f_old ->Get("T_back");


    string fout_name("output_files/SingleMuon_data_slim_july10.root");
    TFile *fout = TFile::Open(fout_name.c_str(), "RECREATE");
    TTree *t_new = t_old->CopyTree("m>130.");
    //TTree *t_new2 = t_old2->CopyTree("m>130.");


    fout->cd();
    t_new->Write();
    //t_new2->Write();
}

