



void SlimTree(){
    TFile *f_old = (TFile*) TFile::Open("analyze/output_files/ElEl_dy_june5.root");
    string fout_name("analyze/output_files/ElEl_dy_slim_june5.root");
    bool do_dy = true;

    TTree *t_old = (TTree *) f_old ->Get("T_data");
    TFile *fout = TFile::Open(fout_name.c_str(), "RECREATE");
    TTree *t_new = t_old->CopyTree("m>130.");

    fout->cd();
    t_new->Write();
    if(do_dy){

        TTree *t_old2 = (TTree *) f_old ->Get("T_back");
        TTree *t_new2 = t_old2->CopyTree("m>130.");
        t_new2->Write();
    }
}

