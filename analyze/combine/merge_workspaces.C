


void merge_workspaces(){

    const TString f1_s("combine/templates/mar18_combined_template_sys.root");
    const TString f2_s("combine/templates/mar18_dilu_templates.root");
    const TString fout_s("combine/templates/mar18_merge.root");
    TFile *f1 = TFile::Open(f1_s, "READ");
    TFile *f2 = TFile::Open(f2_s, "READ");
    TFile *fout = TFile::Open(fout_s, "RECREATE");
    char dirname[40];

    int i_start=0;
    int i_max = 6;



    for(int i=i_start; i<i_max; i++){
        snprintf(dirname, 10, "w%i", i);
        f1->cd(dirname);
        RooWorkspace *w1 = (RooWorkspace *) gDirectory->Get("w");
        f2->cd();
        f2->cd(dirname);
        RooWorkspace *w2 = (RooWorkspace *) gDirectory->Get("w");
        list<RooAbsData *> ds = w2->allData();
        for(auto iter = ds.begin(); iter!= ds.end(); iter++){
            (*iter)->Print();
            w1->import(*(*iter));
        }
        //w1->Print();
        //w2->Print();

        fout->mkdir(dirname);
        fout->cd(dirname);
        w1->Write();
    }
    fout->Print();
    fout->Write();
}

