


void merge_workspaces(){


    const TString f1_s("combine/templates/april9_combined_no_sys.root");
    const TString fout_s("combine/templates/april23_merge.root");
    TFile *f1 = TFile::Open(f1_s, "READ");
    TFile *fout = TFile::Open(fout_s, "RECREATE");
    char dirname[40];

    std::vector<string> fs{"combine/templates/april22_no_sys.root",
        "root://131.225.204.161:1094//store/user/oamram/Condor_outputs/templ_sys_april22/file_0.root",
        "root://131.225.204.161:1094//store/user/oamram/Condor_outputs/templ_sys_april22/file_1.root",
        "root://131.225.204.161:1094//store/user/oamram/Condor_outputs/templ_sys_april22/file_2.root",
        "root://131.225.204.161:1094//store/user/oamram/Condor_outputs/templ_sys_april22/file_3.root",
        "root://131.225.204.161:1094//store/user/oamram/Condor_outputs/templ_sys_april22/file_4.root",
        "root://131.225.204.161:1094//store/user/oamram/Condor_outputs/templ_pdf_sys_april22/file_0.root",
        "root://131.225.204.161:1094//store/user/oamram/Condor_outputs/templ_pdf_sys_april22/file_1.root",
        "root://131.225.204.161:1094//store/user/oamram/Condor_outputs/templ_pdf_sys_april22/file_2.root",
        "root://131.225.204.161:1094//store/user/oamram/Condor_outputs/templ_pdf_sys_april22/file_3.root",
        "root://131.225.204.161:1094//store/user/oamram/Condor_outputs/templ_pdf_sys_april22/file_4.root",
        "root://131.225.204.161:1094//store/user/oamram/Condor_outputs/templ_pdf_sys_april22/file_5.root",
        "root://131.225.204.161:1094//store/user/oamram/Condor_outputs/templ_pdf_sys_april22/file_6.root",
        "root://131.225.204.161:1094//store/user/oamram/Condor_outputs/templ_pdf_sys_april22/file_7.root",
        "root://131.225.204.161:1094//store/user/oamram/Condor_outputs/templ_pdf_sys_april22/file_8.root",
        "root://131.225.204.161:1094//store/user/oamram/Condor_outputs/templ_pdf_sys_april22/file_9.root"

    };

    int i_start=0;
    int i_max = 6;



    for(int i=i_start; i<i_max; i++){
        snprintf(dirname, 10, "w%i", i);
        f1->cd(dirname);
        RooWorkspace *w1 = (RooWorkspace *) gDirectory->Get("w");
        for(auto f_name = fs.begin(); f_name !=fs.end(); f_name++){
            TFile *f2 = TFile::Open(f_name->c_str(), "READ");
            f2->cd();
            f2->cd(dirname);
            RooWorkspace *w2 = (RooWorkspace *) gDirectory->Get("w");
            list<RooAbsData *> ds = w2->allData();
            //doesn't pick up vars, only datasets
            for(auto iter = ds.begin(); iter!= ds.end(); iter++){
                (*iter)->Print();
                w1->import(*(*iter));
            }
            f2->Close();
        }
        //w1->Print();
        //w2->Print();

        fout->mkdir(dirname);
        fout->cd(dirname);
        w1->Write();
    }
    fout->Print();
    fout->Write();
    fout->Close();
}
