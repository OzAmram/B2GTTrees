#include "make_templates.C"
TTree *t_elel_mc_raw, *t_elel_nosig_raw, *t_elel_back_raw;
TTree *t_mumu_mc_raw, *t_mumu_nosig_raw, *t_mumu_back_raw;

void init_condor_files(int year){
    f_elel_mc = (TFile*) TFile::Open("root://131.225.204.161:1094//store/user/oamram/Condor_inputs/ElEl_dy_july1.root");
    f_elel_back = (TFile*) TFile::Open("root://131.225.204.161:1094//store/user/oamram/Condor_inputs/ElEl_comb_back_july1.root");
    f_mumu_mc = (TFile*) TFile::Open("root://131.225.204.161:1094//store/user/oamram/Condor_inputs/MuMu_dy_july15.root");
    f_mumu_back = (TFile*) TFile::Open("root://131.225.204.161:1094//store/user/oamram/Condor_inputs/MuMu_comb_back_july15.root");
    f_mumu_gamgam = (TFile*) TFile::Open("root://131.225.204.161:1094//store/user/oamram/Condor_inputs/MuMu_gamgam_july15.root");
    f_elel_gamgam = (TFile*) TFile::Open("root://131.225.204.161:1094//store/user/oamram/Condor_inputs/ElEl_gamgam_back_june25.root");

    
    t_mumu_mc_raw = (TTree *) f_mumu_mc ->Get("T_data");
    t_mumu_nosig_raw = (TTree *) f_mumu_mc ->Get("T_back");
    t_mumu_back_raw = (TTree *) f_mumu_back ->Get("T_data");
    t_elel_mc_raw = (TTree *) f_elel_mc ->Get("T_data");
    t_elel_nosig_raw = (TTree *) f_elel_mc ->Get("T_back");
    t_elel_back_raw = (TTree *) f_elel_back ->Get("T_data");
    t_mumu_gamgam = (TTree *)f_mumu_gamgam->Get("T_data");
    t_elel_gamgam = (TTree *)f_elel_gamgam->Get("T_data");
}



void make_sys_templates(int nJobs = 1, int iJob =0, int type=0){

    //const TString pdf_fout_name("combine/templates/april18_pdf_test.root");
    const TString pdf_fout_name("output_files/july12_sys.root");
    TFile *pdf_fout = TFile::Open(pdf_fout_name, "RECREATE");

    //init();
    int year = 2016;
    init_condor_files(year);





    printf("Setting up SFs... ");
    setup_all_SFs(year);
    printf("   done \n");


    
    vector<string> sys_labels;
    if (type ==0){
        for(int i =1; i<= 60; i++){
            char name1[20], name2[20];
            sprintf(name1, "_pdf%iUp", i);
            sprintf(name2, "_pdf%iDown", i);
            sys_labels.push_back(string(name1));
            sys_labels.push_back(string(name2));

        }
    }
    else{
      sys_labels =  {"_elScaleStatUp", "_elScaleStatDown", 
          "_elScaleSystUp", "_elScaleSystDown",
          "_elScaleGainUp", "_elScaleGainDown",
          "_elSmearUp", "_elSmearDown", 
       "_muRCUp", "_muRCDown", "_muHLTUp", "_muHLTDown", "_muIDUp", "_muIDDown", "_muISOUp", "_muISODown",  
       "_elHLTUp", "_elHLTDown", "_elIDUp", "_elIDDown", "_elRECOUp", "_elRECODown", 
       "_RENORMUp", "_RENORMDown", "_FACUp", "_FACDown",
       "_PuUp", "_PuDown", "_BTAGUp", "_BTAGDown", 
       "_alphaDenUp", "_alphaDenDown", "_alphaSUp", "_alphaSDown" };


    }
    //sys_labels = {"_alphaDenUp", "_alphaDenDown"};





    char dirname[40];

    int i_start=0;
    int i_max = n_m_bins;


    for(int i=i_start; i<i_max; i++){

        pdf_fout->cd();
        snprintf(dirname, 10, "w%i", i);
        gDirectory->mkdir(dirname);
        gDirectory->cd(dirname);
        w = new RooWorkspace("w", "w");

        m_low = m_bins[i];
        m_high = m_bins[i+1];
        printf("\n \n Start making templates for mass bin %.0f-%.0f \n", m_low, m_high);

        char cut_str[100];
        printf("Cutting templates \n");
        sprintf(cut_str, "(m > %.0f) && (m < %0.f)", m_low - 15., m_high+15.);
        t_elel_mc = t_elel_mc_raw->CopyTree(cut_str);
        t_elel_nosig = t_elel_nosig_raw->CopyTree(cut_str);
        t_elel_back = t_elel_back_raw->CopyTree(cut_str);

        t_mumu_mc = t_mumu_mc_raw->CopyTree(cut_str);
        t_mumu_nosig = t_mumu_nosig_raw->CopyTree(cut_str);
        t_mumu_back = t_mumu_back_raw->CopyTree(cut_str);

        int i_sys =0;
        for(auto iter = sys_labels.begin(); iter !=sys_labels.end(); iter++){
            if(i_sys % nJobs == iJob){
                printf("Making MC templates for sys %s \n", (*iter).c_str());

                Double_t alpha_num = alphas_num[i];
                Double_t alpha_denom = alphas_denom[i];


                if(iter->find("alphaDenUp") != string::npos) alpha_denom = alphas_denom[i] + alpha_denom_unc[i];
                if(iter->find("alphaDenDown") != string::npos) alpha_denom = alphas_denom[i] - alpha_denom_unc[i];

                printf("alpha_denom %.2f \n", alpha_denom);

                make_mc_templates(year, alpha_denom, *iter);
                convert_mc_templates(year, *iter);
            }
            i_sys++;
        }
        pdf_fout->cd();
        gDirectory->cd(dirname);
        w->Write();
        cleanup_mc_templates();
    }


    pdf_fout->Close();
    printf("Templates written to %s \n", pdf_fout_name.Data());

}

