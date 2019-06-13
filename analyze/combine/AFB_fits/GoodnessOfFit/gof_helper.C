



void gof_helper(int label_num=0, int idx=0){

    string label;
    if(label_num == 0) label = string("combined");
    if(label_num == 1) label = string("ee");
    if(label_num == 2) label = string("mumu");

    TFile *_file2 = TFile::Open("higgsCombineTest.GoodnessOfFit.mH120.root");

    Double_t t_obs;
    TTree* data_l = (TTree *) _file2->Get("limit");
    data_l->SetBranchAddress("limit", &t_obs);
    data_l->GetEntry(0);
    _file2->Close();

    TFile *_file1 = TFile::Open("higgsCombineTest.GoodnessOfFit.mH120.123.root");

    TTree* toys = (TTree *) _file1->Get("limit");
    Double_t toy_max =toys->GetMaximum("limit");
    Double_t toy_min =toys->GetMinimum("limit");
    Double_t max = std::max(toy_max, t_obs) + 10.;
    //if is an error in fit can get very large values in toys
    max = std::min(max, 2.*t_obs);
    Double_t min = std::min(toy_min, t_obs) - 5.;
    

    char h_title[80];
    sprintf(h_title, "Goodness of Fit: Mass Bin %i; Test Statitic ", idx);
    TH1D *h_test = new TH1D("h_toys", h_title, 30, min, max);
    toys->Draw("limit>>h_toys");
    int bin_low = h_test->GetXaxis()->FindBin(t_obs);
    int bin_high = h_test->GetXaxis()->FindBin(max);
    float integral = h_test->Integral();
    printf("Data gof is %.0f. p-value is %.3f based on %.0f toys \n", t_obs, h_test->Integral(bin_low, bin_high)/integral, integral);



    char fout_name[100];
    sprintf(fout_name, "GoodnessOfFit/%s_bin%i.png", label.c_str(), idx);
    
    TCanvas *c = new TCanvas("c","",  800,800);
    h_test->Draw("hist");
    Double_t draw_max = h_test->GetMaximum();
    TLine *l = new TLine(t_obs, 0., t_obs, draw_max);
    l->SetLineColor(kRed);
    l->SetLineWidth(2);
    l->Draw("same");
    c->Print(fout_name);


    _file1->Close();
}
