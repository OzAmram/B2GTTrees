


void merge_histos(){
    TFile *f_low = TFile::Open("egammaEffi.txt_EGM2D_runBCDEF_passingRECO_lowEt.root");
    TH2F *h_low = (TH2F *) gDirectory->Get("EGamma_SF2D");



    TFile *f_high = TFile::Open("egammaEffi.txt_EGM2D_runBCDEF_passingRECO.root");
    TH2F *h_high = (TH2F *) gDirectory->Get("EGamma_SF2D");


    Double_t x_bins[] =  {-2.5, -2.0, -1.566, -1.444, -1.0, -0.5, 0., 0.5, 1.0, 1.444, 1.566, 2.0, 2.5};
    int n_x_bins = 12;
    Double_t y_bins[] =  {10., 20., 45., 75., 100., 500.};
    int n_y_bins = 5;

    TH2F *h_new = new TH2F("EGamma_SF2D", "", n_x_bins, x_bins, n_y_bins, y_bins);

    for(int i=1; i<=n_x_bins; i++){
        float cont = h_low->GetBinContent(i,1);
        float err = h_low->GetBinError(i,1);
        h_new->SetBinContent(i,1, cont);
        h_new->SetBinError(i,1, err);
        for(int j=1; j<=n_y_bins-1; j++){
            cont = h_high->GetBinContent(i,j);
            err = h_high->GetBinError(i,j);
            h_new->SetBinContent(i,j+1, cont);
            h_new->SetBinError(i,j+1, err);
        }
    }

    h_low->Print("all");
    h_high->Print("all");
    h_new->Print("all");


    TFile *f_out = new TFile("El_RECO.root", "NEW");
    h_new->Write();

}





