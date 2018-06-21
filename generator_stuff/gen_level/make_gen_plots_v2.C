#include "GenLoader.cc"


void fill_gen_hist(TTree *t, TH1F *h_m, TH1F *h_pt){
    Long64_t size = t->GetEntries();
    GenLoader *g = new GenLoader(t);
    bool got_lep1, got_lep2;
    TLorentzVector lep1, lep2;
    int nEvents=0;
    int nMiss = 0;
    int nExtra = 0;
    int ELECTRON = 11; 
    int MUON = 13;
    int PHOTON = 22;
    int Z=23;
    int GLUON = 21;
    int TAU = 15;
    int PROTON = 2212;
    for(int i =0; i<size; i++){
        got_lep1 = got_lep2 = false;
        g->load(i);
        //printf("event %i \n", i);
        for (int j=0; j<g->fGens->GetEntriesFast(); ++j) {
            TGenParticle *p1 = (TGenParticle*)((*(g->fGens))[j]);
            if(p1->parent < 0) continue;
            TGenParticle *par = (TGenParticle*)((*(g->fGens))[p1->parent]);
            if((abs(p1->pdgId) == 11 || abs(p1->pdgId) == 13  || abs(p1->pdgId) == 15) &&
                    ((par->pdgId == PHOTON || par->pdgId==Z) || (p1->status == 23 && par->pdgId != PROTON))) {
                //printf("ID = %i, parID = %i, parStatus = %i , status = %i \n", p1->pdgId, par->pdgId, par->status, p1->status);
                if(!got_lep1){
                    got_lep1 = true;
                    lep1.SetPtEtaPhiM(p1->pt, p1->eta, p1->phi, p1->mass);
                }
                else if(!got_lep2){
                    got_lep2 = true;
                    lep2.SetPtEtaPhiM(p1->pt, p1->eta, p1->phi, p1->mass);
                }
                else{
                    //printf("Extra lepton! \n");
                    nExtra++;
                }
            }
        }
        if(got_lep1 && got_lep2){
            TLorentzVector cm = lep1 + lep2;
            if(cm.M() > 100. && cm.M() < 200.){
                nEvents++;
                h_m->Fill(cm.M(), g->fGenInfo->weight);
                h_pt->Fill(cm.Pt(), g->fGenInfo->weight);
            }
        }
        else{
            nMiss++;
            //printf("Didn't get leps \n");
        }
    }
    h_m->Scale(1./h_m->Integral());
    h_pt->Scale(1./h_pt->Integral());
    printf("nEvents = %i nMiss  = %i nExtra = %i \n", nEvents, nMiss, nExtra);
    return;
}




void make_ratio_plot(char title[80], TH1F* h1, char h1_label[80], TH1F* h2, char h2_label[80], char ratio_label[80], 
        char axis_label[80], bool logy=false){

    TCanvas *c = new TCanvas("c_m", "Histograms", 200, 10, 900, 700);
    TPad *pad1 = new TPad("pad1", "pad1", 0.,0.3,0.98,1.);
    pad1->SetBottomMargin(0);
    pad1->Draw();
    pad1->cd();
    if(logy) pad1->SetLogy();
    h1->Draw("hist E");
    //m_stack->SetMaximum(65000);
    //h1->SetMinimum(0.1);
    gStyle->SetEndErrorSize(4);
    h2->Draw("hist E same");


    gStyle->SetLegendBorderSize(0);
    TLegend *leg1 = new TLegend(0.15, 0.05, 0.4, 0.3);
    leg1->AddEntry(h1, h1_label, "l");
    leg1->AddEntry(h2, h2_label, "l");
    leg1->Draw();

    //gPad->BuildLegend();
    c->cd();
    TPad *pad2 = new TPad("pad2", "pad2", 0.,0,.98,0.3);
    //pad2->SetTopMargin(0);
    pad2->SetBottomMargin(0.2);
    pad2->SetGridy();
    pad2->Draw();
    pad2->cd();

    auto ratio = (TH1F *) h1->Clone("h_ratio");
    ratio->SetMinimum(0.01);
    ratio->SetMaximum(2.0);
    ratio->Sumw2();
    ratio->SetStats(0);
    ratio->Divide(h2);
    ratio->SetMarkerStyle(21);
    ratio->SetLineColor(kBlack);
    ratio->Draw("ep");
    c->cd();

    ratio->SetTitle("");
    // Y axis ratio plot settings
    ratio->GetYaxis()->SetTitle(ratio_label);
    ratio->GetYaxis()->SetNdivisions(505);
    ratio->GetYaxis()->SetTitleSize(20);
    ratio->GetYaxis()->SetTitleFont(43);
    ratio->GetYaxis()->SetTitleOffset(1.2);
    ratio->GetYaxis()->SetLabelFont(43); // Absolute font size in pixel (precision 3)
    ratio->GetYaxis()->SetLabelSize(15);
    // X axis ratio plot settings
    ratio->GetXaxis()->SetTitle(axis_label);
    ratio->GetXaxis()->SetTitleSize(20);
    ratio->GetXaxis()->SetTitleFont(43);
    ratio->GetXaxis()->SetTitleOffset(3.);
    ratio->GetXaxis()->SetLabelFont(43); // Absolute font size in pixel (precision 3)
    ratio->GetXaxis()->SetLabelSize(20);

    //lumi_sqrtS = "";       // used with iPeriod = 0, e.g. for simulation-only plots (default is an empty string)
    c->Print(title);
    return;
}


void make_gen_plots_v2(){
    gStyle->SetOptStat(0);
    TFile *f_unbinned = TFile::Open("DY_M50_gen_v2.root");
    TTree *t_unbinned = (TTree *)f_unbinned->Get("Events");

    TFile *f_binned = TFile::Open("DY_M100_gen.root");
    TTree *t_binned = (TTree *)f_binned->Get("Events");

    TH1F *h_unbin_m = new TH1F("h_unbin_m", "", 20, 100, 200);
    TH1F *h_unbin_pt = new TH1F("h_unbin_pt", "", 20, 0, 400);

    TH1F *h_bin_m = new TH1F("h_bin_m", "", 20, 100, 200);
    TH1F *h_bin_pt = new TH1F("h_bin_pt", "", 20, 0, 400);
    fill_gen_hist(t_unbinned, h_unbin_m, h_unbin_pt);
    fill_gen_hist(t_binned, h_bin_m, h_bin_pt);

    h_bin_pt->SetLineColor(kRed);
    h_unbin_pt->SetLineColor(kBlue);
    h_bin_pt->SetLineWidth(3);
    h_unbin_pt->SetLineWidth(3);

    h_bin_m->SetLineColor(kRed);
    h_unbin_m->SetLineColor(kBlue);
    h_bin_m->SetLineWidth(3);
    h_unbin_m->SetLineWidth(3);



    make_ratio_plot("gen_m_cmp.pdf", h_bin_m, "Binned MC (DY_M-100to200)",h_unbin_m, "Unbinned MC (DY_M50)", "Binned/Unbinned", "M (GeV)", true);
    make_ratio_plot("gen_pt_cmp.pdf", h_bin_pt, "Binned MC (DY_M-100to200)",h_unbin_pt, "Unbinned MC (DY_M50)", "Binned/Unbinned", "Pt (GeV)", true);
    return;
}
