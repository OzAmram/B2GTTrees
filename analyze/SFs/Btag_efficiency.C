#include "../../Utils/NTupleReader.C"



void Btag_efficiency(int nJobs =1, int iJob = 0, string fin ="")
{

    if(fin == "") fin = string("EOS_files/2016/test.txt");
    NTupleReader nt(fin.c_str(),"output_files/test.root", false);
    nt.year = 2016;

    string fout_name = string("MC_btag_eff.root");
    TFile *fout = TFile::Open(fout_name.c_str(), "RECREATE");

    Double_t Eta_bins[] = {0, 0.9, 1.2, 2.1, 2.4}; 
    Int_t nEta_bins = 4;
    Double_t Pt_bins[] = {0,10,20,30,40,50,60,80,100,120,160,200,250,300,400,500};
    Int_t nPt_bins = 15;
    TH2D *b_num = new TH2D("b_num", "Efficiency for b-tagging b-jets", nPt_bins, Pt_bins, nEta_bins, Eta_bins);
    TH2D *b_denom = new TH2D("b_demon", "Total number of b-jets", nPt_bins, Pt_bins, nEta_bins, Eta_bins);
    TH2D *c_num = new TH2D("c_num", "Efficiency for b-tagging c-jets", nPt_bins, Pt_bins, nEta_bins, Eta_bins);
    TH2D *c_denom = new TH2D("c_denom", "Efficiency for b-tagging c-jets", nPt_bins, Pt_bins, nEta_bins, Eta_bins);
    TH2D *udsg_num = new TH2D("udsg_num", "Efficiency for b-tagging udsg-jets", nPt_bins, Pt_bins, nEta_bins, Eta_bins);
    TH2D *udsg_denom = new TH2D("udsg_denom", "Efficiency for b-tagging udsg-jets", nPt_bins, Pt_bins, nEta_bins, Eta_bins);

    unsigned int nJets=0;
    unsigned int nB=0;
    unsigned int nC=0;
    unsigned int nUDSG=0;



    char out_buff[10000];
    bool print_out = false;

    int B =5;
    int C =4;

    while(nt.getNextFile()){
        for (int i=0; i<nt.tin_nEntries; i++) {
            nt.getEvent(i);
            if(nt.jet_size >= 2){
                for(int j=0; j < nt.jet_size; j++){
                    int flavour = std::abs(nt.jet_hadronflavour[j]);
                    if (flavour == B){

                        nB++;
                        if(nt.jet_btag[j] >= nt.bjet_med_cut) b_num->Fill(nt.jet_Pt[j], std::abs(nt.jet_Eta[j]));
                        b_denom->Fill(nt.jet_Pt[j], std::abs(nt.jet_Eta[j]));
                    }
                    else if (flavour == C){
                        nC++;
                        if(nt.jet_btag[j] >= nt.bjet_med_cut) b_num->Fill(nt.jet_Pt[j], std::abs(nt.jet_Eta[j]));
                        b_denom->Fill(nt.jet_Pt[j], std::abs(nt.jet_Eta[j]));
                    }
                    else if (nt.jet_genPt[j] > 8.){
                        nUDSG++;
                        if(nt.jet_btag[j] >= nt.bjet_med_cut) b_num->Fill(nt.jet_Pt[j], std::abs(nt.jet_Eta[j]));
                        b_denom->Fill(nt.jet_Pt[j], std::abs(nt.jet_Eta[j]));
                    }
                    nJets++;
                }

            }
        }
        printf("moving on to next file, currently %i Jets \n\n", nJets);

    }

    nt.finish();

    printf("There were %i Bs %i Cs and %i UDSGs in %i files.\n",
            nB, nC, nUDSG, nt.fileCount);
    fout->cd();

    TH2D* b_eff = (TH2D *) b_num->Clone("b_eff");
    b_eff->Divide(b_denom);
    b_eff->Write();

    TH2D* c_eff = (TH2D *) c_num->Clone("c_eff");
    c_eff->Divide(c_denom);
    c_eff->Write();

    TH2D* udsg_eff = (TH2D *) udsg_num->Clone("udsg_eff");
    udsg_eff->Divide(udsg_denom);
    udsg_eff->Write();

    printf("Writing output to file at %s \n", fout_name.c_str());
    fout->Close();

    return;
}
