#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/InputTag.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "DataFormats/Common/interface/Handle.h"
#include "SimDataFormats/GeneratorProducts/interface/GenEventInfoProduct.h"
#include "SimDataFormats/PileupSummaryInfo/interface/PileupSummaryInfo.h"
#include "TH1.h"
#include "TH2.h"
#include <iostream>

class  EventCounter : public edm::EDAnalyzer {
public:
  EventCounter(const edm::ParameterSet& iConfig)
  {
    edm::Service<TFileService> fs;
    h_pileup_ = fs->make<TH1D>("pileup",";Pile-up", 100,0,100);
    h_weightsign_ = fs->make<TH1D>("weightsign",";Weight sign", 2, -1,1);
    h_totweight_ = fs->make<TH1D>("totweight",";Sum of Weights", 1,0,1);
    h_totweight_T1tttt_ = fs->make<TH2D>("h_totweight_T1tttt","Total N_Event;M_{#tilde{g}} (GeV);M_{#tilde{#chi}^{0}} (GeV);Total Weight", 201,-12.5,5012.5, 201,-12.5,5012.5);
    h_totweight_T2tt_ = fs->make<TH2D>("h_totweight_T2tt","Total N_Event;M_{#tilde{t}} (GeV);M_{#tilde{#chi}^{0}} (GeV);Total Weight", 401,-2.5,2002.5, 401,-2.5,2002.5);

    edm::EDGetTokenT<int>(mayConsume<int>(edm::InputTag("eventUserData", "puNtrueInt")));
    edm::EDGetTokenT<float>(mayConsume<float>(edm::InputTag("extraVar", "evtGenWeight")));
    edm::EDGetTokenT<float>(mayConsume<float>(edm::InputTag("extraVar", "SUSYStopMass")));
    edm::EDGetTokenT<float>(mayConsume<float>(edm::InputTag("extraVar", "SUSYGluinoMass")));
    edm::EDGetTokenT<float>(mayConsume<float>(edm::InputTag("extraVar", "SUSYLSPMass")));
  }

private:
  TH1D* h_pileup_;
  TH1D* h_weightsign_;
  TH1D* h_totweight_;
  TH2D* h_totweight_T1tttt_;
  TH2D* h_totweight_T2tt_;
  void analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup) {
    if ( iEvent.eventAuxiliary().isRealData() ) {
      h_weightsign_->Fill(0);
      h_totweight_->Fill(0);
    } else {
      // Save Pileup distribution for MC
      edm::Handle<int>  NTrueInt;
      iEvent.getByLabel(edm::InputTag("eventUserData", "puNtrueInt"), NTrueInt);
      h_pileup_->Fill(*NTrueInt);
      
      // Save sign of weights and total weight
      edm::Handle<float> evt_Gen_Weight;
      iEvent.getByLabel(edm::InputTag("extraVar", "evtGenWeight"), evt_Gen_Weight);
      h_weightsign_->Fill((*evt_Gen_Weight>=0)-1);
      h_totweight_->Fill(double(0), *evt_Gen_Weight);

      printf("I am running! \n\n\n");
      
      // Count total weight for SUSY signal events
      edm::Handle<float> SUSY_Stop_Mass;
      edm::Handle<float> SUSY_Gluino_Mass;
      edm::Handle<float> SUSY_LSP_Mass;
      iEvent.getByLabel(edm::InputTag("extraVar", "SUSYStopMass"),   SUSY_Stop_Mass);
      iEvent.getByLabel(edm::InputTag("extraVar", "SUSYGluinoMass"), SUSY_Gluino_Mass);
      iEvent.getByLabel(edm::InputTag("extraVar", "SUSYLSPMass"),    SUSY_LSP_Mass);
      if (*SUSY_Gluino_Mass>0) h_totweight_T1tttt_->Fill(*SUSY_Gluino_Mass, *SUSY_LSP_Mass, *evt_Gen_Weight);
      if (*SUSY_Stop_Mass>0)   h_totweight_T2tt_  ->Fill(*SUSY_Stop_Mass,   *SUSY_LSP_Mass, *evt_Gen_Weight);
    }
  }
};

#include "FWCore/Framework/interface/MakerMacros.h"
DEFINE_FWK_MODULE(EventCounter);
