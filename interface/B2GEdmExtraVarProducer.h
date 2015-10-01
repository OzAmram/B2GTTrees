#include "FWCore/Framework/interface/EDProducer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/InputTag.h"
#include "DataFormats/Common/interface/Handle.h"
#include "DataFormats/Math/interface/LorentzVector.h"

// JEC
#include "CondFormats/JetMETObjects/interface/JetCorrectorParameters.h"
#include "CondFormats/JetMETObjects/interface/FactorizedJetCorrector.h"

class  B2GEdmExtraVarProducer : public edm::EDProducer {
public:
  B2GEdmExtraVarProducer( const edm::ParameterSet & );
  
private:
  void produce(edm::Event&, const edm::EventSetup& );
  void calculate_variables(const edm::Event&, const edm::EventSetup&);
  
  // Input
  bool isData_;
  std::string JEC_location_;
  
  double event_weight_;
  
  std::string filter_label_;
  std::string trigger_label_;
  std::string evt_label_;
  std::string evt_prefix_;
  std::string met_label_;
  std::string met_prefix_;
  std::string gen_label_;
  std::string gen_prefix_;
  std::string electrons_label_;
  std::string electrons_prefix_;
  std::string muons_label_;
  std::string muons_prefix_;
  std::string AK4Jets_label_;
  std::string AK4Jets_prefix_;
  std::string AK8Jets_label_;
  std::string AK8Jets_prefix_;
  std::string AK8Subjets_label_;
  std::string AK8Subjets_prefix_;
  std::string CmsTTSubjets_label_;
  std::string CmsTTSubjets_prefix_;
  std::string AK4JetKeys_label_;
  std::string AK8JetKeys_label_;
  std::string AK8SubjetKeys_label_;
  std::string CmsTTSubjetKeys_label_;
    
  // Handles
  std::map<std::string, edm::Handle<int> > h_bool_;
  std::map<std::string, edm::Handle<int> > h_int_;
  std::map<std::string, edm::Handle<double> > h_double_;
  std::map<std::string, edm::Handle<std::vector<int> > > h_ints_;
  std::map<std::string, edm::Handle<std::vector<float> > > h_floats_;
  std::map<std::string, edm::Handle<std::vector<std::vector<int> > > > h_keys_;
  std::map<std::string, edm::Handle<std::vector<std::string> > > h_strings_;
  
  // Output
  std::vector<std::string> singleB_, singleI_, singleF_, vectorI_, vectorF_;
  std::map<std::string, int> single_bool_;
  std::map<std::string, int> single_int_;
  std::map<std::string, float> single_float_;
  std::map<std::string, std::vector<int> > vector_int_;
  std::map<std::string, std::vector<float> > vector_float_;
  std::vector<std::string> filter_names_, trigger_names_;
  std::map<std::string, size_t > filters_, triggers_;
  
  size_t nfilt_, ntrig_;
  
  // JEC
  FactorizedJetCorrector *AK4_JetCorrector_, *AK8_JetCorrector_;
};


B2GEdmExtraVarProducer::B2GEdmExtraVarProducer(const edm::ParameterSet& iConfig) :
  isData_(iConfig.getUntrackedParameter<bool>("isData", false)),
  JEC_location_(iConfig.getUntrackedParameter<std::string>("JEC_location")),
  event_weight_(iConfig.getUntrackedParameter<double>("event_weight",1)),
  filter_label_(iConfig.getUntrackedParameter<std::string>("filter_label")),
  trigger_label_(iConfig.getUntrackedParameter<std::string>("trigger_label")),
  evt_label_(iConfig.getUntrackedParameter<std::string>("evt_label")),
  evt_prefix_(iConfig.getUntrackedParameter<std::string>("evt_prefix")),
  met_label_(iConfig.getUntrackedParameter<std::string>("met_label")),
  met_prefix_(iConfig.getUntrackedParameter<std::string>("met_prefix")),
  gen_label_(iConfig.getUntrackedParameter<std::string>("gen_label")),
  gen_prefix_(iConfig.getUntrackedParameter<std::string>("gen_prefix")),
  electrons_label_(iConfig.getUntrackedParameter<std::string>("electrons_label")),
  electrons_prefix_(iConfig.getUntrackedParameter<std::string>("electrons_prefix")),
  muons_label_(iConfig.getUntrackedParameter<std::string>("muons_label")),
  muons_prefix_(iConfig.getUntrackedParameter<std::string>("muons_prefix")),
  AK4Jets_label_(iConfig.getUntrackedParameter<std::string>("AK4Jets_label")),
  AK4Jets_prefix_(iConfig.getUntrackedParameter<std::string>("AK4Jets_prefix")),
  AK8Jets_label_(iConfig.getUntrackedParameter<std::string>("AK8Jets_label")),
  AK8Jets_prefix_(iConfig.getUntrackedParameter<std::string>("AK8Jets_prefix")),
  AK8Subjets_label_(iConfig.getUntrackedParameter<std::string>("AK8Subjets_label")),
  AK8Subjets_prefix_(iConfig.getUntrackedParameter<std::string>("AK8Subjets_prefix")),
  CmsTTSubjets_label_(iConfig.getUntrackedParameter<std::string>("CmsTTSubjets_label")),
  CmsTTSubjets_prefix_(iConfig.getUntrackedParameter<std::string>("CmsTTSubjets_prefix")),
  AK4JetKeys_label_(iConfig.getUntrackedParameter<std::string>("AK4JetKeys_label")),
  AK8JetKeys_label_(iConfig.getUntrackedParameter<std::string>("AK8JetKeys_label")),
  AK8SubjetKeys_label_(iConfig.getUntrackedParameter<std::string>("AK8SubjetKeys_label")),
  CmsTTSubjetKeys_label_(iConfig.getUntrackedParameter<std::string>("CmsTTSubjetKeys_label")),
  singleB_(iConfig.getUntrackedParameter<std::vector<std::string> >("singleB")),
  singleI_(iConfig.getUntrackedParameter<std::vector<std::string> >("singleI")),
  singleF_(iConfig.getUntrackedParameter<std::vector<std::string> >("singleF")),
  vectorI_(iConfig.getUntrackedParameter<std::vector<std::string> >("vectorI")),
  vectorF_(iConfig.getUntrackedParameter<std::vector<std::string> >("vectorF"))
{
  for ( auto nameB : singleB_ ) {
    if (nameB.find("Flag_")==0) filter_names_.push_back(nameB);
    if (nameB.find("HLT_")==0) trigger_names_.push_back(nameB);
    size_t f; while ((f=nameB.find("_"))!=std::string::npos) nameB.erase(f,1); // Remove "_" from var name
    produces<bool>(nameB);
  }
  for ( auto nameI : trigger_names_ ) {
    size_t f; while ((f=nameI.find("_"))!=std::string::npos) nameI.erase(f,1); // Remove "_" from var name
    produces<int>(nameI+"prescale");
  }
  for ( auto nameI : singleI_ ) {
    size_t f; while ((f=nameI.find("_"))!=std::string::npos) nameI.erase(f,1); // Remove "_" from var name
    produces<int>(nameI);
  }
  for ( auto nameF : singleF_ ) {
    size_t f; while ((f=nameF.find("_"))!=std::string::npos) nameF.erase(f,1); // Remove "_" from var name
    produces<float>(nameF);
  }
  for ( auto nameVI : vectorI_ ) {
    size_t f; while ((f=nameVI.find("_"))!=std::string::npos) nameVI.erase(f,1); // Remove "_" from var name
    produces<std::vector<int> >(nameVI);
  }
  for ( auto nameVF : vectorF_ ) {
    size_t f; while ((f=nameVF.find("_"))!=std::string::npos) nameVF.erase(f,1); // Remove "_" from var name
    produces<std::vector<float> >(nameVF);
  }
  
  // https://twiki.cern.ch/twiki/bin/view/CMS/JECDataMC#Jet_Energy_Corrections
  // Load the JetCorrectorParameter objects into a vector, IMPORTANT: THE ORDER MATTERS HERE !!!! 
  std::vector<JetCorrectorParameters> AK4_vPar, AK8_vPar;
  AK4_vPar.push_back(*(new JetCorrectorParameters(JEC_location_+"_L1FastJet_AK4PFchs.txt")));
  AK4_vPar.push_back(*(new JetCorrectorParameters(JEC_location_+"_L2Relative_AK4PFchs.txt")));
  AK4_vPar.push_back(*(new JetCorrectorParameters(JEC_location_+"_L3Absolute_AK4PFchs.txt")));
  //if (isData_)
  //  AK4_vPar.push_back(*(new JetCorrectorParameters(JEC_location_+"_L2L3Residual_AK4PFchs.txt")));
  
  AK8_vPar.push_back(*(new JetCorrectorParameters(JEC_location_+"_L1FastJet_AK8PFchs.txt")));
  AK8_vPar.push_back(*(new JetCorrectorParameters(JEC_location_+"_L2Relative_AK8PFchs.txt")));
  AK8_vPar.push_back(*(new JetCorrectorParameters(JEC_location_+"_L3Absolute_AK8PFchs.txt")));
  //if (isData_)
  //  AK8_vPar.push_back(*(new JetCorrectorParameters(JEC_location_+"_L2L3Residual_AK8PFchs.txt")));
  
  AK4_JetCorrector_ = new FactorizedJetCorrector(AK4_vPar);
  AK8_JetCorrector_ = new FactorizedJetCorrector(AK8_vPar);
  
  nfilt_=ntrig_=0;
}

void B2GEdmExtraVarProducer::produce(edm::Event& iEvent, const edm::EventSetup& iSetup) {
  // Initialize containers for each variable
  for ( auto nameB : singleB_ ) single_bool_[nameB] = 0;
  for ( auto nameI : trigger_names_ ) single_int_[nameI+"_prescale"] = -9999;
  for ( auto nameI : singleI_ ) single_int_[nameI] = -9999;
  for ( auto nameF : singleF_ ) single_float_[nameF] = -9999.0;
  /* size of vectors are not known, therefor they have to be initialized in calculate_varibles() */
  
  calculate_variables(iEvent, iSetup);
  
  // Put new variables to the event
  for ( auto nameB : singleB_ ) {
    std::auto_ptr<bool> newB(new bool);
    *newB = single_bool_[nameB];
    size_t f; while ((f=nameB.find("_"))!=std::string::npos) nameB.erase(f,1); // Remove "_" from var name
    iEvent.put(newB,nameB);
  }
  for ( auto nameI : trigger_names_ ) {
    std::auto_ptr<int> newI(new int);
    *newI = single_int_[nameI+"_prescale"];
    size_t f; while ((f=nameI.find("_"))!=std::string::npos) nameI.erase(f,1); // Remove "_" from var name
    iEvent.put(newI,nameI+"prescale");
  }
  for ( auto nameI : singleI_ ) {
    std::auto_ptr<int> newI(new int);
    *newI = single_int_[nameI];
    size_t f; while ((f=nameI.find("_"))!=std::string::npos) nameI.erase(f,1); // Remove "_" from var name
    iEvent.put(newI,nameI);
  }
  for ( auto nameF : singleF_ ) {
    std::auto_ptr<float> newF(new float);
    *newF = single_float_[nameF];
    size_t f; while ((f=nameF.find("_"))!=std::string::npos) nameF.erase(f,1); // Remove "_" from var name
    iEvent.put(newF,nameF);
  }
  for ( auto nameVI : vectorI_ ) {
    std::auto_ptr<std::vector<int> > newVI(new std::vector<int>);
    *newVI = vector_int_[nameVI];
    size_t f; while ((f=nameVI.find("_"))!=std::string::npos) nameVI.erase(f,1); // Remove "_" from var name
    iEvent.put(newVI,nameVI);
  }
  for ( auto nameVF : vectorF_ ) {
    std::auto_ptr<std::vector<float> > newVF(new std::vector<float>);
    *newVF = vector_float_[nameVF];
    size_t f; while ((f=nameVF.find("_"))!=std::string::npos) nameVF.erase(f,1); // Remove "_" from var name
    iEvent.put(newVF,nameVF);
  }
}
