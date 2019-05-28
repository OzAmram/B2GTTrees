#include "Analysis/B2GTTrees/interface/B2GEdmExtraVarProducer.h"
#include "Analysis/B2GTTrees/interface/Razor.h"
#include "Analysis/B2GTTrees/data/GluinoXSec.h"
#include "Analysis/B2GTTrees/data/StopXSec.h"
#include "SimDataFormats/GeneratorProducts/interface/GenEventInfoProduct.h"
#include "SimDataFormats/GeneratorProducts/interface/GenLumiInfoHeader.h"
#include "SimDataFormats/GeneratorProducts/interface/LHERunInfoProduct.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "DataFormats/PatCandidates/interface/MET.h"
#include "DataFormats/PatCandidates/interface/PackedCandidate.h"
#include "DataFormats/PatCandidates/interface/Jet.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/Math/interface/deltaR.h"

void B2GEdmExtraVarProducer::init_tokens_() {
    edm::EDGetTokenT<std::vector<std::string> >(mayConsume<std::vector<std::string>, edm::InRun>(edm::InputTag(trigger_label_, "triggerNameTree")));
    edm::EDGetTokenT<std::vector<std::string> >(mayConsume<std::vector<std::string>, edm::InRun>(edm::InputTag(filter_label_,  "triggerNameTree")));
    edm::EDGetTokenT<std::vector<int> >(consumes<std::vector<int> >(edm::InputTag(trigger_label_, "triggerPrescaleTree")));
    edm::EDGetTokenT<std::vector<float> >(consumes<std::vector<float> >(edm::InputTag(trigger_label_, "triggerBitTree")));
    edm::EDGetTokenT<std::vector<float> >(consumes<std::vector<float> >(edm::InputTag(filter_label_,  "triggerBitTree")));

    edm::EDGetTokenT<double>(consumes<double>(edm::InputTag("fixedGridRhoFastjetAll", "")));
    edm::EDGetTokenT<int>(consumes<int>(edm::InputTag(vtx_label_, vtx_prefix_+"npv")));
    edm::EDGetTokenT<std::vector<int> >(consumes<std::vector<int> >(edm::InputTag(vtx_label_, vtx_prefix_+"ndof")));
    edm::EDGetTokenT<std::vector<float> >(consumes<std::vector<float> >(edm::InputTag(vtx_label_, vtx_prefix_+"z")));
    edm::EDGetTokenT<std::vector<float> >(consumes<std::vector<float> >(edm::InputTag(vtx_label_, vtx_prefix_+"rho")));

    edm::EDGetTokenT<std::vector<float> >(consumes<std::vector<float> >(edm::InputTag(met_label_, met_prefix_+"Pt")));
    edm::EDGetTokenT<std::vector<float> >(consumes<std::vector<float> >(edm::InputTag(met_label_, met_prefix_+"Phi")));

    edm::EDGetTokenT<std::vector<float> >(consumes<std::vector<float> >(edm::InputTag(AK4Jets_label_, AK4Jets_prefix_+"Pt")));
    edm::EDGetTokenT<std::vector<float> >(consumes<std::vector<float> >(edm::InputTag(AK4Jets_label_, AK4Jets_prefix_+"Eta")));
    edm::EDGetTokenT<std::vector<float> >(consumes<std::vector<float> >(edm::InputTag(AK4Jets_label_, AK4Jets_prefix_+"Phi")));
    edm::EDGetTokenT<std::vector<float> >(consumes<std::vector<float> >(edm::InputTag(AK4Jets_label_, AK4Jets_prefix_+"E")));
    edm::EDGetTokenT<std::vector<float> >(consumes<std::vector<float> >(edm::InputTag(AK4Jets_label_, AK4Jets_prefix_+"jecFactor0")));
    edm::EDGetTokenT<std::vector<float> >(consumes<std::vector<float> >(edm::InputTag(AK4Jets_label_, AK4Jets_prefix_+"jetArea")));
    edm::EDGetTokenT<std::vector<float> >(consumes<std::vector<float> >(edm::InputTag(AK4Jets_label_, AK4Jets_prefix_+"SmearedPt")));

    edm::EDGetTokenT<std::vector<float> >(consumes<std::vector<float> >(edm::InputTag(AK8Jets_label_, AK8Jets_prefix_+"Pt")));
    edm::EDGetTokenT<std::vector<float> >(consumes<std::vector<float> >(edm::InputTag(AK8Jets_label_, AK8Jets_prefix_+"Eta")));
    edm::EDGetTokenT<std::vector<float> >(consumes<std::vector<float> >(edm::InputTag(AK8Jets_label_, AK8Jets_prefix_+"Phi")));
    edm::EDGetTokenT<std::vector<float> >(consumes<std::vector<float> >(edm::InputTag(AK8Jets_label_, AK8Jets_prefix_+"E")));
    edm::EDGetTokenT<std::vector<float> >(consumes<std::vector<float> >(edm::InputTag(AK8Jets_label_, AK8Jets_prefix_+"tau2")));
    edm::EDGetTokenT<std::vector<float> >(consumes<std::vector<float> >(edm::InputTag(AK8Jets_label_, AK8Jets_prefix_+"tau3")));
    edm::EDGetTokenT<std::vector<float> >(consumes<std::vector<float> >(edm::InputTag(AK8Jets_label_, AK8Jets_prefix_+"softDropMassPuppi")));
    edm::EDGetTokenT<std::vector<float> >(consumes<std::vector<float> >(edm::InputTag(AK8Jets_label_, AK8Jets_prefix_+"jecFactor0")));
    edm::EDGetTokenT<std::vector<float> >(consumes<std::vector<float> >(edm::InputTag(AK8Jets_label_, AK8Jets_prefix_+"jetArea")));
    edm::EDGetTokenT<std::vector<float> >(consumes<std::vector<float> >(edm::InputTag(AK8Jets_label_, AK8Jets_prefix_+"nSubJets")));
    //edm::EDGetTokenT<std::vector<float> >(consumes<std::vector<float> >(edm::InputTag(AK8Jets_label_, AK8Jets_prefix_+"topSubjetIndex0")));
    //edm::EDGetTokenT<std::vector<float> >(consumes<std::vector<float> >(edm::InputTag(AK8Jets_label_, AK8Jets_prefix_+"topSubjetIndex1")));
    //edm::EDGetTokenT<std::vector<float> >(consumes<std::vector<float> >(edm::InputTag(AK8Jets_label_, AK8Jets_prefix_+"topSubjetIndex2")));
    //edm::EDGetTokenT<std::vector<float> >(consumes<std::vector<float> >(edm::InputTag(AK8Jets_label_, AK8Jets_prefix_+"topSubjetIndex3")));
    edm::EDGetTokenT<std::vector<float> >(consumes<std::vector<float> >(edm::InputTag(AK8Jets_label_, AK8Jets_prefix_+"vSubjetIndex0")));
    edm::EDGetTokenT<std::vector<float> >(consumes<std::vector<float> >(edm::InputTag(AK8Jets_label_, AK8Jets_prefix_+"vSubjetIndex1")));

    edm::EDGetTokenT<std::vector<float> >(consumes<std::vector<float> >(edm::InputTag(AK8Subjets_label_, AK8Subjets_prefix_+"Pt")));
    edm::EDGetTokenT<std::vector<float> >(consumes<std::vector<float> >(edm::InputTag(AK8Subjets_label_, AK8Subjets_prefix_+"Eta")));
    edm::EDGetTokenT<std::vector<float> >(consumes<std::vector<float> >(edm::InputTag(AK8Subjets_label_, AK8Subjets_prefix_+"Phi")));
    edm::EDGetTokenT<std::vector<float> >(consumes<std::vector<float> >(edm::InputTag(AK8Subjets_label_, AK8Subjets_prefix_+"E")));
    edm::EDGetTokenT<std::vector<float> >(consumes<std::vector<float> >(edm::InputTag(AK8Subjets_label_, AK8Subjets_prefix_+"jecFactor0")));
    edm::EDGetTokenT<std::vector<float> >(consumes<std::vector<float> >(edm::InputTag(AK8Subjets_label_, AK8Subjets_prefix_+"jetArea")));
    edm::EDGetTokenT<std::vector<float> >(consumes<std::vector<float> >(edm::InputTag(AK8Subjets_label_, AK8Subjets_prefix_+"CSVv2")));
    edm::EDGetTokenT<std::vector<float> >(consumes<std::vector<float> >(edm::InputTag(AK8Subjets_label_, AK8Subjets_prefix_+"CMVAv2")));

    // Jet ID variables
    edm::EDGetTokenT<std::vector<float> >(consumes<std::vector<float> >(edm::InputTag(AK4Jets_label_, AK4Jets_prefix_+"neutralHadronEnergyFrac")));
    edm::EDGetTokenT<std::vector<float> >(consumes<std::vector<float> >(edm::InputTag(AK4Jets_label_, AK4Jets_prefix_+"neutralEmEnergyFrac")));
    edm::EDGetTokenT<std::vector<float> >(consumes<std::vector<float> >(edm::InputTag(AK4Jets_label_, AK4Jets_prefix_+"chargedHadronEnergyFrac")));
    edm::EDGetTokenT<std::vector<float> >(consumes<std::vector<float> >(edm::InputTag(AK4Jets_label_, AK4Jets_prefix_+"chargedEmEnergyFrac")));
    edm::EDGetTokenT<std::vector<float> >(consumes<std::vector<float> >(edm::InputTag(AK4Jets_label_, AK4Jets_prefix_+"chargedMultiplicity")));
    edm::EDGetTokenT<std::vector<float> >(consumes<std::vector<float> >(edm::InputTag(AK4Jets_label_, AK4Jets_prefix_+"neutralMultiplicity")));
    edm::EDGetTokenT<std::vector<float> >(consumes<std::vector<float> >(edm::InputTag(AK4Jets_label_, AK4Jets_prefix_+"MuonEnergy")));
    edm::EDGetTokenT<std::vector<float> >(consumes<std::vector<float> >(edm::InputTag(AK8Jets_label_, AK8Jets_prefix_+"neutralHadronEnergyFrac")));
    edm::EDGetTokenT<std::vector<float> >(consumes<std::vector<float> >(edm::InputTag(AK8Jets_label_, AK8Jets_prefix_+"neutralEmEnergyFrac")));
    edm::EDGetTokenT<std::vector<float> >(consumes<std::vector<float> >(edm::InputTag(AK8Jets_label_, AK8Jets_prefix_+"chargedHadronEnergyFrac")));
    edm::EDGetTokenT<std::vector<float> >(consumes<std::vector<float> >(edm::InputTag(AK8Jets_label_, AK8Jets_prefix_+"chargedEmEnergyFrac")));
    edm::EDGetTokenT<std::vector<float> >(consumes<std::vector<float> >(edm::InputTag(AK8Jets_label_, AK8Jets_prefix_+"chargedMultiplicity")));
    edm::EDGetTokenT<std::vector<float> >(consumes<std::vector<float> >(edm::InputTag(AK8Jets_label_, AK8Jets_prefix_+"neutralMultiplicity")));
    edm::EDGetTokenT<std::vector<float> >(consumes<std::vector<float> >(edm::InputTag(AK8Jets_label_, AK8Jets_prefix_+"MuonEnergy")));
    edm::EDGetTokenT<std::vector<float> >(consumes<std::vector<float> >(edm::InputTag(AK8Subjets_label_, AK8Subjets_prefix_+"neutralHadronEnergyFrac")));
    edm::EDGetTokenT<std::vector<float> >(consumes<std::vector<float> >(edm::InputTag(AK8Subjets_label_, AK8Subjets_prefix_+"neutralEmEnergyFrac")));
    edm::EDGetTokenT<std::vector<float> >(consumes<std::vector<float> >(edm::InputTag(AK8Subjets_label_, AK8Subjets_prefix_+"chargedHadronEnergyFrac")));
    edm::EDGetTokenT<std::vector<float> >(consumes<std::vector<float> >(edm::InputTag(AK8Subjets_label_, AK8Subjets_prefix_+"chargedEmEnergyFrac")));
    edm::EDGetTokenT<std::vector<float> >(consumes<std::vector<float> >(edm::InputTag(AK8Subjets_label_, AK8Subjets_prefix_+"chargedMultiplicity")));
    edm::EDGetTokenT<std::vector<float> >(consumes<std::vector<float> >(edm::InputTag(AK8Subjets_label_, AK8Subjets_prefix_+"neutralMultiplicity")));
    edm::EDGetTokenT<std::vector<float> >(consumes<std::vector<float> >(edm::InputTag(AK8Subjets_label_, AK8Subjets_prefix_+"MuonEnergy")));

    edm::EDGetTokenT<std::vector<std::vector<int> > >(consumes<std::vector<std::vector<int> > >(edm::InputTag(AK4JetKeys_label_,"")));
    edm::EDGetTokenT<std::vector<std::vector<int> > >(consumes<std::vector<std::vector<int> > >(edm::InputTag(AK8JetKeys_label_,"")));
    edm::EDGetTokenT<std::vector<std::vector<int> > >(consumes<std::vector<std::vector<int> > >(edm::InputTag(AK8SubjetKeys_label_,"")));

    edm::EDGetTokenT<std::vector<float> >(consumes<std::vector<float> >(edm::InputTag(electrons_label_, electrons_prefix_+"Pt")));
    edm::EDGetTokenT<std::vector<float> >(consumes<std::vector<float> >(edm::InputTag(electrons_label_, electrons_prefix_+"Eta")));
    edm::EDGetTokenT<std::vector<float> >(consumes<std::vector<float> >(edm::InputTag(electrons_label_, electrons_prefix_+"Phi")));
    edm::EDGetTokenT<std::vector<float> >(consumes<std::vector<float> >(edm::InputTag(electrons_label_, electrons_prefix_+"E")));
    edm::EDGetTokenT<std::vector<float> >(consumes<std::vector<float> >(edm::InputTag(electrons_label_, electrons_prefix_+"Charge")));
    edm::EDGetTokenT<std::vector<float> >(consumes<std::vector<float> >(edm::InputTag(electrons_label_, electrons_prefix_+"Key")));

    // Electron ID variables
    edm::EDGetTokenT<std::vector<float> >(consumes<std::vector<float> >(edm::InputTag(electrons_label_, electrons_prefix_+"SCEta")));
    edm::EDGetTokenT<std::vector<float> >(consumes<std::vector<float> >(edm::InputTag(electrons_label_, electrons_prefix_+"full5x5siee")));
    edm::EDGetTokenT<std::vector<float> >(consumes<std::vector<float> >(edm::InputTag(electrons_label_, electrons_prefix_+"dEtaInSeed")));
    edm::EDGetTokenT<std::vector<float> >(consumes<std::vector<float> >(edm::InputTag(electrons_label_, electrons_prefix_+"dPhiIn")));
    edm::EDGetTokenT<std::vector<float> >(consumes<std::vector<float> >(edm::InputTag(electrons_label_, electrons_prefix_+"HoE")));
    edm::EDGetTokenT<std::vector<float> >(consumes<std::vector<float> >(edm::InputTag(electrons_label_, electrons_prefix_+"Iso03")));
    edm::EDGetTokenT<std::vector<float> >(consumes<std::vector<float> >(edm::InputTag(electrons_label_, electrons_prefix_+"ooEmooP")));
    edm::EDGetTokenT<std::vector<float> >(consumes<std::vector<float> >(edm::InputTag(electrons_label_, electrons_prefix_+"Dxy")));
    edm::EDGetTokenT<std::vector<float> >(consumes<std::vector<float> >(edm::InputTag(electrons_label_, electrons_prefix_+"Dz")));
    edm::EDGetTokenT<std::vector<float> >(consumes<std::vector<float> >(edm::InputTag(electrons_label_, electrons_prefix_+"missHits")));
    edm::EDGetTokenT<std::vector<float> >(consumes<std::vector<float> >(edm::InputTag(electrons_label_, electrons_prefix_+"hasMatchedConVeto")));

    edm::EDGetTokenT<std::vector<float> >(consumes<std::vector<float> >(edm::InputTag(muons_label_, muons_prefix_+"Pt")));
    edm::EDGetTokenT<std::vector<float> >(consumes<std::vector<float> >(edm::InputTag(muons_label_, muons_prefix_+"Eta")));
    edm::EDGetTokenT<std::vector<float> >(consumes<std::vector<float> >(edm::InputTag(muons_label_, muons_prefix_+"Phi")));
    edm::EDGetTokenT<std::vector<float> >(consumes<std::vector<float> >(edm::InputTag(muons_label_, muons_prefix_+"E")));
    edm::EDGetTokenT<std::vector<float> >(consumes<std::vector<float> >(edm::InputTag(muons_label_, muons_prefix_+"Charge")));
    edm::EDGetTokenT<std::vector<float> >(consumes<std::vector<float> >(edm::InputTag(muons_label_, muons_prefix_+"IsTightMuon")));
    edm::EDGetTokenT<std::vector<float> >(consumes<std::vector<float> >(edm::InputTag(muons_label_, muons_prefix_+"Key")));

    edm::EDGetTokenT<reco::VertexCollection>(consumes<reco::VertexCollection>(edm::InputTag("offlineSlimmedPrimaryVertices")));
    edm::EDGetTokenT<pat::PackedCandidateCollection>(consumes<pat::PackedCandidateCollection>(edm::InputTag("packedPFCandidates")));
    edm::EDGetTokenT<pat::METCollection>(consumes<pat::METCollection>(edm::InputTag("slimmedMETs")));
    //if (isData_) edm::EDGetTokenT<pat::METCollection>(consumes<pat::METCollection>(edm::InputTag("slimmedMETsMuEGClean","","b2gEDMNtuples")));
    //else         edm::EDGetTokenT<pat::METCollection>(consumes<pat::METCollection>(edm::InputTag("slimmedMETsMuClean",  "","b2gEDMNtuples")));
    edm::EDGetTokenT<pat::METCollection>(consumes<pat::METCollection>(edm::InputTag("slimmedMETsPuppi")));
    //edm::EDGetTokenT<pat::JetCollection>(consumes<pat::JetCollection>(edm::InputTag("slimmedJetsAK8")));

    if (!isData_) {
        edm::EDGetTokenT<std::vector<float> >(consumes<std::vector<float> >(edm::InputTag(gen_label_, gen_prefix_+"Pt")));
        edm::EDGetTokenT<std::vector<float> >(consumes<std::vector<float> >(edm::InputTag(gen_label_, gen_prefix_+"Eta")));
        edm::EDGetTokenT<std::vector<float> >(consumes<std::vector<float> >(edm::InputTag(gen_label_, gen_prefix_+"Phi")));
        edm::EDGetTokenT<std::vector<float> >(consumes<std::vector<float> >(edm::InputTag(gen_label_, gen_prefix_+"E")));
        edm::EDGetTokenT<std::vector<float> >(consumes<std::vector<float> >(edm::InputTag(gen_label_, gen_prefix_+"Charge")));
        edm::EDGetTokenT<std::vector<float> >(consumes<std::vector<float> >(edm::InputTag(gen_label_, gen_prefix_+"ID")));
        edm::EDGetTokenT<std::vector<float> >(consumes<std::vector<float> >(edm::InputTag(gen_label_, gen_prefix_+"Status")));
        edm::EDGetTokenT<std::vector<float> >(consumes<std::vector<float> >(edm::InputTag(gen_label_, gen_prefix_+"Mom0ID")));
        edm::EDGetTokenT<std::vector<float> >(consumes<std::vector<float> >(edm::InputTag(gen_label_, gen_prefix_+"Mom0Status")));
        edm::EDGetTokenT<std::vector<float> >(consumes<std::vector<float> >(edm::InputTag(gen_label_, gen_prefix_+"Mom1ID")));
        edm::EDGetTokenT<std::vector<float> >(consumes<std::vector<float> >(edm::InputTag(gen_label_, gen_prefix_+"Mom1Status")));
        edm::EDGetTokenT<std::vector<float> >(consumes<std::vector<float> >(edm::InputTag(gen_label_, gen_prefix_+"Dau0ID")));
        edm::EDGetTokenT<std::vector<float> >(consumes<std::vector<float> >(edm::InputTag(gen_label_, gen_prefix_+"Dau0Status")));
        edm::EDGetTokenT<std::vector<float> >(consumes<std::vector<float> >(edm::InputTag(gen_label_, gen_prefix_+"Dau1ID")));
        edm::EDGetTokenT<std::vector<float> >(consumes<std::vector<float> >(edm::InputTag(gen_label_, gen_prefix_+"Dau1Status")));

        edm::EDGetTokenT<LHERunInfoProduct>(mayConsume<LHERunInfoProduct, edm::InRun>(edm::InputTag(lhe_label_, "")));
        edm::EDGetTokenT<LHEEventProduct>(consumes<LHEEventProduct>(edm::InputTag(lhe_label_, "")));

        edm::EDGetTokenT<GenLumiInfoHeader>(mayConsume<GenLumiInfoHeader, edm::InLumi>(edm::InputTag("generator")));
        edm::EDGetTokenT<GenEventInfoProduct>(consumes<GenEventInfoProduct>(edm::InputTag("generator")));

        edm::EDGetTokenT<reco::GenParticleCollection>(consumes<reco::GenParticleCollection>(edm::InputTag("prunedGenParticles")));
    }
}

void B2GEdmExtraVarProducer::beginRun(edm::Run const& iRun, edm::EventSetup const& iSetup) {

    if (!isData_) {
        // Weights from scale variations, PDFs etc. are stored in the relative product. 
        // Notice that to be used they need to be renormalized to the central event weight
        // at LHE level which may be different from genEvtInfo->weight()

        //int whichWeight = XXX;
        //theWeight *= EvtHandle->weights()[whichWeight].wgt/EvtHandle->originalXWGTUP(); 

        //To know which integer XXX corresponds to which weight you can use:
        edm::Handle<LHERunInfoProduct> lheRunInfo;
        iRun.getByLabel(lhe_label_, lheRunInfo);

        if (lheRunInfo.isValid()) {

            // Check which PDF set was used
            // HEPRUP reference: http://arxiv.org/pdf/hep-ph/0609017.pdf
            // ID reference: https://lhapdf.hepforge.org/pdfsets.html
            lha_pdf_id_ = lheRunInfo->heprup().PDFSUP.first;
            std::cout<<"LHE: LHA PDF ID = "<<lha_pdf_id_<<std::endl;
            std::cout<<"LHE:   --> For more info about the sets, check: https://lhapdf.hepforge.org/pdfsets.html"<<std::endl;

            // Check headers
            std::cout<<"LHE: Weight info in header:"<<std::endl;
            LHERunInfoProduct lheRunInfoProduct = *(lheRunInfo.product());
            typedef std::vector<LHERunInfoProduct::Header>::const_iterator headers_const_iterator;
            size_t iHead = 0;
            for (headers_const_iterator header=lheRunInfoProduct.headers_begin(); header!=lheRunInfoProduct.headers_end(); header++){
                if (header->tag()=="initrwgt") {
                    std::cout<<"LHE: "<<iHead<<" "<<header->tag()<<std::endl;
                    for (auto line : header->lines()) {
                        std::cout<<"LHE: "<<line;
                        // Fix buggy powheg samples
                        if (lha_pdf_id_==-1 && line.find("weight id=\"2001\"")!=std::string::npos) {
                            if (line.find("PDF set = 260001")!=std::string::npos) lha_pdf_id_ = 260000;
                            else if (line.find("PDF set = 260401")!=std::string::npos) lha_pdf_id_ = 260400;
                        }
                    }
                }
                iHead++;
            }

        }
    }

    // ----------------------------
    // - Trigger/Filter names     -
    // ----------------------------

    iRun.getByLabel(edm::InputTag(trigger_label_, "triggerNameTree"),      h_strings_["trigger_names"]);
    iRun.getByLabel(edm::InputTag(filter_label_,  "triggerNameTree"),      h_strings_["filter_names"]);
    bool print_all = false;

    nfilt_=h_strings_["filter_names"]->size();
    filters_.clear();
    for ( auto filter : filter_names_ ) for (size_t i=0; i<nfilt_; ++i) 
        if (h_strings_["filter_names"]->at(i).find(filter)==0) filters_[filter] = i;
    std::cout<<"Filters found: "<<std::endl;
    if (print_all) for (size_t i=0; i<nfilt_; ++i) std::cout<<h_strings_["filter_names"]->at(i)<<std::endl;
    else for ( auto filter : filters_ ) std::cout<<filter.first<<std::endl;

    ntrig_=h_strings_["trigger_names"]->size();
    triggers_.clear();
    for ( auto trig : trigger_names_ ) for (size_t i=0; i<ntrig_; ++i) 
        if (h_strings_["trigger_names"]->at(i).find(trig+"_v")==0) triggers_[trig] = i;
    std::cout<<"Triggers found: "<<std::endl;
    if (print_all) for (size_t i=0; i<ntrig_; ++i) std::cout<<h_strings_["trigger_names"]->at(i)<<std::endl;
    else for ( auto trigger : triggers_ ) std::cout<<trigger.first<<std::endl;
}

void B2GEdmExtraVarProducer::beginLuminosityBlock(edm::LuminosityBlock const& iLumi, edm::EventSetup const& iSetup) {
    if (!isData_) {
        // SUSY signal specific info in 80X
        edm::Handle<GenLumiInfoHeader> genLumiInfo;
        iLumi.getByLabel(edm::InputTag("generator"), genLumiInfo);

        if (genLumiInfo.isValid()) {

            // Check which PDF set was used
            lha_pdf_id_ = 263000;
            std::cout<<"GEN: LHA PDF ID = "<<lha_pdf_id_<<std::endl;
            std::cout<<"GEN:   --> For more info about the sets, check: https://lhapdf.hepforge.org/pdfsets.html"<<std::endl;

            // Print headers
            //std::cout<<"GEN: Weight info in Lumi:"<<std::endl;
            //for (auto wname : genLumiInfo->weightNames()) std::cout<<"GEN:   "<<wname<<std::endl;
            //std::cout<<"GEN: Header info in Lumi:"<<std::endl;
            //for (auto header : genLumiInfo->lheHeaders()) std::cout<<header.first<<"\n"<<header.second<<std::endl;
        }
    }
}

void B2GEdmExtraVarProducer::calculate_variables(edm::Event const& iEvent, edm::EventSetup const& iSetup) {
    // Read variables from EdmNtuple
    iEvent.getByLabel(edm::InputTag(trigger_label_, "triggerBitTree"),       h_floats_["trigger_bits"]);
    iEvent.getByLabel(edm::InputTag(trigger_label_, "triggerPrescaleTree"),  h_ints_["trigger_prescales"]);
    iEvent.getByLabel(edm::InputTag(filter_label_,  "triggerBitTree"),       h_floats_["filter_bits"]);

    iEvent.getByLabel(edm::InputTag("fixedGridRhoFastjetAll", ""),   h_double_["evt_rho"]);
    iEvent.getByLabel(edm::InputTag(vtx_label_, vtx_prefix_+"npv"),  h_int_["vtx_npv"]);
    iEvent.getByLabel(edm::InputTag(vtx_label_, vtx_prefix_+"ndof"), h_ints_["vtx_ndof"]);
    iEvent.getByLabel(edm::InputTag(vtx_label_, vtx_prefix_+"z"),    h_floats_["vtx_z"]);
    iEvent.getByLabel(edm::InputTag(vtx_label_, vtx_prefix_+"rho"),  h_floats_["vtx_rho"]);

    iEvent.getByLabel(edm::InputTag(met_label_, met_prefix_+"Pt"),  h_floats_["met_Pt"]);
    iEvent.getByLabel(edm::InputTag(met_label_, met_prefix_+"Phi"), h_floats_["met_Phi"]);

    iEvent.getByLabel(edm::InputTag(AK4Jets_label_, AK4Jets_prefix_+"Pt"),         h_floats_["AK4_Pt"]);
    iEvent.getByLabel(edm::InputTag(AK4Jets_label_, AK4Jets_prefix_+"Eta"),        h_floats_["AK4_Eta"]);
    iEvent.getByLabel(edm::InputTag(AK4Jets_label_, AK4Jets_prefix_+"Phi"),        h_floats_["AK4_Phi"]);
    iEvent.getByLabel(edm::InputTag(AK4Jets_label_, AK4Jets_prefix_+"E"),          h_floats_["AK4_E"]);
    iEvent.getByLabel(edm::InputTag(AK4Jets_label_, AK4Jets_prefix_+"jecFactor0"), h_floats_["AK4_jecFactor0"]);
    iEvent.getByLabel(edm::InputTag(AK4Jets_label_, AK4Jets_prefix_+"jetArea"),    h_floats_["AK4_jetArea"]);
    iEvent.getByLabel(edm::InputTag(AK4Jets_label_, AK4Jets_prefix_+"SmearedPt"),  h_floats_["AK4_SmearedPt"]);

    iEvent.getByLabel(edm::InputTag(AK8Jets_label_, AK8Jets_prefix_+"Pt"),                h_floats_["AK8_Pt"]);
    iEvent.getByLabel(edm::InputTag(AK8Jets_label_, AK8Jets_prefix_+"Eta"),               h_floats_["AK8_Eta"]);
    iEvent.getByLabel(edm::InputTag(AK8Jets_label_, AK8Jets_prefix_+"Phi"),               h_floats_["AK8_Phi"]);
    iEvent.getByLabel(edm::InputTag(AK8Jets_label_, AK8Jets_prefix_+"E"),                 h_floats_["AK8_E"]);
    iEvent.getByLabel(edm::InputTag(AK8Jets_label_, AK8Jets_prefix_+"tau2"),              h_floats_["AK8_tau2"]);
    iEvent.getByLabel(edm::InputTag(AK8Jets_label_, AK8Jets_prefix_+"tau3"),              h_floats_["AK8_tau3"]);
    iEvent.getByLabel(edm::InputTag(AK8Jets_label_, AK8Jets_prefix_+"softDropMassPuppi"), h_floats_["AK8_softDropMassPuppi"]);
    iEvent.getByLabel(edm::InputTag(AK8Jets_label_, AK8Jets_prefix_+"jecFactor0"),        h_floats_["AK8_jecFactor0"]);
    iEvent.getByLabel(edm::InputTag(AK8Jets_label_, AK8Jets_prefix_+"jetArea"),           h_floats_["AK8_jetArea"]);
    iEvent.getByLabel(edm::InputTag(AK8Jets_label_, AK8Jets_prefix_+"nSubJets"),          h_floats_["AK8_nSubJets"]);
    //iEvent.getByLabel(edm::InputTag(AK8Jets_label_, AK8Jets_prefix_+"topSubjetIndex0"),   h_floats_["AK8_topSubjetIndex0"]);
    //iEvent.getByLabel(edm::InputTag(AK8Jets_label_, AK8Jets_prefix_+"topSubjetIndex1"),   h_floats_["AK8_topSubjetIndex1"]);
    //iEvent.getByLabel(edm::InputTag(AK8Jets_label_, AK8Jets_prefix_+"topSubjetIndex2"),   h_floats_["AK8_topSubjetIndex2"]);
    //iEvent.getByLabel(edm::InputTag(AK8Jets_label_, AK8Jets_prefix_+"topSubjetIndex3"),   h_floats_["AK8_topSubjetIndex3"]);
    iEvent.getByLabel(edm::InputTag(AK8Jets_label_, AK8Jets_prefix_+"vSubjetIndex0"),     h_floats_["AK8_vSubjetIndex0"]);
    iEvent.getByLabel(edm::InputTag(AK8Jets_label_, AK8Jets_prefix_+"vSubjetIndex1"),     h_floats_["AK8_vSubjetIndex1"]);

    iEvent.getByLabel(edm::InputTag(AK8Subjets_label_, AK8Subjets_prefix_+"Pt"),         h_floats_["AK8Sub_Pt"]);
    iEvent.getByLabel(edm::InputTag(AK8Subjets_label_, AK8Subjets_prefix_+"Eta"),        h_floats_["AK8Sub_Eta"]);
    iEvent.getByLabel(edm::InputTag(AK8Subjets_label_, AK8Subjets_prefix_+"Phi"),        h_floats_["AK8Sub_Phi"]);
    iEvent.getByLabel(edm::InputTag(AK8Subjets_label_, AK8Subjets_prefix_+"E"),          h_floats_["AK8Sub_E"]);
    iEvent.getByLabel(edm::InputTag(AK8Subjets_label_, AK8Subjets_prefix_+"jecFactor0"), h_floats_["AK8Sub_jecFactor0"]);
    iEvent.getByLabel(edm::InputTag(AK8Subjets_label_, AK8Subjets_prefix_+"jetArea"),    h_floats_["AK8Sub_jetArea"]);
    iEvent.getByLabel(edm::InputTag(AK8Subjets_label_, AK8Subjets_prefix_+"CSVv2"),      h_floats_["AK8Sub_CSVv2"]);
    iEvent.getByLabel(edm::InputTag(AK8Subjets_label_, AK8Subjets_prefix_+"CMVAv2"),     h_floats_["AK8Sub_CMVAv2"]);

    iEvent.getByLabel(edm::InputTag(AK4JetKeys_label_),      h_keys_["AK4"]);
    iEvent.getByLabel(edm::InputTag(AK8JetKeys_label_),      h_keys_["AK8"]);
    iEvent.getByLabel(edm::InputTag(AK8SubjetKeys_label_),   h_keys_["AK8Sub"]);

    // JET ID
    iEvent.getByLabel(edm::InputTag(AK4Jets_label_, AK4Jets_prefix_+"neutralHadronEnergyFrac"), h_floats_["AK4_neutralHadronEnergyFrac"]);
    iEvent.getByLabel(edm::InputTag(AK4Jets_label_, AK4Jets_prefix_+"neutralEmEnergyFrac"),     h_floats_["AK4_neutralEmEnergyFrac"]);
    iEvent.getByLabel(edm::InputTag(AK4Jets_label_, AK4Jets_prefix_+"chargedHadronEnergyFrac"), h_floats_["AK4_chargedHadronEnergyFrac"]);
    iEvent.getByLabel(edm::InputTag(AK4Jets_label_, AK4Jets_prefix_+"chargedEmEnergyFrac"),     h_floats_["AK4_chargedEmEnergyFrac"]);
    iEvent.getByLabel(edm::InputTag(AK4Jets_label_, AK4Jets_prefix_+"chargedMultiplicity"),     h_floats_["AK4_chargedMultiplicity"]);
    iEvent.getByLabel(edm::InputTag(AK4Jets_label_, AK4Jets_prefix_+"neutralMultiplicity"),     h_floats_["AK4_neutralMultiplicity"]);
    iEvent.getByLabel(edm::InputTag(AK4Jets_label_, AK4Jets_prefix_+"MuonEnergy"),              h_floats_["AK4_MuonEnergy"]);
    iEvent.getByLabel(edm::InputTag(AK8Jets_label_, AK8Jets_prefix_+"neutralHadronEnergyFrac"), h_floats_["AK8_neutralHadronEnergyFrac"]);
    iEvent.getByLabel(edm::InputTag(AK8Jets_label_, AK8Jets_prefix_+"neutralEmEnergyFrac"),     h_floats_["AK8_neutralEmEnergyFrac"]);
    iEvent.getByLabel(edm::InputTag(AK8Jets_label_, AK8Jets_prefix_+"chargedHadronEnergyFrac"), h_floats_["AK8_chargedHadronEnergyFrac"]);
    iEvent.getByLabel(edm::InputTag(AK8Jets_label_, AK8Jets_prefix_+"chargedEmEnergyFrac"),     h_floats_["AK8_chargedEmEnergyFrac"]);
    iEvent.getByLabel(edm::InputTag(AK8Jets_label_, AK8Jets_prefix_+"chargedMultiplicity"),     h_floats_["AK8_chargedMultiplicity"]);
    iEvent.getByLabel(edm::InputTag(AK8Jets_label_, AK8Jets_prefix_+"neutralMultiplicity"),     h_floats_["AK8_neutralMultiplicity"]);
    iEvent.getByLabel(edm::InputTag(AK8Jets_label_, AK8Jets_prefix_+"MuonEnergy"),              h_floats_["AK8_MuonEnergy"]);
    iEvent.getByLabel(edm::InputTag(AK8Subjets_label_, AK8Subjets_prefix_+"neutralHadronEnergyFrac"), h_floats_["AK8Sub_neutralHadronEnergyFrac"]);
    iEvent.getByLabel(edm::InputTag(AK8Subjets_label_, AK8Subjets_prefix_+"neutralEmEnergyFrac"),     h_floats_["AK8Sub_neutralEmEnergyFrac"]);
    iEvent.getByLabel(edm::InputTag(AK8Subjets_label_, AK8Subjets_prefix_+"chargedHadronEnergyFrac"), h_floats_["AK8Sub_chargedHadronEnergyFrac"]);
    iEvent.getByLabel(edm::InputTag(AK8Subjets_label_, AK8Subjets_prefix_+"chargedEmEnergyFrac"),     h_floats_["AK8Sub_chargedEmEnergyFrac"]);
    iEvent.getByLabel(edm::InputTag(AK8Subjets_label_, AK8Subjets_prefix_+"chargedMultiplicity"),     h_floats_["AK8Sub_chargedMultiplicity"]);
    iEvent.getByLabel(edm::InputTag(AK8Subjets_label_, AK8Subjets_prefix_+"neutralMultiplicity"),     h_floats_["AK8Sub_neutralMultiplicity"]);
    iEvent.getByLabel(edm::InputTag(AK8Subjets_label_, AK8Subjets_prefix_+"MuonEnergy"),              h_floats_["AK8Sub_MuonEnergy"]);

    // Leptons
    iEvent.getByLabel(edm::InputTag(electrons_label_, electrons_prefix_+"Pt"),     h_floats_["ele_Pt"]);
    iEvent.getByLabel(edm::InputTag(electrons_label_, electrons_prefix_+"Eta"),    h_floats_["ele_Eta"]);
    iEvent.getByLabel(edm::InputTag(electrons_label_, electrons_prefix_+"Phi"),    h_floats_["ele_Phi"]);
    iEvent.getByLabel(edm::InputTag(electrons_label_, electrons_prefix_+"E"),      h_floats_["ele_E"]);
    iEvent.getByLabel(edm::InputTag(electrons_label_, electrons_prefix_+"Charge"), h_floats_["ele_Charge"]);
    iEvent.getByLabel(edm::InputTag(electrons_label_, electrons_prefix_+"Key"),    h_floats_["ele_Key"]);

    // Electron ID variables
    iEvent.getByLabel(edm::InputTag(electrons_label_, electrons_prefix_+"SCEta"),             h_floats_["ele_SCEta"]);
    iEvent.getByLabel(edm::InputTag(electrons_label_, electrons_prefix_+"full5x5siee"),       h_floats_["ele_full5x5siee"]);
    iEvent.getByLabel(edm::InputTag(electrons_label_, electrons_prefix_+"dEtaInSeed"),        h_floats_["ele_dEtaInSeed"]);
    iEvent.getByLabel(edm::InputTag(electrons_label_, electrons_prefix_+"dPhiIn"),            h_floats_["ele_dPhiIn"]);
    iEvent.getByLabel(edm::InputTag(electrons_label_, electrons_prefix_+"HoE"),               h_floats_["ele_HoE"]);
    iEvent.getByLabel(edm::InputTag(electrons_label_, electrons_prefix_+"Iso03"),             h_floats_["ele_Iso03"]);
    iEvent.getByLabel(edm::InputTag(electrons_label_, electrons_prefix_+"ooEmooP"),           h_floats_["ele_ooEmooP"]);
    iEvent.getByLabel(edm::InputTag(electrons_label_, electrons_prefix_+"Dxy"),               h_floats_["ele_Dxy"]);
    iEvent.getByLabel(edm::InputTag(electrons_label_, electrons_prefix_+"Dz"),                h_floats_["ele_Dz"]);
    iEvent.getByLabel(edm::InputTag(electrons_label_, electrons_prefix_+"missHits"),          h_floats_["ele_missHits"]);
    iEvent.getByLabel(edm::InputTag(electrons_label_, electrons_prefix_+"hasMatchedConVeto"), h_floats_["ele_hasMatchedConVeto"]);

    iEvent.getByLabel(edm::InputTag(muons_label_, muons_prefix_+"Pt"),          h_floats_["mu_Pt"]);
    iEvent.getByLabel(edm::InputTag(muons_label_, muons_prefix_+"Eta"),         h_floats_["mu_Eta"]);
    iEvent.getByLabel(edm::InputTag(muons_label_, muons_prefix_+"Phi"),         h_floats_["mu_Phi"]);
    iEvent.getByLabel(edm::InputTag(muons_label_, muons_prefix_+"E"),           h_floats_["mu_E"]);
    iEvent.getByLabel(edm::InputTag(muons_label_, muons_prefix_+"Charge"),      h_floats_["mu_Charge"]);
    iEvent.getByLabel(edm::InputTag(muons_label_, muons_prefix_+"IsTightMuon"), h_floats_["mu_IsTightMuon"]);
    iEvent.getByLabel(edm::InputTag(muons_label_, muons_prefix_+"Key"),         h_floats_["mu_Key"]);

    if (!isData_) {
        iEvent.getByLabel(edm::InputTag(gen_label_, gen_prefix_+"Pt"),          h_floats_["gen_Pt"]);
        iEvent.getByLabel(edm::InputTag(gen_label_, gen_prefix_+"Eta"),         h_floats_["gen_Eta"]);
        iEvent.getByLabel(edm::InputTag(gen_label_, gen_prefix_+"Phi"),         h_floats_["gen_Phi"]);
        iEvent.getByLabel(edm::InputTag(gen_label_, gen_prefix_+"E"),           h_floats_["gen_E"]);
        iEvent.getByLabel(edm::InputTag(gen_label_, gen_prefix_+"Charge"),      h_floats_["gen_Charge"]);
        iEvent.getByLabel(edm::InputTag(gen_label_, gen_prefix_+"ID"),          h_floats_["gen_ID"]);
        iEvent.getByLabel(edm::InputTag(gen_label_, gen_prefix_+"Status"),      h_floats_["gen_Status"]);
        iEvent.getByLabel(edm::InputTag(gen_label_, gen_prefix_+"Mom0ID"),      h_floats_["gen_Mom0ID"]);
        iEvent.getByLabel(edm::InputTag(gen_label_, gen_prefix_+"Mom0Status"),  h_floats_["gen_Mom0Status"]);
        iEvent.getByLabel(edm::InputTag(gen_label_, gen_prefix_+"Mom1ID"),      h_floats_["gen_Mom1ID"]);
        iEvent.getByLabel(edm::InputTag(gen_label_, gen_prefix_+"Mom1Status"),  h_floats_["gen_Mom1Status"]);
        iEvent.getByLabel(edm::InputTag(gen_label_, gen_prefix_+"Dau0ID"),      h_floats_["gen_Dau0ID"]);
        iEvent.getByLabel(edm::InputTag(gen_label_, gen_prefix_+"Dau0Status"),  h_floats_["gen_Dau0Status"]);
        iEvent.getByLabel(edm::InputTag(gen_label_, gen_prefix_+"Dau1ID"),      h_floats_["gen_Dau1ID"]);
        iEvent.getByLabel(edm::InputTag(gen_label_, gen_prefix_+"Dau1Status"),  h_floats_["gen_Dau1Status"]);
    }

    // ----------------------------
    // -  LHE/Gen/Event weights   -
    // ----------------------------

    // Event weight (xsec/nevent in units of fb), 
    // Usage: Multiply this number by the total int luminosity in units of fb^-1
    single_int_["evt_LHA_PDF_ID"] = -9999;                                               /* evt_LHA_PDF_ID */
    single_float_["evt_XSec"] = isData_ ? -9999 : cross_section_;                        /* evt_Xsec */
    single_float_["evt_Gen_Weight"] = -9999;                                             /* evt_Gen_Weight */
    single_float_["evt_Gen_Ht"] = -9999;                                                 /* evt_Gen_Ht */

    // NLO negative weights, see:
    // https://twiki.cern.ch/twiki/bin/view/CMSPublic/SWGuideDataFormatGeneratorInterface
    // https://twiki.cern.ch/twiki/bin/viewauth/CMS/LHEReaderCMSSW
    //
    //  Comment:
    //    The previous weight were of course not good, due to negative weights
    //    This will not give the correct weight for 1 fb^-1 either
    //    but s = Sum(weight) / Sum( abs(weight) ) fb^-1 instead
    //    In order to get the correct weight for 1 fb^-1, one has to calculate
    //    s on the whole dataset and multiply the currently set weight with 1/s

    // Gen/LHE info
    // For Run II recommendations see:
    // https://indico.cern.ch/event/459797/contributions/1961581/attachments/1181555/1800214/mcaod-Feb15-2016.pdf
    vector_float_["scale_Weights"].clear();
    vector_float_["pdf_Weights"].clear();
    vector_float_["alphas_Weights"].clear();
    if (!isData_) {
        edm::Handle<GenEventInfoProduct> genEvtInfo;
        iEvent.getByLabel("generator", genEvtInfo);

        edm::Handle<LHEEventProduct> lheEvtInfo;
        iEvent.getByLabel(lhe_label_, lheEvtInfo);

        single_int_["evt_LHA_PDF_ID"] = lha_pdf_id_;

        // Generator level HT
        if (lheEvtInfo.isValid()) {
            lhef::HEPEUP lheParticleInfo = lheEvtInfo->hepeup();
            // Get the five vector (Px, Py, Pz, E and M in GeV)
            std::vector<lhef::HEPEUP::FiveVector> allParticles = lheParticleInfo.PUP;
            std::vector<int> statusCodes = lheParticleInfo.ISTUP;
            single_float_["evt_Gen_Ht"] = 0;
            for (unsigned int i = 0; i < allParticles.size(); i++) {
                auto absId = abs(lheParticleInfo.IDUP[i]);
                if (statusCodes[i] == 1 && ( absId < 11 || absId > 16 ) && absId != 22 && !hasAncestor_(i, lheParticleInfo, 6))
                    single_float_["evt_Gen_Ht"] += std::sqrt(std::pow(allParticles[i][0], 2) + std::pow(allParticles[i][1], 2));
            }
        }

        if (lheEvtInfo.isValid() && genEvtInfo.isValid()) {
            // GenHT
            //const lhef::HEPEUP& lheEvent = lheEvtInfo->hepeup();
            //std::vector<lhef::HEPEUP::FiveVector> lheParticles = lheEvent.PUP;
            //single_float_["evt_Gen_Ht"] = 0;
            //for ( size_t idxParticle = 0, numParticles = lheParticles.size(); idxParticle < numParticles; ++idxParticle ) {
            //  int absPdgId = TMath::Abs(lheEvent.IDUP[idxParticle]);
            //  int status = lheEvent.ISTUP[idxParticle];
            //  if ( status == 1 && ((absPdgId >= 1 && absPdgId <= 6) || absPdgId == 21) )   // quarks and gluons
            //    single_float_["evt_Gen_Ht"] += TMath::Sqrt(TMath::Power(lheParticles[idxParticle][0], 2.) + 
            //  					     TMath::Power(lheParticles[idxParticle][1], 2.)); // first entry is px, second py
            //}

            // We only look for the sign of the gen weight but not it's value
            // The xsec weight is already calculated, but certain NLO samples
            // have negative weights and such events need to be subtracted
            // (and consequently the weight also needs to be corrected with a factor)
            single_float_["evt_Gen_Weight"] = genEvtInfo->weight();

            // This weight is used to normalize pdf/scale etc weights
            double lheOrigWeight = lheEvtInfo->originalXWGTUP();

            // Print factors for an event
            //if (nfilt_!=h_strings_["filter_names"]->size()) for (size_t i=0; i<lheEvtInfo->weights().size(); ++i)
            //  std::cout<<"LHE: weights() - index: "<<i<<" id = "<<lheEvtInfo->weights()[i].id<<" wgt = "<<(lheEvtInfo->weights()[i].wgt/lheOrigWeight)<<std::endl;

            // Renormalization/Factorization scale weights
            // These are the first 9 weights for all generators
            // mu_R and mu_F are varied independently (by a factor of 1, 2, 0.5) - check LHE header
            // [0] is the default weight (no variation) - it has worse precision even
            // --> I skip saving it (one can just use 1)
            // --> Also do not save unphysical combinations as recommended
            //    (mu_R = 0.5, mu_F = 2 and mu_R = 2, mu_F = 0.5)
            // Save only: 1,2,3,4,6,8
            if (lheEvtInfo->weights().size()>=9) for (size_t i=0; i<9; ++i) if (i!=0&&i!=5&&i!=7)
                vector_float_["scale_Weights"].push_back(lheEvtInfo->weights()[i].wgt/lheOrigWeight);

            // PDF weights
            // Usually a set of 100 weights (excluding default)
            // Only default PDF variation is saved, but if needed
            // likely others are available depending on the generator
            // index of first weight varies, beware!
            // Sometimes first weight is default=1 weight (has to skip!)
            // Additional info: MC2Hessian conversion will soon be provided
            size_t first = 9;
            // Madgraph nf5 - have to skip 1 weight which is default
            if (lha_pdf_id_ == 263000) first = 10;
            // Madgraph nf4 uses NNPDF30_lo_as_0130_nf_4 (ID=263400)
            // Which is the second 101 pdf set, again has to skip first weight
            if (lha_pdf_id_ == 263400) first = 111;
            // FastSim indices are larger by 1
            if (isFastSim_) ++first;

            //fill vector of set replica weights
            const int nRepWeights = 100;
            const int nEigWeights = 60;

            std::vector<double> inpdfweights(nRepWeights);
            if (lheEvtInfo->weights().size()>=first+100) for (size_t i=0; i<100; ++i){
                int idx = i+first;
                inpdfweights[i] = lheEvtInfo->weights()[idx].wgt;
            }

            std::vector<double> outpdfweights(nEigWeights);
            // do the actual conversion, where the nominal lhe weight is
            //needed as the reference point for the linearization
            pdfweightshelper_.DoMC2Hessian(lheOrigWeight,inpdfweights.data(),outpdfweights.data());
            for (unsigned int iwgt=0; iwgt<nEigWeights; ++iwgt) {
                double wgtval = outpdfweights[iwgt];

                //the is the weight to be used for evaluating uncertainties with hessian weights
                vector_float_["pdf_Weights"].push_back(wgtval/lheOrigWeight);
            }    

            // Alpha_s weights (only for NLO!)
            // A set of two weights for 
            // alpha_s = 0.118 - 0.002 and
            // alpha_s = 0.118 + 0.002
            // is given --> scale result uncertainty by 0.75
            if ( lheEvtInfo->weights().size()>=111 &&
                    ( (lha_pdf_id_ == 260000) || // Powheg 5nf
                      (lha_pdf_id_ == 260400) || // Powheg 4nf 
                      (lha_pdf_id_ == 292000) || // aMC@NLO 5nf
                      (lha_pdf_id_ == 292200)    // aMC@NLO 5nf
                    ) ) {
                if (isFastSim_) {
                    vector_float_["alphas_Weights"].push_back(lheEvtInfo->weights()[110].wgt/lheOrigWeight);
                    vector_float_["alphas_Weights"].push_back(lheEvtInfo->weights()[111].wgt/lheOrigWeight);
                } else {
                    vector_float_["alphas_Weights"].push_back(lheEvtInfo->weights()[109].wgt/lheOrigWeight);
                    vector_float_["alphas_Weights"].push_back(lheEvtInfo->weights()[110].wgt/lheOrigWeight);
                }
            }
        } else if (genEvtInfo.isValid()) {
            // SUSY MC only has genEvtInfo
            // nominal weight is 1
            // --> weight should be weights()[1] = weights()[10] = typically 0.0003*
            single_float_["evt_Gen_Weight"] = genEvtInfo->weight();

            // This weight is used to normalize pdf/scale etc weights
            double genNomWeight = genEvtInfo->weights()[1];

            // Print factors in an event
            //std::cout<<"GEN: weight() = "<<genEvtInfo->weight()<<std::endl;
            //for (size_t i=0; i<genEvtInfo->weights().size(); ++i)
            //  std::cout<<"GEN: weights() - index: "<<i<<" wgt = "<<genEvtInfo->weights()[i]<<" wgtNorm = "<<(genEvtInfo->weights()[i]/genNomWeight)<<std::endl;

            // Renormalization/Factorization scale weights
            // These are the first 9 weights for all generators
            // mu_R and mu_F are varied independently (by a factor of 1, 2, 0.5) - check GEN header
            // [1] is the default weight (no variation)
            // --> I skip saving it (one can just use 1)
            // --> Also do not save unphysical combinations as recommended
            //    (mu_R = 0.5, mu_F = 2 and mu_R = 2, mu_F = 0.5)
            // Save only ids: 2,3,4,5,7,9
            if (genEvtInfo->weights().size()>=10) for (size_t i=1; i<10; ++i) if (i!=1&&i!=6&&i!=8)
                vector_float_["scale_Weights"].push_back(genEvtInfo->weights()[i]/genNomWeight);

            // PDF weights
            // Usually a set of 100 weights (excluding default)
            // Only default PDF variation is saved, but if needed
            // likely others are available depending on the generator
            // index of first weight varies, beware!
            // Sometimes first weight is default=1 weight (has to skip!)
            // Additional info: MC2Hessian conversion will soon be provided
            // NNPDF30_lo_as_0130 (SUSY signals has +1 weight weights()[0])
            size_t first = 11;
            if (genEvtInfo->weights().size()>=first+100) for (size_t i=first; i<first+100; ++i)
                vector_float_["pdf_Weights"].push_back(genEvtInfo->weights()[i]/genNomWeight);

            // Alpha_s weights (only given for NLO!)
        }
    }

    // ----------------------------
    // - Trigger/Filter decisions -
    // ----------------------------
    for ( auto filter : filters_ ) single_int_[filter.first]                              /* Flag_* */
        = h_floats_["filter_bits"]->at(filter.second);
    for ( auto trigger : triggers_ ) {
        single_int_[trigger.first] = h_floats_["trigger_bits"]->at(trigger.second);         /* HLT_* */
        single_int_[trigger.first+"_prescale"]                                              /* HLT_*_prescale */
            = h_ints_["trigger_prescales"]->at(trigger.second);
    }

    // ----------------------------
    // -        Vertices          -
    // ----------------------------

    single_int_["evt_NGoodVtx"] = 0;
    for (int iVtx=0; iVtx<*h_int_["vtx_npv"]; ++iVtx)
        if (h_ints_["vtx_ndof"]->at(iVtx)>=4&&fabs(h_floats_["vtx_z"]->at(iVtx))<24&&fabs(h_floats_["vtx_rho"]->at(iVtx)<2))
            ++single_int_["evt_NGoodVtx"];                                                    /* evt_NGoodVtx */

    // ---------------------
    // - Gen Particle Info -
    // ---------------------

    // Make a list of Generator level objects and save them to vectors
    vector_int_["gen_ID"].clear();
    vector_int_["gen_Status"].clear();
    vector_int_["gen_Mom0ID"].clear();
    vector_int_["gen_Mom0Status"].clear();
    vector_int_["gen_Mom1ID"].clear();
    vector_int_["gen_Mom1Status"].clear();
    vector_int_["gen_Dau0ID"].clear();
    vector_int_["gen_Dau0Status"].clear();
    vector_int_["gen_Dau1ID"].clear();
    vector_int_["gen_Dau1Status"].clear();
    vector_float_["gen_Pt"].clear();
    vector_float_["gen_Phi"].clear();
    vector_float_["gen_Eta"].clear();
    vector_float_["gen_E"].clear();
    vector_float_["gen_Charge"].clear();
    vector_float_["gen_Mass"].clear();

    std::vector<TLorentzVector> gen_top;
    std::vector<size_t > gen_top_index;
    std::vector<TLorentzVector> gen_W_from_top;
    std::vector<TLorentzVector> gen_q_from_top;
    std::vector<bool>           gen_b_from_top;
    std::vector<TLorentzVector> gen_lep_from_W;
    std::vector<TLorentzVector> gen_neu_from_W;
    std::vector<int> gen_top_ID;
    std::vector<int> gen_W_from_top_ID;
    std::vector<int> gen_lep_from_W_ID;
    std::vector<int> gen_neu_from_W_ID;

    bool good_W_matches = true;
    std::vector<TLorentzVector> gen_top_matched_q;
    std::vector<bool>           gen_top_matched_b;
    std::vector<TLorentzVector> gen_top_matched_W;
    std::vector<int> gen_top_matched_W_ID;

    std::vector<int> W_type;
    std::vector<TLorentzVector> gen_top_matched_W_matched_lep;
    std::vector<TLorentzVector> gen_top_matched_W_matched_neu;

    std::map<size_t, size_t > jet_gentop_index;
    std::map<size_t, size_t > ele_genlep_index;
    std::map<size_t, size_t > mu_genlep_index;

    size_t nele = h_floats_["ele_Pt"]->size();
    size_t nmu =  h_floats_["mu_Pt"]->size();
    size_t njet_AK4 = h_floats_["AK4_Pt"]->size();
    size_t njet_AK8 = h_floats_["AK8_Pt"]->size();
    size_t njet_AK8Sub = h_floats_["AK8Sub_Pt"]->size();
    bool useGenParticles = true;

    if (!isData_) {
        // Using GenParticles
        edm::Handle<reco::GenParticleCollection> genParticles;
        iEvent.getByLabel(edm::InputTag("prunedGenParticles"),  genParticles);
        const bool print = false;
        if (print) for(size_t i=0; i<genParticles->size(); ++i) {
            const reco::GenParticle& p = (*genParticles)[i];
            int momId = p.numberOfMothers() ? p.mother()->pdgId() : 0;
            std::cout<<i<<" id="<<p.pdgId()<<" ("<<p.status()<<") mom="<<momId<<", daughters=";
            for(size_t j = 0, n=p.numberOfDaughters(); j<n; ++j) std::cout<<p.daughter(j)->pdgId()<<(j+1<n?", ":"\n");
        }

        size_t ngen =  h_floats_["gen_Pt"]->size();
        for (size_t i=0; i<ngen; ++i) {
            /*
               if (abs(h_floats_["gen_ID"]->at(i))==5||abs(h_floats_["gen_ID"]->at(i))==6||
               (abs(h_floats_["gen_ID"]->at(i))>=11&&abs(h_floats_["gen_ID"]->at(i))<=16)
               ||abs(h_floats_["gen_ID"]->at(i))==24||abs(h_floats_["gen_ID"]->at(i))>1000000) 
               */

            if (abs(h_floats_["gen_ID"]->at(i)) <= 6 || abs(h_floats_["gen_ID"]->at(i)) == 11 ||
                    abs(h_floats_["gen_ID"]->at(i)) ==13 || abs(h_floats_["gen_ID"]->at(i)) == 15 ||
                    (abs(h_floats_["gen_ID"]->at(i)) == 21 && h_floats_["gen_Status"]->at(i) == 21) ||
                    abs(h_floats_["gen_ID"]->at(i)) == 22 ||
                    abs(h_floats_["gen_ID"]->at(i)) ==23 || abs(h_floats_["gen_ID"]->at(i)) == 2212 ||
                    abs(h_floats_["gen_ID"]->at(i))>1000000){

                vector_int_["gen_ID"].push_back(h_floats_["gen_ID"]->at(i));                                /* gen_ID  */
                vector_int_["gen_Status"].push_back(h_floats_["gen_Status"]->at(i));			    /* gen_Status */
                vector_int_["gen_Mom0ID"].push_back(h_floats_["gen_Mom0ID"]->at(i));			    /* gen_Mom0ID */
                vector_int_["gen_Mom0Status"].push_back(h_floats_["gen_Mom0Status"]->at(i));		    /* gen_Mom0Status */
                vector_int_["gen_Mom1ID"].push_back(h_floats_["gen_Mom1ID"]->at(i));			    /* gen_Mom1ID */
                vector_int_["gen_Mom1Status"].push_back(h_floats_["gen_Mom1Status"]->at(i));		    /* gen_Mom1Status */
                vector_int_["gen_Dau0ID"].push_back(h_floats_["gen_Dau0ID"]->at(i));			    /* gen_Dau0ID */
                vector_int_["gen_Dau0Status"].push_back(h_floats_["gen_Dau0Status"]->at(i));		    /* gen_Dau0Status */
                vector_int_["gen_Dau1ID"].push_back(h_floats_["gen_Dau1ID"]->at(i));			    /* gen_Dau1ID */
                vector_int_["gen_Dau1Status"].push_back(h_floats_["gen_Dau1Status"]->at(i));		    /* gen_Dau1Status */
                vector_float_["gen_Pt"].push_back(h_floats_["gen_Pt"]->at(i));				    /* gen_Pt */
                vector_float_["gen_Eta"].push_back(h_floats_["gen_Eta"]->at(i));			    /* gen_Eta */
                vector_float_["gen_Phi"].push_back(h_floats_["gen_Phi"]->at(i));			    /* gen_Phi */
                vector_float_["gen_E"].push_back(h_floats_["gen_E"]->at(i));				    /* gen_E */
                vector_float_["gen_Charge"].push_back(h_floats_["gen_Charge"]->at(i));			    /* gen_Charge */
                TLorentzVector genp; genp.SetPtEtaPhiE(h_floats_["gen_Pt"]->at(i), h_floats_["gen_Eta"]->at(i),
                        h_floats_["gen_Phi"]->at(i), h_floats_["gen_E"]->at(i));
                vector_float_["gen_Mass"].push_back(genp.M());				                    /* gen_Mass */
            }

        }




    } // End !isData

    // ---------------------
    // -      Jets         -
    // ---------------------

    // Read jet correction parameters from DB
    edm::ESHandle<JetCorrectorParametersCollection> JetCorrParColl_AK4, JetCorrParColl_AK8;
    iSetup.get<JetCorrectionsRecord>().get(TString(AK4Jets_label_).Contains("Puppi") ? "AK4PFPuppi" : "AK4PFchs", JetCorrParColl_AK4);
    iSetup.get<JetCorrectionsRecord>().get(TString(AK8Jets_label_).Contains("Puppi") ? "AK8PFPuppi" : "AK8PFchs", JetCorrParColl_AK8);

    // JEC Uncertainty
    // NB: JEC/JER variables are already there in the new versions of B2G ntuples
    /*
       JetCorrectionUncertainty jecUnc_AK4((*JetCorrParColl_AK4)["Uncertainty"]);
       JetCorrectionUncertainty jecUnc_AK8((*JetCorrParColl_AK8)["Uncertainty"]);

    // JER
    // Twiki: https://twiki.cern.ch/twiki/bin/view/CMSPublic/WorkBookJetEnergyResolution#Scale_factors
    // Recipe taken from: https://github.com/blinkseb/cmssw/blob/jer_fix_76x/JetMETCorrections/Modules/plugins/JetResolutionDemo.cc
    JME::JetParameters jetParam;
    JME::JetResolution resolution_AK4 = JME::JetResolution(JER_location_+"_PtResolution_AK4PFchs.txt");
    JME::JetResolution resolution_AK8 = JME::JetResolution(JER_location_+"_PtResolution_AK8PFchs.txt");
    JME::JetResolutionScaleFactor res_sf_AK4 = JME::JetResolutionScaleFactor(JER_location_+"_SF_AK4PFchs.txt");
    JME::JetResolutionScaleFactor res_sf_AK8 = JME::JetResolutionScaleFactor(JER_location_+"_SF_AK8PFchs.txt");

    vector_float_[AK4Jets_prefix_+"_jecUncertainty"].assign(njet_AK4,-9999);
    vector_float_[AK4Jets_prefix_+"_PtResolution"].assign(njet_AK4,-9999);
    vector_float_[AK4Jets_prefix_+"_JERSF"].assign(njet_AK4,-9999);
    vector_float_[AK4Jets_prefix_+"_JERSFDown"].assign(njet_AK4,-9999);
    vector_float_[AK4Jets_prefix_+"_JERSFUp"].assign(njet_AK4,-9999);
    for (size_t iJet=0; iJet<njet_AK4; ++iJet) {
    jecUnc_AK4.setJetPt(h_floats_["AK4_Pt"]->at(iJet));
    jecUnc_AK4.setJetEta(h_floats_["AK4_Eta"]->at(iJet));
    jetParam.setJetPt(h_floats_["AK4_Pt"]->at(iJet)).setJetEta(h_floats_["AK4_Eta"]->at(iJet)).setRho(*h_double_["evt_rho"]);
    vector_float_[AK4Jets_prefix_+"_jecUncertainty"][iJet] = jecUnc_AK4.getUncertainty(true);
    vector_float_[AK4Jets_prefix_+"_PtResolution"][iJet] = resolution_AK4.getResolution(jetParam);
    vector_float_[AK4Jets_prefix_+"_JERSF"][iJet]     = res_sf_AK4.getScaleFactor(jetParam);
    vector_float_[AK4Jets_prefix_+"_JERSFDown"][iJet] = res_sf_AK4.getScaleFactor(jetParam, Variation::UP);
    vector_float_[AK4Jets_prefix_+"_JERSFUp"][iJet]   = res_sf_AK4.getScaleFactor(jetParam, Variation::DOWN);
    }
    vector_float_[AK8Jets_prefix_+"_jecUncertainty"].assign(njet_AK8,-9999);
    vector_float_[AK8Jets_prefix_+"_PtResolution"].assign(njet_AK8,-9999);
    vector_float_[AK8Jets_prefix_+"_JERSF"].assign(njet_AK8,-9999);
    vector_float_[AK8Jets_prefix_+"_JERSFDown"].assign(njet_AK8,-9999);
    vector_float_[AK8Jets_prefix_+"_JERSFUp"].assign(njet_AK8,-9999);
    for (size_t iJet=0; iJet<njet_AK8; ++iJet) {
    jecUnc_AK8.setJetPt(h_floats_["AK8_Pt"]->at(iJet));
    jecUnc_AK8.setJetEta(h_floats_["AK8_Eta"]->at(iJet));
    jetParam.setJetPt(h_floats_["AK8_Pt"]->at(iJet)).setJetEta(h_floats_["AK8_Eta"]->at(iJet)).setRho(*h_double_["evt_rho"]);
    vector_float_[AK8Jets_prefix_+"_jecUncertainty"][iJet] = jecUnc_AK8.getUncertainty(true);
    vector_float_[AK8Jets_prefix_+"_PtResolution"][iJet] = resolution_AK8.getResolution(jetParam);
    vector_float_[AK8Jets_prefix_+"_JERSF"][iJet]     = res_sf_AK8.getScaleFactor(jetParam);
    vector_float_[AK8Jets_prefix_+"_JERSFDown"][iJet] = res_sf_AK8.getScaleFactor(jetParam, Variation::UP);
    vector_float_[AK8Jets_prefix_+"_JERSFUp"][iJet]   = res_sf_AK8.getScaleFactor(jetParam, Variation::DOWN);
    }
    */

    // Jet Correctors (Recalculate after lepton cleaning, or do on-the-fly correction)
    std::vector<JetCorrectorParameters> AK4_vPar;
    AK4_vPar.push_back((*JetCorrParColl_AK4)["L1FastJet"]);
    AK4_vPar.push_back((*JetCorrParColl_AK4)["L2Relative"]);
    AK4_vPar.push_back((*JetCorrParColl_AK4)["L3Absolute"]);
    if (isData_) AK4_vPar.push_back((*JetCorrParColl_AK4)["L2L3Residual"]);
    FactorizedJetCorrector AK4_JetCorrector(AK4_vPar);
    std::vector<JetCorrectorParameters> AK8_vPar;
    AK8_vPar.push_back((*JetCorrParColl_AK8)["L1FastJet"]);
    AK8_vPar.push_back((*JetCorrParColl_AK8)["L2Relative"]);
    AK8_vPar.push_back((*JetCorrParColl_AK8)["L3Absolute"]);
    if (isData_) AK8_vPar.push_back((*JetCorrParColl_AK8)["L2L3Residual"]);
    FactorizedJetCorrector AK8_JetCorrector(AK8_vPar);

    // GEN infos
    vector_int_[AK8Jets_prefix_+"_HasNearGenTop"].assign(njet_AK8,-9999);
    vector_int_[AK8Jets_prefix_+"_NearGenTopIsHadronic"].assign(njet_AK8,-9999);
    vector_int_[AK8Jets_prefix_+"_NearGenWIsHadronic"].assign(njet_AK8,-9999);
    vector_int_[AK8Jets_prefix_+"_NearGenWToENu"].assign(njet_AK8,-9999);
    vector_int_[AK8Jets_prefix_+"_NearGenWToMuNu"].assign(njet_AK8,-9999);
    vector_int_[AK8Jets_prefix_+"_NearGenWToTauNu"].assign(njet_AK8,-9999);
    vector_float_[AK8Jets_prefix_+"_DRNearGenTop"].assign(njet_AK8,9999);
    vector_float_[AK8Jets_prefix_+"_DRNearGenWFromTop"].assign(njet_AK8,9999);
    vector_float_[AK8Jets_prefix_+"_DRNearGenBFromTop"].assign(njet_AK8,9999);
    vector_float_[AK8Jets_prefix_+"_DRNearGenLepFromSLTop"].assign(njet_AK8,9999);
    vector_float_[AK8Jets_prefix_+"_DRNearGenNuFromSLTop"].assign(njet_AK8,9999);
    vector_float_[AK8Jets_prefix_+"_PtNearGenTop"].assign(njet_AK8,-9999);
    vector_float_[AK8Jets_prefix_+"_PtNearGenBFromTop"].assign(njet_AK8,-9999);
    vector_float_[AK8Jets_prefix_+"_PtNearGenWFromTop"].assign(njet_AK8,-9999);
    vector_float_[AK8Jets_prefix_+"_PtNearGenLepFromSLTop"].assign(njet_AK8,-9999);
    vector_float_[AK8Jets_prefix_+"_PtNearGenNuFromSLTop"].assign(njet_AK8,-9999);


    /*
       Jet ID
    https://twiki.cern.ch/twiki/bin/view/CMS/JetID13TeVRun2016?rev=4

    For |eta|<=2.7 Apply
    looseJetID = (NHF<0.99 && NEMF<0.99 && NumConst>1) && ((abs(eta)<=2.4 && CHF>0 && CHM>0 && CEMF<0.99) || abs(eta)>2.4) && abs(eta)<=2.7
    tightJetID = (NHF<0.90 && NEMF<0.90 && NumConst>1) && ((abs(eta)<=2.4 && CHF>0 && CHM>0 && CEMF<0.99) || abs(eta)>2.4) && abs(eta)<=2.7
    tightLepVetoJetID = (NHF<0.90 && NEMF<0.90 && NumConst>1 && MUF<0.8) && ((abs(eta)<=2.4 && CHF>0 && CHM>0 && CEMF<0.90) || abs(eta)>2.4) && abs(eta)<=2.7

    For 2.7<|eta|<= 3.0 Apply
    looseJetID = (NHF<0.98 && NEMF>0.01 && NumNeutralParticle>2 && abs(eta)>2.7 && abs(eta)<=3.0 )
    tightJetID = (NHF<0.98 && NEMF>0.01 && NumNeutralParticle>2 && abs(eta)>2.7 && abs(eta)<=3.0 )

    For |eta|> 3.0 Apply
    looseJetID = (NEMF<0.90 && NumNeutralParticle>10 && abs(eta)>3.0 )
    tightJetID = (NEMF<0.90 && NumNeutralParticle>10 && abs(eta)>3.0 ) 
*/

    vector_int_[AK4Jets_prefix_+"_looseJetID"].clear();
    vector_int_[AK4Jets_prefix_+"_tightJetID"].clear();
    vector_int_[AK4Jets_prefix_+"_tightLepVetoJetID"].clear();
    for (size_t i=0; i<njet_AK4; ++i) {
        float eta  = h_floats_["AK4_Eta"]->at(i);
        float NHF  = h_floats_["AK4_neutralHadronEnergyFrac"]->at(i);
        float NEMF = h_floats_["AK4_neutralEmEnergyFrac"]->at(i);
        float MUF = h_floats_["AK4_MuonEnergy"]->at(i) / h_floats_["AK4_E"]->at(i);
        float CHF  = h_floats_["AK4_chargedHadronEnergyFrac"]->at(i);
        float CEMF = h_floats_["AK4_chargedEmEnergyFrac"]->at(i);
        int CHM  = h_floats_["AK4_chargedMultiplicity"]->at(i);
        int NumNeutralParticle   = h_floats_["AK4_neutralMultiplicity"]->at(i);
        int NumConst = CHM + NumNeutralParticle;
        bool looseJetID = 0, tightJetID = 0, tightLepVetoJetID = 0;
        if (std::abs(eta)<=2.7) {
            looseJetID = (NHF<0.99 && NEMF<0.99 && NumConst>1) && ((std::abs(eta)<=2.4 && CHF>0 && CHM>0 && CEMF<0.99) || std::abs(eta)>2.4);
            tightJetID = (NHF<0.90 && NEMF<0.90 && NumConst>1) && ((std::abs(eta)<=2.4 && CHF>0 && CHM>0 && CEMF<0.99) || std::abs(eta)>2.4);
            tightLepVetoJetID = (NHF<0.90 && NEMF<0.90 && NumConst>1 && MUF<0.8) && ((std::abs(eta)<=2.4 && CHF>0 && CHM>0 && CEMF<0.90) || std::abs(eta)>2.4);
        } else if (std::abs(eta)>2.7&&std::abs(eta)<=3.0) {
            looseJetID = (NHF<0.98 && NEMF>0.01 && NumNeutralParticle>2);
            tightJetID = (NHF<0.98 && NEMF>0.01 && NumNeutralParticle>2);
        } else {
            looseJetID = (NEMF<0.90 && NumNeutralParticle>10);
            tightJetID = (NEMF<0.90 && NumNeutralParticle>10);
        }
        vector_int_[AK4Jets_prefix_+"_looseJetID"].push_back(looseJetID);                       /* jetAK4_looseJetID  */
        vector_int_[AK4Jets_prefix_+"_tightJetID"].push_back(tightJetID);                       /* jetAK4_tightJetID  */
        vector_int_[AK4Jets_prefix_+"_tightLepVetoJetID"].push_back(tightLepVetoJetID);         /* jetAK4_tightLepVetoJetID  */
    }

    vector_int_[AK8Jets_prefix_+"_looseJetID"].clear();
    vector_int_[AK8Jets_prefix_+"_tightJetID"].clear();
    vector_int_[AK8Jets_prefix_+"_tightLepVetoJetID"].clear();
    vector_float_[AK8Jets_prefix_+"_maxSubjetCSVv2"].clear();
    vector_float_[AK8Jets_prefix_+"_maxSubjetCMVAv2"].clear();
    for (size_t i=0; i<njet_AK8; ++i) {
        float eta  = h_floats_["AK8_Eta"]->at(i);
        float NHF  = h_floats_["AK8_neutralHadronEnergyFrac"]->at(i);
        float NEMF = h_floats_["AK8_neutralEmEnergyFrac"]->at(i);
        float MUF = h_floats_["AK8_MuonEnergy"]->at(i) / h_floats_["AK8_E"]->at(i);
        float CHF  = h_floats_["AK8_chargedHadronEnergyFrac"]->at(i);
        float CEMF = h_floats_["AK8_chargedEmEnergyFrac"]->at(i);
        int CHM  = h_floats_["AK8_chargedMultiplicity"]->at(i);
        int NumNeutralParticle   = h_floats_["AK8_neutralMultiplicity"]->at(i);
        int NumConst = CHM + NumNeutralParticle;
        bool looseJetID = 0, tightJetID = 0, tightLepVetoJetID = 0;
        if (std::abs(eta)<=2.7) {
            looseJetID = (NHF<0.99 && NEMF<0.99 && NumConst>1) && ((std::abs(eta)<=2.4 && CHF>0 && CHM>0 && CEMF<0.99) || std::abs(eta)>2.4);
            tightJetID = (NHF<0.90 && NEMF<0.90 && NumConst>1) && ((std::abs(eta)<=2.4 && CHF>0 && CHM>0 && CEMF<0.99) || std::abs(eta)>2.4);
            tightLepVetoJetID = (NHF<0.90 && NEMF<0.90 && NumConst>1 && MUF<0.8) && ((std::abs(eta)<=2.4 && CHF>0 && CHM>0 && CEMF<0.90) || std::abs(eta)>2.4);
        } else if (std::abs(eta)>2.7&&std::abs(eta)<=3.0) {
            looseJetID = (NHF<0.98 && NEMF>0.01 && NumNeutralParticle>2);
            tightJetID = (NHF<0.98 && NEMF>0.01 && NumNeutralParticle>2);
        } else {
            looseJetID = (NEMF<0.90 && NumNeutralParticle>10);
            tightJetID = (NEMF<0.90 && NumNeutralParticle>10);
        }
        vector_int_[AK8Jets_prefix_+"_looseJetID"].push_back(looseJetID);                       /* jetAK8_looseJetID  */
        vector_int_[AK8Jets_prefix_+"_tightJetID"].push_back(tightJetID);                       /* jetAK8_tightJetID  */
        vector_int_[AK8Jets_prefix_+"_tightLepVetoJetID"].push_back(tightLepVetoJetID);         /* jetAK8_tightLepVetoJetID  */

        // Subjet btag info
        int i_sj0 = h_floats_["AK8_vSubjetIndex0"]->at(i), i_sj1 = h_floats_["AK8_vSubjetIndex1"]->at(i);
        float maxCSVv2 = -9999, maxCMVAv2 = -9999;
        if (i_sj0 != -1) if (h_floats_["AK8Sub_CSVv2"]->at(i_sj0) > maxCSVv2) maxCSVv2 = h_floats_["AK8Sub_CSVv2"]->at(i_sj0);
        if (i_sj1 != -1) if (h_floats_["AK8Sub_CSVv2"]->at(i_sj1) > maxCSVv2) maxCSVv2 = h_floats_["AK8Sub_CSVv2"]->at(i_sj1);
        if (i_sj0 != -1) if (h_floats_["AK8Sub_CMVAv2"]->at(i_sj0) > maxCMVAv2) maxCMVAv2 = h_floats_["AK8Sub_CMVAv2"]->at(i_sj0);
        if (i_sj1 != -1) if (h_floats_["AK8Sub_CMVAv2"]->at(i_sj1) > maxCMVAv2) maxCMVAv2 = h_floats_["AK8Sub_CMVAv2"]->at(i_sj1);
        vector_float_[AK8Jets_prefix_+"_maxSubjetCSVv2"].push_back(maxCSVv2);
        vector_float_[AK8Jets_prefix_+"_maxSubjetCMVAv2"].push_back(maxCMVAv2);
    }

    vector_int_[AK8Subjets_prefix_+"_looseJetID"].clear();
    vector_int_[AK8Subjets_prefix_+"_tightJetID"].clear();
    vector_int_[AK8Subjets_prefix_+"_tightLepVetoJetID"].clear();
    for (size_t i=0; i<njet_AK8Sub; ++i) {
        float eta  = h_floats_["AK8Sub_Eta"]->at(i);
        float NHF  = h_floats_["AK8Sub_neutralHadronEnergyFrac"]->at(i);
        float NEMF = h_floats_["AK8Sub_neutralEmEnergyFrac"]->at(i);
        float MUF = h_floats_["AK8Sub_MuonEnergy"]->at(i) / h_floats_["AK8Sub_E"]->at(i);
        float CHF  = h_floats_["AK8Sub_chargedHadronEnergyFrac"]->at(i);
        float CEMF = h_floats_["AK8Sub_chargedEmEnergyFrac"]->at(i);
        int CHM  = h_floats_["AK8Sub_chargedMultiplicity"]->at(i);
        int NumNeutralParticle   = h_floats_["AK8Sub_neutralMultiplicity"]->at(i);
        int NumConst = CHM + NumNeutralParticle;
        bool looseJetID = 0, tightJetID = 0, tightLepVetoJetID = 0;
        if (std::abs(eta)<=2.7) {
            looseJetID = (NHF<0.99 && NEMF<0.99 && NumConst>1) && ((std::abs(eta)<=2.4 && CHF>0 && CHM>0 && CEMF<0.99) || std::abs(eta)>2.4);
            tightJetID = (NHF<0.90 && NEMF<0.90 && NumConst>1) && ((std::abs(eta)<=2.4 && CHF>0 && CHM>0 && CEMF<0.99) || std::abs(eta)>2.4);
            tightLepVetoJetID = (NHF<0.90 && NEMF<0.90 && NumConst>1 && MUF<0.8) && ((std::abs(eta)<=2.4 && CHF>0 && CHM>0 && CEMF<0.90) || std::abs(eta)>2.4);
        } else if (std::abs(eta)>2.7&&std::abs(eta)<=3.0) {
            looseJetID = (NHF<0.98 && NEMF>0.01 && NumNeutralParticle>2);
            tightJetID = (NHF<0.98 && NEMF>0.01 && NumNeutralParticle>2);
        } else {
            looseJetID = (NEMF<0.90 && NumNeutralParticle>10);
            tightJetID = (NEMF<0.90 && NumNeutralParticle>10);
        }
        vector_int_[AK8Subjets_prefix_+"_looseJetID"].push_back(looseJetID);                       /* subjetAK8_looseJetID  */
        vector_int_[AK8Subjets_prefix_+"_tightJetID"].push_back(tightJetID);                       /* subjetAK8_tightJetID  */
        vector_int_[AK8Subjets_prefix_+"_tightLepVetoJetID"].push_back(tightLepVetoJetID);         /* subjetAK8_tightLepVetoJetID  */
    }



    // ---------------------
    // -        MET        -
    // ---------------------

    // Uncertainties
    //edm::Handle<pat::METCollection> mets_MuCleanOnly;
    //iEvent.getByLabel(edm::InputTag("slimmedMETs"), mets_MuCleanOnly);
    //const pat::MET &met_MuCleanOnly = mets_MuCleanOnly->front();
    //edm::Handle<pat::METCollection> mets;
    //if (isData_) iEvent.getByLabel(edm::InputTag("slimmedMETsMuEGClean","","b2gEDMNtuples"), mets);
    //else         iEvent.getByLabel(edm::InputTag("slimmedMETsMuClean",  "","b2gEDMNtuples"), mets);
    //const pat::MET &met = mets->front();
    //edm::Handle<pat::METCollection> puppimets;
    //iEvent.getByLabel(edm::InputTag("slimmedMETsPuppi"), puppimets);
    //const pat::MET &puppimet = puppimets->front();
    //
    //vector_float_["metsyst_MuCleanOnly_Pt"].clear();
    //vector_float_["metsyst_MuCleanOnly_Phi"].clear();
    //vector_float_["metsyst_Pt"].clear();
    //vector_float_["metsyst_Phi"].clear();
    //vector_float_["puppimetsyst_Pt"].clear();
    //vector_float_["puppimetsyst_Phi"].clear();
    //for (int shift=0; shift<pat::MET::METUncertainty::METUncertaintySize; ++shift)
    //  if (shift != pat::MET::METUncertainty::NoShift) {
    //    float met_MuCleanOnly_shiftedPt  = met_MuCleanOnly.shiftedPt ((pat::MET::METUncertainty)shift, pat::MET::METCorrectionLevel::Type1);
    //    float met_MuCleanOnly_shiftedPhi = met_MuCleanOnly.shiftedPhi((pat::MET::METUncertainty)shift, pat::MET::METCorrectionLevel::Type1);
    //    float met_shiftedPt  = met.shiftedPt ((pat::MET::METUncertainty)shift, pat::MET::METCorrectionLevel::Type1);
    //    float met_shiftedPhi = met.shiftedPhi((pat::MET::METUncertainty)shift, pat::MET::METCorrectionLevel::Type1);
    //    float puppimet_shiftedPt  = puppimet.shiftedPt ((pat::MET::METUncertainty)shift, pat::MET::METCorrectionLevel::Type1);
    //    float puppimet_shiftedPhi = puppimet.shiftedPhi((pat::MET::METUncertainty)shift, pat::MET::METCorrectionLevel::Type1);
    //    vector_float_["metsyst_MuCleanOnly_Pt"].push_back(met_MuCleanOnly_shiftedPt);
    //    vector_float_["metsyst_MuCleanOnly_Phi"].push_back(met_MuCleanOnly_shiftedPhi);
    //    vector_float_["metsyst_Pt"].push_back(met_shiftedPt);
    //    vector_float_["metsyst_Phi"].push_back(met_shiftedPhi);
    //    vector_float_["puppimetsyst_Pt"].push_back(puppimet_shiftedPt);
    //    vector_float_["puppimetsyst_Phi"].push_back(puppimet_shiftedPhi);
    //  }


    // ---------------------
    // - Lepton Selection  -
    // ---------------------


    // ---------------------
    // -- Isolated tracks --
    // ---------------------

    edm::Handle<pat::PackedCandidateCollection> packedPFCands;
    iEvent.getByLabel(edm::InputTag("packedPFCandidates"), packedPFCands);
    const pat::PackedCandidateCollection& pfCands = *packedPFCands.product();

    single_int_["evt_NIsoTrk"] = 0;
    for (size_t i=0, n=pfCands.size(); i<n; ++i) {
        if (pfCands[i].charge()==0) continue;
        if (pfCands[i].pt()<5) continue;
        if (std::abs(pfCands[i].dz())>=0.1) continue;

        // Calculate track isolation
        float iso_trk = 0;
        for (size_t j=0; j<n; ++j) {
            if (j==i) continue;
            if (pfCands[j].charge()==0) continue;
            float dR = reco::deltaR(pfCands[i].eta(), pfCands[i].phi(), pfCands[j].eta(), pfCands[j].phi());
            if (dR >= 0.3) continue;
            if (std::abs(pfCands[j].dz())>=0.1) continue;
            iso_trk += pfCands[j].pt();
        }
        iso_trk /= pfCands[i].pt();
        if (std::abs(pfCands[i].pdgId())==11||std::abs(pfCands[i].pdgId())==13) {
            if (iso_trk<0.2) ++single_int_["evt_NIsoTrk"];
        } else {
            if (iso_trk<0.1 && pfCands[i].pt()>=10) ++single_int_["evt_NIsoTrk"];
        }
    }

    // ---------------------
    // -- Razor variables --
    // ---------------------

    // Jet selection for jet combiner
    std::vector<TLorentzVector> jets_AK4, jets_AK4_smear;
    for (size_t i=0; i<njet_AK4; ++i) {
        // Cut in MINIAOD: pt>20
        // 2016/10/14: JetID cut added, pt lowered from 40 to 30, |eta| lowered to 2.4
        if (std::abs(h_floats_["AK4_Eta"]->at(i)) < 2.4 && vector_int_[AK4Jets_prefix_+"_looseJetID"][i]) {
            TLorentzVector jl;
            jl.SetPtEtaPhiE(h_floats_["AK4_Pt"]->at(i), h_floats_["AK4_Eta"]->at(i),
                    h_floats_["AK4_Phi"]->at(i), h_floats_["AK4_E"]->at(i));
            float pt = h_floats_["AK4_Pt"]->at(i), smeared_pt = h_floats_["AK4_SmearedPt"]->at(i);
            if (pt >= 30) jets_AK4.push_back(jl);
            // Use also JER Smeared jets in order to propagate the uncertainty
            if (smeared_pt >= 30) jets_AK4_smear.push_back(jl*(smeared_pt/pt));
        }
    }

    TVector3 metl;
    metl.SetPtEtaPhi(h_floats_["met_Pt"]->at(0), 0, h_floats_["met_Phi"]->at(0));


}


// checks if a particle has a special mother. Treats anti-particles as particles
bool B2GEdmExtraVarProducer::hasAncestor_(int index, const lhef::HEPEUP& info, int searchId) {
    if (index < 2 || index > info.NUP) return false;
    else if (std::abs(info.IDUP[index]) == searchId) return true;
    else {
        auto mothers = info.MOTHUP[index];
        return
            (index != mothers.first-1 && hasAncestor_(mothers.first-1, info, searchId)) ||
            (index != mothers.second-1 && hasAncestor_(mothers.second-1, info, searchId));
    }
}

#include "FWCore/Framework/interface/MakerMacros.h"
DEFINE_FWK_MODULE(B2GEdmExtraVarProducer);
