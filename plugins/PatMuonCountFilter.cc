
#include "FWCore/Framework/interface/MakerMacros.h"
#include "DataFormats/PatCandidates/interface/Muon.h"
#include "CommonTools/UtilAlgos/interface/ObjectCountFilter.h"
#include "CommonTools/UtilAlgos/interface/StringCutObjectSelector.h"

typedef ObjectCountFilter<pat::MuonCollection, StringCutObjectSelector<pat::Muon> >::type PatMuonCountFilter;

DEFINE_FWK_MODULE( PatMuonCountFilter );
