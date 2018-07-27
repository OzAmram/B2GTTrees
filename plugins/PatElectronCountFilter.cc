#include "FWCore/Framework/interface/MakerMacros.h"
#include "DataFormats/PatCandidates/interface/Electron.h"
#include "CommonTools/UtilAlgos/interface/ObjectCountFilter.h"
#include "CommonTools/UtilAlgos/interface/StringCutObjectSelector.h"

typedef ObjectCountFilter<pat::ElectronCollection, StringCutObjectSelector<pat::Electron> >::type PatElectronCountFilter;

DEFINE_FWK_MODULE( PatElectronCountFilter );
