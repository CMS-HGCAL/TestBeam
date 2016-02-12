#include "HGCal/DataFormats/interface/HGCalTBRecHit.h"
#include "HGCal/DataFormats/interface/HGCalTBRecHitCollections.h"
#include "HGCal/DataFormats/interface/HGCalTBDataFrameContainers.h"

#include "DataFormats/Common/interface/RefProd.h" 
#include "DataFormats/Common/interface/Wrapper.h"
#include "DataFormats/Common/interface/RefToBase.h"
#include "DataFormats/Common/interface/Holder.h"
#include <vector>

namespace DataFormats_HGCTB {
  struct dictionary {
    HGCalTBRecHit _aRecHit;
    std::vector<HGCalTBRecHit> _HGCTBRHitVect;
    edm::SortedCollection<HGCalTBRecHit> _theHGCTBRsc;
    edm::Wrapper< HGCalTBRecHitCollection > _HGCTBeeRHitProd;

    SKIROC2DigiCollection _SR2DC;
    edm::Wrapper<SKIROC2DigiCollection> _theSR2DC;
    
    
  };
}

