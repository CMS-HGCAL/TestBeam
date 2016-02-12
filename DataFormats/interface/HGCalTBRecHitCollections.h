#ifndef DATAFORMATS_HGCRECHIT_HGCRECHITCOLLECTION_H
#define DATAFORMATS_HGCRECHIT_HGCRECHITCOLLECTION_H

#include "DataFormats/Common/interface/SortedCollection.h"
#include "HGCal/DataFormats/interface/HGCalTBRecHit.h"
#include "DataFormats/Common/interface/Ref.h"
#include "DataFormats/Common/interface/RefVector.h"


typedef edm::SortedCollection<HGCalTBRecHit> HGCalTBRecHitCollection;
typedef edm::Ref<HGCalTBRecHitCollection> HGCalTBRecHitRef;
typedef edm::RefVector<HGCalTBRecHitCollection> HGCalTBRecHitRefs;
typedef edm::RefProd<HGCalTBRecHitCollection> HGCalTBRecHitsRef;

#endif
