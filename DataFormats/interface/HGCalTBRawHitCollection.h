#ifndef DATAFORMATS_HGCALTBRAWHITCOLLECTION_H
#define DATAFORMATS_HGCALTBRAWHITCOLLECTION_H

#include "DataFormats/Common/interface/EDCollection.h"
#include "HGCal/DataFormats/interface/HGCalTBRawHit.h"
#include "DataFormats/Common/interface/Ref.h"
#include "DataFormats/Common/interface/RefVector.h"

typedef edm::EDCollection<HGCalTBRawHit> HGCalTBRawHitCollection;
typedef edm::Ref<HGCalTBRawHitCollection> HGCalTBRawHitRef;
typedef edm::RefVector<HGCalTBRawHitCollection> HGCalTBRawHitRefs;
typedef edm::RefProd<HGCalTBRawHitCollection> HGCalTBRawHitRefProd;


#endif
