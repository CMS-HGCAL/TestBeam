#ifndef DATAFORMATS_HGCCALOTRACK_HGCCALOTRACKCOLLECTION_H
#define DATAFORMATS_HGCCALOTRACK_HGCCALOTRACKCOLLECTION_H

#include "DataFormats/Common/interface/EDCollection.h"
#include "HGCal/DataFormats/interface/HGCalTBCaloTrack.h"
#include "DataFormats/Common/interface/Ref.h"
#include "DataFormats/Common/interface/RefVector.h"

// ~ copy paste from "HGCal/DataFormats/interface/HGCalTBRecHitCollection.h"

namespace reco
{
  typedef edm::EDCollection<HGCalTBCaloTrack> HGCalTBCaloTrackCollection;
  typedef edm::Ref<HGCalTBCaloTrackCollection> HGCalTBCaloTrackRef;
  typedef edm::RefVector<HGCalTBCaloTrackCollection> HGCalTBCaloTrackRefs;
  typedef edm::RefProd<HGCalTBCaloTrackCollection> HGCalTBCaloTracksRef;
}
#endif
