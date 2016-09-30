#ifndef DATAFORMATS_HGCCLUSTER_HGCCLUSTERCOLLECTION_H
#define DATAFORMATS_HGCCLUSTER_HGCCLUSTERCOLLECTION_H

#include "DataFormats/Common/interface/EDCollection.h"
#include "HGCal/DataFormats/interface/HGCalTBCluster.h"
#include "DataFormats/Common/interface/Ref.h"
#include "DataFormats/Common/interface/RefVector.h"

// ~ copy paste from "HGCal/DataFormats/interface/HGCalTBRecHitCollection.h"

namespace reco
{
  typedef edm::EDCollection<HGCalTBCluster> HGCalTBClusterCollection;
  typedef edm::Ref<HGCalTBClusterCollection> HGCalTBClusterRef;
  typedef edm::RefVector<HGCalTBClusterCollection> HGCalTBClusterRefs;
  typedef edm::RefProd<HGCalTBClusterCollection> HGCalTBClustersRef;
}
#endif
