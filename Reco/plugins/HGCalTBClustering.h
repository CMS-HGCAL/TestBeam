#ifndef RECO_HGCALTBCLUSTERING
#define RECO_HGCALTBCLUSTERING 1

#include "HGCal/DataFormats/interface/HGCalTBClusterCollection.h"
#include "HGCal/DataFormats/interface/HGCalTBRecHitCollections.h"
#include "HGCal/DataFormats/interface/HGCalTBDetId.h"
#include "iostream"
#include "vector"

struct HGCalTBClusteringParameterSetting{
  int maxTransverse;
  bool useDistanceInsteadCellID;
  bool useEnergyToWeightPosition;
  HGCalTBClusteringParameterSetting() : maxTransverse(1) ,
    useDistanceInsteadCellID(false),
    useEnergyToWeightPosition(false)
  {;}
};

class HGCalTBClustering
{
 public:
  HGCalTBClustering(){;}
  ~HGCalTBClustering(){;}
  void Run(HGCalTBRecHitCollection hitcol, std::vector<reco::HGCalTBCluster> &outClusterColl);
  void RunSimple(HGCalTBRecHitCollection hitcol, reco::HGCalTBCluster &cluster);
  inline void SetHGCalTBClusteringParameterSetting(HGCalTBClusteringParameterSetting val){settings=val;}

 private:
  void BuildCluster(HGCalTBRecHitCollection hitcol,
		    std::vector<HGCalTBDetId> &temp,
		    std::vector<HGCalTBDetId> &clusterDetIDs);

  //std::vector<HGCalTBDetId> _clusterDetIDList;
    
  HGCalTBClusteringParameterSetting settings;
};

#endif
