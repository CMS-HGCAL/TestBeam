#ifndef HGCAL_COMMONMODE
#define HGCAL_COMMONMODE

#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "HGCal/DataFormats/interface/HGCalTBRawHitCollection.h"
#include "HGCal/DataFormats/interface/HGCalTBDetId.h"
#include "HGCal/DataFormats/interface/HGCalTBCommonModeNoise.h"
#include "HGCal/Geometry/interface/HGCalTBCellVertices.h"
#include "HGCal/CondObjects/interface/HGCalElectronicsMap.h"
#include "HGCal/DataFormats/interface/HGCalTBElectronicsId.h"

  
class CommonMode{
 public:
  CommonMode( HGCalElectronicsMap &map, bool useMedian=true, bool cmPerChip=true );
  ~CommonMode(){}
  void Evaluate( edm::Handle<HGCalTBRawHitCollection> hits );
  std::map<int,commonModeNoise> CommonModeNoiseMap() const {return _cmMap;} 
 private:
  void EvaluateMedianPerChip( edm::Handle<HGCalTBRawHitCollection> hits );
  void EvaluateMedianPerLayer( edm::Handle<HGCalTBRawHitCollection> hits );
  void EvaluateAveragePerChip( edm::Handle<HGCalTBRawHitCollection> hits );
  void EvaluateAveragePerLayer( edm::Handle<HGCalTBRawHitCollection> hits );
  HGCalElectronicsMap _emap;
  bool _useMedian;
  bool _cmPerChip;
    
  std::map<int,commonModeNoise> _cmMap;
  
};

#endif
