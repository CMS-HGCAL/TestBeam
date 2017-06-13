#ifndef HGCAL_COMMONMODE
#define HGCAL_COMMONMODE

#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "HGCal/DataFormats/interface/HGCalTBRawHitCollection.h"
#include "HGCal/DataFormats/interface/HGCalTBDetId.h"
#include "HGCal/Geometry/interface/HGCalTBCellVertices.h"
#include "HGCal/CondObjects/interface/HGCalElectronicsMap.h"
#include "HGCal/DataFormats/interface/HGCalTBElectronicsId.h"

struct commonModeNoise{
  commonModeNoise(){
    for(int it=0;it<NUMBER_OF_TIME_SAMPLES; it++){
      fullHG[it] = halfHG[it] = mouseBiteHG[it] = outerHG[it] = 0.;
      fullLG[it] = halfLG[it] = mouseBiteLG[it] = outerLG[it] = 0.;
      fullCounter[it] = halfCounter[it] = mouseBiteCounter[it] = outerCounter[it] = 0.;
    }
  }
  float fullHG[NUMBER_OF_TIME_SAMPLES],halfHG[NUMBER_OF_TIME_SAMPLES],mouseBiteHG[NUMBER_OF_TIME_SAMPLES],outerHG[NUMBER_OF_TIME_SAMPLES];
  float fullLG[NUMBER_OF_TIME_SAMPLES],halfLG[NUMBER_OF_TIME_SAMPLES],mouseBiteLG[NUMBER_OF_TIME_SAMPLES],outerLG[NUMBER_OF_TIME_SAMPLES];
  int fullCounter[NUMBER_OF_TIME_SAMPLES],halfCounter[NUMBER_OF_TIME_SAMPLES],mouseBiteCounter[NUMBER_OF_TIME_SAMPLES],outerCounter[NUMBER_OF_TIME_SAMPLES];
};
  
class CommonMode{
 public:
  CommonMode( HGCalElectronicsMap &map, bool useMedian=true, bool cmPerChip=true, float threshold=100 );
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
  float _threshold;

  std::map<int,commonModeNoise> _cmMap;
  
};

#endif
