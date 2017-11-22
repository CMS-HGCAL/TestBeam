#ifndef HGCalTBCommonModeNoise_H
#define HGCalTBCommonModeNoise_H
#include <map>
#include "HGCal/DataFormats/interface/HGCalTBRawHitCollection.h"

struct commonModeNoise{
  commonModeNoise(){
    for(int it=0;it<NUMBER_OF_TIME_SAMPLES; it++){
      fullHG[it] = halfHG[it] = mouseBiteHG[it] = outerHG[it] = mergedHG[it] = 0.;
      fullLG[it] = halfLG[it] = mouseBiteLG[it] = outerLG[it] = mergedLG[it] = 0.;
      fullCounter[it] = halfCounter[it] = mouseBiteCounter[it] = outerCounter[it] = 0.;
    }
  }
  float fullHG[NUMBER_OF_TIME_SAMPLES],halfHG[NUMBER_OF_TIME_SAMPLES],mouseBiteHG[NUMBER_OF_TIME_SAMPLES],outerHG[NUMBER_OF_TIME_SAMPLES],mergedHG[NUMBER_OF_TIME_SAMPLES];
  float fullLG[NUMBER_OF_TIME_SAMPLES],halfLG[NUMBER_OF_TIME_SAMPLES],mouseBiteLG[NUMBER_OF_TIME_SAMPLES],outerLG[NUMBER_OF_TIME_SAMPLES],mergedLG[NUMBER_OF_TIME_SAMPLES];
  int fullCounter[NUMBER_OF_TIME_SAMPLES],halfCounter[NUMBER_OF_TIME_SAMPLES],mouseBiteCounter[NUMBER_OF_TIME_SAMPLES],outerCounter[NUMBER_OF_TIME_SAMPLES];
};


#endif