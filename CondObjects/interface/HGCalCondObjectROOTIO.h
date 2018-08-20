#ifndef HGCAL_CONDOBJECTS_HGCALCONDOBJECTROOTIO_H
#define HGCAL_CONDOBJECTS_HGCALCONDOBJECTROOTIO_H 1

#include <string>
#include "HGCal/CondObjects/interface/HGCalCondObjectContainer.h"
#include "HGCal/CondObjects/interface/HGCalTBADCConversionsMap.h"
#include "HGCal/CondObjects/interface/HGCalTBDetectorLayout.h"
#include "HGCal/CondObjects/interface/HGCalElectronicsMap.h"

class HGCalCondObjectROOTIO
{
 public:
 HGCalCondObjectROOTIO(HGCalTBDetectorLayout layout, HGCalElectronicsMap emap) : m_layout(layout), m_emap(emap){;}
  bool loadADCConversion(const std::string& filename, HGCalTBADCConversionsMap& adcConvMap);
 private:
  HGCalTBDetectorLayout m_layout;
  HGCalElectronicsMap m_emap;
};


#endif
