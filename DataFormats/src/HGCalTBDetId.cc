#include "HGCal/DataFormats/interface/HGCalTBDetId.h"
#include <iostream>

//#define DebugLog

HGCalTBDetId::HGCalTBDetId() : DetId() {
}

HGCalTBDetId::HGCalTBDetId(uint32_t rawid) : DetId(rawid) {
}

HGCalTBDetId::HGCalTBDetId(int lay, int sen, int ix, int iy, bool isCalib) : DetId(0) {  
  uint32_t rawid=0;
  rawid |= ((ix   & kHGCalTBXMask)        << kHGCalTBXOffset);
  rawid |= ((iy   & kHGCalTBYMask)        << kHGCalTBYOffset);
  rawid |= ((sen    & kHGCalTBSensorMask)      << kHGCalTBSensorOffset);  
  rawid |= ((lay    & kHGCalLayerMask)       << kHGCalLayerOffset);
  if (isCalib) rawid |= (kHGCalTBCalibFlagMask << kHGCalTBCalibFlagOffset);
  id_ = rawid;
}

HGCalTBDetId::HGCalTBDetId(const DetId& gen) {
  id_ = gen.rawId();
}

HGCalTBDetId& HGCalTBDetId::operator=(const DetId& gen) {
  id_ = gen.rawId();
  return (*this);
}

std::ostream& operator<<(std::ostream& s,const HGCalTBDetId& id) {
  return s << "HGC(TB)" 
	   << " layer=" << id.layer()  
	   << " sensor=" << id.sensor() << " x,y=" << id.ix() << "," << id.iy() << (id.isCalib()?("CALIB"):(""));
}
