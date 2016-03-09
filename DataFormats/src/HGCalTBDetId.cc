#include "HGCal/DataFormats/interface/HGCalTBDetId.h"
#include <iostream>

//#define DebugLog

HGCalTBDetId::HGCalTBDetId() : DetId()
{
}

HGCalTBDetId::HGCalTBDetId(uint32_t rawid) : DetId(rawid)
{
}

HGCalTBDetId::HGCalTBDetId(int lay, int sen_iu, int sen_iv, int iu, int iv, int cellType) : DetId(0)
{
	uint32_t rawid = 0;
	rawid |= ((iu   & (kHGCalTBXSignMask | kHGCalTBXMask))        << kHGCalTBXOffset);
	rawid |= ((iv   & (kHGCalTBVSignMask | kHGCalTBVMask))        << kHGCalTBVOffset);
	rawid |= ((sen_iu  & (kHGCalTBSensorXSignMask | kHGCalTBSensorXMask))      << kHGCalTBSensorXOffset);
	rawid |= ((sen_iv  & (kHGCalTBSensorVSignMask || kHGCalTBSensorVMask))      << kHGCalTBSensorVOffset);
	rawid |= ((lay    & kHGCalLayerMask)       << kHGCalLayerOffset);
	rawid |= ((cellType & kHGCalTBCellTypeMask) << kHGCalTBCellTypeOffset);
	id_ = rawid;
}

HGCalTBDetId::HGCalTBDetId(const DetId& gen)
{
	id_ = gen.rawId();
}

HGCalTBDetId& HGCalTBDetId::operator=(const DetId& gen)
{
	id_ = gen.rawId();
	return (*this);
}

std::ostream& operator<<(std::ostream& s, const HGCalTBDetId& id)
{
	s << "HGC(TB)"
	  << " layer=" << id.layer()
	  << " sensor(iu,iv)=" << id.sensorIU() << "," << id.sensorIV() << " iu,iv=" << id.iu() << "," << id.iv();
	switch (id.cellType()) {
	case HGCalTBDetId::kCellTypeCalibInner:
		s << "Calib(Inner)";
		break;
	case HGCalTBDetId::kCellTypeCalibOuter:
		s << "Calib(Outer)";
		break;
	default:
		break;
	}

	return s;
}
