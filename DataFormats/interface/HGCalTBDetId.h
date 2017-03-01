#ifndef DataFormats_ForwardDetId_HGCalTBDetId_H
#define DataFormats_ForwardDetId_HGCalTBDetId_H 1

#include <iosfwd>
#include "DataFormats/DetId/interface/DetId.h"
#include "DataFormats/ForwardDetId/interface/ForwardSubdetector.h"
#include <iostream>

/*
   On a sensor, the indexes are X and V.  X is horizontal (in most diagrams) and V increases along the hexagon faces towards the upper right.

 */
class HGCalTBDetId : public DetId
{

public:
	static const int kHGCalTBXOffset         = 0;
	static const int kHGCalTBXMask           = 0x1F;
	static const int kHGCalTBXSignMask       = 0x20;
	static const int kHGCalTBVOffset         = 6;
	static const int kHGCalTBVMask           = 0x1F;
	static const int kHGCalTBVSignMask       = 0x20;
	static const int kHGCalTBSensorXOffset   = 12;
	static const int kHGCalTBSensorXMask     = 0x7;
	static const int kHGCalTBSensorXSignMask = 0x8;
	static const int kHGCalTBSensorVOffset   = 16;
	static const int kHGCalTBSensorVMask     = 0x7;
	static const int kHGCalTBSensorVSignMask = 0x8;
	static const int kHGCalTBCellTypeOffset  = 20;
	static const int kHGCalTBCellTypeMask    = 0xF;
	static const int kHGCalLayerOffset       = 24;
	static const int kHGCalLayerMask         = 0x7F;
	static const int kHGCalZsideOffset       = 31;
	static const int kHGCalZsideMask         = 0x1;

	static const int kCellTypeStandard      =   0;
	static const int kCellTypeCalibInner    =   1;
	static const int kCellTypeCalibOuter    =   2;

	/** Create a null cellid*/
	HGCalTBDetId();
	/** Create cellid from raw id (0=invalid tower id) */
	HGCalTBDetId(uint32_t rawid);
	/** Constructor from layer, sensor_iu, sensor_iv, iu, iv, calibr cell numbers */
	HGCalTBDetId(int lay, int sensor_iu, int sensor_iv, int iu, int iv, int cellType);
	/** Constructor from a generic cell id */
	HGCalTBDetId(const DetId& id);
	/** Assignment from a generic cell id */
	HGCalTBDetId& operator=(const DetId& id);

	/// get the absolute value of the cell #'s
	int iu() const
	{
		uint32_t v = id_;
		return (v & kHGCalTBXSignMask) ? ((v & kHGCalTBXMask) - kHGCalTBXSignMask) : (v & kHGCalTBXMask);
	}
	int iv() const
	{
		uint32_t v = id_ >> kHGCalTBVOffset;
		return (v & kHGCalTBVSignMask) ? ((v & kHGCalTBVMask) - kHGCalTBVSignMask) : (v & kHGCalTBVMask);
	}

	/// get the sensor #
	int sensorIU() const
	{
		uint32_t v = id_ >> kHGCalTBSensorXOffset;
		return (v & kHGCalTBSensorXSignMask) ? (-(v & kHGCalTBSensorXMask)) : (v & kHGCalTBSensorXMask);
	}
	int sensorIV() const
	{
		uint32_t v = id_ >> kHGCalTBSensorVOffset;
		return (v & kHGCalTBSensorVSignMask) ? (-(v & kHGCalTBSensorVMask)) : (v & kHGCalTBSensorVMask);
	}

	/// cell type
	int cellType() const
	{
		return ( (id_ >> kHGCalTBCellTypeOffset)&kHGCalTBCellTypeMask);
	}

	/// get the layer #
	int layer() const
	{
		return (id_ >> kHGCalLayerOffset)&kHGCalLayerMask;
	}

	/// get the z-side of the cell (1/-1)
	int zside() const
	{
		return ((id_ >> kHGCalZsideOffset) & kHGCalZsideMask ? -1 : 1);
	}

	/// consistency check : no bits left => no overhead
	bool isHGCal()   const
	{
		return true;
	}
	bool isForward() const
	{
		return true;
	}

};

std::ostream& operator<<(std::ostream&, const HGCalTBDetId& id);

#endif
