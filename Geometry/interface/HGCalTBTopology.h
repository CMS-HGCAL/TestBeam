#ifndef HGCAL_GEOMETRY_HGCALTBTOPOLOGY_H
#define HGCAL_GEOMETRY_HGCALTBTOPOLOGY_H 1
#include <vector>
#include "HGCal/DataFormats/interface/HGCalTBDetId.h"
/** \class HGCalTBTopology
  *
  * Reference: https://indico.cern.ch/event/456955/
  *
  * $Date: $
  * $Revision: $
  * \author J. Mans - Minnesota
  */
class HGCalTBTopology
{
public:
	// valid sensorSizes are 128 and 256
	bool iu_iv_valid(int layer, int sensor_iu, int sensor_iv, int iu, int iv, int sensorSize) const;
	inline bool isValidDetId(HGCalTBDetId detId, int sensorSize) const
	{
		return iu_iv_valid(detId.layer(), detId.sensorIU(), detId.sensorIV(), detId.iu(), detId.iv(), sensorSize);
	};
	double Cell_Area(int cell_type) const;//returns area in cm*cm
	std::vector<HGCalTBDetId> validDetIds(unsigned int layer, unsigned int sensorsize);

};

#endif
