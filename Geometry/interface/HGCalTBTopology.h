#ifndef HGCAL_GEOMETRY_HGCALTBTOPOLOGY_H
#define HGCAL_GEOMETRY_HGCALTBTOPOLOGY_H 1

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
	double Cell_Area(int cell_type) const;//returns area in cm*cm
};

#endif
