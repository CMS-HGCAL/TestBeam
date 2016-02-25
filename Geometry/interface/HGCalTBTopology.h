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
	bool ix_iv_valid(int ix, int iv, int sensorSize) const;
};


#endif
