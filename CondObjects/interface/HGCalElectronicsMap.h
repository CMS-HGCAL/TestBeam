#ifndef HGCAL_CONDOBJECTS_HGCALELECTRONICSMAP_H
#define HGCAL_CONDOBJECTS_HGCALELECTRONICSMAP_H 1

#include "DataFormats/DetId/interface/DetId.h"
#include <vector>

/** \class HGCalElectronicsMap
  *
  * $Date: $
  * $Revision: $
  * \author J. Mans - Minnesota
  *
  * \brief provides the conversion between electronics Id to DetId
  *
  *
  */
class HGCalElectronicsMap
{
public:
	bool existsDetId(DetId did) const;
	bool existsEId(uint32_t eid) const;

	DetId eid2detId(uint32_t eid) const;
	uint32_t detId2eid(DetId did) const;

	void insert(uint32_t, DetId did);

	size_t size() const
	{
		return m_map.size();
	}
	uint32_t eidAt(size_t i) const;
	DetId didAt(size_t i) const;

	struct MapEntry {
		uint32_t eid;
		uint32_t detid;
		bool operator<(const MapEntry&) const;
	};
private:
	std::vector<MapEntry> m_map; // ordered for eid2detid
};


#endif
