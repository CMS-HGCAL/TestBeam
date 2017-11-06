#include "HGCal/CondObjects/interface/HGCalElectronicsMap.h"
#include <algorithm>

class DetIdMatch
{
public:
	DetIdMatch(uint32_t matchTo) : m_id(matchTo) { }
	bool operator()(const HGCalElectronicsMap::MapEntry& a)
	{
		return a.detid == m_id;
	}
private:
	uint32_t m_id;
};

bool  HGCalElectronicsMap::MapEntry::operator<(const HGCalElectronicsMap::MapEntry& a) const
{
	return eid < a.eid;
}

bool HGCalElectronicsMap::existsDetId(DetId did) const
{
	DetIdMatch x(did.rawId());
	return std::find_if(m_map.begin(), m_map.end(), x) != m_map.end();
}

bool HGCalElectronicsMap::existsEId(uint32_t eid) const
{
	MapEntry me;
	me.eid = eid;

	std::vector<HGCalElectronicsMap::MapEntry>::const_iterator where = std::lower_bound(m_map.begin(), m_map.end(), me);
	return (where != m_map.end() && where->eid == eid);
}

DetId HGCalElectronicsMap::eid2detId(uint32_t eid) const
{
	MapEntry me;
	me.eid = eid;

	std::vector<HGCalElectronicsMap::MapEntry>::const_iterator where = std::lower_bound(m_map.begin(), m_map.end(), me);
	return (where == m_map.end() || where->eid != eid) ? (DetId(0)) : (DetId(where->detid));
}

uint32_t HGCalElectronicsMap::detId2eid(DetId did) const
{
	DetIdMatch x(did.rawId());
	std::vector<MapEntry>::const_iterator i = std::find_if(m_map.begin(), m_map.end(), x);
	if (i == m_map.end()) return 0;
	else return i->eid;
}

void HGCalElectronicsMap::insert(uint32_t eid, DetId did)
{
	MapEntry me;
	me.eid = eid;
	me.detid = did.rawId();

	std::vector<HGCalElectronicsMap::MapEntry>::const_iterator where = std::lower_bound(m_map.begin(), m_map.end(), me);

	if (where != m_map.end() && where->eid == eid) return; // duplicate
	m_map.insert(where, me);
}

uint32_t HGCalElectronicsMap::eidAt(size_t i) const
{
	return (i >= m_map.size()) ? (0) : (m_map[i].eid);
}
DetId HGCalElectronicsMap::didAt(size_t i) const
{
	return (i >= m_map.size()) ? (DetId(0)) : (DetId(m_map[i].detid));
}
