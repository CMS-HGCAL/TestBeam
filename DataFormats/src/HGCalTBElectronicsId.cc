#include "HGCal/DataFormats/interface/HGCalTBElectronicsId.h"
#include <iostream>


HGCalTBElectronicsId::HGCalTBElectronicsId(int iskiroc, int ichan) : m_id(0)
{
	m_id |= (iskiroc & kISkiRocMask) << kISkiRocOffset;
	m_id |= (ichan & kIChanMask);
}

std::ostream& operator<<(std::ostream& s, const HGCalTBElectronicsId& id)
{
	return s << "HGCalTB skiroc " << id.iskiroc() << "." << id.ichan();
}
