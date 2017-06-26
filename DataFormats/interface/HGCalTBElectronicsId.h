#ifndef HGCAL_DATAFORMATS_HGCALTBELECTRONICSID_H
#define HGCAL_DATAFORMATS_HGCALTBELECTRONICSID_H 1

#include <iosfwd>
#include <stdint.h>

/** \class HGCalTBElectronicsId
  *
  * [6:0]  ICHAN
  * [10:7] ISKIROC
  *
  * $Date: $
  * $Revision: $
  * \author J. Mans - Minnesota
  */
class HGCalTBElectronicsId
{
public:
	static const int kIChanOffset   = 0;
	static const int kIChanMask     = 0x7F;
	static const int kISkiRocOffset = 7;
	static const int kISkiRocMask   = 0x1FF;

	HGCalTBElectronicsId(int iskiroc, int ichan);
	HGCalTBElectronicsId() : m_id(0) { }
	HGCalTBElectronicsId(uint32_t rawId) : m_id(rawId) { }

	int ichan() const
	{
		return (m_id & kIChanMask);
	}
	int iskiroc() const
	{
		return (m_id >> kISkiRocOffset)&kISkiRocMask;
	}

	uint32_t rawId() const
	{
		return m_id;
	}

	bool null() const
	{
		return m_id == 0;
	}

	operator uint32_t() const
	{
		return m_id;
	}
private:
	uint32_t m_id;
};

std::ostream& operator<<(std::ostream&, const HGCalTBElectronicsId& id);


#endif
