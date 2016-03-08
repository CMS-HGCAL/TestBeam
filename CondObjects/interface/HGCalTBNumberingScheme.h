#ifndef HGCAL_CONDOBJECTS_HGCALTBNUMBERINGSCHEME_H
#define HGCAL_CONDOBJECTS_HGCALTBNUMBERINGSCHEME_H 1

#include "HGCal/CondObjects/interface/HGCalCondObjectContainer.h"

/** \class HGCalTBNumberingScheme
  *
  * Simple-minded numbering scheme appropriate for a testbeam geometry.
  * Scheme version 0 assumes up to 28 planes of single 128-cell sensors.
  *
  * $Date: $
  * $Revision: $
  * \author J. Mans - Minnesota
  */
class HGCalTBNumberingScheme : public HGCalCondObjectNumberingScheme
{
public:
	virtual size_t rangeFor(uint64_t scheme) const;
	virtual size_t denseIndexFor(uint32_t rawDetId, uint64_t scheme) const;
	static const HGCalCondObjectNumberingScheme* scheme()
	{
		return &the_scheme;
	}
private:
	HGCalTBNumberingScheme() { }
	static HGCalTBNumberingScheme the_scheme;
};

#endif
