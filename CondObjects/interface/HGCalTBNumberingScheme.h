#ifndef HGCAL_CONDOBJECTS_HGCALTBNUMBERINGSCHEME_H
#define HGCAL_CONDOBJECTS_HGCALTBNUMBERINGSCHEME_H 1

#include "HGCal/CondObjects/interface/HGCalCondObjectContainer.h"

/** \class HGCalTBNumberingScheme
  *  
  * $Date: $
  * $Revision: $
  * \author J. Mans - Minnesota
  */
class HGCalTBNumberingScheme : public HGCalCondObjectNumberingScheme {
public:
  virtual size_t rangeFor(uint64_t scheme) const;
  virtual size_t denseIndexFor(uint32_t rawDetId, uint64_t scheme) const;
};

#endif
