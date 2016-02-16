#ifndef HGCAL_CONDOBJECTS_HGCALCONDOBJECTTEXTIO_H
#define HGCAL_CONDOBJECTS_HGCALCONDOBJECTTEXTIO_H 1

#include <string>
#include "HGCal/CondObjects/interface/HGCalCondObjectContainer.h"

/** \class HGCalCondObjectTextIO
  *  
  * $Date: $
  * $Revision: $
  * \author J. Mans - Minnesota
  */
class HGCalCondObjectTextIO {
public:
  HGCalCondObjectTextIO(const HGCalCondObjectNumberingScheme* scheme) : p_scheme(scheme) { }
  bool load(const std::string& filename, HGCalCondObjectContainer<double>&);
  bool store(const std::string& filename, const HGCalCondObjectContainer<double>&);
private:
  const HGCalCondObjectNumberingScheme* p_scheme;
};


#endif
