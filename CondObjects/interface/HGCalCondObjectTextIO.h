#ifndef HGCAL_CONDOBJECTS_HGCALCONDOBJECTTEXTIO_H
#define HGCAL_CONDOBJECTS_HGCALCONDOBJECTTEXTIO_H 1

#include <string>
#include "HGCal/CondObjects/interface/HGCalCondObjectContainer.h"
#include "HGCal/CondObjects/interface/HGCalElectronicsMap.h"

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

  bool load(const std::string& filename, HGCalElectronicsMap&);
  bool store(const std::string& filename, const HGCalElectronicsMap&);
private:
  const HGCalCondObjectNumberingScheme* p_scheme;
};


#endif
