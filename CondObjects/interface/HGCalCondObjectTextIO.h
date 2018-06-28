#ifndef HGCAL_CONDOBJECTS_HGCALCONDOBJECTTEXTIO_H
#define HGCAL_CONDOBJECTS_HGCALCONDOBJECTTEXTIO_H 1

#include <string>
#include "HGCal/CondObjects/interface/HGCalCondObjectContainer.h"
#include "HGCal/CondObjects/interface/HGCalElectronicsMap.h"
#include "HGCal/CondObjects/interface/HGCalTBDetectorLayout.h"
#include "HGCal/CondObjects/interface/HGCalTBADCConversionsMap.h"
#include "HGCal/CondObjects/interface/HGCalTBADCConversionsMap_perChannel.h"

/** \class HGCalCondObjectTextIO
  *
  * $Date: $
  * $Revision: $
  * \author J. Mans - Minnesota
  *
  * \todo load and store are HGCalCondObjectContainer not templated.... so now takes only floats!!!
  */
class HGCalCondObjectTextIO
{
public:
	HGCalCondObjectTextIO(const HGCalCondObjectNumberingScheme* scheme) : p_scheme(scheme) { }
	bool load(const std::string& filename, HGCalCondObjectContainer<float>&); ///<load conditions from file
	bool store(const std::string& filename, const HGCalCondObjectContainer<float>&); ///< saves condition to file

	bool load(const std::string& filename, HGCalElectronicsMap&);
	bool store(const std::string& filename, const HGCalElectronicsMap&);

	bool load(const std::string& filename, HGCalTBDetectorLayout& hgcalGeom);
	bool load(const std::string& filename, HGCalTBADCConversionsMap& adcConvMap);
	bool load(const std::string& filename, HGCalTBADCConversionsMap_perChannel& adcConvMap);
private:
	const HGCalCondObjectNumberingScheme* p_scheme;
};


#endif
