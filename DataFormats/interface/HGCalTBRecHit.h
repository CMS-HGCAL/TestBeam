#ifndef DATAFORMATS_HGCALTBRECHIT_H
#define DATAFORMATS_HGCALTBRECHIT_H 1

#include "DataFormats/CaloRecHit/interface/CaloRecHit.h"
#include "HGCal/DataFormats/interface/HGCalTBDetId.h"
#include <vector>

/** \class HGCalTBRecHit
 *
 * \author Jeremy Mans
 */

class HGCalTBRecHit : public CaloRecHit
{
public:
	typedef DetId key_type;

	HGCalTBRecHit();
	// by default a recHit is greated with no flag
	//	HGCalTBRecHit(const DetId& id, float energyLow, float energyHigh, float time, uint32_t flags = 0); // when constructing from digis using 2 gains for the ADC
	HGCalTBRecHit(const DetId& id, float energyLow, float energyHigh, float time, double LG2HG, double gainThrE, uint32_t flags = 0); // when constructing from digis using 2 gains for the ADC
	/// get the id
	HGCalTBDetId id() const
	{
		return HGCalTBDetId(detid());
	}
	/////  bool isRecovered() const;
	float _energyLow, _energyHigh, _energy;

	float energyLow() const
	{
		return _energyLow;
	};

	float energyHigh() const
	{
		return _energyHigh;
	};

	float energy() const
	{
	  return _energy;
	};

};

std::ostream& operator<<(std::ostream& s, const HGCalTBRecHit& hit);

#endif
