#ifndef DATAFORMATS_HGCALTBRECHIT_H
#define DATAFORMATS_HGCALTBRECHIT_H 1

#include "DataFormats/CaloRecHit/interface/CaloRecHit.h"
#include "HGCal/DataFormats/interface/HGCalTBDetId.h"
#include <vector>

/** \class HGCRecHit
 *
 * \author Valeri Andreev
 */

class HGCalTBRecHit : public CaloRecHit
{
public:
	typedef DetId key_type;

	HGCalTBRecHit();
	// by default a recHit is greated with no flag
	HGCalTBRecHit(const DetId& id, float energy, float time, uint32_t flags = 0);
	/// get the id
	HGCalTBDetId id() const
	{
		return HGCalTBDetId(detid());
	}
	/////  bool isRecovered() const;
};

std::ostream& operator<<(std::ostream& s, const HGCalTBRecHit& hit);

#endif
