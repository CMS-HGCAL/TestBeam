#include "HGCal/DataFormats/interface/HGCalTBRecHit.h"
#include <cassert>
#include <math.h>

HGCalTBRecHit::HGCalTBRecHit() : CaloRecHit()
{
}


HGCalTBRecHit::HGCalTBRecHit(const DetId& id, float energyLow, float energyHigh, float time, uint32_t flags) :
	CaloRecHit(id, energyHigh, time, flags),
	_energyLow(energyLow),
	_energyHigh(energyHigh)
{

	///\todo set the default recHit energy to the highGain values unless saturated
	setEnergy( ( energyLow < _lowGainSaturationThreshold ) ? energyLow : energyHigh);
}

std::ostream& operator<<(std::ostream& s, const HGCalTBRecHit& hit)
{
	return s << hit.id() << ": " << hit.energy() << " GeV, " << hit.time() << " ns";

}
