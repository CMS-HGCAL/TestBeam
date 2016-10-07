#include "HGCal/DataFormats/interface/HGCalTBRecHit.h"
#include <cassert>
#include <math.h>

HGCalTBRecHit::HGCalTBRecHit() : CaloRecHit()
{
}


HGCalTBRecHit::HGCalTBRecHit(const DetId& id, float energy, float energyLow, float energyHigh, float time, uint32_t flags) :
	CaloRecHit(id, energy, time, flags),
	_energyLow(energyLow),
	_energyHigh(energyHigh)
{

}

std::ostream& operator<<(std::ostream& s, const HGCalTBRecHit& hit)
{
	return s << hit.id() << ": " << hit.energy() << " GeV, " << hit.time() << " ns";

}
