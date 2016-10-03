#include "HGCal/DataFormats/interface/HGCalTBRecHit.h"
#include <cassert>
#include <math.h>

HGCalTBRecHit::HGCalTBRecHit() : CaloRecHit()
{
}


HGCalTBRecHit::HGCalTBRecHit(const DetId& id, float energyLow, float energyHigh, float time, double LG2HG, double gainThrE, uint32_t flags) :
	CaloRecHit(id, energyHigh, time, flags),
	_energyLow(energyLow),
	_energyHigh(energyHigh)
{
  if(_energyHigh < gainThrE) _energy = _energyHigh;
  else _energy = _energyLow * LG2HG;

	///\todo set the default recHit energy to the highGain values unless saturated
}

std::ostream& operator<<(std::ostream& s, const HGCalTBRecHit& hit)
{
  return s << hit.id() << ": " << hit.energy() << " GeV, " << hit.time() << " ns";

}
