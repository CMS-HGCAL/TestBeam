#include "HGCal/DataFormats/interface/HGCalTBRecHit.h"
#include <cassert>
#include <math.h>

HGCalTBRecHit::HGCalTBRecHit() : CaloRecHit() {
}

HGCalTBRecHit::HGCalTBRecHit(const DetId& id, float energy, float time, uint32_t flags) :
  CaloRecHit(id,energy,time,flags) {
}

std::ostream& operator<<(std::ostream& s, const HGCalTBRecHit& hit) {
  return s << hit.id() << ": " << hit.energy() << " GeV, " << hit.time() << " ns";

}
