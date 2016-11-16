#include "HGCal/DataFormats/interface/HGCalTBRecHit.h"
#include <iostream>
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
  cartesian_coordinates = new double[3];
}

void HGCalTBRecHit::setCartesianCoordinates(double x, double y, double z){
  cartesian_coordinates[0] = x;
  cartesian_coordinates[1] = y;
  cartesian_coordinates[2] = z;
}

double HGCalTBRecHit::getCartesianCoordinate(int index){
  assert(index<=2);
  return cartesian_coordinates[index];
}

std::ostream& operator<<(std::ostream& s, const HGCalTBRecHit& hit)
{
	return s << hit.id() << ": " << hit.energy() << " GeV, " << hit.time() << " ns";

}
