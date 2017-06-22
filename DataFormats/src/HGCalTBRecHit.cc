#include "HGCal/DataFormats/interface/HGCalTBRecHit.h"
#include "DataFormats/CaloRecHit/interface/CaloRecHit.h"

#include <iostream>
#include <cassert>
#include <math.h>

HGCalTBRecHit::HGCalTBRecHit() : CaloRecHit()
{
}


HGCalTBRecHit::HGCalTBRecHit(const DetId& id, float energy, float energyLow, float energyHigh, float energyTot, float time, uint32_t flags) :
  CaloRecHit(id, energy, time, flags),
  _energyLow(energyLow),
  _energyHigh(energyHigh),
  _energyTot(energyTot),
  cellCenter_x(0),
  cellCenter_y(0)
{
}


void HGCalTBRecHit::setTime(float time){
  _time = time;
  return;
}

void HGCalTBRecHit::setCellCenterCoordinate(float x, float y) {
  cellCenter_x = x;
  cellCenter_y = y;
}

float HGCalTBRecHit::getCellCenterCartesianCoordinate(int index) {
  switch(index) {
    case 0:
      return cellCenter_x;
    case 1:
      return cellCenter_y;
    default:
      return -999;
  }
}


std::ostream& operator<<(std::ostream& s, const HGCalTBRecHit& hit)
{
  return s << hit.id() << ": " << hit.energy() << " GeV, " << hit._time << " ns";

}
