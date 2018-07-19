#include "HGCal/DataFormats/interface/HGCalTBRecHit.h"
#include "DataFormats/CaloRecHit/interface/CaloRecHit.h"

#include <iostream>
#include <cassert>
#include <math.h>

HGCalTBRecHit::HGCalTBRecHit() : CaloRecHit()
{
}


HGCalTBRecHit::HGCalTBRecHit(const DetId& id, Float16_t energy, Float16_t energyLow, Float16_t energyHigh, unsigned int short energyTot, Float16_t time, uint32_t flags) :
  CaloRecHit(id, energy, time, flags),
  _energyLow(energyLow),
  _energyHigh(energyHigh),
  _energyTot(energyTot),
  cellCenter_x(0),
  cellCenter_y(0),
  m_underSaturationHG(false),
  m_underSaturationLG(false)
{
  _toaRise = -1;
  _toaFall = -1;
  _timeMaxHG = -1;
  _timeMaxLG = -1;
  
  _energyTSLow = std::make_pair(-1, -1);
  _energyTSHigh = std::make_pair(-1, -1);
    
  _energy_HGExcl = -1;    //
}



void HGCalTBRecHit::setCellCenterCoordinate(Float16_t x, Float16_t y) {
  cellCenter_x = x;
  cellCenter_y = y;
}

Float16_t HGCalTBRecHit::getCellCenterCartesianCoordinate(int index) {
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
