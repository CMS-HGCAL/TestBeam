#include "HGCal/DataFormats/interface/HGCalTBRawHit.h"

std::ostream& operator<<(std::ostream& s, HGCalTBRawHit& hit)
{
  s << "Channel det id : " << hit.detid() << "\n";
  s << "High gain ADC => " ;
  for (size_t i = 0; i < NUMBER_OF_TIME_SAMPLES; i++)
    s << hit.highGainADC(i) << " ";
  s << "Low gain ADC => " ;
  for (size_t i = 0; i < NUMBER_OF_TIME_SAMPLES; i++)
    s << hit.lowGainADC(i) << " ";
  s << std::endl;
  return s;
}
