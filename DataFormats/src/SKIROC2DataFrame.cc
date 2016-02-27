#include "HGCal/DataFormats/interface/SKIROC2DataFrame.h"

void SKIROC2DataFrame::copyContent(const SKIROC2DataFrame& src) {

}
void SKIROC2DataFrame::setSample(edm::DataFrame::size_type isample, int adc, int tdc) {
  m_data[HEADER_WORDS+WORDS_PER_SAMPLE*isample]=uint16_t(adc);
  m_data[HEADER_WORDS+WORDS_PER_SAMPLE*isample+1]=uint16_t(tdc);
}

std::ostream& operator<<(std::ostream& s, const SKIROC2DataFrame& ski) {
  s << ski.detid() << std::endl;
  for (int i=0; i<ski.samples(); i++)
    s << "    " << ski[i].adc() << " " << ski[i].tdc() << std::endl;      
  return s;
}
