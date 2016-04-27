#include "HGCal/DataFormats/interface/SKIROC2DataFrame.h"

void SKIROC2DataFrame::copyContent(const SKIROC2DataFrame& src)
{

}
void SKIROC2DataFrame::setSample(edm::DataFrame::size_type isample, int adcLow, int adcHigh, int tdc)
{
	unsigned int pos = HEADER_WORDS + WORDS_PER_SAMPLE * isample;
	m_data[ pos + Sample::ADCLOW_SHIFT  ] = uint16_t(adcLow);
	m_data[ pos + Sample::ADCHIGH_SHIFT ] = uint16_t(adcHigh);
	m_data[ pos + Sample::TDC_SHIFT     ] = uint16_t(tdc);
}

std::ostream& operator<<(std::ostream& s, const SKIROC2DataFrame& ski)
{
	s << ski.detid() << std::endl;
	for (int i = 0; i < ski.samples(); i++)
		s << "    " << ski[i].adcLow()
		  << " " << ski[i].adcHigh()
		  << " " << ski[i].tdc()
		  << " " << ski.samples()  << std::endl;
	return s;
}
