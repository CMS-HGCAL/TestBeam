#include "HGCal/CondObjects/interface/HGCalCondObjectTextIO.h"
#include <HGCal/CondObjects/interface/HGCalTBADCConversionsMap.h>
#include <algorithm>

ASIC_ADC_Conversions HGCalTBADCConversionsMap::getASICConversions(uint32_t moduleId, uint32_t asicId)
{
  ASIC_ADC_Conversions adcConv(moduleId,asicId);
  return (*std::find(_vec.begin(),_vec.end(),adcConv));
}

float HGCalTBADCConversionsMap::adc_to_MIP(uint32_t moduleId, uint32_t asicId)
{
  ASIC_ADC_Conversions adcConv(moduleId,asicId);
  std::vector<ASIC_ADC_Conversions>::iterator it=std::find(_vec.begin(),_vec.end(),adcConv);
  return it!=_vec.end() ? (*it).adc_to_MIP() : 0;
}

float HGCalTBADCConversionsMap::lowGain_highGain_transition(uint32_t moduleId, uint32_t asicId)
{
  ASIC_ADC_Conversions adcConv(moduleId,asicId);
  std::vector<ASIC_ADC_Conversions>::iterator it=std::find(_vec.begin(),_vec.end(),adcConv);
  return it!=_vec.end() ? (*it).lowGain_highGain_transition() : 0;
}

float HGCalTBADCConversionsMap::lowGain_to_highGain(uint32_t moduleId, uint32_t asicId)
{
  ASIC_ADC_Conversions adcConv(moduleId,asicId);
  std::vector<ASIC_ADC_Conversions>::iterator it=std::find(_vec.begin(),_vec.end(),adcConv);
  return it!=_vec.end() ? (*it).lowGain_to_highGain() : 0;
}

float HGCalTBADCConversionsMap::TOT_lowGain_transition(uint32_t moduleId, uint32_t asicId)
{
  ASIC_ADC_Conversions adcConv(moduleId,asicId);
  std::vector<ASIC_ADC_Conversions>::iterator it=std::find(_vec.begin(),_vec.end(),adcConv);
  return it!=_vec.end() ? (*it).TOT_lowGain_transition() : 0;
}

float HGCalTBADCConversionsMap::TOT_to_lowGain(uint32_t moduleId, uint32_t asicId)
{
  ASIC_ADC_Conversions adcConv(moduleId,asicId);
  std::vector<ASIC_ADC_Conversions>::iterator it=std::find(_vec.begin(),_vec.end(),adcConv);
  return it!=_vec.end() ? (*it).TOT_to_lowGain() : 0;
}

std::ostream& operator<<(std::ostream& s, ASIC_ADC_Conversions& adc){
  s << adc.moduleId() << " "
    << adc.asicId() << " "
    << adc.adc_to_MIP() << " "
    << adc.lowGain_highGain_transition() << " "
    << adc.lowGain_to_highGain() << " "
    << adc.TOT_lowGain_transition() << " "
    << adc.TOT_to_lowGain() << "\n";
  return s;
}


std::ostream& operator<<(std::ostream& s, HGCalTBADCConversionsMap& a_map)
{
  s << "moduleId "
    << "asicId "
    << "adc_to_MIP "
    << "lowGain_highGain_transition "
    << "lowGain_to_highGain "
    << "TOT_lowGain_transition "
    << "TOT_to_lowGain \n";
  for( std::vector<ASIC_ADC_Conversions>::iterator it=a_map.getEntries().begin(); it!=a_map.getEntries().end(); ++it )
    s << (*it) ;
  return s;
}
