#include "HGCal/CondObjects/interface/HGCalCondObjectTextIO.h"
#include <HGCal/CondObjects/interface/HGCalTBADCConversionsMap_perChannel.h>
#include <algorithm>

ASIC_ADC_Conversions_perChannel HGCalTBADCConversionsMap_perChannel::getASICConversions(uint32_t moduleId, uint32_t asicId, uint32_t channelId)
{
  ASIC_ADC_Conversions_perChannel adcConv(moduleId,asicId, channelId);
  return (*std::find(_vec.begin(),_vec.end(),adcConv));
}

float HGCalTBADCConversionsMap_perChannel::adc_to_MIP(uint32_t moduleId, uint32_t asicId, uint32_t channelId)
{
  ASIC_ADC_Conversions_perChannel adcConv(moduleId,asicId, channelId);
  std::vector<ASIC_ADC_Conversions_perChannel>::iterator it=std::find(_vec.begin(),_vec.end(),adcConv);
  return it!=_vec.end() ? (*it).adc_to_MIP() : -1;
}

float HGCalTBADCConversionsMap_perChannel::lowGain_highGain_transition(uint32_t moduleId, uint32_t asicId, uint32_t channelId)
{
  ASIC_ADC_Conversions_perChannel adcConv(moduleId,asicId, channelId);
  std::vector<ASIC_ADC_Conversions_perChannel>::iterator it=std::find(_vec.begin(),_vec.end(),adcConv);
  return it!=_vec.end() ? (*it).lowGain_highGain_transition() : -1;
}

float HGCalTBADCConversionsMap_perChannel::lowGain_to_highGain(uint32_t moduleId, uint32_t asicId, uint32_t channelId)
{
  ASIC_ADC_Conversions_perChannel adcConv(moduleId,asicId, channelId);
  std::vector<ASIC_ADC_Conversions_perChannel>::iterator it=std::find(_vec.begin(),_vec.end(),adcConv);
  return it!=_vec.end() ? (*it).lowGain_to_highGain() : -1;
}

float HGCalTBADCConversionsMap_perChannel::TOT_lowGain_transition(uint32_t moduleId, uint32_t asicId, uint32_t channelId)
{
  ASIC_ADC_Conversions_perChannel adcConv(moduleId,asicId, channelId);
  std::vector<ASIC_ADC_Conversions_perChannel>::iterator it=std::find(_vec.begin(),_vec.end(),adcConv);
  return it!=_vec.end() ? (*it).TOT_lowGain_transition() : -1;
}

float HGCalTBADCConversionsMap_perChannel::TOT_to_lowGain(uint32_t moduleId, uint32_t asicId, uint32_t channelId)
{
  ASIC_ADC_Conversions_perChannel adcConv(moduleId,asicId, channelId);
  std::vector<ASIC_ADC_Conversions_perChannel>::iterator it=std::find(_vec.begin(),_vec.end(),adcConv);
  return it!=_vec.end() ? (*it).TOT_to_lowGain() : -1;
}

float HGCalTBADCConversionsMap_perChannel::TOT_offset(uint32_t moduleId, uint32_t asicId, uint32_t channelId)
{
  ASIC_ADC_Conversions_perChannel adcConv(moduleId,asicId, channelId);
  std::vector<ASIC_ADC_Conversions_perChannel>::iterator it=std::find(_vec.begin(),_vec.end(),adcConv);
  return it!=_vec.end() ? (*it).TOT_offset() : -1;
}

int HGCalTBADCConversionsMap_perChannel::fully_calibrated(uint32_t moduleId, uint32_t asicId, uint32_t channelId)
{
  ASIC_ADC_Conversions_perChannel adcConv(moduleId,asicId, channelId);
  std::vector<ASIC_ADC_Conversions_perChannel>::iterator it=std::find(_vec.begin(),_vec.end(),adcConv);
  return it!=_vec.end() ? (*it).fully_calibrated() : -1;
}
std::ostream& operator<<(std::ostream& s, ASIC_ADC_Conversions_perChannel& adc){
  s << adc.moduleId() << " "
    << adc.asicId() << " "
    << adc.channelId() << " "
    << adc.adc_to_MIP() << " "
    << adc.lowGain_highGain_transition() << " "
    << adc.lowGain_to_highGain() << " "
    << adc.TOT_lowGain_transition() << " "
    << adc.TOT_to_lowGain() << " "
    << adc.TOT_offset() << " "
    << adc.fully_calibrated() << "\n";
  return s;
}


std::ostream& operator<<(std::ostream& s, HGCalTBADCConversionsMap_perChannel& a_map)
{
  s << "moduleId "
    << "asicId "
    << "channelId "
    << "adc_to_MIP "
    << "lowGain_highGain_transition "
    << "lowGain_to_highGain "
    << "TOT_lowGain_transition "
    << "TOT_to_lowGain "
    << "TOT_offset "
    << "fully_calibrated \n";
  for( std::vector<ASIC_ADC_Conversions_perChannel>::iterator it=a_map.getEntries().begin(); it!=a_map.getEntries().end(); ++it )
    s << (*it) ;
  return s;
}
