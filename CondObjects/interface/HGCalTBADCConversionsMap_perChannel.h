#ifndef HGCAL_TB_ADC_CONVERSION_MAP_PERCHANNEL_H
#define HGCAL_TB_ADC_CONVERSION_MAP_PERCHANNEL_H

#include <iostream>
#include <vector>

class ASIC_ADC_Conversions_perChannel{
 public:
  ASIC_ADC_Conversions_perChannel(){;}
  ASIC_ADC_Conversions_perChannel(uint32_t moduleId,uint32_t asicId,uint32_t channelId)
    {
      _moduleId=moduleId;
      _asicId=asicId;
      _channelId=channelId;
    }
  ASIC_ADC_Conversions_perChannel(uint32_t moduleId,uint32_t asicId,uint32_t channelId,
		       float adc_to_MIP,
		       float lowGain_highGain_transition, float lowGain_to_highGain,
		       float TOT_lowGain_transition, float TOT_to_lowGain, float TOT_offset, int fully_calibrated)
    {
      _moduleId=moduleId;
      _asicId=asicId;
      _channelId=channelId;
      _adc_to_MIP=adc_to_MIP;
      _lowGain_highGain_transition=lowGain_highGain_transition;
      _lowGain_to_highGain=lowGain_to_highGain;
      _TOT_lowGain_transition=TOT_lowGain_transition;
      _TOT_to_lowGain=TOT_to_lowGain;
      _TOT_offset=TOT_offset;
      _fully_calibrated=fully_calibrated;
    }
  ~ASIC_ADC_Conversions_perChannel(){;}
  uint32_t moduleId() const {return _moduleId;}
  uint32_t asicId() const {return _asicId;}
  uint32_t channelId() const {return _channelId;}
  float adc_to_MIP() const {return _adc_to_MIP;}
  float lowGain_highGain_transition() const {return _lowGain_highGain_transition;}
  float lowGain_to_highGain() const {return _lowGain_to_highGain;}
  float TOT_lowGain_transition() const {return _TOT_lowGain_transition;}
  float TOT_to_lowGain() const {return _TOT_to_lowGain;}
  float TOT_offset() const {return _TOT_offset;}
  float fully_calibrated() const {return _fully_calibrated;}
  
 private:
  uint32_t _moduleId;
  uint32_t _asicId;
  uint32_t _channelId;
  float _adc_to_MIP;
  float _lowGain_highGain_transition;
  float _lowGain_to_highGain;
  float _TOT_lowGain_transition;
  float _TOT_to_lowGain;
  float _TOT_offset;
  int _fully_calibrated;
};


std::ostream& operator<<(std::ostream&, ASIC_ADC_Conversions_perChannel&);

inline bool operator==(const ASIC_ADC_Conversions_perChannel& rhs0, const ASIC_ADC_Conversions_perChannel& rhs1) { 
  return rhs0.moduleId() == rhs1.moduleId() && rhs0.asicId() == rhs1.asicId() && rhs0.channelId() == rhs1.channelId();
} 

class HGCalTBADCConversionsMap_perChannel{
 public:
  HGCalTBADCConversionsMap_perChannel(){;}
  ~HGCalTBADCConversionsMap_perChannel(){;}
    
  void addEntry( ASIC_ADC_Conversions_perChannel conv){ _vec.push_back(conv); }
  std::vector<ASIC_ADC_Conversions_perChannel> &getEntries(){ return _vec; }
  ASIC_ADC_Conversions_perChannel getASICConversions(uint32_t moduleId, uint32_t asicId, uint32_t channelId);
  float adc_to_MIP(uint32_t moduleId, uint32_t asicId, uint32_t channelId);
  float lowGain_highGain_transition(uint32_t moduleId, uint32_t asicId, uint32_t channelId);
  float lowGain_to_highGain(uint32_t moduleId, uint32_t asicId, uint32_t channelId);
  float TOT_lowGain_transition(uint32_t moduleId, uint32_t asicId, uint32_t channelId);
  float TOT_to_lowGain(uint32_t moduleId, uint32_t asicId, uint32_t channelId);
  float TOT_offset(uint32_t moduleId, uint32_t asicId, uint32_t channelId);
  int fully_calibrated(uint32_t moduleId, uint32_t asicId, uint32_t channelId);

 private:
  std::vector<ASIC_ADC_Conversions_perChannel> _vec;
};

std::ostream& operator<<(std::ostream&, HGCalTBADCConversionsMap_perChannel&);

#endif
