#ifndef HGCAL_TB_ADC_CONVERSION_MAP_H
#define HGCAL_TB_ADC_CONVERSION_MAP_H

#include <iostream>
#include <vector>

class ASIC_ADC_Conversions{
 public:
  ASIC_ADC_Conversions(){;}
  ASIC_ADC_Conversions(uint32_t moduleId,uint32_t asicId)
    {
      _moduleId=moduleId;
      _asicId=asicId;
    }
  ASIC_ADC_Conversions(uint32_t moduleId,uint32_t asicId,
		       float adc_to_MIP,
		       float lowGain_highGain_transition, float lowGain_to_highGain,
		       float TOT_lowGain_transition, float TOT_to_lowGain)
    {
      _moduleId=moduleId;
      _asicId=asicId;
      _adc_to_MIP=adc_to_MIP;
      _lowGain_highGain_transition=lowGain_highGain_transition;
      _lowGain_to_highGain=lowGain_to_highGain;
      _TOT_lowGain_transition=TOT_lowGain_transition;
      _TOT_to_lowGain=TOT_to_lowGain;
    }
  ~ASIC_ADC_Conversions(){;}
  uint32_t moduleId() const {return _moduleId;}
  uint32_t asicId() const {return _asicId;}
  float adc_to_MIP() const {return _adc_to_MIP;}
  float lowGain_highGain_transition() const {return _lowGain_highGain_transition;}
  float lowGain_to_highGain() const {return _lowGain_to_highGain;}
  float TOT_lowGain_transition() const {return _TOT_lowGain_transition;}
  float TOT_to_lowGain() const {return _TOT_to_lowGain;}
  
 private:
  uint32_t _moduleId;
  uint32_t _asicId;
  float _adc_to_MIP;
  float _lowGain_highGain_transition;
  float _lowGain_to_highGain;
  float _TOT_lowGain_transition;
  float _TOT_to_lowGain;
};

std::ostream& operator<<(std::ostream&, ASIC_ADC_Conversions&);

inline bool operator==(const ASIC_ADC_Conversions& rhs0, const ASIC_ADC_Conversions& rhs1) { 
  return rhs0.moduleId() == rhs1.moduleId() && rhs0.asicId() == rhs1.asicId() ;
} 

class HGCalTBADCConversionsMap{
 public:
  HGCalTBADCConversionsMap(){;}
  ~HGCalTBADCConversionsMap(){;}
    
  void addEntry( ASIC_ADC_Conversions conv){ _vec.push_back(conv); }
  std::vector<ASIC_ADC_Conversions> &getEntries(){ return _vec; }
  ASIC_ADC_Conversions getASICConversions(uint32_t moduleId, uint32_t asicId);
  float adc_to_MIP(uint32_t moduleId, uint32_t asicId);
  float lowGain_highGain_transition(uint32_t moduleId, uint32_t asicId);
  float lowGain_to_highGain(uint32_t moduleId, uint32_t asicId);
  float TOT_lowGain_transition(uint32_t moduleId, uint32_t asicId);
  float TOT_to_lowGain(uint32_t moduleId, uint32_t asicId);

 private:
  std::vector<ASIC_ADC_Conversions> _vec;
};

std::ostream& operator<<(std::ostream&, HGCalTBADCConversionsMap&);

#endif
