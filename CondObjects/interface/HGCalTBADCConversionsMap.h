#ifndef HGCAL_TB_ADC_CONVERSION_MAP_H
#define HGCAL_TB_ADC_CONVERSION_MAP_H

#include <HGCal/DataFormats/interface/HGCalTBDetId.h>
#include <iostream>
#include <vector>

struct ADCConversionParameters{
  float hgtomip;
  float lgtomip;
  float totcoeff,totped,totthr,totnorm,totpower;
  float hgsat,lgsat;  
};

class ADCConversions{
 public:
  ADCConversions(HGCalTBDetId id, float hgToMip=0, float lgToMip=0, float totCoeff=0, float totPed=0, float totThr=0, float totNorm=0, float totPower=0, int hgSat=1600, int lgSat=1100)
    {
      m_detid=id;
      m_params.hgtomip=hgToMip;
      m_params.lgtomip=lgToMip;
      m_params.totcoeff=totCoeff;
      m_params.totped=totPed;
      m_params.totthr=totThr;
      m_params.totnorm=totNorm;
      m_params.totpower=totPower;
      m_params.hgsat=hgSat;
      m_params.lgsat=lgSat;
    }
  bool getCalibEnergy(float hg, float lg, int tot, float &energy);
  const ADCConversionParameters &getParameters() const {return m_params;}
  const HGCalTBDetId getDetId() const {return m_detid;}
  bool operator==(const ADCConversions& adcConv) const
  {
    return m_detid==adcConv.getDetId();
  } 
  friend std::ostream&operator<<(std::ostream&, const ADCConversions &adcConv);
  
 private:
  float totShape(float nmip);
 private:
  HGCalTBDetId m_detid;
  ADCConversionParameters m_params;
};


class HGCalTBADCConversionsMap{
 public:
  HGCalTBADCConversionsMap(){;}
  ~HGCalTBADCConversionsMap(){;}
  
  const std::vector<ADCConversions> &getADCConversionsMap(){ return m_vec; }
  void addEntry( ADCConversions conv);
  ADCConversions getADCConversions(HGCalTBDetId detid);
  bool getCalibratedEnergy(HGCalTBDetId id, float hg, float lg, int tot, float& energy);

 private:
  std::vector<ADCConversions> m_vec;
};

#endif
