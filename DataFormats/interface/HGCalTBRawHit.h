#ifndef HGCALTBRAWHIT_H_INCLUDED
#define HGCALTBRAWHIT_H_INCLUDED 1

#include <vector>
#include "HGCal/DataFormats/interface/HGCalTBDetId.h"

const static int NUMBER_OF_TIME_SAMPLES = 11;

class HGCalTBRawHit
{
 public:

  HGCalTBRawHit( ){;}
 HGCalTBRawHit( unsigned int rawid, unsigned int skiroc, unsigned int channel, std::vector<float> &adcHigh, std::vector<float> &adcLow) :
  m_rawid(rawid),
    m_skiroc(skiroc),
    m_channel(channel),
    m_adcHigh(adcHigh),
    m_adcLow(adcLow)
    {;}
 HGCalTBRawHit( unsigned int rawid, unsigned int skiroc, unsigned int channel, std::vector<float> &adcHigh, std::vector<float> &adcLow, float toaRise, float toaFall, float totSlow, float totFast ) :
  m_rawid(rawid),
    m_skiroc(skiroc),
    m_channel(channel),
    m_adcHigh(adcHigh),
    m_adcLow(adcLow),
    m_toaRise(toaRise),
    m_toaFall(toaFall),
    m_totSlow(totSlow),
    m_totFast(totFast)
    {;}

  HGCalTBDetId detid() const {return HGCalTBDetId(m_rawid);}
  void setRawId(unsigned int rawid){m_rawid=rawid;}
  
  unsigned int skiroc() const{return m_skiroc;}
  unsigned int channel() const{return m_channel;}
  
  void setHighGainADCs(std::vector<float> vec){m_adcHigh=vec;}
  void setLowGainADCs(std::vector<float> vec){m_adcLow=vec;}
  float highGainADC(int timeSample){return m_adcHigh.at(timeSample);}
  float lowGainADC(int timeSample){return m_adcLow.at(timeSample);}

  void setToaRise(float val){m_toaRise=val;}
  void setToaFall(float val){m_toaFall=val;}
  void setTotFast(float val){m_totFast=val;}
  void setTotSlow(float val){m_totSlow=val;}
  float toaRise(){return m_toaRise;}
  float toaFall(){return m_toaFall;}
  float totFast(){return m_totFast;}
  float totSlow(){return m_totSlow;}

 private:
  unsigned int m_rawid; //for some reasons (I don't know) root does not allow saving HGCalTBDetId because of some dictionary issue. 
  // because of non-connected channels for which we assign dummy det id, we must keep skiroc and channel here
  unsigned int m_skiroc;
  unsigned int m_channel;
  std::vector<float> m_adcHigh;
  std::vector<float> m_adcLow;
  float m_toaRise;
  float m_toaFall;
  float m_totSlow;
  float m_totFast;
};

std::ostream& operator<<(std::ostream&, HGCalTBRawHit&);

#endif
