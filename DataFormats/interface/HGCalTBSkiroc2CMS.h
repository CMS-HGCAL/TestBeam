#ifndef HGCALTBSKIROC2CMS_H_INCLUDED
#define HGCALTBSKIROC2CMS_H_INCLUDED 1

#include <stdint.h>
#include <array>
#include "HGCal/DataFormats/interface/HGCalTBDetId.h"

static const int MASK_ADC = 0x0FFF;
static const int MASK_ROLL = 0x1FFF;
static const int MASK_GTS_MSB = 0x2FFF;
static const int MASK_GTS_LSB = 0x1FFF;
static const int MASK_ID = 0xFF;
static const int MASK_HEAD = 0xF000;

static const size_t NUMBER_OF_SCA = 15;
static const size_t NUMBER_OF_CHANNELS = 64;
static const size_t ADCLOW_SHIFT = 0;
static const size_t ADCHIGH_SHIFT = 64;
static const size_t SCA_SHIFT = 128;
static const size_t SKIROC_DATA_SIZE = 1924; //number of 16 bits words

class HGCalTBSkiroc2CMS
{
 public:

 HGCalTBSkiroc2CMS( const std::array<uint16_t,SKIROC_DATA_SIZE>& data, const std::array<HGCalTBDetId,NUMBER_OF_CHANNELS> &ids ) :
  m_data(data),
  m_id(ids)
  {;}
  uint16_t gray_to_brady(const uint16_t gray) const;
    
  uint16_t ADCLow( int chan, int sca ) const { return (sca>0 && sca<14) ? gray_to_brady( m_data.at(chan+ADCLOW_SHIFT+SCA_SHIFT*sca) & MASK_ADC ) : 10000 ; }
  uint16_t ADCHigh( int chan, int sca ) const {return (sca>0 && sca<14) ? gray_to_brady( m_data.at(chan+ADCHIGH_SHIFT+SCA_SHIFT*sca) & MASK_ADC ) : 10000;}
  uint16_t TOTSlow( int chan ) const {return gray_to_brady( m_data.at(chan+ADCLOW_SHIFT+SCA_SHIFT*15) & MASK_ADC );}
  uint16_t TOTFast( int chan ) const {return gray_to_brady( m_data.at(chan+ADCHIGH_SHIFT+SCA_SHIFT*15) & MASK_ADC );}
  uint16_t TOAFall( int chan ) const {return gray_to_brady( m_data.at(chan+ADCLOW_SHIFT+SCA_SHIFT*14) & MASK_ADC );}
  uint16_t TOARise( int chan ) const {return gray_to_brady( m_data.at(chan+ADCHIGH_SHIFT+SCA_SHIFT*14) & MASK_ADC );}
  bool TOAHitFall(int chan) const { return ((m_data.at(chan+ADCLOW_SHIFT)&~MASK_ADC)>>4*3) ;}
  bool TOAHitRise(int chan) const { return ((m_data.at(chan+ADCHIGH_SHIFT)&~MASK_ADC)>>4*3) ;}
  uint16_t rollMask() const { return (m_data.at(SKIROC_DATA_SIZE-4)&MASK_ROLL); }
    
 private:
  uint32_t globalTS_MSB() const { return (m_data.at(SKIROC_DATA_SIZE-3) & MASK_GTS_MSB ); }
  uint32_t globalTS_LSB() const { return ( (m_data.at(SKIROC_DATA_SIZE-2) & MASK_GTS_LSB)>>1 ); }
 public:
  uint32_t globalTS()const{ return gray_to_brady( globalTS_MSB() | globalTS_LSB()); }
  int skirocId()const{ return (m_data.at(SKIROC_DATA_SIZE-1) & MASK_ID); }
  bool check();
  

  HGCalTBDetId detid( int chan ) const
  {
    return HGCalTBDetId(m_id.at(chan));
  }

 private:
  std::array<uint16_t,SKIROC_DATA_SIZE> m_data;
  std::array<HGCalTBDetId,NUMBER_OF_CHANNELS> m_id;
};

std::ostream& operator<<(std::ostream&, const HGCalTBSkiroc2CMS&);

#endif // SKIROC2DATAFRAME_H_INCLUDED
