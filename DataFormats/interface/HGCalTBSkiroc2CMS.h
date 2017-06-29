#ifndef HGCALTBSKIROC2CMS_H_INCLUDED
#define HGCALTBSKIROC2CMS_H_INCLUDED 1

#include <vector>
#include "HGCal/DataFormats/interface/HGCalTBDetId.h"
#include "HGCal/Geometry/interface/HGCalTBGeometryParameters.h"

static const int MASK_ADC = 0x0FFF;
static const int MASK_ROLL = 0x1FFF;
static const int MASK_GTS_MSB = 0x2FFF;
static const int MASK_GTS_LSB = 0x1FFF;
static const int MASK_ID = 0xFF;
static const int MASK_HEAD = 0xF000;

static const int NUMBER_OF_SCA = 13;
static const int ADCLOW_SHIFT = 0;
static const int ADCHIGH_SHIFT = 64;
static const int SCA_SHIFT = 128;
static const int SKIROC_DATA_SIZE = 1924; //number of 16 bits words

class HGCalTBSkiroc2CMS
{
 public:

 HGCalTBSkiroc2CMS( const std::vector<uint16_t> data, const std::vector<HGCalTBDetId> &ids ) :
  m_data(data),
  m_id(ids)
  {;}
  uint16_t gray_to_brady(const uint16_t gray) const;
    
  uint16_t ADCLow( int chan, int sca ) const {chan=HGCAL_TB_GEOMETRY::N_CHANNELS_PER_SKIROC-1-chan; sca=NUMBER_OF_SCA-1-sca; return (sca>=0 && sca<NUMBER_OF_SCA) ? gray_to_brady( m_data.at(chan+ADCLOW_SHIFT+SCA_SHIFT*sca) & MASK_ADC ) : 10000;}
  uint16_t ADCHigh( int chan, int sca ) const {chan=HGCAL_TB_GEOMETRY::N_CHANNELS_PER_SKIROC-1-chan; sca=NUMBER_OF_SCA-1-sca; return (sca>=0 && sca<NUMBER_OF_SCA) ? gray_to_brady( m_data.at(chan+ADCHIGH_SHIFT+SCA_SHIFT*sca) & MASK_ADC ) : 10000;}
  uint16_t TOTFast( int chan ) const {chan=HGCAL_TB_GEOMETRY::N_CHANNELS_PER_SKIROC-1-chan; return gray_to_brady( m_data.at(chan+ADCLOW_SHIFT+SCA_SHIFT*(NUMBER_OF_SCA+1)) & MASK_ADC );}
  uint16_t TOTSlow( int chan ) const {chan=HGCAL_TB_GEOMETRY::N_CHANNELS_PER_SKIROC-1-chan; return gray_to_brady( m_data.at(chan+ADCHIGH_SHIFT+SCA_SHIFT*(NUMBER_OF_SCA+1)) & MASK_ADC );}
  uint16_t TOAFall( int chan ) const {chan=HGCAL_TB_GEOMETRY::N_CHANNELS_PER_SKIROC-1-chan; return gray_to_brady( m_data.at(chan+ADCLOW_SHIFT+SCA_SHIFT*NUMBER_OF_SCA) & MASK_ADC );}
  uint16_t TOARise( int chan ) const {chan=HGCAL_TB_GEOMETRY::N_CHANNELS_PER_SKIROC-1-chan; return gray_to_brady( m_data.at(chan+ADCHIGH_SHIFT+SCA_SHIFT*NUMBER_OF_SCA) & MASK_ADC );}
  bool TOAHitFall(int chan) const {chan=HGCAL_TB_GEOMETRY::N_CHANNELS_PER_SKIROC-1-chan; return ((m_data.at(chan+ADCLOW_SHIFT)&~MASK_ADC)>>4*3)&0x1 ;}
  bool TOAHitRise(int chan) const {chan=HGCAL_TB_GEOMETRY::N_CHANNELS_PER_SKIROC-1-chan; return ((m_data.at(chan+ADCHIGH_SHIFT)&~MASK_ADC)>>4*3)&0x1 ;}
  uint16_t rollMask() const { return (m_data.at(SKIROC_DATA_SIZE-4)&MASK_ROLL); }
  std::vector<int> rollPositions() const ;

 private:
  uint32_t globalTS_MSB() const { return (m_data.at(SKIROC_DATA_SIZE-3) & MASK_GTS_MSB ); }
  uint32_t globalTS_LSB() const { return ( (m_data.at(SKIROC_DATA_SIZE-2) & MASK_GTS_LSB)>>1 ); }
  
 public:
  uint32_t globalTS()const{ return gray_to_brady( globalTS_MSB()<<12 | globalTS_LSB()); }
  int skirocId()const{ return (m_data.at(SKIROC_DATA_SIZE-1) & MASK_ID); }
  bool check();

  HGCalTBDetId detid( int chan ) const
  {
    return HGCalTBDetId(m_id.at(chan));
  }

 private:
  std::vector<uint16_t> m_data;
  std::vector<HGCalTBDetId> m_id;
};

std::ostream& operator<<(std::ostream&, HGCalTBSkiroc2CMS&);

#endif
