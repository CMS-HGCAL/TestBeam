#include "HGCal/DataFormats/interface/HGCalTBSkiroc2CMS.h"

uint16_t HGCalTBSkiroc2CMS::gray_to_brady(const uint16_t gray) const
{
  uint16_t result = gray & (1 << 11);
  result |= (gray ^ (result >> 1)) & (1 << 10);
  result |= (gray ^ (result >> 1)) & (1 << 9);
  result |= (gray ^ (result >> 1)) & (1 << 8);
  result |= (gray ^ (result >> 1)) & (1 << 7);
  result |= (gray ^ (result >> 1)) & (1 << 6);
  result |= (gray ^ (result >> 1)) & (1 << 5);
  result |= (gray ^ (result >> 1)) & (1 << 4);
  result |= (gray ^ (result >> 1)) & (1 << 3);
  result |= (gray ^ (result >> 1)) & (1 << 2);
  result |= (gray ^ (result >> 1)) & (1 << 1);
  result |= (gray ^ (result >> 1)) & (1 << 0);
  return result;

}


bool HGCalTBSkiroc2CMS::check()
{
  for( size_t j=0; j<NUMBER_OF_CHANNELS; j++ ){
    uint16_t head=(m_data.at(j)&MASK_HEAD)>>4*3;
    if(head!=8&&head!=9){
      std::cout << "ISSUE : we expected 8(1000) or 9(1001) for the adc header and I find " << head << std::endl;
      return false;
    }
    for( size_t k=0; k<NUMBER_OF_SCA-1; k++){
      if( ((m_data.at(j+SCA_SHIFT*k)&MASK_HEAD)>>4*3)!=head ){
	std::cout << "\n We have a major issue (LG)-> " << head << " should be the same as " << ((m_data.at(j+SCA_SHIFT*k)&MASK_HEAD)>>4*3) << std::endl;
	return false;
      }
    }
    head=(m_data.at(j+NUMBER_OF_CHANNELS)&MASK_HEAD)>>4*3;
    if(head!=8&&head!=9){
      std::cout << "ISSUE : we expected 8(1000) or 9(1001) for the adc header and I find " << head << std::endl;
      return false;
    }
    for( size_t k=0; k<NUMBER_OF_SCA-1; k++){
      if( ((m_data.at(j+SCA_SHIFT*k+NUMBER_OF_CHANNELS)&MASK_HEAD)>>4*3)!=head ){
	std::cout << "\n We have a major issue (HG)-> " << head << " should be the same as " << ((m_data.at(j+SCA_SHIFT*k+NUMBER_OF_CHANNELS)&MASK_HEAD)>>4*3) << std::endl;
	return false;
      }
    }
  }
  return true;
}

std::ostream& operator<<(std::ostream& s, const HGCalTBSkiroc2CMS& ski)
{
  for (size_t i = 0; i < NUMBER_OF_CHANNELS; i++){
    s << "\n Channel det id : " << ski.detid(i) << "\n High gain ADC => " << ski.TOAHitRise(i) ;
    for( size_t j=0; j<NUMBER_OF_SCA-2; j++)
      s << " " << ski.ADCHigh(i,j) ;
    s << " " << ski.TOARise(i) ;
    s << " " << ski.TOTFast(i) ;
    s << "\n Low gain ADC => " << ski.TOAHitFall(i) ;
    for( size_t j=0; j<NUMBER_OF_SCA-2; j++)
      s << " " << ski.ADCLow(i,j) ;
    s << " " << ski.TOAFall(i) ;
    s << " " << ski.TOTSlow(i) ;
  }
  s << std::endl;
  return s;
}


  
