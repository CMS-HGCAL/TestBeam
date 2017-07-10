#include "HGCal/DataFormats/interface/HGCalTBSkiroc2CMS.h"
#include <bitset>
#include <cstring>
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
  for( size_t j=0; j<HGCAL_TB_GEOMETRY::N_CHANNELS_PER_SKIROC; j++ ){
    uint16_t head=(m_data.at(j)&MASK_HEAD)>>4*3;
    if(head!=8&&head!=9){
      std::cout << "ISSUE : we expected 8(1000) or 9(1001) for the adc header and I find " << head << std::endl;
      return false;
    }
    for( size_t k=0; k<NUMBER_OF_SCA+1; k++){
      if( ((m_data.at(j+SCA_SHIFT*k)&MASK_HEAD)>>4*3)!=head ){
	std::cout << "\n We have a major issue (LG)-> " << head << " should be the same as " << ((m_data.at(j+SCA_SHIFT*k)&MASK_HEAD)>>4*3) << std::endl;
	return false;
      }
    }
    head=(m_data.at(j+HGCAL_TB_GEOMETRY::N_CHANNELS_PER_SKIROC)&MASK_HEAD)>>4*3;
    if(head!=8&&head!=9){
      std::cout << "ISSUE : we expected 8(1000) or 9(1001) for the adc header and I find " << head << std::endl;
      return false;
    }
    for( size_t k=0; k<NUMBER_OF_SCA+1; k++){
      if( ((m_data.at(j+SCA_SHIFT*k+HGCAL_TB_GEOMETRY::N_CHANNELS_PER_SKIROC)&MASK_HEAD)>>4*3)!=head ){
	std::cout << "\n We have a major issue (HG)-> " << head << " should be the same as " << ((m_data.at(j+SCA_SHIFT*k+HGCAL_TB_GEOMETRY::N_CHANNELS_PER_SKIROC)&MASK_HEAD)>>4*3) << std::endl;
	return false;
      }
    }
  }
  return true;
}

std::vector<int> HGCalTBSkiroc2CMS::rollPositions() const
{
  std::bitset<NUMBER_OF_SCA> bitstmp=rollMask();
  std::bitset<NUMBER_OF_SCA> bits;
  for( size_t i=0; i<NUMBER_OF_SCA; i++ )
    bits[i]=bitstmp[12-i];
  std::vector<int> rollpositions(NUMBER_OF_SCA);
  if(bits.test(0)&&bits.test(12)){
    rollpositions[0]=12;
    for(size_t i=1; i<NUMBER_OF_SCA; i++)
      rollpositions[i]=(i-1);
  }
  else{
    int pos_trk1 = -1;
    for(size_t i=0; i<NUMBER_OF_SCA; i++)
      if(bits.test(i)){
	pos_trk1 = i;
	break;
      }
    for(size_t i=0; i<NUMBER_OF_SCA; i++)
      if( (int)i <= pos_trk1 + 1 )
	rollpositions[i]=i + 12 - (pos_trk1 + 1);
      else
	rollpositions[i]=i - 1 - (pos_trk1 + 1);
  }
  return rollpositions;
}

std::ostream& operator<<(std::ostream& s, HGCalTBSkiroc2CMS& ski)
{
  std::vector<int> rollpositions=ski.rollPositions();
  std::cout << "rollMask = " << ski.rollMask() << "\t";
  for(size_t i=0; i<13; i++)
    std::cout << rollpositions[i] << " ";
  std::cout << "\n";
  for (size_t i = 0; i < HGCAL_TB_GEOMETRY::N_CHANNELS_PER_SKIROC; i++){
    s << "\n Channel det id : " << ski.detid(i) << "\n High gain ADC => " << ski.TOAHitRise(i) ;
    for( size_t j=0; j<NUMBER_OF_SCA; j++)
      s << " " << ski.ADCHigh(i,j) ;
    s << " " << ski.TOARise(i) ;
    s << " " << ski.TOTFast(i) ;
    s << "\n Low gain ADC => " << ski.TOAHitFall(i) ;
    for( size_t j=0; j<NUMBER_OF_SCA; j++)
      s << " " << ski.ADCLow(i,j) ;
    s << " " << ski.TOAFall(i) ;
    s << " " << ski.TOTSlow(i) ;
  }
  s << std::endl;
  return s;
}
