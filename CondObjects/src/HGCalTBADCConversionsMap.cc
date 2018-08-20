#include "HGCal/CondObjects/interface/HGCalCondObjectTextIO.h"
#include <HGCal/CondObjects/interface/HGCalTBADCConversionsMap.h>
#include <algorithm>


float ADCConversions::totShape(float nmip)
{
  return nmip<m_params.totthr ? 0 : m_params.totcoeff*nmip+m_params.totped-m_params.totnorm/(std::pow(nmip,m_params.totpower)-m_params.totthr);
}

bool ADCConversions::getCalibEnergy(float hg, float lg, int tot, float &energy)
{
  bool correctCalib=false;
  energy=0; //MIP unit
  if( lg>m_params.lgsat && tot>10){ //LG saturated -> use ToT
    float nmin(0),nmax(2000);
    int index=0;
    float diff(2000);
    float centralValue(0);
    while(fabs(diff)>1e-2&&index<100){
      centralValue=nmin+fabs(nmax-nmin)/2;
      float totApprox=totShape(centralValue);
      diff=tot-totApprox;
      if( diff>0 )
	nmin=centralValue;
      else
	nmax=centralValue;
      index++;
    }
    if(centralValue>0&&index<1000){
      energy=centralValue;
      correctCalib=true;
    }
  }
  else if( hg>m_params.hgsat ){
    correctCalib=true;
    energy=lg*m_params.lgtomip;
  }
  else{
    correctCalib=true;
    energy=hg*m_params.hgtomip;
  }
  if(!correctCalib)
    energy=lg*m_params.lgtomip;
  return correctCalib;
}


std::ostream&operator<<(std::ostream& s, const ADCConversions& adc){
  ADCConversionParameters params=adc.getParameters();
  s << adc.getDetId()
    << params.hgtomip << " " << params.hgsat << " "
    << params.lgtomip << " " << params.lgsat << " "
    << params.totcoeff << " " << params.totthr << " "
    << params.totped << " " << params.totnorm << " " << params.totpower;
  return s;
}

void HGCalTBADCConversionsMap::addEntry(ADCConversions adcConv)
{
  if( std::find(m_vec.begin(),m_vec.end(),adcConv)==m_vec.end() )
    m_vec.push_back(adcConv);
  else{
    std::cout << "Error when loading ADC conversion parameters: same detid was found more than one time -> exit(1)" << std::endl;
    exit(1);
  }
}

ADCConversions HGCalTBADCConversionsMap::getADCConversions(HGCalTBDetId detid)
{
  ADCConversions tmp(detid);
  std::vector<ADCConversions>::iterator it=std::find(m_vec.begin(),m_vec.end(),tmp);
  if( it==m_vec.end() ){
    std::cout << "Error the detid : " << detid << " could not be found in the HGCalTBADCCalibrationMap -> exit(1)" << std::endl;
    exit(1);
  }
  else return (*it);
}

bool HGCalTBADCConversionsMap::getCalibratedEnergy(HGCalTBDetId detid, float hg, float lg, int tot, float &energy)
{
  return getADCConversions(detid).getCalibEnergy(hg,lg,tot,energy);
}
