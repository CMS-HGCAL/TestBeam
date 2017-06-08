#include "HGCal/Reco/plugins/HGCalTBRecHitProducer.h"
#include <iostream>

const static size_t N_SKIROC_PER_HEXA = 4;

HGCalTBRecHitProducer::HGCalTBRecHitProducer(const edm::ParameterSet& cfg) : 
  m_outputCollectionName(cfg.getParameter<std::string>("OutputCollectionName")),
  m_electronicMap(cfg.getUntrackedParameter<std::string>("ElectronicsMap","HGCal/CondObjects/data/map_CERN_Hexaboard_OneLayers_May2017.txt")),
  m_commonModeThreshold(cfg.getUntrackedParameter<double>("CommonModeThreshold",100)),
  m_highGainADCSaturation(cfg.getUntrackedParameter<double>("HighGainADCSaturation",1800)),
  m_lowGainADCSaturation(cfg.getUntrackedParameter<double>("LowGainADCSaturation",1800)),
  m_keepOnlyTimeSample3(cfg.getUntrackedParameter<bool>("KeepOnlyTimeSample3",true))
{

  m_HGCalTBRawHitCollection = consumes<HGCalTBRawHitCollection>(cfg.getParameter<edm::InputTag>("InputCollection"));

  produces <HGCalTBRecHitCollection>(m_outputCollectionName);
  std::vector<double> v0(1,10.);
  m_LG2HG_value = cfg.getUntrackedParameter<std::vector<double> >("LG2HG",v0);
  std::vector<double> v1(1,10.);
  m_TOT2LG_value = cfg.getUntrackedParameter<std::vector<double> >("TOT2LG",v1);

  std::cout << cfg.dump() << std::endl;
}

void HGCalTBRecHitProducer::beginJob()
{
  HGCalCondObjectTextIO io(0);
  edm::FileInPath fip(m_electronicMap);
  if (!io.load(fip.fullPath(), essource_.emap_)) {
    throw cms::Exception("Unable to load electronics map");
  };
}

void HGCalTBRecHitProducer::produce(edm::Event& event, const edm::EventSetup& iSetup)
{
  std::auto_ptr<HGCalTBRecHitCollection> rechits(new HGCalTBRecHitCollection);

  edm::Handle<HGCalTBRawHitCollection> rawhits;
  event.getByToken(m_HGCalTBRawHitCollection, rawhits);
  commonModeNoise cm[NUMBER_OF_TIME_SAMPLES][N_SKIROC_PER_HEXA];
  for( auto rawhit : *rawhits ){
    if( !essource_.emap_.existsDetId(rawhit.detid()) ) continue;
    HGCalTBElectronicsId eid( essource_.emap_.detId2eid(rawhit.detid()) );
    int iboard=(eid.iskiroc()-1)/N_SKIROC_PER_HEXA;
    int iski=(N_SKIROC_PER_HEXA-(eid.iskiroc()-1))%N_SKIROC_PER_HEXA;//from 0 to 3
    iski+=iboard*N_SKIROC_PER_HEXA;
    for( int it=0; it<NUMBER_OF_TIME_SAMPLES; it++ ){
      if( rawhit.highGainADC(it)>m_commonModeThreshold ) continue;
      float highGain = rawhit.highGainADC(it);
      float lowGain = rawhit.lowGainADC(it);
      switch ( rawhit.detid().cellType() ){
      case 0 : cm[it][iski].fullHG += highGain; cm[it][iski].fullLG += lowGain; cm[it][iski].fullCounter++; break;
      case 2 : cm[it][iski].halfHG += highGain; cm[it][iski].halfLG += lowGain; cm[it][iski].halfCounter++; break;
      case 3 : cm[it][iski].mouseBiteHG += highGain; cm[it][iski].mouseBiteLG += lowGain; cm[it][iski].mouseBiteCounter++; break;
      case 4 : cm[it][iski].outerHG += highGain; cm[it][iski].outerLG += lowGain; cm[it][iski].outerCounter++; break;
      }
    }
  }
  for( auto rawhit : *rawhits ){
    if( !essource_.emap_.existsDetId(rawhit.detid()) ) continue;
    HGCalTBElectronicsId eid( essource_.emap_.detId2eid(rawhit.detid()) );
    int iboard=(eid.iskiroc()-1)/N_SKIROC_PER_HEXA;
    int iski=(N_SKIROC_PER_HEXA-(eid.iskiroc()-1))%N_SKIROC_PER_HEXA;//from 0 to 3
    iski+=iboard*N_SKIROC_PER_HEXA;
    if(m_keepOnlyTimeSample3){
      float highGain,lowGain;
      float subHG[NUMBER_OF_TIME_SAMPLES],subLG[NUMBER_OF_TIME_SAMPLES];
      for( int it=0; it<NUMBER_OF_TIME_SAMPLES; it++ ){
      	subHG[it]=0;
      	subLG[it]=0;
      }
      switch ( rawhit.detid().cellType() ){
      case 0 : 
      	for( int it=0; it<NUMBER_OF_TIME_SAMPLES; it++ ){
      	  subHG[it]=cm[it][iski].fullCounter>0 ? cm[it][iski].fullHG/cm[it][iski].fullCounter : 0; 
      	  subLG[it]=cm[it][iski].fullCounter>0 ? cm[it][iski].fullLG/cm[it][iski].fullCounter : 0; 
      	}
      	break;
      case 2 : 
      	for( int it=0; it<NUMBER_OF_TIME_SAMPLES; it++ ){
      	  subHG[it]=cm[it][iski].halfCounter>0 ? cm[it][iski].halfHG/cm[it][iski].halfCounter : 0; 
      	  subLG[it]=cm[it][iski].halfCounter>0 ? cm[it][iski].halfLG/cm[it][iski].halfCounter : 0; 
      	}
      	break;
      case 3 : 
      	for( int it=0; it<NUMBER_OF_TIME_SAMPLES; it++ ){
      	  subHG[it]=cm[it][iski].mouseBiteCounter>0 ? cm[it][iski].mouseBiteHG/cm[it][iski].mouseBiteCounter : 0; 
      	  subLG[it]=cm[it][iski].mouseBiteCounter>0 ? cm[it][iski].mouseBiteLG/cm[it][iski].mouseBiteCounter : 0; 
      	}
      	break;
      case 4 : for( int it=0; it<NUMBER_OF_TIME_SAMPLES; it++ ){
       	  subHG[it]=cm[it][iski].outerCounter>0 ? cm[it][iski].outerHG/cm[it][iski].outerCounter : 0; 
       	  subLG[it]=cm[it][iski].outerCounter>0 ? cm[it][iski].outerLG/cm[it][iski].outerCounter : 0; 
       	}
       	break;
      }
      //this is a just try to isolate hits with signal
      float en2=rawhit.highGainADC(2)-subHG[2];
      float en3=rawhit.highGainADC(3)-subHG[3];
      float en4=rawhit.highGainADC(4)-subHG[4];
      float en5=rawhit.highGainADC(5)-subHG[5];
      float en6=rawhit.highGainADC(6)-subHG[6];
      if( en2<en3 && en3>en5 && en4>en5 && (en5>en6||(en5<0&&en6<0)) ){
	highGain=rawhit.highGainADC(3)-subHG[3];
	lowGain=rawhit.lowGainADC(3)-subLG[3];
      }
      else{
       	highGain=-500;
       	lowGain=-500;
      }
      float energy = (highGain<m_highGainADCSaturation) ? highGain : lowGain*m_LG2HG_value.at(iboard);
      float time = rawhit.toaRise();
      rechits->push_back( HGCalTBRecHit(rawhit.detid(), energy, lowGain, highGain, time) );
    }
    else{
      std::cout << "Should run with m_keepOnlyTimeSample3 sets to true, other method not yet implemented -> exit" << std::endl;
      exit(1);
    }
  }
  event.put(rechits, m_outputCollectionName);
}

DEFINE_FWK_MODULE(HGCalTBRecHitProducer);
