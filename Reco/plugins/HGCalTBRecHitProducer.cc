#include "HGCal/Reco/plugins/HGCalTBRecHitProducer.h"
#include "HGCal/Reco/interface/CommonMode.h"
#include <iostream>

const static size_t N_SKIROC_PER_HEXA = 4;

HGCalTBRecHitProducer::HGCalTBRecHitProducer(const edm::ParameterSet& cfg) : 
  m_outputCollectionName(cfg.getParameter<std::string>("OutputCollectionName")),
  m_electronicMap(cfg.getUntrackedParameter<std::string>("ElectronicsMap","HGCal/CondObjects/data/map_CERN_Hexaboard_OneLayers_May2017.txt")),
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

  CommonMode cm(essource_.emap_); //default is common mode per chip using the median
  cm.Evaluate( rawhits );
  std::map<int,commonModeNoise> cmMap=cm.CommonModeNoiseMap();
  for( auto rawhit : *rawhits ){
    HGCalTBElectronicsId eid( essource_.emap_.detId2eid(rawhit.detid().rawId()) );
    if( !essource_.emap_.existsEId(eid) ) continue;
    int iboard=(eid.iskiroc()-1)/N_SKIROC_PER_HEXA;
    int iski=eid.iskiroc();
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
      	  subHG[it]=cmMap[iski].fullHG[it]; 
      	  subLG[it]=cmMap[iski].fullLG[it]; 
      	}
      	break;
      case 2 : 
      	for( int it=0; it<NUMBER_OF_TIME_SAMPLES; it++ ){
      	  subHG[it]=cmMap[iski].halfHG[it]; 
      	  subLG[it]=cmMap[iski].halfLG[it]; 
      	}
      	break;
      case 3 : 
      	for( int it=0; it<NUMBER_OF_TIME_SAMPLES; it++ ){
      	  subHG[it]=cmMap[iski].mouseBiteHG[it]; 
      	  subLG[it]=cmMap[iski].mouseBiteLG[it]; 
      	}
      	break;
      case 4 : for( int it=0; it<NUMBER_OF_TIME_SAMPLES; it++ ){
      	  subHG[it]=cmMap[iski].outerHG[it]; 
      	  subLG[it]=cmMap[iski].outerLG[it]; 
       	}
       	break;
      }
      //this is a just try to isolate hits with signal
      float en2=rawhit.highGainADC(2)-subHG[2];
      float en3=rawhit.highGainADC(3)-subHG[3];
      float en4=rawhit.highGainADC(4)-subHG[4];
      float en6=rawhit.highGainADC(6)-subHG[6];
      if( en2<en3 && en3>en6 && en4>en6 ){
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
