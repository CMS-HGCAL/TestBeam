#include "HGCal/RawToDigi/plugins/HGCalTBRawHitProducer.h"
#include <iostream>


HGCalTBRawHitProducer::HGCalTBRawHitProducer(const edm::ParameterSet& cfg) : 
  m_outputCollectionName(cfg.getParameter<std::string>("OutputCollectionName"))
{
  m_HGCalTBSkiroc2CMSCollection = consumes<HGCalTBSkiroc2CMSCollection>(cfg.getParameter<edm::InputTag>("InputCollection"));
  produces <HGCalTBRawHitCollection>(m_outputCollectionName);
  std::cout << cfg.dump() << std::endl;
}

void HGCalTBRawHitProducer::beginJob()
{
}

void HGCalTBRawHitProducer::produce(edm::Event& event, const edm::EventSetup& iSetup)
{

  std::auto_ptr<HGCalTBRawHitCollection> hits(new HGCalTBRawHitCollection);

  edm::Handle<HGCalTBSkiroc2CMSCollection> skirocs;
  event.getByToken(m_HGCalTBSkiroc2CMSCollection, skirocs);
  for( size_t iski=0; iski<skirocs->size(); iski++ ){
    for( int ichan=0; ichan<NUMBER_OF_CHANNELS; ichan++ ){
      HGCalTBSkiroc2CMS skiroc=skirocs->at(iski);
      unsigned int rawid=skiroc.detid(ichan).rawId();
      std::vector<float> adchigh;
      std::vector<float> adclow;
      std::vector<int> rollpositions=skiroc.rollPositions();
      for( int it=0; it<NUMBER_OF_TIME_SAMPLES; it++ ){
      	adchigh.push_back( skiroc.ADCHigh(ichan,rollpositions[it]) );
      	adclow.push_back( skiroc.ADCLow(ichan,rollpositions[it]) );
      }
      HGCalTBRawHit hit( rawid, iski, ichan, adchigh, adclow);
      hits->push_back(hit);
    }
  }
  event.put(hits, m_outputCollectionName);
}

DEFINE_FWK_MODULE(HGCalTBRawHitProducer);
