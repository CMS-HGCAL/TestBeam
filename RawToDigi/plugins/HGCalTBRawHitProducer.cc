#include "HGCal/RawToDigi/plugins/HGCalTBRawHitProducer.h"
#include <iostream>


HGCalTBRawHitProducer::HGCalTBRawHitProducer(const edm::ParameterSet& cfg) : 
  m_outputCollectionName(cfg.getParameter<std::string>("OutputCollectionName")),
  m_electronicMap(cfg.getUntrackedParameter<std::string>("ElectronicMap","HGCal/CondObjects/data/map_CERN_Hexaboard_OneLayers_May2017.txt")),
  m_pedestalLow_filename(cfg.getParameter<std::string>("pedestalLow")),
  m_pedestalHigh_filename(cfg.getParameter<std::string>("pedestalHigh"))
{
  m_HGCalTBSkiroc2CMSCollection = consumes<HGCalTBSkiroc2CMSCollection>(cfg.getParameter<edm::InputTag>("InputCollection"));
  produces <HGCalTBRawHitCollection>(m_outputCollectionName);

  HGCalCondObjectTextIO io(0);
  edm::FileInPath fip(m_electronicMap);
  if (!io.load(fip.fullPath(), essource_.emap_)) {
    throw cms::Exception("Unable to load electronics map");
  };

  std::cout << cfg.dump() << std::endl;
}

void HGCalTBRawHitProducer::beginJob()
{
  //maybe not the correct place to start pedestal subtraction, but OK this is a try
  FILE* file;
  char buffer[300];
  file = fopen (m_pedestalHigh_filename.c_str() , "r");
  if (file == NULL){ perror ("Error opening pedestal high gain"); exit(1); }
  else{
    while ( ! feof (file) ){
       if ( fgets (buffer , 300 , file) == NULL ) break;
       pedestalChannel ped;
       const char* index = buffer;
       int hexaboard,skiroc,channel,ptr,nval;
       nval=sscanf( index, "%d %d %d %n",&hexaboard,&skiroc,&channel,&ptr );
       if( nval==3 ){
	 HGCalTBElectronicsId eid(4-skiroc+1,channel);
	 if (!essource_.emap_.existsEId(eid.rawId()))
	   ped.id = HGCalTBDetId(-1);
	 else
	   ped.id = essource_.emap_.eid2detId(eid);
	 //std::cout << hexaboard << " " << skiroc << " " << channel << " " ;
	 index+=ptr;
       }else continue;
       for( int ii=0; ii<NUMBER_OF_SCA; ii++ ){
	 float mean,rms;
	 nval = sscanf( index, "%f %f %n",&mean,&rms,&ptr );
	 if( nval==2 ){
	   ped.pedHGMean[ii]=mean;
	   ped.pedHGRMS[ii]=rms;
	   //std::cout << mean << " " << rms << " " ;
	   index+=ptr;
	 }else continue;
       }
       m_pedMap.insert( std::pair<int,pedestalChannel>(10000*hexaboard+100*skiroc+channel,ped) );
       //std::cout << std::endl;
    }
    fclose (file);
  }
  
  file = fopen (m_pedestalLow_filename.c_str() , "r");
  if (file == NULL){ perror ("Error opening pedestal low gain"); exit(1); }
  else{
    while ( ! feof (file) ){
       if ( fgets (buffer , 300 , file) == NULL ) break;
       pedestalChannel ped;
       const char* index = buffer;
       int hexaboard,skiroc,channel,ptr,nval,key;
       nval=sscanf( index, "%d %d %d %n",&hexaboard,&skiroc,&channel,&ptr );
       if( nval==3 ){
	 key=10000*hexaboard+100*skiroc+channel;
	 index+=ptr;
       }else continue;
       for( int ii=0; ii<NUMBER_OF_SCA; ii++ ){
	 float mean,rms;
	 nval = sscanf( index, "%f %f %n",&mean,&rms,&ptr );
	 if( nval==2 ){
	   m_pedMap[key].pedLGMean[ii]=mean;
	   m_pedMap[key].pedLGRMS[ii]=rms;
	   index+=ptr;
	 }else continue;
       }
    }
    fclose (file);
  }
  for( std::map<int,pedestalChannel>::iterator it=m_pedMap.begin(); it!=m_pedMap.end(); ++it ){
    std::cout << it->first/10000 << " " << (it->first%10000)/100 << " " << it->first%100 ;
    for( int ii=0; ii<NUMBER_OF_SCA; ii++ )
      std::cout << " " << it->second.pedHGMean[ii] << " " << it->second.pedHGRMS[ii] ;
    std::cout << std::endl;
    for( int ii=0; ii<NUMBER_OF_SCA; ii++ )
      std::cout << " " << it->second.pedLGMean[ii] << " " << it->second.pedLGRMS[ii] ;
    std::cout << std::endl;
  }
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
      for( int it=0; it<NUMBER_OF_TIME_SAMPLES; it++ ){//the two last time samples are probably not good as the on track sca, may be we should not keep them?
      	adchigh.push_back( skiroc.ADCHigh(ichan,rollpositions[it])-m_pedMap[10000*iski/4+100*iski+ichan].pedHGMean[it] );
      	adclow.push_back( skiroc.ADCLow(ichan,rollpositions[it])-m_pedMap[10000*iski/4+100*iski+ichan].pedLGMean[it] );
      }
      HGCalTBRawHit hit( rawid, adchigh, adclow);
      hits->push_back(hit);
    }
  }
  event.put(hits, m_outputCollectionName);
}

DEFINE_FWK_MODULE(HGCalTBRawHitProducer);
