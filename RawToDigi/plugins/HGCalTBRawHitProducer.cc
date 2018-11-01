#include "HGCal/RawToDigi/plugins/HGCalTBRawHitProducer.h"
#include "HGCal/Geometry/interface/HGCalTBGeometryParameters.h"
#include <iostream>

HGCalTBRawHitProducer::HGCalTBRawHitProducer(const edm::ParameterSet& cfg) : 
  m_electronicMap(cfg.getUntrackedParameter<std::string>("ElectronicMap","HGCal/CondObjects/data/map_CERN_Hexaboard_OneLayers_May2017.txt")),
  m_outputCollectionName(cfg.getParameter<std::string>("OutputCollectionName")),
  m_subtractPedestal(cfg.getUntrackedParameter<bool>("SubtractPedestal",false)),
  m_maskNoisyChannels(cfg.getUntrackedParameter<bool>("MaskNoisyChannels",false)),
  m_pedestalHigh_filename(cfg.getUntrackedParameter<std::string>("HighGainPedestalFileName",std::string("pedestalHG.txt"))),
  m_pedestalLow_filename(cfg.getUntrackedParameter<std::string>("LowGainPedestalFileName",std::string("pedestalLG.txt"))),
  m_channelsToMask_filename(cfg.getUntrackedParameter<std::string>("ChannelsToMaskFileName","HGCal/CondObjects/data/noisyChannels.txt")),
  m_underSaturationADC(cfg.getUntrackedParameter<int>("UnderSaturationADC",4)),
  m_minTimeSampleForSaturation(cfg.getUntrackedParameter<int>("MinTimeSampleForSaturation",2)),
  m_maxTimeSampleForSaturation(cfg.getUntrackedParameter<int>("MaxTimeSampleForSaturation",4))
{
  m_HGCalTBSkiroc2CMSCollection = consumes<HGCalTBSkiroc2CMSCollection>(cfg.getParameter<edm::InputTag>("InputCollection"));
  produces <HGCalTBRawHitCollection>(m_outputCollectionName);
  std::cout << cfg.dump() << std::endl;
}

void HGCalTBRawHitProducer::beginJob()
{
  HGCalCondObjectTextIO io(0);
  edm::FileInPath fip(m_electronicMap);
  if (!io.load(fip.fullPath(), essource_.emap_)) {
    throw cms::Exception("Unable to load electronics map");
  };
  
  if( m_subtractPedestal ){
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
	  int skiId=HGCAL_TB_GEOMETRY::N_SKIROC_PER_HEXA*hexaboard+(HGCAL_TB_GEOMETRY::N_SKIROC_PER_HEXA-skiroc)%HGCAL_TB_GEOMETRY::N_SKIROC_PER_HEXA+1;
	  HGCalTBElectronicsId eid(skiId,channel);      
	  if (!essource_.emap_.existsEId(eid.rawId()))
	    ped.id = HGCalTBDetId(-1);
	  else
	    ped.id = essource_.emap_.eid2detId(eid);
	  index+=ptr;
	}else continue;
	for( unsigned int ii=0; ii<NUMBER_OF_SCA; ii++ ){
	  float mean,rms;
	  nval = sscanf( index, "%f %f %n",&mean,&rms,&ptr );
	  if( nval==2 ){
	    ped.pedHGMean[ii]=mean;
	    ped.pedHGRMS[ii]=rms;
	    index+=ptr;
	  }else continue;
	}
	m_pedMap.insert( std::pair<int,pedestalChannel>(10000*hexaboard+100*skiroc+channel,ped) );
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
	for( unsigned int ii=0; ii<NUMBER_OF_SCA; ii++ ){
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
  }
  
  if( m_maskNoisyChannels ){
    FILE* file;
    char buffer[300];
    //edm::FileInPath fip();
    file = fopen (m_channelsToMask_filename.c_str() , "r");
    if (file == NULL){ perror ("Error opening noisy channels file"); exit(1); }
    else{
      while ( ! feof (file) ){
	if ( fgets (buffer , 300 , file) == NULL ) break;
	const char* index = buffer;
	int layer,skiroc,channel,ptr,nval;
	nval=sscanf( index, "%d %d %d %n",&layer,&skiroc,&channel,&ptr );
	int skiId=HGCAL_TB_GEOMETRY::N_SKIROC_PER_HEXA*layer+(HGCAL_TB_GEOMETRY::N_SKIROC_PER_HEXA-skiroc)%HGCAL_TB_GEOMETRY::N_SKIROC_PER_HEXA+1;
	if( nval==3 ){
	  HGCalTBElectronicsId eid(skiId,channel);      
	  if (essource_.emap_.existsEId(eid.rawId()))
	    m_noisyChannels.push_back(eid.rawId());
	}else continue;
      }
      fclose (file);
    }
  }
}

void HGCalTBRawHitProducer::produce(edm::Event& event, const edm::EventSetup& iSetup)
{

  std::auto_ptr<HGCalTBRawHitCollection> hits(new HGCalTBRawHitCollection);

  edm::Handle<HGCalTBSkiroc2CMSCollection> skirocs;
  event.getByToken(m_HGCalTBSkiroc2CMSCollection, skirocs);

  if( !skirocs->size() ){
    event.put(hits, m_outputCollectionName);
    return;
  }
  
  for( size_t iski=0; iski<skirocs->size(); iski++ ){
    HGCalTBSkiroc2CMS skiroc=skirocs->at(iski);
    //if( !skiroc.check(1) ){
    //  std::cout << "!!! No hit will be created for this chip !!!" << std::endl;
    //  //continue;
    //}
    std::vector<int> rollpositions=skiroc.rollPositions();
    for( size_t ichan=0; ichan<HGCAL_TB_GEOMETRY::N_CHANNELS_PER_SKIROC; ichan++ ){
      int iboard=iski/HGCAL_TB_GEOMETRY::N_SKIROC_PER_HEXA;
      int iskiroc=iski%HGCAL_TB_GEOMETRY::N_SKIROC_PER_HEXA;
      int skiId=HGCAL_TB_GEOMETRY::N_SKIROC_PER_HEXA*iboard+(HGCAL_TB_GEOMETRY::N_SKIROC_PER_HEXA-iskiroc)%HGCAL_TB_GEOMETRY::N_SKIROC_PER_HEXA+1;
      HGCalTBElectronicsId eid(skiId,ichan);
      if( !essource_.emap_.existsEId(eid.rawId()) || std::find(m_noisyChannels.begin(),m_noisyChannels.end(),eid.rawId())!=m_noisyChannels.end() )
	continue;
      unsigned int rawid=skiroc.detid(ichan).rawId();
      std::vector<float> adchigh(NUMBER_OF_SCA,0);
      std::vector<float> adclow(NUMBER_OF_SCA,0);
      bool hgSat(false),lgSat(false);
      for( int it=0; it<NUMBER_OF_SCA; it++ ){
	if( it>=m_minTimeSampleForSaturation && it<=m_maxTimeSampleForSaturation && skiroc.ADCHigh(ichan,it)==m_underSaturationADC ) hgSat=true;
	if( it>=m_minTimeSampleForSaturation && it<=m_maxTimeSampleForSaturation && skiroc.ADCLow(ichan,it)==m_underSaturationADC ) lgSat=true;
	if( m_subtractPedestal ){
	  adchigh.at( rollpositions[it] ) = skiroc.ADCHigh(ichan,it)>4 ? skiroc.ADCHigh(ichan,it)-m_pedMap[10000*iboard+100*iskiroc+ichan].pedHGMean[it] : -10000;
	  adclow.at( rollpositions[it] ) = skiroc.ADCLow(ichan,it)>4 ? skiroc.ADCLow(ichan,it)-m_pedMap[10000*iboard+100*iskiroc+ichan].pedLGMean[it] : -10000;
	}
	else{
	  adchigh.at( rollpositions[it] ) = skiroc.ADCHigh(ichan,it) ;
	  adclow.at( rollpositions[it] ) = skiroc.ADCLow(ichan,it) ;
	}
      }
      for( int it=0; it<NUMBER_OF_SCA-NUMBER_OF_TIME_SAMPLES; it++ ){
      	adchigh.pop_back();
	adclow.pop_back();
      }
      HGCalTBRawHit hit(rawid, iski, ichan, adchigh, adclow,
			skiroc.TOARise(ichan), skiroc.TOAFall(ichan),
			skiroc.TOTSlow(ichan), skiroc.TOTFast(ichan));
      if(hgSat) hit.setUnderSaturationForHighGain();
      if(lgSat) hit.setUnderSaturationForLowGain();
      hits->push_back(hit);
    }
  }
  event.put(hits, m_outputCollectionName);
}

DEFINE_FWK_MODULE(HGCalTBRawHitProducer);
