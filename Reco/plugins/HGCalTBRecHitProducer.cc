#include "HGCal/Reco/plugins/HGCalTBRecHitProducer.h"
#include <iostream>

const static size_t N_SKIROC_PER_HEXA = 4;

HGCalTBRecHitProducer::HGCalTBRecHitProducer(const edm::ParameterSet& cfg) : 
  m_pedestalHigh_filename(cfg.getParameter<std::string>("HighGainPedestalFileName")),
  m_pedestalLow_filename(cfg.getParameter<std::string>("LowGainPedestalFileName")),
  m_outputCollectionName(cfg.getParameter<std::string>("OutputCollectionName")),
  m_electronicMap(cfg.getUntrackedParameter<std::string>("ElectronicsMap","HGCal/CondObjects/data/map_CERN_Hexaboard_OneLayers_May2017.txt")),
  m_commonModeThreshold(cfg.getUntrackedParameter<double>("CommonModeThreshold",3)),
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

  HGCalCondObjectTextIO io(0);
  edm::FileInPath fip(m_electronicMap);
  if (!io.load(fip.fullPath(), essource_.emap_)) {
    throw cms::Exception("Unable to load electronics map");
  };

  std::cout << cfg.dump() << std::endl;
}

void HGCalTBRecHitProducer::beginJob()
{
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
	HGCalTBElectronicsId eid( (4-skiroc)%4+1 );
	if (!essource_.emap_.existsEId(eid.rawId()))
	  ped.id = HGCalTBDetId(-1);
	else
	  ped.id = essource_.emap_.eid2detId(eid);
	index+=ptr;
      }else continue;
      for( unsigned int ii=0; ii<NUMBER_OF_TIME_SAMPLES-0; ii++ ){
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
      for( unsigned int ii=0; ii<NUMBER_OF_TIME_SAMPLES-0; ii++ ){
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

void HGCalTBRecHitProducer::produce(edm::Event& event, const edm::EventSetup& iSetup)
{
  std::auto_ptr<HGCalTBRecHitCollection> rechits(new HGCalTBRecHitCollection);

  edm::Handle<HGCalTBRawHitCollection> rawhits;
  event.getByToken(m_HGCalTBRawHitCollection, rawhits);
  commonModeNoise cm[NUMBER_OF_TIME_SAMPLES][4];
  for( auto rawhit : *rawhits ){
    if( !essource_.emap_.existsDetId(rawhit.detid()) ) continue;
    HGCalTBElectronicsId eid( essource_.emap_.detId2eid(rawhit.detid()) );
    int iboard=(eid.iskiroc()-1)/N_SKIROC_PER_HEXA;
    int iski=eid.iskiroc()-1;
    int ichan=eid.ichan();
    if(m_keepOnlyTimeSample3){
      if( fabs(rawhit.highGainADC(3)-m_pedMap[10000*iboard+100*iski+ichan].pedHGMean[3])>m_commonModeThreshold*m_pedMap[10000*iboard+100*iski+ichan].pedHGRMS[3] ) continue;
      float highGain = rawhit.highGainADC(3)-m_pedMap[10000*iboard+100*iski+ichan].pedHGMean[3];
      float lowGain = rawhit.lowGainADC(3)-m_pedMap[10000*iboard+100*iski+ichan].pedLGMean[3];
      switch ( rawhit.detid().cellType() ){
      case 0 : cm[3][iski].fullHG += highGain; cm[3][iski].fullLG += lowGain; cm[3][iski].fullCounter++; break;
      case 2 : cm[3][iski].fullHG += highGain; cm[3][iski].fullLG += lowGain; cm[3][iski].fullCounter++; break;
      //case 2 : cm[3][iski].halfHG += highGain; cm[3][iski].halfLG += lowGain; cm[3][iski].halfCounter++; break;
      case 3 : cm[3][iski].mouseBiteHG += highGain; cm[3][iski].mouseBiteLG += lowGain; cm[3][iski].mouseBiteCounter++; break;
      case 4 : cm[3][iski].outerHG += highGain; cm[3][iski].outerLG += lowGain; cm[3][iski].outerCounter++; break;
      }
    }
    else{
      std::cout << "Should run with m_keepOnlyTimeSample3 sets to true, other method not yet implemented -> exit" << std::endl;
      exit(1);
    }
  }
  for( auto rawhit : *rawhits ){
    if( !essource_.emap_.existsDetId(rawhit.detid()) ) continue;
    HGCalTBElectronicsId eid( essource_.emap_.detId2eid(rawhit.detid()) );
    int iboard=(eid.iskiroc()-1)/N_SKIROC_PER_HEXA;
    int iski=eid.iskiroc()-1;
    int ichan=eid.ichan();
    if(m_keepOnlyTimeSample3){
      float highGain,lowGain;
      float subHG(0),subLG(0);
      switch ( rawhit.detid().cellType() ){
      case 0 : subHG=cm[3][iski].fullHG/cm[3][iski].fullCounter; subLG=cm[3][iski].fullLG/cm[3][iski].fullCounter; break;
      case 2 : subHG=cm[3][iski].fullHG/cm[3][iski].fullCounter; subLG=cm[3][iski].fullLG/cm[3][iski].fullCounter; break;
      //case 2 : subHG=cm[3][iski].halfHG/cm[3][iski].halfCounter; subLG=cm[3][iski].halfLG/cm[3][iski].halfCounter; break;
      case 3 : subHG=cm[3][iski].mouseBiteHG/cm[3][iski].mouseBiteCounter; subLG=cm[3][iski].mouseBiteLG/cm[3][iski].mouseBiteCounter; break;
      case 4 : subHG=cm[3][iski].outerHG/cm[3][iski].outerCounter; subLG=cm[3][iski].outerLG/cm[3][iski].outerCounter; break;
      }
      highGain=rawhit.highGainADC(3)-m_pedMap[10000*iboard+100*iski+ichan].pedHGMean[3]-subHG;
      lowGain=rawhit.lowGainADC(3)-m_pedMap[10000*iboard+100*iski+ichan].pedLGMean[3]-subLG;

      float energy = (highGain<m_highGainADCSaturation) ? highGain : lowGain*m_LG2HG_value.at(iboard);
      float time = rawhit.toaRise();

      HGCalTBRecHit recHit(rawhit.detid(), energy, lowGain, highGain, time);
      
      CellCenterXY = TheCell.GetCellCentreCoordinatesForPlots((recHit.id()).layer(), (recHit.id()).sensorIU(), (recHit.id()).sensorIV(), (recHit.id()).iu(), (recHit.id()).iv(), 128); 
      recHit.setCellCenterCoordinate(CellCenterXY.first, CellCenterXY.second);

      rechits->push_back( recHit );
    }
    else{
      std::cout << "Should run with m_keepOnlyTimeSample3 sets to true, other method not yet implemented -> exit" << std::endl;
      exit(1);
    }
  }
  event.put(rechits, m_outputCollectionName);
}

DEFINE_FWK_MODULE(HGCalTBRecHitProducer);
