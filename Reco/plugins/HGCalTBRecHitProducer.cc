#include "HGCal/Reco/plugins/HGCalTBRecHitProducer.h"
#include "HGCal/Reco/interface/PulseFitter.h"
#include "HGCal/Reco/interface/CommonMode.h"
#include "HGCal/Geometry/interface/HGCalTBGeometryParameters.h"
//source before compilation: source /cvmfs/cms.cern.ch/slc7_amd64_gcc630/external/tensorflow-c/1.1.0/etc/profile.d/init.sh;

#include <iostream>


const static int SENSORSIZE = 128;

HGCalTBRecHitProducer::HGCalTBRecHitProducer(const edm::ParameterSet& cfg) : 
  m_CommonModeNoiseCollectionName(cfg.getUntrackedParameter<std::string>("CommonModeNoiseCollectionName", "HGCALTBCOMMONMODENOISEMAP")),
  m_outputCollectionName(cfg.getParameter<std::string>("OutputCollectionName")),
  m_electronicMap(cfg.getUntrackedParameter<std::string>("ElectronicsMap","HGCal/CondObjects/data/map_CERN_Hexaboard_28Layers_AllFlipped.txt")),
  m_detectorLayoutFile(cfg.getUntrackedParameter<std::string>("DetectorLayout","HGCal/CondObjects/data/layerGeom_oct2017_h2_17layers.txt")),
  m_adcCalibrationsFile(cfg.getUntrackedParameter<std::string>("ADCCalibrations","HGCal/CondObjects/data/hgcal_calibration.txt")),
  m_timeSample3ADCCut(cfg.getUntrackedParameter<double>("TimeSample3ADCCut",15)),
  m_maskNoisyChannels(cfg.getUntrackedParameter<bool>("MaskNoisyChannels",false)),
  m_channelsToMask_filename(cfg.getUntrackedParameter<std::string>("ChannelsToMaskFileName","HGCal/CondObjects/data/noisyChannels.txt")),
  m_NHexaBoards(cfg.getUntrackedParameter<int>("NHexaBoards", 10)) { 
  
  m_HGCalTBRawHitCollection = consumes<HGCalTBRawHitCollection>(cfg.getParameter<edm::InputTag>("InputCollection"));
  RunDataToken= consumes<RunData>(cfg.getParameter<edm::InputTag>("RUNDATA"));
  HGCalTBGlobalTimestampsToken= consumes<HGCalTBGlobalTimestamps>(cfg.getParameter<edm::InputTag>("GlobalTimestampCollectionName"));
  investigatePulseShape = cfg.getUntrackedParameter<bool>("investigatePulseShape", false);
  std::cout << cfg.dump() << std::endl;

  produces <std::map<int, commonModeNoise> >(m_CommonModeNoiseCollectionName);
  produces <HGCalTBRecHitCollection>(m_outputCollectionName);
}

void HGCalTBRecHitProducer::beginJob()
{
  HGCalCondObjectTextIO io(0);
  edm::FileInPath fip(m_electronicMap);
  if (!io.load(fip.fullPath(), essource_.emap_)) {
    throw cms::Exception("Unable to load electronics map");
  };

  edm::Service<TFileService> fs;


  std::ostringstream os( std::ostringstream::ate );
  for(int ib = 0; ib<m_NHexaBoards; ib++) {
    for( size_t iski=0; iski<HGCAL_TB_GEOMETRY::N_SKIROC_PER_HEXA; iski++ ){
      os.str("");os<<"HexaBoard"<<ib<<"_Skiroc"<<iski;
      TFileDirectory dir = fs->mkdir( os.str().c_str() );
      for( size_t ichan=0; ichan<HGCAL_TB_GEOMETRY::N_CHANNELS_PER_SKIROC; ichan++ ){
        if ((ichan % 2) == 1) continue;

        int key = ib * 10000 + iski * 100 + ichan;

        os.str("");os<<"Channel"<<ichan<<"__LGShape";
        shapesLG[key] = dir.make<TH2F>(os.str().c_str(),os.str().c_str(), 100, -75, 225, 150, -0.75, 1.);
        os.str("");os<<"Channel"<<ichan<<"__HGShape";
        shapesHG[key] = dir.make<TH2F>(os.str().c_str(),os.str().c_str(), 100, -75, 225, 150, -0.75, 1.);
        
      }
    }
  }

  if( m_maskNoisyChannels ){
    FILE* file;
    char buffer[300];
    //edm::FileInPath fip();
    file = fopen (m_channelsToMask_filename.c_str() , "r");
    if (file == NULL){
      perror ("Error opening noisy channels file"); exit(1); 
    } else{
    
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
        } else continue;
      }
    }
    fclose (file);
  }

  fip=edm::FileInPath(m_detectorLayoutFile);
  if (!io.load(fip.fullPath(), essource_.layout_)) {
    throw cms::Exception("Unable to load detector layout file");
  };
  for( auto layer : essource_.layout_.layers() )
    layer.print();
  
  fip=edm::FileInPath(m_adcCalibrationsFile);
  if (!io.load(fip.fullPath(), essource_.adccalibmap_)) {
    throw cms::Exception("Unable to load ADC conversions map");
  };

}

void HGCalTBRecHitProducer::produce(edm::Event& event, const edm::EventSetup& iSetup)
{
  edm::Handle<RunData> rd;
  //get the relevant event information
  event.getByToken(RunDataToken, rd);


  //edm::Handle<HGCalTBGlobalTimestamps> gts;
  ////event.getByToken(HGCalTBGlobalTimestampsToken, gts);

  std::unique_ptr<HGCalTBRecHitCollection> rechits(new HGCalTBRecHitCollection);

  edm::Handle<HGCalTBRawHitCollection> rawhits;
  event.getByToken(m_HGCalTBRawHitCollection, rawhits);

  //compute the common mode noise 
  CommonMode cm(essource_.emap_, true, true, -1); //default is common mode per chip using the median
  cm.Evaluate( rawhits );

  //..and at it to the edm event
  std::unique_ptr<std::map<int, commonModeNoise> > cmMap(new std::map<int, commonModeNoise>);
  (*cmMap)=cm.CommonModeNoiseMap();
  

  std::vector<std::pair<double, double> > CellXY;


  PulseFitter fitter(0,150);
  PulseFitterResult fitresultLG;
  PulseFitterResult fitresultHG;

  for( auto rawhit : *rawhits ){


    HGCalTBElectronicsId eid( essource_.emap_.detId2eid(rawhit.detid().rawId()) );
    if( !essource_.emap_.existsEId(eid.rawId()) || std::find(m_noisyChannels.begin(),m_noisyChannels.end(),eid.rawId())!=m_noisyChannels.end() )
      continue;
    
    int iski=rawhit.skiroc();
    int iboard=iski/HGCAL_TB_GEOMETRY::N_SKIROC_PER_HEXA;
    int ichannel=eid.ichan();
    int key = iboard * 10000 + (iski % 4) * 100 + ichannel;

    //uint32_t globalTimestamp = gts->skiroc_to_timestamps.at(iski);



    std::vector<double> sampleHG, sampleLG, sampleT;

    float highGain(0), lowGain(0), totGain(0),toaRise(0), toaFall(0), energy(0), _time(-1);

    float subHG[NUMBER_OF_TIME_SAMPLES],subLG[NUMBER_OF_TIME_SAMPLES];

    totGain = rawhit.totSlow();
    toaRise = rawhit.toaRise();
    toaFall = rawhit.toaFall();

    switch ( rawhit.detid().cellType() ){
    default :
      for( int it=0; it<NUMBER_OF_TIME_SAMPLES; it++ ){
      	subHG[it]=0;
      	subLG[it]=0;
      }
    case 0 : 
      for( int it=0; it<NUMBER_OF_TIME_SAMPLES; it++ ){
      	subHG[it]=cmMap->at(iski).fullHG[it]; 
      	subLG[it]=cmMap->at(iski).fullLG[it]; 
      }
      break;
    case 2 : 
      for( int it=0; it<NUMBER_OF_TIME_SAMPLES; it++ ){
      	subHG[it]=cmMap->at(iski).halfHG[it]; 
      	subLG[it]=cmMap->at(iski).halfLG[it]; 
      }
      break;
    case 3 : 
      for( int it=0; it<NUMBER_OF_TIME_SAMPLES; it++ ){
      	subHG[it]=cmMap->at(iski).mouseBiteHG[it]; 
      	subLG[it]=cmMap->at(iski).mouseBiteLG[it]; 
      }
      break;
    case 4 : for( int it=0; it<NUMBER_OF_TIME_SAMPLES; it++ ){
      	subHG[it]=cmMap->at(iski).outerHG[it]; 
      	subLG[it]=cmMap->at(iski).outerLG[it]; 
      }
      break;
    case 5 : for( int it=0; it<NUMBER_OF_TIME_SAMPLES; it++ ){
        subHG[it]=cmMap->at(iski).mergedHG[it]; 
        subLG[it]=cmMap->at(iski).mergedLG[it]; 
      }
      break;
    }

    
    for( int it=0; it<NUMBER_OF_TIME_SAMPLES; it++ ){
      sampleHG.push_back(rawhit.highGainADC(it)-subHG[it]);
      sampleLG.push_back((rawhit.lowGainADC(it)-subLG[it]));
      if( highGain<sampleHG[it] && it>1 && it<5 ) highGain=sampleHG[it];
      if( lowGain<sampleLG[it] && it>1 && it<5 ) lowGain=sampleLG[it];
      sampleT.push_back(25*it+12.5);
    }

    //this is a just try to isolate hits with signal
    float en1=sampleHG[1];
    float en2=sampleHG[2];
    float en3=sampleHG[3];
    float en4=sampleHG[4];
    float en6=sampleHG[6];
    if( en1<en3 && en3>en6 && (en4>en6||en2>en6) && en3>m_timeSample3ADCCut){
      HGCalTBRecHit recHit(rawhit.detid(), energy, lowGain, highGain, totGain, _time);
      if( rawhit.isUnderSaturationForHighGain() ) recHit.setUnderSaturationForHighGain();
      if( rawhit.isUnderSaturationForLowGain() ) recHit.setUnderSaturationForLowGain();
      HGCalTBDetId detid = rawhit.detid();
      CellCentreXY = TheCell.GetCellCentreCoordinatesForPlots(detid.layer(), detid.sensorIU(), detid.sensorIV(), detid.iu(), detid.iv(), SENSORSIZE );
      double iux = (CellCentreXY.first < 0 ) ? (CellCentreXY.first + HGCAL_TB_GEOMETRY::DELTA) : (CellCentreXY.first - HGCAL_TB_GEOMETRY::DELTA);
      double iuy = (CellCentreXY.second < 0 ) ? (CellCentreXY.second + HGCAL_TB_GEOMETRY::DELTA) : (CellCentreXY.second - HGCAL_TB_GEOMETRY::DELTA);
      recHit.setCellCenterCoordinate(iux, iuy);
        
      HGCalTBLayer layer= essource_.layout_.at(rawhit.detid().layer()-1);
      int moduleId= layer.at( recHit.id().sensorIU(),recHit.id().sensorIV() ).moduleID();
      iski = rawhit.skiroc()%4;
      ASIC_ADC_Conversions adcConv=essource_.adccalibmap_.getASICConversions(moduleId,iski);

      //fit in any case
      fitter.run(sampleT, sampleLG, fitresultLG);
      fitter.run(sampleT, sampleHG, fitresultHG);
      
      recHit.setEnergyTOT(totGain);

      if (fitresultLG.status==0) {
        recHit.setEnergyLow(fitresultLG.amplitude);
        recHit.setTimeMaxLG(fitresultLG.tmax - fitresultLG.trise);
      }
      if (fitresultHG.status==0) {
        recHit.setEnergyHigh(fitresultHG.amplitude);
        recHit.setTimeMaxHG(fitresultHG.tmax - fitresultHG.trise);
      }


      //switches between the gains for the rechit
      if( rawhit.lowGainADC(3) > adcConv.TOT_lowGain_transition() ){
        energy = totGain * adcConv.TOT_to_lowGain() * adcConv.lowGain_to_highGain();
        recHit.setFlag(HGCalTBRecHit::kLowGainSaturated);
        recHit.setFlag(HGCalTBRecHit::kGood);
      } else {

        

        if( rawhit.highGainADC(3) > adcConv.lowGain_highGain_transition() ){
          recHit.setFlag(HGCalTBRecHit::kHighGainSaturated);
          if( fitresultLG.status==0 ){
            energy = fitresultLG.amplitude * adcConv.lowGain_to_highGain();
            _time = fitresultLG.tmax - fitresultLG.trise;
            recHit.setFlag(HGCalTBRecHit::kGood);

            if (investigatePulseShape) {
              for( int it=0; it<NUMBER_OF_TIME_SAMPLES; it++) {
                shapesLG[key]->Fill(25*it+12.5-(fitresultLG.tmax - fitresultLG.trise), sampleLG[it]/fitresultLG.amplitude);
              }
            }
          }
        } else{
          if( fitresultHG.status==0 ){
            energy = fitresultHG.amplitude;
            _time = fitresultHG.tmax - fitresultHG.trise;
            recHit.setFlag(HGCalTBRecHit::kGood);

            if (investigatePulseShape) {
              for( int it=0; it<NUMBER_OF_TIME_SAMPLES; it++) {
                shapesHG[key]->Fill(25*it+12.5-(fitresultHG.tmax - fitresultHG.trise), sampleHG[it]/fitresultHG.amplitude);
              }
            }
          }
        }
      }

      //std::cout<<"Setting energy and time of: "<<detid.layer()<<"  "<<detid.iu()+7<<"  "<<detid.iv()+7<<"  "<<energy*adcConv.adc_to_MIP()<<"  "<<_time<<std::endl;
      recHit.setEnergy(energy*adcConv.adc_to_MIP());
      recHit.setTime(_time);


      recHit.setEnergyTSLow(sampleLG[2]>sampleLG[3] ? sampleLG[2] : sampleLG[3]);
      recHit.setEnergyTSHigh(sampleHG[2]>sampleHG[3] ? sampleHG[2] : sampleHG[3]);
      recHit.setToaRise(toaRise);
      recHit.setToaFall(toaFall);

      rechits->push_back(recHit);
    }

  }
  //event.put(std::move(cmMap), m_CommonModeNoiseCollectionName);
  event.put(std::move(rechits), m_outputCollectionName);
  #ifdef DEBUG
    eventCounter++;
  #endif
}

void HGCalTBRecHitProducer::endJob() {
}

DEFINE_FWK_MODULE(HGCalTBRecHitProducer);
