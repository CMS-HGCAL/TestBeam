#include "HGCal/Reco/plugins/HGCalTBRecHitProducer.h"


const static int SENSORSIZE = 128;

HGCalTBRecHitProducer::HGCalTBRecHitProducer(const edm::ParameterSet& cfg) : 
  m_outputCollectionName(cfg.getParameter<std::string>("OutputCollectionName")),
  m_electronicMap(cfg.getUntrackedParameter<std::string>("ElectronicsMap","HGCal/CondObjects/data/map_CERN_Hexaboard_28Layers_AllFlipped.txt")),
  m_maskNoisyChannels(cfg.getUntrackedParameter<bool>("MaskNoisyChannels",false)),
  m_channelsToMask_filename(cfg.getUntrackedParameter<std::string>("ChannelsToMaskFileName","HGCal/CondObjects/data/noisyChannels.txt")),
  m_NHexaBoards(cfg.getUntrackedParameter<int>("NHexaBoards", 10)), 
  m_timeSample3ADCCut(cfg.getUntrackedParameter<double>("TimeSample3ADCCut",15))
{
  m_HGCalTBRawHitCollection = consumes<HGCalTBRawHitCollection>(cfg.getParameter<edm::InputTag>("InputCollection"));

  performPulseFit = cfg.getUntrackedParameter<bool>("performPulseFit", true);
  performAveraging = cfg.getUntrackedParameter<bool>("performAveraging", false);

  produces <HGCalTBRecHitCollection>(m_outputCollectionName);
  std::vector<double> v0(1,10.);
  m_ADC_per_MIP = cfg.getUntrackedParameter<std::vector<double> >("ADC_per_MIP",v0);
  std::vector<double> v1(1,10.);
  m_highGainADCSaturation = cfg.getUntrackedParameter<std::vector<double> >("HighGainADCSaturation",v1);
  std::vector<double> v2(1,10.);
  m_lowGainADCSaturation = cfg.getUntrackedParameter<std::vector<double> >("LowGainADCSaturation",v2);
  std::vector<double> v3(1,10.);
  m_LG2HG_value = cfg.getUntrackedParameter<std::vector<double> >("LG2HG",v3);
  std::vector<double> v4(1,10.);
  m_TOT2LG_value = cfg.getUntrackedParameter<std::vector<double> >("TOT2LG",v4);

  investigatePulseShape = cfg.getUntrackedParameter<bool>("investigatePulseShape", false);
  std::cout << cfg.dump() << std::endl;
}

void HGCalTBRecHitProducer::beginJob()
{
  HGCalCondObjectTextIO io(0);
  edm::FileInPath fip(m_electronicMap);
  if (!io.load(fip.fullPath(), essource_.emap_)) {
    throw cms::Exception("Unable to load electronics map");
  };
  #ifdef DEBUG
    eventCounter=1;
  #endif


  //usesResource("TFileService");
  edm::Service<TFileService> fs;

  std::ostringstream os( std::ostringstream::ate );
  for(int ib = 0; ib<m_NHexaBoards; ib++) {
    for( size_t iski=0; iski<HGCAL_TB_GEOMETRY::N_SKIROC_PER_HEXA; iski++ ){
      os.str("");os<<"HexaBoard"<<ib<<"_Skiroc"<<iski;
      TFileDirectory dir = fs->mkdir( os.str().c_str() );
      for( size_t ichan=0; ichan<HGCAL_TB_GEOMETRY::N_CHANNELS_PER_SKIROC; ichan++ ){
      if ((ichan % 2) == 1) continue;

      int key = ib * 10000 + iski * 100 + ichan;

      //std::cout<<"Creating key: "<<key<<std::endl;
      os.str("");os<<"Channel"<<ichan<<"__LGShape";
      shapesLG[key] = dir.make<TH2F>(os.str().c_str(),os.str().c_str(), 100, 0, 225, 100, -0.2, 1.);
      os.str("");os<<"Channel"<<ichan<<"__HGShape";
      shapesHG[key] = dir.make<TH2F>(os.str().c_str(),os.str().c_str(), 100, 0, 225, 100, -0.2, 1.);
      
      os.str("");os<<"Channel"<<ichan<<"__ToARiseVsTMaxLG";
      ToARisevsTMaxLG[key] = dir.make<TH2F>(os.str().c_str(),os.str().c_str(), 100, 50., 150., 100, 4., 3500.);
      os.str("");os<<"Channel"<<ichan<<"__ToARiseVsTMaxHG";
      ToARisevsTMaxHG[key] = dir.make<TH2F>(os.str().c_str(),os.str().c_str(), 100, 50., 150., 100, 4., 3500.);

      os.str("");os<<"Channel"<<ichan<<"__ToAFallVsTMaxLG";
      ToAFallvsTMaxLG[key] = dir.make<TH2F>(os.str().c_str(),os.str().c_str(), 100, 50., 150., 100, 4., 3500.);
      os.str("");os<<"Channel"<<ichan<<"__ToAFallVsTMaxHG";
      ToAFallvsTMaxHG[key] = dir.make<TH2F>(os.str().c_str(),os.str().c_str(), 100, 50., 150., 100, 4., 3500.);
      
      
      os.str("");os<<"Channel"<<ichan<<"__TMaxHGVsTMaxLG";
      TMaxHGvsTMaxLG[key] = dir.make<TH2F>(os.str().c_str(),os.str().c_str(), 100, 50., 150., 100, 50., 150.);
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
}

void HGCalTBRecHitProducer::produce(edm::Event& event, const edm::EventSetup& iSetup)
{
  std::auto_ptr<HGCalTBRecHitCollection> rechits(new HGCalTBRecHitCollection);

  edm::Handle<HGCalTBRawHitCollection> rawhits;
  event.getByToken(m_HGCalTBRawHitCollection, rawhits);

  CommonMode cm(essource_.emap_, true, true, -1); //default is common mode per chip using the median
  cm.Evaluate( rawhits );
  std::map<int,commonModeNoise> cmMap=cm.CommonModeNoiseMap();

  std::vector<std::pair<double, double> > CellXY;


  for( auto rawhit : *rawhits ){
    HGCalTBElectronicsId eid( essource_.emap_.detId2eid(rawhit.detid().rawId()) );
    if( !essource_.emap_.existsEId(eid) ) continue;
    int skiID=rawhit.skiroc();
    int iboard=(skiID-1)/HGCAL_TB_GEOMETRY::N_SKIROC_PER_HEXA;
    int ichannel=eid.ichan();
    int key = iboard * 10000 + ((skiID-1) % 4) * 100 + ichannel;

    if( !essource_.emap_.existsEId(eid.rawId()) || std::find(m_noisyChannels.begin(),m_noisyChannels.end(),eid.rawId())!=m_noisyChannels.end() )
      continue;

    std::vector<double> sampleHG, sampleLG, sampleT;

    float highGain(0), lowGain(0), totGain(0), toaRise(0), toaFall(0);

    int hgStatus = -1;
    int lgStatus = -1;
    float timeHG = 0.;
    float timeLG = 0.;
    float subHG[NUMBER_OF_TIME_SAMPLES],subLG[NUMBER_OF_TIME_SAMPLES];

    totGain = rawhit.totSlow();
    toaRise = rawhit.toaRise();
    toaFall = rawhit.toaFall();


    for( int it=0; it<NUMBER_OF_TIME_SAMPLES; it++ ){
      subHG[it]=0;
      subLG[it]=0;
    }
    switch ( rawhit.detid().cellType() ){
    case 0 : 
      for( int it=0; it<NUMBER_OF_TIME_SAMPLES; it++ ){
      	subHG[it]=cmMap[skiID].fullHG[it]; 
      	subLG[it]=cmMap[skiID].fullLG[it]; 
      }
      break;
    case 2 : 
      	for( int it=0; it<NUMBER_OF_TIME_SAMPLES; it++ ){
      	  subHG[it]=cmMap[skiID].halfHG[it]; 
      	  subLG[it]=cmMap[skiID].halfLG[it]; 
      	}
      	break;
    case 3 : 
      for( int it=0; it<NUMBER_OF_TIME_SAMPLES; it++ ){
      	subHG[it]=cmMap[skiID].mouseBiteHG[it]; 
      	subLG[it]=cmMap[skiID].mouseBiteLG[it]; 
      }
      break;
      case 4 : for( int it=0; it<NUMBER_OF_TIME_SAMPLES; it++ ){
      	subHG[it]=cmMap[skiID].outerHG[it]; 
      	subLG[it]=cmMap[skiID].outerLG[it]; 
      }
       	break;
    }

    
    for( int it=0; it<NUMBER_OF_TIME_SAMPLES; it++ ){
      sampleHG.push_back(rawhit.highGainADC(it)-subHG[it]);
      sampleLG.push_back((rawhit.lowGainADC(it)-subLG[it]) * m_LG2HG_value.at(skiID-1));
      sampleT.push_back(25*it+12.5);
    }
    
    if (performPulseFit) {    
      //pulse fitting
      PulseFitter fitter(0,150);

      //first HG
      //this is a just try to isolate hits with signal
      float en1=rawhit.highGainADC(1)-subHG[1];
      float en3=rawhit.highGainADC(3)-subHG[3];
      float en4=rawhit.highGainADC(4)-subHG[4];
      float en6=rawhit.highGainADC(6)-subHG[6];
      
      PulseFitterResult fithg;
      if( en1<en3 && en3>en6 && en4>en6 && en3>m_timeSample3ADCCut){
        fitter.run(sampleT, sampleHG, fithg);
        
        highGain = fithg.amplitude;
        timeHG = fithg.tmax - fithg.trise;
        hgStatus = fithg.status;
      }
      if(hgStatus != 0) {
        highGain=0;
        timeHG=-1;
      } else if (investigatePulseShape) {
        int key = iboard * 10000 + (skiID % 4) * 100 + ichannel;
        //std::cout<<"Filling HG key: "<<key<<std::endl;
        for( int it=0; it<NUMBER_OF_TIME_SAMPLES; it++) 
          shapesHG[key]->Fill(25*it+12.5-(fithg.tmax - fithg.trise), (rawhit.highGainADC(it)-subHG[it])/fithg.amplitude);

        ToARisevsTMaxHG[key]->Fill(fithg.tmax, toaRise);
        ToAFallvsTMaxHG[key]->Fill(fithg.tmax, toaFall);

      }

      //second LG
      //this is a just try to isolate hits with signal
      en1=(rawhit.lowGainADC(1)-subLG[1]);
      en3=(rawhit.lowGainADC(3)-subLG[3]);
      en4=(rawhit.lowGainADC(4)-subLG[4]);
      en6=(rawhit.lowGainADC(6)-subLG[6]);

      PulseFitterResult fitlg;
      if( en1<en3 && en3>en6 && en4>en6 && en3>m_timeSample3ADCCut){
        fitter.run(sampleT, sampleLG, fitlg);
        
        lowGain = fitlg.amplitude;
        timeLG = fitlg.tmax - fitlg.trise;
        lgStatus = fitlg.status;
      }

      if(lgStatus != 0) {
        lowGain=0;
        timeLG=-1;
      } else if (investigatePulseShape) {
        for( int it=0; it<NUMBER_OF_TIME_SAMPLES; it++) 
          shapesLG[key]->Fill(25*it+12.5-(fitlg.tmax - fitlg.trise), (rawhit.lowGainADC(it)-subLG[it])*m_LG2HG_value.at(skiID-1)/fitlg.amplitude);
        
        ToARisevsTMaxLG[key]->Fill(fitlg.tmax, toaRise);
        ToAFallvsTMaxLG[key]->Fill(fitlg.tmax, toaFall);
      }

      if (lgStatus != 0 && hgStatus != 0 && investigatePulseShape)
        TMaxHGvsTMaxLG[key]->Fill(fithg.tmax, fitlg.tmax);

    } else if (performAveraging) {
      //averaging of TS2-TS5
      float en2=rawhit.highGainADC(2)-subHG[2];
      float en3=rawhit.highGainADC(3)-subHG[3];
      float en4=rawhit.highGainADC(4)-subHG[4];
      float en5=rawhit.highGainADC(5)-subHG[5];
      hgStatus=0;
      timeHG = 50.;
      highGain = (en2+en3+en4+en5)/4.;

      en2=rawhit.lowGainADC(2)-subLG[2];
      en3=rawhit.lowGainADC(3)-subLG[3];
      en4=rawhit.lowGainADC(4)-subLG[4];
      en5=rawhit.lowGainADC(5)-subLG[5];
      lgStatus=0;
      timeLG = 50.;
      lowGain = (en2+en3+en4+en5)/4.;
    } else  {
      //simply use TS3
      hgStatus=0;
      timeHG = 50.;
      highGain = rawhit.highGainADC(3)-subHG[3];

      lgStatus=0;
      timeLG = 50.;
      lowGain = rawhit.lowGainADC(3)-subLG[3];

    }
    
    float energy = -1;
    float time = -1.;
    HGCalTBRecHit recHit(rawhit.detid(), energy, lowGain, highGain, totGain, time);
    if(rawhit.highGainADC(3) < m_highGainADCSaturation.at(skiID-1) && hgStatus == 0){
      energy = highGain;
      time = timeHG;
      recHit.setFlag(HGCalTBRecHit::kGood);
    }     
    else if(rawhit.lowGainADC(3)-subHG[3] < m_lowGainADCSaturation.at(skiID-1) && lgStatus == 0){
      energy = lowGain;
      time = timeLG;
      recHit.setFlag(HGCalTBRecHit::kHighGainSaturated);
      recHit.setFlag(HGCalTBRecHit::kGood);
    }
    else if(totGain > 10 && toaRise > 0){
      //std::cout<<"Setting totGain: "<<totGain<<std::endl;
      energy = totGain * m_TOT2LG_value.at(skiID-1) * m_LG2HG_value.at(skiID-1);
      recHit.setFlag(HGCalTBRecHit::kLowGainSaturated);
      recHit.setFlag(HGCalTBRecHit::kGood);
    }
    else {
      energy = 0;
    }
    if( rawhit.isUnderSaturationForHighGain() ) recHit.setUnderSaturationForHighGain();
    if( rawhit.isUnderSaturationForLowGain() ) recHit.setUnderSaturationForLowGain();

    energy = energy/m_ADC_per_MIP.at(skiID-1);
    recHit.setEnergy(energy);
    recHit.setTime(time);
    
    HGCalTBDetId detid = rawhit.detid();
    CellCentreXY = TheCell.GetCellCentreCoordinatesForPlots(detid.layer(), detid.sensorIU(), detid.sensorIV(), detid.iu(), detid.iv(), SENSORSIZE );
    double iux = (CellCentreXY.first < 0 ) ? (CellCentreXY.first + HGCAL_TB_GEOMETRY::DELTA) : (CellCentreXY.first - HGCAL_TB_GEOMETRY::DELTA);
    double iuy = (CellCentreXY.second < 0 ) ? (CellCentreXY.second + HGCAL_TB_GEOMETRY::DELTA) : (CellCentreXY.second - HGCAL_TB_GEOMETRY::DELTA);
    recHit.setCellCenterCoordinate(iux, iuy);
    
    #ifdef DEBUG
      if ((energy > 1000. ) && (iboard==0 || iboard==1)) {
        std::cout<< "Event: "<<eventCounter<<std::endl;
        std::cout<<"iboard: "<<iboard<<"  skiID: "<<skiID<<"   ichannel: "<<ichannel<<"   LG3: "<<rawhit.lowGainADC(3)<<"  - "<<subLG[3]<<"    HG3: "<<rawhit.highGainADC(3)<<"  - "<<subHG[3]<<std::endl;
      }
    #endif
    
    rechits->push_back(recHit);
  }
  event.put(rechits, m_outputCollectionName);
  #ifdef DEBUG
    eventCounter++;
  #endif
}

DEFINE_FWK_MODULE(HGCalTBRecHitProducer);
