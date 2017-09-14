#include "HGCal/Reco/plugins/HGCalTBRecHitProducer.h"
#include "HGCal/Reco/interface/PulseFitter.h"
#include "HGCal/Reco/interface/CommonMode.h"
#include "HGCal/Geometry/interface/HGCalTBGeometryParameters.h"

#include <iostream>

//#define DEBUG

const static int SENSORSIZE = 128;

HGCalTBRecHitProducer::HGCalTBRecHitProducer(const edm::ParameterSet& cfg) : 
  m_outputCollectionName(cfg.getParameter<std::string>("OutputCollectionName")),
  m_electronicMap(cfg.getUntrackedParameter<std::string>("ElectronicsMap","HGCal/CondObjects/data/map_CERN_Hexaboard_28Layers_AllFlipped.txt")),
  m_timeSample3ADCCut(cfg.getUntrackedParameter<double>("TimeSample3ADCCut",15))
{
  m_HGCalTBRawHitCollection = consumes<HGCalTBRawHitCollection>(cfg.getParameter<edm::InputTag>("InputCollection"));

  produces <HGCalTBRecHitCollection>(m_outputCollectionName);
  std::vector<double> v0(1,10.);
  m_highGainADCSaturation = cfg.getUntrackedParameter<std::vector<double> >("HighGainADCSaturation",v0);
  std::vector<double> v1(1,10.);
  m_lowGainADCSaturation = cfg.getUntrackedParameter<std::vector<double> >("LowGainADCSaturation",v1);
  std::vector<double> v2(1,10.);
  m_LG2HG_value = cfg.getUntrackedParameter<std::vector<double> >("LG2HG",v2);
  std::vector<double> v3(1,10.);
  m_TOT2LG_value = cfg.getUntrackedParameter<std::vector<double> >("TOT2LG",v3);

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
}

void HGCalTBRecHitProducer::produce(edm::Event& event, const edm::EventSetup& iSetup)
{
  std::auto_ptr<HGCalTBRecHitCollection> rechits(new HGCalTBRecHitCollection);

  edm::Handle<HGCalTBRawHitCollection> rawhits;
  event.getByToken(m_HGCalTBRawHitCollection, rawhits);

  CommonMode cm(essource_.emap_); //default is common mode per chip using the median
  cm.Evaluate( rawhits );
  std::map<int,commonModeNoise> cmMap=cm.CommonModeNoiseMap();

  std::vector<std::pair<double, double> > CellXY;

  #ifdef DEBUG
    std::cout<<"++++++++++++++++++++++  Event: "<<eventCounter<<"   +++++++++++++++++++++"<<std::endl;
  #endif

  for( auto rawhit : *rawhits ){
    HGCalTBElectronicsId eid( essource_.emap_.detId2eid(rawhit.detid().rawId()) );
    if( !essource_.emap_.existsEId(eid) ) continue;
    int iboard=rawhit.detid().layer()-1;
    int iski=rawhit.skiroc();
    int ichannel=rawhit.channel();

    std::vector<double> sampleHG, sampleLG, sampleT;

    float highGain, lowGain, totGain;
    int hgStatus = -1;
    int lgStatus = -1;
    float timeHG = 0.;
    float timeLG = 0.;
    float subHG[NUMBER_OF_TIME_SAMPLES],subLG[NUMBER_OF_TIME_SAMPLES];

    totGain = rawhit.totSlow();

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
    for( int it=0; it<NUMBER_OF_TIME_SAMPLES; it++ ){
      sampleHG.push_back(rawhit.highGainADC(it)-subHG[it]);
      sampleLG.push_back(rawhit.lowGainADC(it)-subLG[it]);
      sampleT.push_back(25*it+12.5);
    }
    //this is a just try to isolate hits with signal
    float en2=rawhit.highGainADC(2)-subHG[2];
    float en3=rawhit.highGainADC(3)-subHG[3];
    float en4=rawhit.highGainADC(4)-subHG[4];
    float en6=rawhit.highGainADC(6)-subHG[6];
    
    if( en2<en3 && en3>en6 && en4>en6 && en3>m_timeSample3ADCCut){
      PulseFitter fitter(0,150);
      PulseFitterResult fithg;
      fitter.run(sampleT, sampleHG, fithg);
      PulseFitterResult fitlg;
      fitter.run(sampleT, sampleLG, fitlg);
      
      highGain = fithg.amplitude;
      lowGain = fitlg.amplitude;
      timeHG = fithg.tmax - fithg.trise;
      timeLG = fitlg.tmax - fitlg.trise;
      hgStatus = fithg.status;
      lgStatus = fitlg.status;
    }
    if(hgStatus != 0)
      highGain=0;
    if(lgStatus != 0)
      lowGain=0;
    
    float energy = -1;
    float time = -1.;
    HGCalTBRecHit recHit(rawhit.detid(), energy, lowGain, highGain, totGain, time);
    if(rawhit.highGainADC(3) < m_highGainADCSaturation.at(iboard) && hgStatus == 0){
      energy = highGain;
      time = timeHG;
      recHit.setFlag(HGCalTBRecHit::kGood);
    }     
    else if(rawhit.lowGainADC(3)-subHG[3] < m_lowGainADCSaturation.at(iboard) && lgStatus == 0){
      energy = lowGain * m_LG2HG_value.at(iboard);
      time = timeLG;
      recHit.setFlag(HGCalTBRecHit::kHighGainSaturated);
      recHit.setFlag(HGCalTBRecHit::kGood);
    }
    else if(totGain > 50){
      energy = totGain * m_TOT2LG_value.at(iboard) * m_LG2HG_value.at(iboard);
      recHit.setFlag(HGCalTBRecHit::kLowGainSaturated);
      recHit.setFlag(HGCalTBRecHit::kGood);
    }
    else {
      energy = 0;
    }
    if( rawhit.isUnderSaturationForHighGain() ) recHit.setUnderSaturationForHighGain();
    if( rawhit.isUnderSaturationForLowGain() ) recHit.setUnderSaturationForLowGain();

    recHit.setEnergy(energy);
    recHit.setTime(time);
    recHit.setChannel(ichannel);
    recHit.setSkiroc(iski);

    HGCalTBDetId detid = rawhit.detid();
    CellCentreXY = TheCell.GetCellCentreCoordinatesForPlots(detid.layer(), detid.sensorIU(), detid.sensorIV(), detid.iu(), detid.iv(), SENSORSIZE );
    double iux = (CellCentreXY.first < 0 ) ? (CellCentreXY.first + HGCAL_TB_GEOMETRY::DELTA) : (CellCentreXY.first - HGCAL_TB_GEOMETRY::DELTA);
    double iuy = (CellCentreXY.second < 0 ) ? (CellCentreXY.second + HGCAL_TB_GEOMETRY::DELTA) : (CellCentreXY.second - HGCAL_TB_GEOMETRY::DELTA);
    recHit.setCellCenterCoordinate(iux, iuy);

    #ifdef DEBUG
      if (totGain>1000) {
        std::cout<<"Board: "<<iboard+1<<"  and skiroc "<<iski;
        std::cout<<"  Channel: "<<ichannel<<std::endl;
        std::cout<<"subHG[3] = "<<subHG[3]<<"   subLG[3] = "<<subLG[3]<<std::endl;
        std::cout<<"energy: "<<energy<<"   lowGain (pulse): "<<lowGain<<"   sample3 LG: "<<sampleLG[3]<<"  highGain (pulse): "<<highGain<<"   sample3 HG: " <<sampleHG[3]<<"  totGain: "<<totGain<<std::endl<<std::endl;
      }
    #endif

    rechits->push_back(recHit);
  }
  event.put(rechits, m_outputCollectionName);
  #ifdef DEBUG
    std::cout<<std::endl<<std::endl;
    eventCounter++;
  #endif
}

DEFINE_FWK_MODULE(HGCalTBRecHitProducer);
