#include "HGCal/Reco/plugins/HGCalTBRecHitProducer.h"
#include "HGCal/Reco/interface/PulseFitter.h"
#include "HGCal/Reco/interface/CommonMode.h"
#include <iostream>

const static size_t N_SKIROC_PER_HEXA = 4;
static const double deltaCellBoundary = 0.00001;
const static int SENSORSIZE = 128;


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

  std::vector<std::pair<double, double> > CellXY;

  for( auto rawhit : *rawhits ){
    HGCalTBElectronicsId eid( essource_.emap_.detId2eid(rawhit.detid().rawId()) );
    if( !essource_.emap_.existsEId(eid) ) continue;
    int iboard=(eid.iskiroc()-1)/N_SKIROC_PER_HEXA;
    int iski=eid.iskiroc();

    std::vector<double> sampleHG, sampleLG, sampleT;

    float highGain = -500, lowGain = -500, totGain;
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
      sampleT.push_back(25*it);
    }
    //this is a just try to isolate hits with signal
    float en2=rawhit.highGainADC(2)-subHG[2];
    float en3=rawhit.highGainADC(3)-subHG[3];
    float en4=rawhit.highGainADC(4)-subHG[4];
    float en6=rawhit.highGainADC(6)-subHG[6];
    
    if(!m_keepOnlyTimeSample3 && 
       (en2<en3 && en3>en6 && en4>en6 && en3>20)){
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
    if(m_keepOnlyTimeSample3 && 
       (en2<en3 && en3>en6 && en4>en6)){
      highGain=rawhit.highGainADC(3)-subHG[3];
      lowGain=rawhit.lowGainADC(3)-subLG[3];
      hgStatus = 3;
      lgStatus = 3;
    }
    if(hgStatus == -1 && lgStatus == -1){
      highGain=-500;
      lowGain=-500;
    }
    
    float energy = -1;
    float time = -1.;
    HGCalTBRecHit recHit(rawhit.detid(), energy, lowGain, highGain, totGain, time);
    if(highGain < m_highGainADCSaturation && hgStatus == 0){
      energy = highGain;
      time = timeHG;
      recHit.setFlag(HGCalTBRecHit::kGood);
    }     
    else if(lowGain < m_lowGainADCSaturation && hgStatus == 0 && lgStatus == 0){
      energy = lowGain * m_LG2HG_value.at(iboard);
      time = timeLG;
      recHit.setFlag(HGCalTBRecHit::kHighGainSaturated);
      recHit.setFlag(HGCalTBRecHit::kGood);
    }
    else if(hgStatus == 0 && lgStatus == 0 && totGain > 50){
      energy = totGain * m_TOT2LG_value.at(iboard) * m_LG2HG_value.at(iboard);
      recHit.setFlag(HGCalTBRecHit::kLowGainSaturated);
      recHit.setFlag(HGCalTBRecHit::kGood);
    }
    else if(hgStatus == 3 && lgStatus == 3){
      energy = (highGain<m_highGainADCSaturation) ? highGain : lowGain * m_LG2HG_value.at(iboard);
      recHit.setFlag(HGCalTBRecHit::kThirdSample);
    }
    else {
      energy = -500;
    }


    recHit.setEnergy(energy);
    recHit.setTime(time);

    HGCalTBDetId detid = rawhit.detid();
    CellCentreXY = TheCell.GetCellCentreCoordinatesForPlots(detid.layer(), detid.sensorIU(), detid.sensorIV(), detid.iu(), detid.iv(), SENSORSIZE );
    double iux = (CellCentreXY.first < 0 ) ? (CellCentreXY.first + deltaCellBoundary) : (CellCentreXY.first - deltaCellBoundary);
    double iuy = (CellCentreXY.second < 0 ) ? (CellCentreXY.second + deltaCellBoundary) : (CellCentreXY.second - deltaCellBoundary);
    recHit.setCellCenterCoordinate(iux, iuy);

    rechits->push_back(recHit);
  }
  event.put(rechits, m_outputCollectionName);
}

DEFINE_FWK_MODULE(HGCalTBRecHitProducer);
