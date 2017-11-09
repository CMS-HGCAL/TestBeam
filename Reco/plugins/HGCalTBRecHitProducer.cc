#include "HGCal/Reco/plugins/HGCalTBRecHitProducer.h"
#include "HGCal/Reco/interface/PulseFitter.h"
#include "HGCal/Reco/interface/CommonMode.h"
#include "HGCal/Geometry/interface/HGCalTBGeometryParameters.h"

#include <iostream>

const static int SENSORSIZE = 128;

HGCalTBRecHitProducer::HGCalTBRecHitProducer(const edm::ParameterSet& cfg) : 
  m_outputCollectionName(cfg.getParameter<std::string>("OutputCollectionName")),
  m_electronicMap(cfg.getUntrackedParameter<std::string>("ElectronicsMap","HGCal/CondObjects/data/map_CERN_Hexaboard_28Layers_AllFlipped.txt")),
  m_detectorLayoutFile(cfg.getUntrackedParameter<std::string>("DetectorLayout","HGCal/CondObjects/data/layerGeom_oct2017_h2_17layers.txt")),
  m_adcCalibrationsFile(cfg.getUntrackedParameter<std::string>("ADCCalibrations","HGCal/CondObjects/data/hgcal_calibration.txt")),
  m_timeSample3ADCCut(cfg.getUntrackedParameter<double>("TimeSample3ADCCut",15))
{
  m_HGCalTBRawHitCollection = consumes<HGCalTBRawHitCollection>(cfg.getParameter<edm::InputTag>("InputCollection"));
  produces <HGCalTBRecHitCollection>(m_outputCollectionName);
  std::cout << cfg.dump() << std::endl;
}

void HGCalTBRecHitProducer::beginJob()
{
  HGCalCondObjectTextIO io(0);
  edm::FileInPath fip(m_electronicMap);
  if (!io.load(fip.fullPath(), essource_.emap_)) {
    throw cms::Exception("Unable to load electronics map");
  };
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
  //  std::cout << essource_.adccalibmap_ << std::endl;
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

  PulseFitter fitter(0,150);
  PulseFitterResult fitresult;
  for( auto rawhit : *rawhits ){
    HGCalTBElectronicsId eid( essource_.emap_.detId2eid(rawhit.detid().rawId()) );
    if( !essource_.emap_.existsEId(eid) ) continue;
    int iski=eid.iskiroc();

    std::vector<double> sampleHG, sampleLG, sampleT;

    float highGain(0), lowGain(0), totGain(0), energy(0), time(-1);
    float subHG[NUMBER_OF_TIME_SAMPLES],subLG[NUMBER_OF_TIME_SAMPLES];

    totGain = rawhit.totSlow();

    switch ( rawhit.detid().cellType() ){
    default :
      for( int it=0; it<NUMBER_OF_TIME_SAMPLES; it++ ){
	subHG[it]=0;
	subLG[it]=0;
      }
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
    case 5 : for( int it=0; it<NUMBER_OF_TIME_SAMPLES; it++ ){
	subHG[it]=cmMap[iski].mergedHG[it]; 
	subLG[it]=cmMap[iski].mergedLG[it]; 
      }
      break;
    }
    for( int it=0; it<NUMBER_OF_TIME_SAMPLES; it++ ){
      sampleHG.push_back(rawhit.highGainADC(it)-subHG[it]);
      sampleLG.push_back(rawhit.lowGainADC(it)-subLG[it]);
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
    HGCalTBRecHit recHit(rawhit.detid(), energy, lowGain, highGain, totGain, time);
    if( rawhit.isUnderSaturationForHighGain() ) recHit.setUnderSaturationForHighGain();
    if( rawhit.isUnderSaturationForLowGain() ) recHit.setUnderSaturationForLowGain();
    HGCalTBDetId detid = rawhit.detid();
    CellCentreXY = TheCell.GetCellCentreCoordinatesForPlots(detid.layer(), detid.sensorIU(), detid.sensorIV(), detid.iu(), detid.iv(), SENSORSIZE );
    double iux = (CellCentreXY.first < 0 ) ? (CellCentreXY.first + HGCAL_TB_GEOMETRY::DELTA) : (CellCentreXY.first - HGCAL_TB_GEOMETRY::DELTA);
    double iuy = (CellCentreXY.second < 0 ) ? (CellCentreXY.second + HGCAL_TB_GEOMETRY::DELTA) : (CellCentreXY.second - HGCAL_TB_GEOMETRY::DELTA);
    recHit.setCellCenterCoordinate(iux, iuy);

    if( en1<en3 && en3>en6 && (en4>en6||en2>en6) && en3>m_timeSample3ADCCut){

      HGCalTBLayer layer= essource_.layout_.at(rawhit.detid().layer()-1);
      int moduleId= layer.at( recHit.id().sensorIU(),recHit.id().sensorIV() ).moduleID();
      iski = rawhit.skiroc()%4;
      ASIC_ADC_Conversions adcConv=essource_.adccalibmap_.getASICConversions(moduleId,iski);
      
      if( rawhit.lowGainADC(3) > adcConv.TOT_lowGain_transition() ){
	energy = totGain * adcConv.TOT_to_lowGain() * adcConv.lowGain_to_highGain();
	recHit.setFlag(HGCalTBRecHit::kLowGainSaturated);
	recHit.setFlag(HGCalTBRecHit::kGood);
      }
      else if( rawhit.highGainADC(3) > adcConv.lowGain_highGain_transition() ){
	fitter.run(sampleT, sampleLG, fitresult);
	recHit.setFlag(HGCalTBRecHit::kHighGainSaturated);
	if( fitresult.status==0 ){
	  energy = fitresult.amplitude * adcConv.lowGain_to_highGain();
	  time = fitresult.tmax - fitresult.trise;
	  recHit.setFlag(HGCalTBRecHit::kGood);
	}
      }
      else{
	fitter.run(sampleT, sampleHG, fitresult);
	if( fitresult.status==0 ){
	  energy = fitresult.amplitude;
	  time = fitresult.tmax - fitresult.trise;
	  recHit.setFlag(HGCalTBRecHit::kGood);
	}
      }
      recHit.setEnergy(energy*adcConv.adc_to_MIP());
      recHit.setTime(time);
    }
    rechits->push_back(recHit);
  }
  event.put(rechits, m_outputCollectionName);
}

DEFINE_FWK_MODULE(HGCalTBRecHitProducer);
