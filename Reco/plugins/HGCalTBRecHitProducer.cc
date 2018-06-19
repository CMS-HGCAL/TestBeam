#include "HGCal/Reco/plugins/HGCalTBRecHitProducer.h"
#include "HGCal/Reco/interface/PulseFitter.h"
#include "HGCal/Reco/interface/CommonMode.h"
#include "HGCal/Geometry/interface/HGCalTBGeometryParameters.h"

#include <iostream>
#include <ctime>

const static int SENSORSIZE = 128;

HGCalTBRecHitProducer::HGCalTBRecHitProducer(const edm::ParameterSet& cfg) : 
  m_outputCollectionName(cfg.getParameter<std::string>("OutputCollectionName")),
  m_electronicMap(cfg.getUntrackedParameter<std::string>("ElectronicsMap","HGCal/CondObjects/data/map_CERN_Hexaboard_28Layers_AllFlipped.txt")),
  m_detectorLayoutFile(cfg.getUntrackedParameter<std::string>("DetectorLayout","HGCal/CondObjects/data/layerGeom_oct2017_h2_17layers.txt")),
  m_adcCalibrationsFile(cfg.getUntrackedParameter<std::string>("ADCCalibrations","HGCal/CondObjects/data/hgcal_calibration.txt")),
  m_expectedMaxTimeSample(cfg.getUntrackedParameter<int>("ExpectedMaxTimeSample",3)),
  m_maxADCCut(cfg.getUntrackedParameter<double>("MaxADCCut",15))
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
  PulseFitter fitter(0);
  PulseFitterResult fitresultLG;
  PulseFitterResult fitresultHG;
  for( auto rawhit : *rawhits ){
    HGCalTBElectronicsId eid( essource_.emap_.detId2eid(rawhit.detid().rawId()) );
    if( !essource_.emap_.existsEId(eid) ) continue;
    int iski=rawhit.skiroc();

    std::vector<double> sampleHG, sampleLG, sampleT;

    float highgain(0), lowgain(0), totgain(0), energy(0), time(-1);
    float subHG[NUMBER_OF_TIME_SAMPLES],subLG[NUMBER_OF_TIME_SAMPLES];

    totgain = rawhit.totSlow();

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
      sampleT.push_back(25*it+12.5);
    }
    // if( rawhit.isUnderSaturationForHighGain() ) recHit.setUnderSaturationForHighGain();
    // if( rawhit.isUnderSaturationForLowGain() ) recHit.setUnderSaturationForLowGain();
    HGCalTBDetId detid = rawhit.detid();
    HGCalTBLayer layer= essource_.layout_.at(detid.layer()-1);

    float max_minus=rawhit.highGainADC(m_expectedMaxTimeSample-2)-subHG[m_expectedMaxTimeSample-2];
    float themax=rawhit.highGainADC(m_expectedMaxTimeSample)-subHG[m_expectedMaxTimeSample];
    float max_plus=rawhit.highGainADC(m_expectedMaxTimeSample+1)-subHG[m_expectedMaxTimeSample+1];
    float undershoot=rawhit.highGainADC(m_expectedMaxTimeSample+3)-subHG[m_expectedMaxTimeSample+3];
    if( themax>500||(max_minus<themax && themax>undershoot && max_plus>undershoot && themax>m_maxADCCut) ){
      //int moduleId= layer.at( detid.sensorIU(),detid.sensorIV() ).moduleID();
      iski = rawhit.skiroc()%4;
      //ASIC_ADC_Conversions adcConv=essource_.adccalibmap_.getASICConversions(moduleId,iski);
      fitter.run(sampleT, sampleLG, fitresultLG);
      fitter.run(sampleT, sampleHG, fitresultHG);

      lowgain=fitresultLG.amplitude;
      highgain=fitresultHG.amplitude;
      time=fitresultLG.tmax;

      HGCalTBRecHit recHit(rawhit.detid(), energy, lowgain, highgain, totgain, time);
      CellCentreXY = TheCell.GetCellCentreCoordinatesForPlots(detid.layer(), detid.sensorIU(), detid.sensorIV(), detid.iu(), detid.iv(), SENSORSIZE );
      double iux = (CellCentreXY.first < 0 ) ? (CellCentreXY.first + HGCAL_TB_GEOMETRY::DELTA) : (CellCentreXY.first - HGCAL_TB_GEOMETRY::DELTA);
      double iuy = (CellCentreXY.second < 0 ) ? (CellCentreXY.second + HGCAL_TB_GEOMETRY::DELTA) : (CellCentreXY.second - HGCAL_TB_GEOMETRY::DELTA);
      recHit.setPosition( math::XYZPoint(iux,iuy,layer.z()) );
      // if( rawhit.lowGainADC(3) > adcConv.TOT_lowGain_transition() ){
      // 	energy = totgain * adcConv.TOT_to_lowGain() * adcConv.lowGain_to_highGain();
      // 	recHit.setFlag(HGCalTBRecHit::kLowGainSaturated);
      // 	recHit.setFlag(HGCalTBRecHit::kGood);
      // }
      // else if( rawhit.highGainADC(3) > adcConv.lowGain_highGain_transition() ){
      // 	recHit.setFlag(HGCalTBRecHit::kHighGainSaturated);
      // 	if( fitresultLG.status==0 ){
      // 	  energy = lowgain * adcConv.lowGain_to_highGain();
      // 	  recHit.setFlag(HGCalTBRecHit::kGood);
      // 	}
      // }
      // else{
      // 	if( fitresultHG.status==0 ){
      // 	  energy = highgain;
      // 	  recHit.setFlag(HGCalTBRecHit::kGood);
      // 	}
      // }
      // recHit.setEnergy(energy*adcConv.adc_to_MIP());
      recHit.setTime(time);
      rechits->push_back(recHit);
    }
  }
  event.put(rechits, m_outputCollectionName);
}

DEFINE_FWK_MODULE(HGCalTBRecHitProducer);
