#include "HGCal/Reco/plugins/HGCalTBRecHitProducer.h"
#include "HGCal/Reco/interface/PulseFitter.h"
#include "HGCal/Reco/interface/CommonMode.h"
#include "HGCal/Geometry/interface/HGCalTBGeometryParameters.h"
#include "HGCal/CondObjects/interface/HGCalCondObjectROOTIO.h"

#include <iostream>
#include <ctime>

const static int SENSORSIZE = 128;

HGCalTBRecHitProducer::HGCalTBRecHitProducer(const edm::ParameterSet& cfg) : 
  m_outputCollectionName(cfg.getParameter<std::string>("OutputCollectionName")),
  m_electronicMap(cfg.getUntrackedParameter<std::string>("ElectronicsMap","HGCal/CondObjects/data/map_CERN_Hexaboard_28Layers_AllFlipped.txt")),
  m_detectorLayoutFile(cfg.getUntrackedParameter<std::string>("DetectorLayout","HGCal/CondObjects/data/layerGeom_oct2017_h2_17layers.txt")),
  m_adcCalibrationsFile(cfg.getUntrackedParameter<std::string>("ADCCalibrationsFile","ADCCalibrations.root")),
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
  if (!io.load(fip.fullPath(), m_emap)) {
    throw cms::Exception("Unable to load electronics map");
  };
  fip=edm::FileInPath(m_detectorLayoutFile);
  if (!io.load(fip.fullPath(), m_layout)) {
    throw cms::Exception("Unable to load detector layout file");
  };
  for( auto layer : m_layout.layers() )
    layer.print();
  
  HGCalCondObjectROOTIO rio(m_layout,m_emap);
  if (!rio.loadADCConversion(m_adcCalibrationsFile, m_adccalibmap)) {
    throw cms::Exception("Unable to load electronics map");
  };
  
  std::cout << "number of ADC calibrations = " << m_adccalibmap.getADCConversionsMap().size() << std::endl;
  std::cout << "number of electronic channels  = " << m_emap.size() << std::endl;

}

void HGCalTBRecHitProducer::produce(edm::Event& event, const edm::EventSetup& iSetup)
{
  std::auto_ptr<HGCalTBRecHitCollection> rechits(new HGCalTBRecHitCollection);

  edm::Handle<HGCalTBRawHitCollection> rawhits;
  event.getByToken(m_HGCalTBRawHitCollection, rawhits);

  CommonMode cm(m_emap); //default is common mode per chip using the median
  cm.Evaluate( rawhits );
  std::map<int,commonModeNoise> cmMap=cm.CommonModeNoiseMap();

  std::vector<std::pair<double, double> > CellXY;
  PulseFitter fitter(0);
  PulseFitterResult fitresultLG;
  PulseFitterResult fitresultHG;
  
  for( auto rawhit : *rawhits ){
    HGCalTBElectronicsId eid( m_emap.detId2eid(rawhit.detid().rawId()) );
    if( !m_emap.existsEId(eid) ) continue;
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
      sampleT.push_back(25*it);
    }
    HGCalTBDetId detid = rawhit.detid();
    HGCalTBLayer layer= m_layout.at(detid.layer()-1);

    float max_minus=rawhit.highGainADC(m_expectedMaxTimeSample-2)-subHG[m_expectedMaxTimeSample-2];
    float themax=rawhit.highGainADC(m_expectedMaxTimeSample)-subHG[m_expectedMaxTimeSample];
    float max_plus=rawhit.highGainADC(m_expectedMaxTimeSample+1)-subHG[m_expectedMaxTimeSample+1];
    float undershoot=rawhit.highGainADC(m_expectedMaxTimeSample+3)-subHG[m_expectedMaxTimeSample+3];
    if( themax>500||(max_minus<themax && themax>undershoot && max_plus>undershoot && themax>m_maxADCCut) ){
      fitter.run(sampleT, sampleLG, fitresultLG);
      fitter.run(sampleT, sampleHG, fitresultHG);

      lowgain=fitresultLG.amplitude;
      highgain=fitresultHG.amplitude;
      time=fitresultLG.tmax;
      
      bool correctCalib=m_adccalibmap.getCalibratedEnergy(detid,highgain,lowgain,totgain,energy);
      HGCalTBRecHit recHit(detid, energy, lowgain, highgain, totgain, time);

      if(!correctCalib)
      	std::cout << "Calibration for detid: " << detid << " is not correct" << std::endl;
      else
      	recHit.setFlag(HGCalTBRecHit::kGood);
      if( highgain>1600 )
      	recHit.setFlag(HGCalTBRecHit::kHighGainSaturated);
      if( lowgain>1100 )
      	recHit.setFlag(HGCalTBRecHit::kLowGainSaturated);

      CellCentreXY = TheCell.GetCellCentreCoordinatesForPlots(detid.layer(), detid.sensorIU(), detid.sensorIV(), detid.iu(), detid.iv(), SENSORSIZE );
      double iux = (CellCentreXY.first < 0 ) ? (CellCentreXY.first + HGCAL_TB_GEOMETRY::DELTA) : (CellCentreXY.first - HGCAL_TB_GEOMETRY::DELTA);
      double iuy = (CellCentreXY.second < 0 ) ? (CellCentreXY.second + HGCAL_TB_GEOMETRY::DELTA) : (CellCentreXY.second - HGCAL_TB_GEOMETRY::DELTA);
      recHit.setPosition( math::XYZPoint(iux,iuy,layer.z()) );
      recHit.setTime(time);
      rechits->push_back(recHit);
    }
  }
  event.put(rechits, m_outputCollectionName);
}

DEFINE_FWK_MODULE(HGCalTBRecHitProducer);
