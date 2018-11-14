#include "HGCal/Reco/plugins/HGCalTBRecHitProducer.h"
#include "HGCal/Reco/interface/PulseFitter.h"
#include "HGCal/Reco/interface/CommonMode.h"
#include "HGCal/Geometry/interface/HGCalTBGeometryParameters.h"

#include <iostream>
#include <ctime>

const static int SENSORSIZE = 128;

HGCalTBRecHitProducer::HGCalTBRecHitProducer(const edm::ParameterSet& cfg) :
  m_outputCollectionName(cfg.getParameter<std::string>("OutputCollectionName")),
  m_electronicMap(cfg.getUntrackedParameter<std::string>("ElectronicsMap", "HGCal/CondObjects/data/map_CERN_Hexaboard_28Layers_AllFlipped.txt")),
  m_detectorLayoutFile(cfg.getUntrackedParameter<std::string>("DetectorLayout", "HGCal/CondObjects/data/layerGeom_oct2017_h2_17layers.txt")),
  m_adcCalibrationsFile(cfg.getUntrackedParameter<std::string>("ADCCalibrations", "HGCal/CondObjects/data/hgcal_calibration.txt")),
  m_calibrationPerChannel(cfg.getUntrackedParameter<bool>("calibrationPerChannel", false)),
  m_expectedMaxTimeSample(cfg.getUntrackedParameter<int>("ExpectedMaxTimeSample", 3)),
  m_maxADCCut(cfg.getUntrackedParameter<double>("MaxADCCut", 15)),
  m_subtractCommonMode(cfg.getUntrackedParameter<bool>("subtractCommonMode", true)),
  m_TSForCommonModeNoiseSubtraction(cfg.getUntrackedParameter<int>("TSForCommonModeNoiseSubtraction", -1)), //-1: use all TS
  m_preselectionMethod(cfg.getUntrackedParameter<std::string>("preselectionMethod", "TB2018"))
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
  fip = edm::FileInPath(m_detectorLayoutFile);
  if (!io.load(fip.fullPath(), essource_.layout_)) {
    throw cms::Exception("Unable to load detector layout file");
  };
  for ( auto layer : essource_.layout_.layers() )
    layer.print();


  if (!m_calibrationPerChannel) {
    fip = edm::FileInPath(m_adcCalibrationsFile);
    if (!io.load(fip.fullPath(), essource_.adccalibmap_)) {
      throw cms::Exception("Unable to load ADC conversions map");
    };

  }
  else {
    fip = edm::FileInPath(m_adcCalibrationsFile);
    if (!io.load(fip.fullPath(), essource_.adccalibmap_perchannel_)) {
      throw cms::Exception("Unable to load ADC conversions map");
    };
  }

  if (m_preselectionMethod == "TB2017") _preselectionMethod = TB2017;
  else if (m_preselectionMethod == "TB2018") _preselectionMethod = TB2018;
  else _preselectionMethod = NONE;
  //  std::cout << essource_.adccalibmap_ << std::endl;

  edm::Service<TFileService> fs;


  std::ostringstream os( std::ostringstream::ate );
  for (int ib = 0; ib < HGCAL_TB_GEOMETRY::NUMBER_OF_HEXABOARD; ib++) {
    for ( size_t iski = 0; iski < HGCAL_TB_GEOMETRY::N_SKIROC_PER_HEXA; iski++ ) {
      os.str(""); os << "HexaBoard" << ib << "_Skiroc" << iski;
      TFileDirectory dir = fs->mkdir( os.str().c_str() );
      for ( size_t ichan = 0; ichan < HGCAL_TB_GEOMETRY::N_CHANNELS_PER_SKIROC; ichan++ ) {
        if ((ichan % 2) == 1) continue;

        int key = ib * 10000 + iski * 100 + ichan;

        os.str(""); os << "Channel" << ichan << "__LGShape";
        shapesLG[key] = dir.make<TH2F>(os.str().c_str(), os.str().c_str(), 100, -75, 225, 150, -0.75, 1.);
        os.str(""); os << "Channel" << ichan << "__HGShape";
        shapesHG[key] = dir.make<TH2F>(os.str().c_str(), os.str().c_str(), 100, -75, 225, 150, -0.75, 1.);

      }
    }
  }

}

void HGCalTBRecHitProducer::produce(edm::Event& event, const edm::EventSetup& iSetup)
{
  std::unique_ptr<HGCalTBRecHitCollection> rechits(new HGCalTBRecHitCollection);

  edm::Handle<HGCalTBRawHitCollection> rawhits;
  event.getByToken(m_HGCalTBRawHitCollection, rawhits);

  CommonMode cm(essource_.emap_); //default is common mode per chip using the median, here: median per layer
  if (m_subtractCommonMode) cm.Evaluate( rawhits );
  std::map<int, commonModeNoise> cmMap = cm.CommonModeNoiseMap();
  

  std::vector<std::pair<double, double> > CellXY;
  PulseFitter fitter(0);
  for ( auto rawhit : *rawhits ) {
    PulseFitterResult fitresultLG;
    PulseFitterResult fitresultHG;
    HGCalTBElectronicsId eid( essource_.emap_.detId2eid(rawhit.detid().rawId()) );
    if ( !essource_.emap_.existsEId(eid) ) continue;
    int iski = rawhit.skiroc();
    int iboard = iski / HGCAL_TB_GEOMETRY::N_SKIROC_PER_HEXA;
    int ichannel = eid.ichan();
    int key = iboard * 10000 + (iski % 4) * 100 + ichannel;

    std::vector<Float16_t> sampleHG, sampleLG, sampleT;

    Float16_t highgain(0), lowgain(0),  energy(-1), time(-1);
    unsigned int short totgain(0);
    Float16_t subHG[NUMBER_OF_TIME_SAMPLES], subLG[NUMBER_OF_TIME_SAMPLES];

    totgain = rawhit.totSlow();
    unsigned int short toaRise = rawhit.toaRise();
    unsigned int short toaFall = rawhit.toaFall();

    switch ( rawhit.detid().cellType() ) {
    default :
      for ( int it = 0; it < NUMBER_OF_TIME_SAMPLES; it++ ) {
        subHG[it] = 0;
        subLG[it] = 0;
      }
    case 0 :
      for ( int it = 0; it < NUMBER_OF_TIME_SAMPLES; it++ ) {
        subHG[it] = m_subtractCommonMode ? cmMap[iski].fullHG[(m_TSForCommonModeNoiseSubtraction==-1) ? it : m_TSForCommonModeNoiseSubtraction] : 0.;
        subLG[it] = m_subtractCommonMode ? cmMap[iski].fullLG[(m_TSForCommonModeNoiseSubtraction==-1) ? it : m_TSForCommonModeNoiseSubtraction] : 0.;
      }
      break;
    case 2 :
      for ( int it = 0; it < NUMBER_OF_TIME_SAMPLES; it++ ) {
        subHG[it] = m_subtractCommonMode ? cmMap[iski].halfHG[(m_TSForCommonModeNoiseSubtraction==-1) ? it : m_TSForCommonModeNoiseSubtraction] : 0.;
        subLG[it] = m_subtractCommonMode ? cmMap[iski].halfLG[(m_TSForCommonModeNoiseSubtraction==-1) ? it : m_TSForCommonModeNoiseSubtraction] : 0.;
      }
      break;
    case 3 :
      for ( int it = 0; it < NUMBER_OF_TIME_SAMPLES; it++ ) {
        subHG[it] = m_subtractCommonMode ? cmMap[iski].mouseBiteHG[(m_TSForCommonModeNoiseSubtraction==-1) ? it : m_TSForCommonModeNoiseSubtraction] : 0.;
        subLG[it] = m_subtractCommonMode ? cmMap[iski].mouseBiteLG[(m_TSForCommonModeNoiseSubtraction==-1) ? it : m_TSForCommonModeNoiseSubtraction] : 0.;
      }
      break;
    case 4 : for ( int it = 0; it < NUMBER_OF_TIME_SAMPLES; it++ ) {
        subHG[it] = m_subtractCommonMode ? cmMap[iski].outerHG[(m_TSForCommonModeNoiseSubtraction==-1) ? it : m_TSForCommonModeNoiseSubtraction] : 0.;
        subLG[it] = m_subtractCommonMode ? cmMap[iski].outerLG[(m_TSForCommonModeNoiseSubtraction==-1) ? it : m_TSForCommonModeNoiseSubtraction] : 0.;
      }
      break;
    case 5 : for ( int it = 0; it < NUMBER_OF_TIME_SAMPLES; it++ ) {
        subHG[it] = m_subtractCommonMode ? cmMap[iski].mergedHG[(m_TSForCommonModeNoiseSubtraction==-1) ? it : m_TSForCommonModeNoiseSubtraction] : 0.;
        subLG[it] = m_subtractCommonMode ? cmMap[iski].mergedLG[(m_TSForCommonModeNoiseSubtraction==-1) ? it : m_TSForCommonModeNoiseSubtraction] : 0.;
      }
      break;
    }
    for ( int it = 0; it < NUMBER_OF_TIME_SAMPLES; it++ ) {
      sampleHG.push_back(rawhit.highGainADC(it) - subHG[it]);
      sampleLG.push_back(rawhit.lowGainADC(it) - subLG[it]);
      sampleT.push_back(25 * it + 12.5);
    }
    // if( rawhit.isUnderSaturationForHighGain() ) recHit.setUnderSaturationForHighGain();
    // if( rawhit.isUnderSaturationForLowGain() ) recHit.setUnderSaturationForLowGain();
    HGCalTBDetId detid = rawhit.detid();
    HGCalTBLayer layer = essource_.layout_.at(detid.layer() - 1);

    bool passPreselection = false;

    if (_preselectionMethod == TB2018) {
      Float16_t max_minus = rawhit.highGainADC(m_expectedMaxTimeSample - 2) - subHG[m_expectedMaxTimeSample - 2];
      Float16_t themax = rawhit.highGainADC(m_expectedMaxTimeSample) - subHG[m_expectedMaxTimeSample];
      Float16_t max_plus = rawhit.highGainADC(m_expectedMaxTimeSample + 1) - subHG[m_expectedMaxTimeSample + 1];
      Float16_t undershoot = rawhit.highGainADC(m_expectedMaxTimeSample + 3) - subHG[m_expectedMaxTimeSample + 3];
      passPreselection = ( themax > 500 || (max_minus < themax && themax > undershoot && max_plus > undershoot && themax > m_maxADCCut) );
    } else if (_preselectionMethod == TB2017) {
      float en1 = sampleHG[1];
      float en2 = sampleHG[2];
      float en3 = sampleHG[3];
      float en4 = sampleHG[4];
      float en6 = sampleHG[6];
      passPreselection = ( en1 < en3 && en3 > en6 && (en4 > en6 || en2 > en6) && en3 > m_maxADCCut);
    } else {
      passPreselection = true;
    }


    if (passPreselection) {
      int moduleId = layer.at( detid.sensorIU(), detid.sensorIV() ).moduleID();
      iski = rawhit.skiroc() % 4;

      fitter.run(sampleT, sampleLG, fitresultLG);
      fitter.run(sampleT, sampleHG, fitresultHG);
      lowgain = fitresultLG.amplitude;
      highgain = fitresultHG.amplitude;
      time = -1;
      if ( fitresultLG.status == 0 )for ( int it = 0; it < NUMBER_OF_TIME_SAMPLES; it++) shapesLG[key]->Fill(25 * it + 12.5 - (fitresultLG.tmax - fitresultLG.trise), sampleLG[it] / fitresultLG.amplitude);
      if ( fitresultHG.status == 0 )for ( int it = 0; it < NUMBER_OF_TIME_SAMPLES; it++) shapesHG[key]->Fill(25 * it + 12.5 - (fitresultHG.tmax - fitresultHG.trise), sampleHG[it] / fitresultHG.amplitude);


      HGCalTBRecHit recHit(rawhit.detid(), energy, lowgain, highgain, totgain, time);
      CellCentreXY = TheCell.GetCellCentreCoordinatesForPlots(detid.layer(), detid.sensorIU(), detid.sensorIV(), detid.iu(), detid.iv(), SENSORSIZE );
      Float16_t iux = (CellCentreXY.first < 0 ) ? (CellCentreXY.first + HGCAL_TB_GEOMETRY::DELTA) : (CellCentreXY.first - HGCAL_TB_GEOMETRY::DELTA);
      Float16_t iuy = (CellCentreXY.second < 0 ) ? (CellCentreXY.second + HGCAL_TB_GEOMETRY::DELTA) : (CellCentreXY.second - HGCAL_TB_GEOMETRY::DELTA);
      recHit.setCellCenterCoordinate(iux, iuy);

      recHit.setTimeMaxLG(fitresultLG.tmax);
      recHit.setTimeMaxHG(fitresultHG.tmax);
      recHit.setEnergyTSLow(sampleLG[m_expectedMaxTimeSample - 1], sampleLG[m_expectedMaxTimeSample]);
      recHit.setEnergyTSHigh(sampleHG[m_expectedMaxTimeSample - 1], sampleHG[m_expectedMaxTimeSample]);
      recHit.setToaRise(toaRise);
      recHit.setToaFall(toaFall);

      //copy the noise flag
      recHit.setNoiseFlag(rawhit.getNoiseFlag());

      //energy conversion default
      if (!m_calibrationPerChannel) {
        ASIC_ADC_Conversions adcConv = essource_.adccalibmap_.getASICConversions(moduleId, iski);
        if ( rawhit.lowGainADC(3) > adcConv.TOT_lowGain_transition() ) {
          energy = totgain * adcConv.TOT_to_lowGain() * adcConv.lowGain_to_highGain();
          recHit.setFlag(HGCalTBRecHit::kLowGainSaturated);
          recHit.setFlag(HGCalTBRecHit::kGood);
        }
        else if ( rawhit.highGainADC(3) > adcConv.lowGain_highGain_transition() ) {
          recHit.setFlag(HGCalTBRecHit::kHighGainSaturated);
          if ( fitresultLG.status == 0 ) {
            energy = lowgain * adcConv.lowGain_to_highGain();
            recHit.setFlag(HGCalTBRecHit::kGood);
          }
        }
        else {
          if ( fitresultHG.status == 0 ) {
            energy = highgain;
            recHit.setFlag(HGCalTBRecHit::kGood);
          }
        }
        energy *= adcConv.adc_to_MIP();

      } else {

        ASIC_ADC_Conversions_perChannel adcConv = essource_.adccalibmap_perchannel_.getASICConversions(moduleId, iski, ichannel);

        if (fitresultLG.status != 0) recHit.setFlag(HGCalTBRecHit::kLGFitFailed);
        if (fitresultHG.status != 0) recHit.setFlag(HGCalTBRecHit::kHGFitFailed);

        //totgain = 1/provided constant from the calibration
        Float16_t energy_TOT_contrib = adcConv.TOT_to_lowGain() * (totgain - adcConv.TOT_offset()) * adcConv.lowGain_to_highGain() * adcConv.adc_to_MIP();

        Float16_t energy_LG_contrib = 0;
        if (((fitresultLG.status == 0) && (lowgain < adcConv.TOT_lowGain_transition())) || ((fitresultLG.status != 0) && (rawhit.lowGainADC(m_expectedMaxTimeSample) < adcConv.TOT_lowGain_transition()))) {
          energy_LG_contrib = lowgain * adcConv.lowGain_to_highGain() *  adcConv.adc_to_MIP();
        } else recHit.setFlag(HGCalTBRecHit::kLowGainSaturated);

        Float16_t energy_HG_contrib = 0;
        if (((fitresultHG.status == 0) && (highgain < adcConv.lowGain_highGain_transition())) || ((fitresultHG.status != 0) && (rawhit.highGainADC(m_expectedMaxTimeSample) < adcConv.lowGain_highGain_transition()))) {
          energy_HG_contrib = highgain *  adcConv.adc_to_MIP();
        } else recHit.setFlag(HGCalTBRecHit::kHighGainSaturated);

        //gain switching
        if (energy_HG_contrib > 0) {
          energy = energy_HG_contrib;
          recHit.setEnergy_HGExcl(energy_LG_contrib);
        } else if (energy_LG_contrib > 0) {
          energy = energy_LG_contrib;
          recHit.setEnergy_HGExcl(energy_LG_contrib);
        } else if (energy_TOT_contrib > 0) {
          energy = energy_TOT_contrib;
          recHit.setEnergy_HGExcl(energy_TOT_contrib);
        }



        if (adcConv.fully_calibrated() == 1) recHit.setFlag(HGCalTBRecHit::kFullyCalibrated);
      }


      if (energy < 0) continue;
      recHit.setEnergy(energy);

      rechits->push_back(recHit);
    }
  }
  event.put(std::move(rechits), m_outputCollectionName);
}

DEFINE_FWK_MODULE(HGCalTBRecHitProducer);