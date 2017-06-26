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
  m_keepOnlyTimeSample3(cfg.getUntrackedParameter<bool>("KeepOnlyTimeSample3",true)),
  m_performParabolicFit(cfg.getUntrackedParameter<bool>("performParabolicFit",true))   //checked secondly
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

    float highGain, lowGain, totGain;
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
    

    
    if(!m_keepOnlyTimeSample3 && !m_performParabolicFit && 
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

    }
    else if(m_keepOnlyTimeSample3 && 
      (en2<en3 && en3>en6 && en4>en6)){
      
      highGain=rawhit.highGainADC(3)-subHG[3];
      lowGain=rawhit.lowGainADC(3)-subLG[3];
    } 

    else if (m_performParabolicFit) {
      int N_entries = 3;    //number of data points for the parabolic fit
      int index_of_first = 2;

      float HG_tmp[N_entries]; float LG_tmp[N_entries];
      float _x[N_entries];
      
      for (int i=0; i<N_entries; i++) {
        int time_stamp = index_of_first+i;
        _x[i] = sampleT[time_stamp];
            
        HG_tmp[i] = sampleHG[time_stamp];
        LG_tmp[i] = sampleLG[time_stamp];

      }

      //energy reconstruction from parabolic function a*x^2 + b*x + c
      double _aHG = ( (HG_tmp[2]-HG_tmp[1])/(_x[2]-_x[1]) - (HG_tmp[1]-HG_tmp[0])/(_x[1]-_x[0]) )/( (pow(_x[2],2)-pow(_x[1],2))/(_x[2]-_x[1]) - (pow(_x[1],2)-pow(_x[0],2))/(_x[1]-_x[0]) );
      double _aLG = ( (LG_tmp[2]-LG_tmp[1])/(_x[2]-_x[1]) - (LG_tmp[1]-LG_tmp[0])/(_x[1]-_x[0]) )/( (pow(_x[2],2)-pow(_x[1],2))/(_x[2]-_x[1]) - (pow(_x[1],2)-pow(_x[0],2))/(_x[1]-_x[0]) );
      double _bHG = (HG_tmp[2]-HG_tmp[1]+_aHG*(pow(_x[1],2)-pow(_x[2],2)))/(_x[2]-_x[1]);
      double _bLG = (LG_tmp[2]-LG_tmp[1]+_aLG*(pow(_x[1],2)-pow(_x[2],2)))/(_x[2]-_x[1]);
      double _cHG = HG_tmp[0]-_aHG*pow(_x[0],2)-_bHG*_x[0];
      double _cLG = LG_tmp[0]-_aLG*pow(_x[0],2)-_bLG*_x[0];

      
      timeHG = (_aHG < 0) ? -_bHG/(2*_aHG) : 0;   //require maximum <--> a<0, unit is ns
      timeLG = (_aLG < 0) ? -_bLG/(2*_aLG) : 0;
      
      highGain = (timeHG<=125. && timeHG>=25.) ? _aHG*pow(timeHG,2)+_bHG*timeHG+_cHG : -1.;
      lowGain = (timeLG<=125. && timeLG>=25.) ? _aLG*pow(timeLG,2)+_bLG*timeLG+_cLG : -1.;
      
    }

    else{
      highGain=-1.;
      lowGain=-1.;
    }

    float energy = -1;
    float time = -1.;
    
    HGCalTBRecHit recHit(rawhit.detid(), energy, lowGain, highGain, totGain, time);
    if(highGain < m_highGainADCSaturation && highGain!=-1.){
      energy = highGain;
      time = timeHG;
      recHit.setFlag(HGCalTBRecHit::kHighGainSaturated); //std::cout<<"High gain saturated   ";
    }     
    else if(lowGain < m_lowGainADCSaturation && lowGain!=-1.){
      energy = lowGain * m_LG2HG_value.at(iboard);
      time = timeLG;
      recHit.setFlag(HGCalTBRecHit::kLowGainSaturated); //std::cout<<"Low gain saturated   ";
    }
    else{
      energy = totGain * m_TOT2LG_value.at(iboard) * m_LG2HG_value.at(iboard);
      recHit.setFlag(HGCalTBRecHit::kTotGainSaturated); //std::cout<<"Tot gain saturated   ";
    }


    if(m_keepOnlyTimeSample3) recHit.setFlag(HGCalTBRecHit::kThirdSample);
    else if(m_performParabolicFit) recHit.setFlag(HGCalTBRecHit::kParabolic);

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