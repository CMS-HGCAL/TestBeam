#include "HGCal/Reco/plugins/HGCalTBRecHitProducer.h"
#include <iostream>

const static size_t N_SKIROC_PER_HEXA = 4;

HGCalTBRecHitProducer::HGCalTBRecHitProducer(const edm::ParameterSet& cfg) : 
  m_outputCollectionName(cfg.getParameter<std::string>("OutputCollectionName")),
  m_electronicMap(cfg.getUntrackedParameter<std::string>("ElectronicsMap","HGCal/CondObjects/data/map_CERN_Hexaboard_OneLayers_May2017.txt")),
  m_commonModeThreshold(cfg.getUntrackedParameter<double>("CommonModeThreshold",100)),
  m_highGainADCSaturation(cfg.getUntrackedParameter<double>("HighGainADCSaturation",1800)),
  m_lowGainADCSaturation(cfg.getUntrackedParameter<double>("LowGainADCSaturation",1800)),
  m_keepOnlyTimeSample3(cfg.getUntrackedParameter<bool>("KeepOnlyTimeSample3",true)),   //checked first
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

  commonModeNoise cm[NUMBER_OF_TIME_SAMPLES][N_SKIROC_PER_HEXA];

  for( auto rawhit : *rawhits ){
    if( !essource_.emap_.existsDetId(rawhit.detid()) ) continue;
    HGCalTBElectronicsId eid( essource_.emap_.detId2eid(rawhit.detid()) );
    int iboard=(eid.iskiroc()-1)/N_SKIROC_PER_HEXA;

    int iski=(N_SKIROC_PER_HEXA-(eid.iskiroc()-1))%N_SKIROC_PER_HEXA;//from 0 to 3
    iski+=iboard*N_SKIROC_PER_HEXA;
    for( int it=0; it<NUMBER_OF_TIME_SAMPLES; it++ ){
      if( rawhit.highGainADC(it)>m_commonModeThreshold ) continue;
      float highGain = rawhit.highGainADC(it);
      float lowGain = rawhit.lowGainADC(it);
      switch ( rawhit.detid().cellType() ){
      case 0 : cm[it][iski].fullHG += highGain; cm[it][iski].fullLG += lowGain; cm[it][iski].fullCounter++; break;
      case 2 : cm[it][iski].halfHG += highGain; cm[it][iski].halfLG += lowGain; cm[it][iski].halfCounter++; break;
      case 3 : cm[it][iski].mouseBiteHG += highGain; cm[it][iski].mouseBiteLG += lowGain; cm[it][iski].mouseBiteCounter++; break;
      case 4 : cm[it][iski].outerHG += highGain; cm[it][iski].outerLG += lowGain; cm[it][iski].outerCounter++; break;
      }
    }

  }

  for( auto rawhit : *rawhits ){

    if( !essource_.emap_.existsDetId(rawhit.detid()) ) continue;
    HGCalTBElectronicsId eid( essource_.emap_.detId2eid(rawhit.detid()) );
    int iboard=(eid.iskiroc()-1)/N_SKIROC_PER_HEXA;

    int iski=(N_SKIROC_PER_HEXA-(eid.iskiroc()-1))%N_SKIROC_PER_HEXA;//from 0 to 3
    iski+=iboard*N_SKIROC_PER_HEXA;

    if(m_keepOnlyTimeSample3){
      float highGain,lowGain;
      float subHG[NUMBER_OF_TIME_SAMPLES],subLG[NUMBER_OF_TIME_SAMPLES];
      for( int it=0; it<NUMBER_OF_TIME_SAMPLES; it++ ){
      	subHG[it]=0;
      	subLG[it]=0;
      }
      switch ( rawhit.detid().cellType() ){
      case 0 : 
      	for( int it=0; it<NUMBER_OF_TIME_SAMPLES; it++ ){
      	  subHG[it]=cm[it][iski].fullCounter>0 ? cm[it][iski].fullHG/cm[it][iski].fullCounter : 0; 
      	  subLG[it]=cm[it][iski].fullCounter>0 ? cm[it][iski].fullLG/cm[it][iski].fullCounter : 0; 
      	}
      	break;
      case 2 : 
      	for( int it=0; it<NUMBER_OF_TIME_SAMPLES; it++ ){
      	  subHG[it]=cm[it][iski].halfCounter>0 ? cm[it][iski].halfHG/cm[it][iski].halfCounter : 0; 
      	  subLG[it]=cm[it][iski].halfCounter>0 ? cm[it][iski].halfLG/cm[it][iski].halfCounter : 0; 
      	}
      	break;
      case 3 : 
      	for( int it=0; it<NUMBER_OF_TIME_SAMPLES; it++ ){
      	  subHG[it]=cm[it][iski].mouseBiteCounter>0 ? cm[it][iski].mouseBiteHG/cm[it][iski].mouseBiteCounter : 0; 
      	  subLG[it]=cm[it][iski].mouseBiteCounter>0 ? cm[it][iski].mouseBiteLG/cm[it][iski].mouseBiteCounter : 0; 
      	}
      	break;
      case 4 : for( int it=0; it<NUMBER_OF_TIME_SAMPLES; it++ ){
       	  subHG[it]=cm[it][iski].outerCounter>0 ? cm[it][iski].outerHG/cm[it][iski].outerCounter : 0; 
       	  subLG[it]=cm[it][iski].outerCounter>0 ? cm[it][iski].outerLG/cm[it][iski].outerCounter : 0; 
       	}
       	break;
      }
      //this is a just try to isolate hits with signal
      float en2=rawhit.highGainADC(2)-subHG[2];
      float en3=rawhit.highGainADC(3)-subHG[3];
      float en4=rawhit.highGainADC(4)-subHG[4];
      float en5=rawhit.highGainADC(5)-subHG[5];
      float en6=rawhit.highGainADC(6)-subHG[6];
      if( en2<en3 && en3>en5 && en4>en5 && (en5>en6||(en5<0&&en6<0)) ){
    	highGain=rawhit.highGainADC(3)-subHG[3];
    	lowGain=rawhit.lowGainADC(3)-subLG[3];
      }
      else{
       	highGain=-500;
       	lowGain=-500;
      }
      float energy = (highGain<m_highGainADCSaturation) ? highGain : lowGain*m_LG2HG_value.at(iboard);
      energy = (energy > 0.) ? energy : 0.;
      float time = rawhit.toaRise();

      HGCalTBRecHit recHit(rawhit.detid(), energy, lowGain, highGain, time);
      
      CellCenterXY = TheCell.GetCellCentreCoordinatesForPlots((recHit.id()).layer(), (recHit.id()).sensorIU(), (recHit.id()).sensorIV(), (recHit.id()).iu(), (recHit.id()).iv(), 128); 
      recHit.setCellCenterCoordinate(CellCenterXY.first, CellCenterXY.second);

      rechits->push_back( recHit );

    } else if (m_performParabolicFit) {
      int N_entries = 3;    //number of data points for the parabolic fit
      int index_of_first = 2;

      float highGain[N_entries]; float lowGain[N_entries];
      float subHG[N_entries]; float subLG[N_entries];
      float _x[N_entries];
      
      for (int i=0; i<N_entries; i++) {
        int time_stamp = index_of_first+i;
        _x[i] = time_stamp;
      
        switch ( rawhit.detid().cellType() ){
          case 0 : subHG[i]=cm[time_stamp][iski].fullHG/cm[time_stamp][iski].fullCounter; subLG[i]=cm[time_stamp][iski].fullLG/cm[time_stamp][iski].fullCounter; break;
          //case 2 : subHG[i]=cm[time_stamp][iski].fullHG/cm[time_stamp][iski].fullCounter; subLG[i]=cm[time_stamp][iski].fullLG/cm[time_stamp][iski].fullCounter; break;
          case 2 : subHG[i]=cm[time_stamp][iski].halfHG/cm[time_stamp][iski].halfCounter; subLG[i]=cm[time_stamp][iski].halfLG/cm[time_stamp][iski].halfCounter; break;
          case 3 : subHG[i]=cm[time_stamp][iski].mouseBiteHG/cm[time_stamp][iski].mouseBiteCounter; subLG[i]=cm[time_stamp][iski].mouseBiteLG/cm[time_stamp][iski].mouseBiteCounter; break;
          case 4 : subHG[i]=cm[time_stamp][iski].outerHG/cm[time_stamp][iski].outerCounter; subLG[i]=cm[time_stamp][iski].outerLG/cm[time_stamp][iski].outerCounter; break;
        }
      
        highGain[i] = rawhit.highGainADC(time_stamp)-subHG[i];
        lowGain[i] = rawhit.lowGainADC(time_stamp)-subLG[i];

      }

      //energy reconstruction from parabolic function a*x^2 + b*x + c
      double _aHG = ( (highGain[2]-highGain[1])/(_x[2]-_x[1]) - (highGain[1]-highGain[0])/(_x[1]-_x[0]) )/( (pow(_x[2],2)-pow(_x[1],2))/(_x[2]-_x[1]) - (pow(_x[1],2)-pow(_x[0],2))/(_x[1]-_x[0]) );
      double _aLG = ( (lowGain[2]-lowGain[1])/(_x[2]-_x[1]) - (lowGain[1]-lowGain[0])/(_x[1]-_x[0]) )/( (pow(_x[2],2)-pow(_x[1],2))/(_x[2]-_x[1]) - (pow(_x[1],2)-pow(_x[0],2))/(_x[1]-_x[0]) );
      double _bHG = (highGain[2]-highGain[1]+_aHG*(pow(_x[1],2)-pow(_x[2],2)))/(_x[2]-_x[1]);
      double _bLG = (lowGain[2]-lowGain[1]+_aLG*(pow(_x[1],2)-pow(_x[2],2)))/(_x[2]-_x[1]);
      double _cHG = highGain[0]-_aHG*pow(_x[0],2)-_bHG*_x[0];
      double _cLG = lowGain[0]-_aLG*pow(_x[0],2)-_bLG*_x[0];

      
      double _x_max_HG = (_aHG < 0) ? -_bHG/(2*_aHG) : 0;   //require maximum <--> a<0
      double _x_max_LG = (_aLG < 0) ? -_bLG/(2*_aLG) : 0;
      
      double energyHG = (_x_max_HG<=4 && _x_max_HG>=2) ? _aHG*pow(_x_max_HG,2)+_bHG*_x_max_HG+_cHG : 0;
      double energyLG = (_x_max_LG<=4 && _x_max_LG>=2) ? _aLG*pow(_x_max_LG,2)+_bLG*_x_max_LG+_cLG : 0;
      
      float energyTOT = rawhit.totFast();
      
      float time = rawhit.toaRise();

      HGCalTBRecHit recHit(rawhit.detid(), energyTOT, energyLG, energyHG, time);
      
      CellCenterXY = TheCell.GetCellCentreCoordinatesForPlots((recHit.id()).layer(), (recHit.id()).sensorIU(), (recHit.id()).sensorIV(), (recHit.id()).iu(), (recHit.id()).iv(), 128); 
      recHit.setCellCenterCoordinate(CellCenterXY.first, CellCenterXY.second);

      rechits->push_back( recHit );
      

    } else{   //todo: perform parabolic fit to determine the maximum
      std::cout << "Should run with m_keepOnlyTimeSample3 sets to true, other method not yet implemented -> exit" << std::endl;
      exit(1);
    }
  }
  event.put(rechits, m_outputCollectionName);
}

DEFINE_FWK_MODULE(HGCalTBRecHitProducer);
