#include <iostream>
#include "TH1F.h"
#include "TH2F.h"
#include "TH2Poly.h"
#include "TTree.h"
#include <fstream>
#include <sstream>
// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/one/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ServiceRegistry/interface/Service.h"

#include "CommonTools/UtilAlgos/interface/TFileService.h"

#include "HGCal/DataFormats/interface/HGCalTBRawHitCollection.h"
#include "HGCal/DataFormats/interface/HGCalTBDetId.h"
#include "HGCal/DataFormats/interface/HGCalTBElectronicsId.h"
#include "HGCal/CondObjects/interface/HGCalElectronicsMap.h"
#include "HGCal/CondObjects/interface/HGCalCondObjectTextIO.h"
#include "HGCal/Geometry/interface/HGCalTBGeometryParameters.h"
#include "HGCal/Reco/interface/CommonMode.h"
#include "HGCal/Reco/interface/PulseFitter.h"

#include <iomanip>
#include <set>

class PulseShapePlotter : public edm::one::EDAnalyzer<edm::one::SharedResources>
{
public:
  explicit PulseShapePlotter(const edm::ParameterSet&);
  ~PulseShapePlotter();
  static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);
private:
  virtual void beginJob() override;
  void analyze(const edm::Event& , const edm::EventSetup&) override;
  virtual void endJob() override;

  std::string m_electronicMap;
  struct {
    HGCalElectronicsMap emap_;
  } essource_;
  double m_commonModeThreshold; //currently not use (need to implement the "average" option in CommonMode.cc)
  int m_expectedMaxTimeSample;
  
  edm::EDGetTokenT<HGCalTBRawHitCollection> m_HGCalTBRawHitCollection;
  TTree *m_tree;
  int m_evtID;
  int m_channelID;
  int m_skirocID;
  int m_layerID;
  float m_hgADC;
  float m_hgTmax;
  float m_hgChi2;
  float m_hgErrorADC;
  float m_hgErrorTmax;
  int m_hgStatus;
  int m_hgNCalls;
  float m_lgADC;
  float m_lgTmax;
  float m_lgChi2;
  float m_lgErrorADC;
  float m_lgErrorTmax;
  int m_lgStatus;
  int m_lgNCalls;
  int m_totSlow;
  int m_toaRise;
  int m_toaFall;
};

PulseShapePlotter::PulseShapePlotter(const edm::ParameterSet& iConfig) :
  m_electronicMap(iConfig.getUntrackedParameter<std::string>("ElectronicMap","HGCal/CondObjects/data/map_CERN_Hexaboard_OneLayers_May2017.txt")),
  m_commonModeThreshold(iConfig.getUntrackedParameter<double>("CommonModeThreshold",100)),
  m_expectedMaxTimeSample(iConfig.getUntrackedParameter<int>("ExpectedMaxTimeSample",3))
{
  m_HGCalTBRawHitCollection = consumes<HGCalTBRawHitCollection>(iConfig.getParameter<edm::InputTag>("InputCollection"));
  m_evtID=0;
  usesResource("TFileService");
  edm::Service<TFileService> fs;
  m_tree=fs->make<TTree>("tree","Pulse shape fitter results");
  m_tree->Branch("eventID",&m_evtID);
  m_tree->Branch("skirocID",&m_skirocID);
  m_tree->Branch("boardID",&m_layerID);
  m_tree->Branch("channelID",&m_channelID);
  m_tree->Branch("HighGainADC",&m_hgADC);
  m_tree->Branch("HighGainTmax",&m_hgTmax);
  m_tree->Branch("HighGainChi2",&m_hgChi2);
  m_tree->Branch("HighGainErrorADC",&m_hgErrorADC);
  m_tree->Branch("HighGainErrorTmax",&m_hgErrorTmax);
  m_tree->Branch("HighGainStatus",&m_hgStatus);
  m_tree->Branch("HighGainNCalls",&m_hgNCalls);

  m_tree->Branch("LowGainADC",&m_lgADC);
  m_tree->Branch("LowGainTmax",&m_lgTmax);
  m_tree->Branch("LowGainChi2",&m_lgChi2);
  m_tree->Branch("LowGainErrorADC",&m_lgErrorADC);
  m_tree->Branch("LowGainErrorTmax",&m_lgErrorTmax);
  m_tree->Branch("LowGainStatus",&m_lgStatus);
  m_tree->Branch("LowGainNCalls",&m_lgNCalls);
  m_tree->Branch("TotSlow",&m_totSlow);
  m_tree->Branch("ToaRise",&m_toaRise);
  m_tree->Branch("ToaFall",&m_toaFall);
  std::cout << iConfig.dump() << std::endl;
}


PulseShapePlotter::~PulseShapePlotter()
{

}

void PulseShapePlotter::beginJob()
{
  HGCalCondObjectTextIO io(0);
  edm::FileInPath fip(m_electronicMap);
  if (!io.load(fip.fullPath(), essource_.emap_)) {
    throw cms::Exception("Unable to load electronics map");
  };
}

void PulseShapePlotter::analyze(const edm::Event& event, const edm::EventSetup& setup)
{
  usesResource("TFileService");
  edm::Service<TFileService> fs;
  std::map<int,TH1F*>  hMapHG,hMapLG;
  std::ostringstream os( std::ostringstream::ate );
  os.str("");
  os << "Event" << event.id().event();
  TFileDirectory dir = fs->mkdir( os.str().c_str() );
  for(size_t ib = 0; ib<HGCAL_TB_GEOMETRY::NUMBER_OF_HEXABOARD; ib++) {
    for(size_t is = 0; is<HGCAL_TB_GEOMETRY::N_SKIROC_PER_HEXA; is++) {
      os.str("");
      os << "Hexaboard" << ib << "_Skiroc" << is ;
      TFileDirectory subdir = dir.mkdir( os.str().c_str() );
      for( size_t ic=0; ic<HGCAL_TB_GEOMETRY::N_CHANNELS_PER_SKIROC; ic++ ){
	int skiId=HGCAL_TB_GEOMETRY::N_SKIROC_PER_HEXA*ib+(HGCAL_TB_GEOMETRY::N_SKIROC_PER_HEXA-is)%HGCAL_TB_GEOMETRY::N_SKIROC_PER_HEXA+1;
	HGCalTBElectronicsId eid(skiId,ic);      
	if (!essource_.emap_.existsEId(eid.rawId())) continue;
	os.str("");
	os<<"HighGain_Channel"<<ic;
	TH1F *hHG=subdir.make<TH1F>(os.str().c_str(),os.str().c_str(),NUMBER_OF_TIME_SAMPLES,0,NUMBER_OF_TIME_SAMPLES*25);
	hMapHG.insert( std::pair<int,TH1F*>(1000*ib+100*is+ic,hHG) );
	os.str("");
	os<<"LowGain_Channel"<<ic;
	TH1F* hLG=subdir.make<TH1F>(os.str().c_str(),os.str().c_str(),NUMBER_OF_TIME_SAMPLES,0,NUMBER_OF_TIME_SAMPLES*25);
	hMapLG.insert( std::pair<int,TH1F*>(1000*ib+100*is+ic,hLG) );
      }
    }
  }
  
  edm::Handle<HGCalTBRawHitCollection> hits;
  event.getByToken(m_HGCalTBRawHitCollection, hits);
  
  CommonMode cm(essource_.emap_); //default is common mode per chip using the median
  cm.Evaluate( hits );
  std::map<int,commonModeNoise> cmMap=cm.CommonModeNoiseMap();
  PulseFitter fitter(0,200);
  for( auto hit : *hits ){
    HGCalTBElectronicsId eid( essource_.emap_.detId2eid(hit.detid().rawId()) );
    if( essource_.emap_.existsEId(eid) ){
      int iski=hit.skiroc();
      float highGain,lowGain;
      float subHG[NUMBER_OF_TIME_SAMPLES],subLG[NUMBER_OF_TIME_SAMPLES];
      for( int it=0; it<NUMBER_OF_TIME_SAMPLES; it++ ){
      	subHG[it]=0;
      	subLG[it]=0;
      }
      switch ( hit.detid().cellType() ){
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
      int iboard=hit.skiroc()/HGCAL_TB_GEOMETRY::N_SKIROC_PER_HEXA;
      int ichan=hit.channel();
      std::vector<double> hg,lg,time;
      for( int it=0; it<NUMBER_OF_TIME_SAMPLES; it++ ){
	highGain=hit.highGainADC(it)-subHG[it];
	lowGain=hit.lowGainADC(it)-subLG[it];
	hg.push_back(highGain);
	lg.push_back(lowGain);
	time.push_back(25*it+12.5);
	hMapHG[1000*iboard+100*(iski%HGCAL_TB_GEOMETRY::N_SKIROC_PER_HEXA)+ichan]->Fill(25*it,highGain);
	hMapLG[1000*iboard+100*(iski%HGCAL_TB_GEOMETRY::N_SKIROC_PER_HEXA)+ichan]->Fill(25*it,lowGain);
      }
      float max_minus=hit.highGainADC(m_expectedMaxTimeSample-1)-subHG[m_expectedMaxTimeSample-1];
      float themax=hit.highGainADC(m_expectedMaxTimeSample)-subHG[m_expectedMaxTimeSample];
      float max_plus=hit.highGainADC(m_expectedMaxTimeSample+1)-subHG[m_expectedMaxTimeSample+1];
      float undershoot=hit.highGainADC(m_expectedMaxTimeSample+3)-subHG[m_expectedMaxTimeSample+3];
      if( themax>500||(max_minus<themax && themax>undershoot && max_plus>undershoot && themax>20) ){
	// std::cout << iboard << " " << iski%HGCAL_TB_GEOMETRY::N_SKIROC_PER_HEXA << " " << ichan << "\t" << max_minus << " " << themax << " " << max_plus << " " << undershoot << std::endl;
	PulseFitterResult fithg;
	fitter.run( time,hg,fithg,8. );
	PulseFitterResult fitlg;
	fitter.run( time,lg,fitlg,2. );
	// std::cout << "\t" << fithg.amplitude << " " << fithg.tmax << " " << fithg.chi2 << std::endl;
	// std::cout << "\t" << fitlg.amplitude << " " << fitlg.tmax << " " << fitlg.chi2 << std::endl;
	// getchar();
	m_channelID=ichan;
	m_skirocID=iski%HGCAL_TB_GEOMETRY::N_SKIROC_PER_HEXA;
	m_layerID=iboard;
	m_hgADC=fithg.amplitude;
	m_hgTmax=fithg.tmax;
	m_hgChi2=fithg.chi2;
	m_hgErrorADC=fithg.erroramplitude;
	m_hgErrorTmax=fithg.errortmax;
	m_hgStatus=fithg.status;
	m_hgNCalls=fithg.ncalls;
	m_lgADC=fitlg.amplitude;
	m_lgTmax=fitlg.tmax;
	m_lgChi2=fitlg.chi2;
	m_lgErrorADC=fitlg.erroramplitude;
	m_lgErrorTmax=fitlg.errortmax;
	m_lgStatus=fitlg.status;
	m_lgNCalls=fitlg.ncalls;
	m_totSlow=hit.totSlow();
	m_toaRise=hit.toaRise();
	m_toaFall=hit.toaFall();
	m_tree->Fill();
      }
    }
  }
  m_evtID++;
}

void PulseShapePlotter::endJob()
{
}

void PulseShapePlotter::fillDescriptions(edm::ConfigurationDescriptions& descriptions)
{
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

DEFINE_FWK_MODULE(PulseShapePlotter);
