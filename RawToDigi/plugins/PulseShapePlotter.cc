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
#include "HGCal/DataFormats/interface/HGCalTBRawHit.h"
#include "HGCal/DataFormats/interface/HGCalTBDetId.h"
#include "HGCal/DataFormats/interface/HGCalTBElectronicsId.h"
#include "HGCal/CondObjects/interface/HGCalElectronicsMap.h"
#include "HGCal/CondObjects/interface/HGCalCondObjectTextIO.h"
#include "HGCal/Geometry/interface/HGCalTBGeometryParameters.h"
#include "HGCal/Reco/interface/CommonMode.h"
#include "HGCal/Reco/interface/PulseFitter.h"

#include <iomanip>
#include <set>
#include <boost/thread/thread.hpp>
//#include <boost/mutex.hpp>

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
  static void fittingPulseShapes(HGCalTBRawHit* hit,float *subLG, float* subHG);

  std::string m_electronicMap;
  std::string m_detectorLayoutFile;
  boost::mutex m_mutex;
  struct {
    HGCalElectronicsMap emap_;
    HGCalTBDetectorLayout layout_;
  } essource_;
  double m_commonModeThreshold; //currently not use (need to implement the "average" option in CommonMode.cc)
  int m_expectedMaxTimeSample;
  int m_maxADCCut;
  bool m_savePulseShapes;
  
  edm::EDGetTokenT<HGCalTBRawHitCollection> m_HGCalTBRawHitCollection;
  TTree *m_tree;
  int m_evtID;
  std::vector<int> m_channelID;
  std::vector<int> m_skirocID;
  std::vector<int> m_layerID;
  std::vector<int> m_moduleID;
  std::vector<float> m_hgADC;
  std::vector<float> m_hgTmax;
  std::vector<float> m_hgChi2;
  std::vector<float> m_hgErrorADC;
  std::vector<float> m_hgErrorTmax;
  std::vector<int> m_hgStatus;
  std::vector<int> m_hgNCalls;
  std::vector<float> m_lgADC;
  std::vector<float> m_lgTmax;
  std::vector<float> m_lgChi2;
  std::vector<float> m_lgErrorADC;
  std::vector<float> m_lgErrorTmax;
  std::vector<int> m_lgStatus;
  std::vector<int> m_lgNCalls;
  std::vector<int> m_totSlow;
  std::vector<int> m_toaRise;
  std::vector<int> m_toaFall;
};

PulseShapePlotter::PulseShapePlotter(const edm::ParameterSet& iConfig) :
  m_electronicMap(iConfig.getUntrackedParameter<std::string>("ElectronicMap","elecmap.txt")),
  m_detectorLayoutFile(iConfig.getUntrackedParameter<std::string>("DetectorLayout","geom.txt")),
  m_commonModeThreshold(iConfig.getUntrackedParameter<double>("CommonModeThreshold",100)),
  m_expectedMaxTimeSample(iConfig.getUntrackedParameter<int>("ExpectedMaxTimeSample",3)),
  m_maxADCCut(iConfig.getUntrackedParameter<double>("MaxADCCut",15)),
  m_savePulseShapes(iConfig.getUntrackedParameter<bool>("SavePulseShapes",false))
{
  m_HGCalTBRawHitCollection = consumes<HGCalTBRawHitCollection>(iConfig.getParameter<edm::InputTag>("InputCollection"));
  m_evtID=0;
  usesResource("TFileService");
  edm::Service<TFileService> fs;
  m_tree=fs->make<TTree>("tree","Pulse shape fitter results");
  m_tree->Branch("eventID",&m_evtID);
  m_tree->Branch("skirocID",&m_skirocID);
  m_tree->Branch("boardID",&m_layerID);
  m_tree->Branch("moduleID",&m_moduleID);
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
  fip=edm::FileInPath(m_detectorLayoutFile);
  if (!io.load(fip.fullPath(), essource_.layout_)) {
    throw cms::Exception("Unable to load detector layout file");
  };
  for( auto layer : essource_.layout_.layers() )
    layer.print();
}

void PulseShapePlotter::fittingPulseShapes(HGCalTBRawHit* hit,float *subLG, float* subHG)
{
}

void PulseShapePlotter::analyze(const edm::Event& event, const edm::EventSetup& setup)
{
  m_channelID.clear();
  m_skirocID.clear();
  m_layerID.clear();
  m_moduleID.clear();
  m_hgADC.clear();
  m_hgTmax.clear();
  m_hgChi2.clear();
  m_hgErrorADC.clear();
  m_hgErrorTmax.clear();
  m_hgStatus.clear();
  m_hgNCalls.clear();
  m_lgADC.clear();
  m_lgTmax.clear();
  m_lgChi2.clear();
  m_lgErrorADC.clear();
  m_lgErrorTmax.clear();
  m_lgStatus.clear();
  m_lgNCalls.clear();
  m_totSlow.clear();
  m_toaRise.clear();
  m_toaFall.clear();
  
  std::map<int,TH1F*>  hMapHG,hMapLG;
  if( m_savePulseShapes ){
    usesResource("TFileService");
    edm::Service<TFileService> fs;
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
	  TH1F *hHG=subdir.make<TH1F>(os.str().c_str(),os.str().c_str(),NUMBER_OF_TIME_SAMPLES+2,-25,(1+NUMBER_OF_TIME_SAMPLES)*25);
	  hMapHG.insert( std::pair<int,TH1F*>(1000*ib+100*is+ic,hHG) );
	  os.str("");
	  os<<"LowGain_Channel"<<ic;
	  TH1F* hLG=subdir.make<TH1F>(os.str().c_str(),os.str().c_str(),NUMBER_OF_TIME_SAMPLES+2,-25,(1+NUMBER_OF_TIME_SAMPLES)*25);
	  hMapLG.insert( std::pair<int,TH1F*>(1000*ib+100*is+ic,hLG) );
	}
      }
    }
  }
  
  edm::Handle<HGCalTBRawHitCollection> hits;
  event.getByToken(m_HGCalTBRawHitCollection, hits);
  
  CommonMode cm(essource_.emap_); //default is common mode per chip using the median
  cm.Evaluate( hits );
  std::map<int,commonModeNoise> cmMap=cm.CommonModeNoiseMap();
  PulseFitter fitter(0);
  for( auto hit : *hits ){
    HGCalTBElectronicsId eid( essource_.emap_.detId2eid(hit.detid().rawId()) );
    if( essource_.emap_.existsEId(eid) ){
      int iski=hit.skiroc();
      float highGain,lowGain;
      float subHG[NUMBER_OF_TIME_SAMPLES];
      float subLG[NUMBER_OF_TIME_SAMPLES];
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
      // double *hg[NUMBER_OF_TIME_SAMPLES];
      // double *lg[NUMBER_OF_TIME_SAMPLES];
      // double *time[NUMBER_OF_TIME_SAMPLES];
      std::vector<double> hg(NUMBER_OF_TIME_SAMPLES),lg(NUMBER_OF_TIME_SAMPLES),time(NUMBER_OF_TIME_SAMPLES);
      for( int it=0; it<NUMBER_OF_TIME_SAMPLES; it++ ){
	highGain=hit.highGainADC(it)-subHG[it];
	lowGain=hit.lowGainADC(it)-subLG[it];
	hg[it]=highGain;
	lg[it]=lowGain;
	time[it]=25*it+12.5;
	if( m_savePulseShapes ){
	  hMapHG[1000*iboard+100*(iski%HGCAL_TB_GEOMETRY::N_SKIROC_PER_HEXA)+ichan]->Fill(time[it],highGain);
	  hMapLG[1000*iboard+100*(iski%HGCAL_TB_GEOMETRY::N_SKIROC_PER_HEXA)+ichan]->Fill(time[it],lowGain);
	}
      }
      float max_minus=hg[m_expectedMaxTimeSample-2];
      float themax=hg[m_expectedMaxTimeSample];
      float max_plus=hg[m_expectedMaxTimeSample+1];
      float undershoot=hg[m_expectedMaxTimeSample+3];
      if( themax>500||(max_minus<themax && themax>undershoot && max_plus>undershoot && themax>20) ){
	// std::cout << iboard << " " << iski%HGCAL_TB_GEOMETRY::N_SKIROC_PER_HEXA << " " << ichan << "\t" << max_minus << " " << themax << " " << max_plus << " " << undershoot << std::endl;
	PulseFitterResult fithg;
	fitter.run( time,hg,fithg,5. );
	PulseFitterResult fitlg;
	fitter.run( time,lg,fitlg,2. );
	// std::cout << "\t" << fithg.amplitude << " " << fithg.tmax << " " << fithg.chi2 << std::endl;
	// std::cout << "\t" << fitlg.amplitude << " " << fitlg.tmax << " " << fitlg.chi2 << std::endl;
	// getchar();
	m_channelID.push_back(ichan);
	m_skirocID.push_back(iski%HGCAL_TB_GEOMETRY::N_SKIROC_PER_HEXA);
	m_layerID.push_back(iboard);
	m_moduleID.push_back(essource_.layout_.at(iboard).at(0).moduleID());
	m_hgADC.push_back(fithg.amplitude);
	m_hgTmax.push_back(fithg.tmax);
	m_hgChi2.push_back(fithg.chi2);
	m_hgErrorADC.push_back(fithg.erroramplitude);
	m_hgErrorTmax.push_back(fithg.errortmax);
	m_hgStatus.push_back(fithg.status);
	m_hgNCalls.push_back(fithg.ncalls);
	m_lgADC.push_back(fitlg.amplitude);
	m_lgTmax.push_back(fitlg.tmax);
	m_lgChi2.push_back(fitlg.chi2);
	m_lgErrorADC.push_back(fitlg.erroramplitude);
	m_lgErrorTmax.push_back(fitlg.errortmax);
	m_lgStatus.push_back(fitlg.status);
	m_lgNCalls.push_back(fitlg.ncalls);
	m_totSlow.push_back(hit.totSlow());
	m_toaRise.push_back(hit.toaRise());
	m_toaFall.push_back(hit.toaFall());
      }
    }
  }
  m_tree->Fill();
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
