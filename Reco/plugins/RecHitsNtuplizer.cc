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
#include "HGCal/CondObjects/interface/HGCalElectronicsMap.h"
#include "HGCal/CondObjects/interface/HGCalCondObjectTextIO.h"
#include "HGCal/DataFormats/interface/HGCalTBElectronicsId.h"
#include "HGCal/Geometry/interface/HGCalTBCellVertices.h"
#include "HGCal/Reco/interface/PulseFitter.h"

#include <iomanip>
#include <set>
#include <vector>

#ifdef __MAKECINT__
#pragma link C++ class vector<float>+;
#endif


const static size_t N_HEXABOARDS = 1;
const static size_t N_SKIROC_PER_HEXA = 4;
const static size_t N_CHANNELS_PER_SKIROC = 64;
static const double deltaCellBoundary = 0.00001;
const static int SENSORSIZE = 128;

class RecHitsNtuplizer : public edm::one::EDAnalyzer<edm::one::SharedResources>
{
public:
  explicit RecHitsNtuplizer(const edm::ParameterSet&);
  ~RecHitsNtuplizer();
  static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);

  void initTree();

private:
  virtual void beginJob() override;
  void analyze(const edm::Event& , const edm::EventSetup&) override;
  virtual void endJob() override;

  std::string m_electronicMap;
  struct {
    HGCalElectronicsMap emap_;
  } essource_;
  double m_commonModeThreshold; //number of sigmas from ped mean

  std::pair<double, double> CellCentreXY;
  HGCalTBCellVertices TheCell;

  std::vector<double> m_LG2HG_value;
  std::vector<double> m_TOT2LG_value;

  int m_highGainADCSaturation;
  int m_lowGainADCSaturation;
  int m_evtID;
  edm::EDGetTokenT<HGCalTBRawHitCollection> m_HGCalTBRawHitCollection;

  struct commonModeNoise{
    commonModeNoise():fullHG(0),halfHG(0),mouseBiteHG(0),outerHG(0),fullLG(0),halfLG(0),mouseBiteLG(0),outerLG(0),fullCounter(0),halfCounter(0),mouseBiteCounter(0),outerCounter(0){;}
    float fullHG,halfHG,mouseBiteHG,outerHG;
    float fullLG,halfLG,mouseBiteLG,outerLG;
    int fullCounter,halfCounter,mouseBiteCounter,outerCounter;
  };

  TTree* tree;
  edm::Service<TFileService> fs;
  std::vector<double> tsHG;
  std::vector<double> tsLG;
  Float_t posX;
  Float_t posY;
  Float_t amp;
  Float_t ampHG;
  Float_t ampLG;
  Float_t ampTOT;
  Float_t startTHG;
  Float_t startTLG;
  Float_t toaR;
  Float_t toaF;
  Int_t hgSat;
  Int_t lgSat;
  Int_t evtGOOD;
  Int_t nHexB;
  Int_t nSki;
  Int_t nCha;
  Int_t nEvt;
  Int_t nRun;
};

RecHitsNtuplizer::RecHitsNtuplizer(const edm::ParameterSet& iConfig) :
  m_electronicMap(iConfig.getUntrackedParameter<std::string>("ElectronicMap","HGCal/CondObjects/data/map_CERN_Hexaboard_OneLayers_May2017.txt")),
  m_commonModeThreshold(iConfig.getUntrackedParameter<double>("CommonModeThreshold",100)),
  m_highGainADCSaturation(iConfig.getUntrackedParameter<double>("HighGainADCSaturation",1800)),
  m_lowGainADCSaturation(iConfig.getUntrackedParameter<double>("LowGainADCSaturation",1800))
{
  m_HGCalTBRawHitCollection = consumes<HGCalTBRawHitCollection>(iConfig.getParameter<edm::InputTag>("InputCollection"));

  std::vector<double> v0(1,10.);
  m_LG2HG_value = iConfig.getUntrackedParameter<std::vector<double> >("LG2HG",v0);
  std::vector<double> v1(1,10.);
  m_TOT2LG_value = iConfig.getUntrackedParameter<std::vector<double> >("TOT2LG",v1);

  m_evtID=0;
  std::cout << iConfig.dump() << std::endl;

  usesResource("TFileService");
  tree = new TTree("hgcRH","");
 
  tree->Branch("tsHG",&tsHG);
  tree->Branch("tsLG",&tsLG);
  tree->Branch("posX",&posX,"posX/F");
  tree->Branch("posY",&posY,"posY/F");
  tree->Branch("amp",&amp,"amp/F");
  tree->Branch("ampHG",&ampHG,"ampHG/F");
  tree->Branch("ampLG",&ampLG,"ampLG/F");
  tree->Branch("ampTOT",&ampLG,"ampTOT/F");
  tree->Branch("startTHG",&startTHG,"startTHG/F");
  tree->Branch("startTLG",&startTLG,"startTLG/F");
  tree->Branch("toaR",&toaR,"toaR/F");
  tree->Branch("toaF",&toaF,"toaF/F");
  tree->Branch("hgSat",&hgSat,"hgSat/I");
  tree->Branch("lgSat",&lgSat,"lgSat/I");
  tree->Branch("evtGOOD",&evtGOOD,"evtGOOD/I");
  tree->Branch("nHexB",&nHexB,"nHexB/I");
  tree->Branch("nSki",&nSki,"nSki/I");
  tree->Branch("nCha",&nCha,"nCha/I");
  tree->Branch("nEvt",&nEvt,"nEvt/I");
  tree->Branch("nRun",&nRun,"nRun/I");


}


RecHitsNtuplizer::~RecHitsNtuplizer()
{

}

void RecHitsNtuplizer::beginJob()
{
  HGCalCondObjectTextIO io(0);
  edm::FileInPath fip(m_electronicMap);
  if (!io.load(fip.fullPath(), essource_.emap_)) {
    throw cms::Exception("Unable to load electronics map");
  };
}

void RecHitsNtuplizer::initTree()
{
  tsHG.clear();
  tsLG.clear();
  posX = -1.;
  posY = -1.;
  amp = -1.;
  ampHG = -1.;
  ampLG = -1.;
  ampTOT = -1.;
  startTHG = -1.;
  startTLG = -1.;
  toaR = -1.;
  toaF = -1.;
  hgSat = -1.;
  lgSat = -1.;
  evtGOOD = -1.;
  nHexB = -1.;
  nSki = -1.;
  nCha = -1.;
  nEvt = -1.;
  nRun = -1.;

}


void RecHitsNtuplizer::analyze(const edm::Event& event, const edm::EventSetup& setup)
{
  initTree();


  nEvt = event.id().event();
  nRun = event.id().run();
  
  edm::Handle<HGCalTBRawHitCollection> hits;
  event.getByToken(m_HGCalTBRawHitCollection, hits);
  
  commonModeNoise cm[NUMBER_OF_TIME_SAMPLES-0][4];
  for( auto hit : *hits ){
    int iski=hit.skiroc();
    if( !essource_.emap_.existsDetId(hit.detid()) ) continue;
    for( size_t it=0; it<NUMBER_OF_TIME_SAMPLES; it++ ){
      if( hit.highGainADC(it)>m_commonModeThreshold ) continue;
      float highGain = hit.highGainADC(it);
      float lowGain = hit.lowGainADC(it);
      switch ( hit.detid().cellType() ){
      case 0 : cm[it][iski].fullHG += highGain; cm[it][iski].fullLG += lowGain; cm[it][iski].fullCounter++; break;
      case 2 : cm[it][iski].halfHG += highGain; cm[it][iski].halfLG += lowGain; cm[it][iski].halfCounter++; break;
      case 3 : cm[it][iski].mouseBiteHG += highGain; cm[it][iski].mouseBiteLG += lowGain; cm[it][iski].mouseBiteCounter++; break;
      case 4 : cm[it][iski].outerHG += highGain; cm[it][iski].outerLG += lowGain; cm[it][iski].outerCounter++; break;
      }
    }
  }
  
  for( auto hit : *hits ){
    toaR = hit.toaRise();
    toaF = hit.toaFall();

    if(toaR == 4.) continue;

    int iboard = hit.skiroc()/N_SKIROC_PER_HEXA;
    int iski = hit.skiroc();
    int ichan = hit.channel();

    nHexB = 1000*iboard;
    nSki = 100*iski;
    nCha = ichan;


    if( essource_.emap_.existsDetId(hit.detid()) ){
      float highGain, lowGain;
      int hgStatus = -1;
      int lgStatus = -1;
      std::vector<double> sampleT;
      sampleT.clear();
      tsHG.clear();
      tsLG.clear();
      ampTOT = hit.totSlow();

      float subHG[NUMBER_OF_TIME_SAMPLES],subLG[NUMBER_OF_TIME_SAMPLES];
      for( int it=0; it<NUMBER_OF_TIME_SAMPLES; it++ ){
      	subHG[it]=0;
      	subLG[it]=0;
      }
      switch ( hit.detid().cellType() ){
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
      for( int it=0; it<NUMBER_OF_TIME_SAMPLES; it++ ){
	highGain = hit.highGainADC(it)-subHG[it];
	lowGain = hit.lowGainADC(it)-subLG[it];
	
	tsHG.push_back(highGain);
	tsLG.push_back(lowGain);
	sampleT.push_back(it*25);
      }
      //this is a just try to isolate hits with signal 
      float en2=hit.highGainADC(2)-subHG[2];
      float en3=hit.highGainADC(3)-subHG[3];
      float en4=hit.highGainADC(4)-subHG[4];
      float en6=hit.highGainADC(6)-subHG[6];


      if(en2<en3 && en3>en6 && en4>en6 && en3>20){
	PulseFitter fitter(0,150);
	PulseFitterResult fithg;
	fitter.run(sampleT, tsHG, fithg);
	PulseFitterResult fitlg;
	fitter.run(sampleT, tsLG, fitlg);
	
	ampHG = fithg.amplitude;
	ampLG = fitlg.amplitude;
	startTHG = fithg.tmax - fithg.trise;
	startTLG = fitlg.tmax - fitlg.trise;	

	hgStatus = fithg.status;
	lgStatus = fitlg.status;

	if(ampHG < m_highGainADCSaturation && hgStatus == 0){
	  amp = ampHG;
	  hgSat = 0;
	  evtGOOD = 1;
	}
	else if(ampLG < m_lowGainADCSaturation && hgStatus == 0 && lgStatus == 0){
	  amp = ampLG * m_LG2HG_value.at(iboard);
	  hgSat = 1;
	  lgSat = 0;
	  evtGOOD = 1;
	}
	else if(hgStatus == 0 && lgStatus == 0 && ampTOT > 50){
	  amp = ampTOT * m_TOT2LG_value.at(iboard) * m_LG2HG_value.at(iboard);
	  hgSat = 1;
	  lgSat = 1;
	  evtGOOD = 1;
	}
	else{
	  amp = -1;
	  evtGOOD = 0;
	}
      }

      HGCalTBDetId detid = hit.detid();
      CellCentreXY = TheCell.GetCellCentreCoordinatesForPlots(detid.layer(), detid.sensorIU(), detid.sensorIV(), detid.iu(), detid.iv(), SENSORSIZE );
      posX = (CellCentreXY.first < 0 ) ? (CellCentreXY.first + deltaCellBoundary) : (CellCentreXY.first - deltaCellBoundary);
      posY = (CellCentreXY.second < 0 ) ? (CellCentreXY.second + deltaCellBoundary) : (CellCentreXY.second - deltaCellBoundary);      
    }
    if(startTHG > 0. && startTLG > 0.) tree->Fill();
  }
  m_evtID++;
}

void RecHitsNtuplizer::endJob()
{
}

void RecHitsNtuplizer::fillDescriptions(edm::ConfigurationDescriptions& descriptions)
{
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

DEFINE_FWK_MODULE(RecHitsNtuplizer);
