#include <iostream>
#include "TH1F.h"
#include "TH2F.h"
#include "TH2Poly.h"
#include <fstream>
#include <sstream>
#include <algorithm>
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
#include "HGCal/Geometry/interface/HGCalTBTopology.h"
#include "HGCal/Geometry/interface/HGCalTBGeometryParameters.h"
#include "HGCal/Reco/interface/CommonMode.h"
#include <iomanip>
#include <set>

#define MAXVERTICES 6
static const double delta = 0.00001;//Add/subtract delta = 0.00001 to x,y of a cell centre so the TH2Poly::Fill doesnt have a problem at the edges where the centre of a half-hex cell passes through the sennsor boundary line.

class RawHitPlotter : public edm::one::EDAnalyzer<edm::one::SharedResources>
{
public:
  explicit RawHitPlotter(const edm::ParameterSet&);
  ~RawHitPlotter();
  static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);
private:
  virtual void beginJob() override;
  void analyze(const edm::Event& , const edm::EventSetup&) override;
  virtual void endJob() override;
  void InitTH2Poly(TH2Poly& poly, int layerID, int sensorIU, int sensorIV);

  std::string m_electronicMap;

  struct {
    HGCalElectronicsMap emap_;
  } essource_;
  int m_sensorsize;
  bool m_eventPlotter;
  bool m_subtractCommonMode;
  double m_commonModeThreshold; //currently not use (need to implement the "average" option in CommonMode.cc)

  int m_evtID;
  uint16_t m_numberOfBoards;
  std::map<int,float> m_meanHGMap;
  std::map<int,float> m_meanLGMap;
  std::map<int,float> m_rmsHGMap;
  std::map<int,float> m_rmsLGMap;
  std::map<int,int> m_counterMap;

  edm::EDGetTokenT<HGCalTBRawHitCollection> m_HGCalTBRawHitCollection;

  HGCalTBTopology IsCellValid;
  HGCalTBCellVertices TheCell;
  std::vector<std::pair<double, double>> CellXY;
  std::pair<double, double> CellCentreXY;
  std::set< std::pair<int,HGCalTBDetId> > setOfConnectedDetId;
  
};

RawHitPlotter::RawHitPlotter(const edm::ParameterSet& iConfig) :
  m_electronicMap(iConfig.getUntrackedParameter<std::string>("ElectronicMap","HGCal/CondObjects/data/map_CERN_Hexaboard_OneLayers_May2017.txt")),
  m_sensorsize(iConfig.getUntrackedParameter<int>("SensorSize",128)),
  m_eventPlotter(iConfig.getUntrackedParameter<bool>("EventPlotter",false)),
  m_subtractCommonMode(iConfig.getUntrackedParameter<bool>("SubtractCommonMode",false))
{
  usesResource("TFileService");
  edm::Service<TFileService> fs;

  m_HGCalTBRawHitCollection = consumes<HGCalTBRawHitCollection>(iConfig.getParameter<edm::InputTag>("InputCollection"));

  m_evtID=0;
  
  std::cout << iConfig.dump() << std::endl;
}


RawHitPlotter::~RawHitPlotter()
{

}

void RawHitPlotter::beginJob()
{
  HGCalCondObjectTextIO io(0);
  edm::FileInPath fip(m_electronicMap);
  if (!io.load(fip.fullPath(), essource_.emap_)) {
    throw cms::Exception("Unable to load electronics map");
  };  
}

void RawHitPlotter::analyze(const edm::Event& event, const edm::EventSetup& setup)
{
  usesResource("TFileService");
  edm::Service<TFileService> fs;

  edm::Handle<HGCalTBRawHitCollection> hits;
  event.getByToken(m_HGCalTBRawHitCollection, hits);
  if( !hits->size() )
    return;

  std::map<int,TH2Poly*>  polyMap;
  if( m_eventPlotter ){
    std::ostringstream os( std::ostringstream::ate );
    os << "Event" << event.id().event();
    TFileDirectory dir = fs->mkdir( os.str().c_str() );
    for(size_t ib = 0; ib<HGCAL_TB_GEOMETRY::NUMBER_OF_HEXABOARD; ib++) {
      for( size_t it=0; it<NUMBER_OF_TIME_SAMPLES; it++ ){
	TH2Poly *h=dir.make<TH2Poly>();
	os.str("");
	os<<"HexaBoard"<<ib<<"_TimeSample"<<it;
	h->SetName(os.str().c_str());
	h->SetTitle(os.str().c_str());
	InitTH2Poly(*h, (int)ib, 0, 0);
	polyMap.insert( std::pair<int,TH2Poly*>(100*ib+it,h) );
      }
    }
  }
  
  CommonMode cm(essource_.emap_); //default is common mode per chip using the median
  cm.Evaluate( hits );
  std::map<int,commonModeNoise> cmMap=cm.CommonModeNoiseMap();

  for( auto hit : *hits ){
    for( size_t it=0; it<NUMBER_OF_TIME_SAMPLES; it++ ){
      float highGain,lowGain;
      HGCalTBElectronicsId eid( essource_.emap_.detId2eid(hit.detid().rawId()) );
      if( !essource_.emap_.existsEId(eid) ) continue;
      if( m_subtractCommonMode ){
  	int iski = eid.iskiroc();
	if( cmMap[iski].fullHG[it]==4 ) continue;
  	float subHG(0),subLG(0);
  	switch ( hit.detid().cellType() ){
  	case 0 : subHG=cmMap[iski].fullHG[it]; subLG=cmMap[iski].fullLG[it]; break;
  	case 2 : subHG=cmMap[iski].halfHG[it]; subLG=cmMap[iski].halfLG[it]; break;
  	case 3 : subHG=cmMap[iski].mouseBiteHG[it]; subLG=cmMap[iski].mouseBiteLG[it]; break;
  	case 4 : subHG=cmMap[iski].outerHG[it]; subLG=cmMap[iski].outerLG[it]; break;
  	}
  	highGain=hit.highGainADC(it)-subHG;
  	lowGain=hit.lowGainADC(it)-subLG;
      }
      else{
  	highGain=hit.highGainADC(it);
  	lowGain=hit.lowGainADC(it);
      }
      int iboard=hit.skiroc()/HGCAL_TB_GEOMETRY::N_SKIROC_PER_HEXA;
      int iski=hit.skiroc();
      int ichan=hit.channel();
      std::pair<int,HGCalTBDetId> p( iboard*1000+(iski%HGCAL_TB_GEOMETRY::N_SKIROC_PER_HEXA)*100+ichan,hit.detid() );
      setOfConnectedDetId.insert(p);
      uint32_t key=iboard*100000+(iski%HGCAL_TB_GEOMETRY::N_SKIROC_PER_HEXA)*10000+ichan*100+it;
      if( m_meanHGMap.find(key)==m_meanHGMap.end() ){
	m_meanHGMap[key]=highGain;
	m_meanLGMap[key]=lowGain;
	m_rmsHGMap[key]=highGain*highGain;
	m_rmsLGMap[key]=lowGain*lowGain;
	m_counterMap[key]=1;
      }
      else{
	m_meanHGMap[key]+=highGain;
	m_meanLGMap[key]+=lowGain;
	m_rmsHGMap[key]+=highGain*highGain;
	m_rmsLGMap[key]+=lowGain*lowGain;
	m_counterMap[key]+=1;
      }
      if(!m_eventPlotter||!IsCellValid.iu_iv_valid(hit.detid().layer(),hit.detid().sensorIU(),hit.detid().sensorIV(),hit.detid().iu(),hit.detid().iv(),m_sensorsize))  continue;
      CellCentreXY=TheCell.GetCellCentreCoordinatesForPlots(hit.detid().layer(),hit.detid().sensorIU(),hit.detid().sensorIV(),hit.detid().iu(),hit.detid().iv(),m_sensorsize);
      double iux = (CellCentreXY.first < 0 ) ? (CellCentreXY.first + delta) : (CellCentreXY.first - delta) ;
      double iuy = (CellCentreXY.second < 0 ) ? (CellCentreXY.second + delta) : (CellCentreXY.second - delta);
      polyMap[ 100*iboard+it ]->Fill(iux/2 , iuy, highGain);
    }
  }
}

void RawHitPlotter::InitTH2Poly(TH2Poly& poly, int layerID, int sensorIU, int sensorIV)
{
  double HexX[MAXVERTICES] = {0.};
  double HexY[MAXVERTICES] = {0.};
  for(int iv = -7; iv < 8; iv++) {
    for(int iu = -7; iu < 8; iu++) {
      if(!IsCellValid.iu_iv_valid(layerID, sensorIU, sensorIV, iu, iv, m_sensorsize)) continue;
      CellXY = TheCell.GetCellCoordinatesForPlots(layerID, sensorIU, sensorIV, iu, iv, m_sensorsize);
      assert(CellXY.size() == 4 || CellXY.size() == 6);
      unsigned int iVertex = 0;
      for(std::vector<std::pair<double, double>>::const_iterator it = CellXY.begin(); it != CellXY.end(); it++) {
	HexX[iVertex] =  it->first;
	HexY[iVertex] =  it->second;
	++iVertex;
      }
      //Somehow cloning of the TH2Poly was not working. Need to look at it. Currently physically booked another one.
      poly.AddBin(CellXY.size(), HexX, HexY);
    }//loop over iu
  }//loop over iv
}


void RawHitPlotter::endJob()
{
  usesResource("TFileService");
  edm::Service<TFileService> fs;
  std::map<int,TH2Poly*>  hgMeanMap;
  std::map<int,TH2Poly*>  lgMeanMap;
  std::map<int,TH2Poly*>  hgRMSMap;
  std::map<int,TH2Poly*>  lgRMSMap;
  std::ostringstream os( std::ostringstream::ate );
  TH2Poly *h;
  for(size_t ib = 0; ib<HGCAL_TB_GEOMETRY::NUMBER_OF_HEXABOARD; ib++) {
    os.str("");
    os << "HexaBoard" << ib ;
    TFileDirectory dir = fs->mkdir( os.str().c_str() );
    TFileDirectory hgpdir = dir.mkdir( "HighGainPedestal" );
    TFileDirectory lgpdir = dir.mkdir( "LowGainPedestal" );
    TFileDirectory hgndir = dir.mkdir( "HighGainNoise" );
    TFileDirectory lgndir = dir.mkdir( "LowGainNoise" );
    for( size_t it=0; it<NUMBER_OF_TIME_SAMPLES; it++ ){
      h=hgpdir.make<TH2Poly>();
      os.str("");
      os<<"TS"<<it;
      h->SetName(os.str().c_str());
      h->SetTitle(os.str().c_str());
      InitTH2Poly(*h, (int)ib, 0, 0);
      hgMeanMap.insert( std::pair<int,TH2Poly*>(100*ib+it,h) );

      h=lgpdir.make<TH2Poly>();
      os.str("");
      os<<"TS"<<it;
      h->SetName(os.str().c_str());
      h->SetTitle(os.str().c_str());
      InitTH2Poly(*h, (int)ib, 0, 0);
      lgMeanMap.insert( std::pair<int,TH2Poly*>(100*ib+it,h) );

      h=hgndir.make<TH2Poly>();
      os.str("");
      os<<"TS"<<it;
      h->SetName(os.str().c_str());
      h->SetTitle(os.str().c_str());
      InitTH2Poly(*h, (int)ib, 0, 0);
      hgRMSMap.insert( std::pair<int,TH2Poly*>(100*ib+it,h) );

      h=lgndir.make<TH2Poly>();
      os.str("");
      os<<"TS"<<it;
      h->SetName(os.str().c_str());
      h->SetTitle(os.str().c_str());
      InitTH2Poly(*h, (int)ib, 0, 0);
      lgRMSMap.insert( std::pair<int,TH2Poly*>(100*ib+it,h) );
    }
  }
  for( std::set< std::pair<int,HGCalTBDetId> >::iterator it=setOfConnectedDetId.begin(); it!=setOfConnectedDetId.end(); ++it ){
    int iboard=(*it).first/1000;
    int iski=((*it).first%1000)/100;
    int ichan=(*it).first%100;
    HGCalTBDetId detid=(*it).second;
    CellCentreXY = TheCell.GetCellCentreCoordinatesForPlots( detid.layer(), detid.sensorIU(), detid.sensorIV(), detid.iu(), detid.iv(), m_sensorsize );
    double iux = (CellCentreXY.first < 0 ) ? (CellCentreXY.first + HGCAL_TB_GEOMETRY::DELTA) : (CellCentreXY.first - HGCAL_TB_GEOMETRY::DELTA) ;
    double iuy = (CellCentreXY.second < 0 ) ? (CellCentreXY.second + HGCAL_TB_GEOMETRY::DELTA) : (CellCentreXY.second - HGCAL_TB_GEOMETRY::DELTA);
    for( size_t it=0; it<NUMBER_OF_TIME_SAMPLES; it++ ){
      int key=iboard*100000+iski*10000+ichan*100+it;
      float hgMean=m_meanHGMap[key]/m_counterMap[key];
      float lgMean=m_meanLGMap[key]/m_counterMap[key];
      float hgRMS=std::sqrt(m_rmsHGMap[key]/m_counterMap[key]-m_meanHGMap[key]/m_counterMap[key]*m_meanHGMap[key]/m_counterMap[key]);
      float lgRMS=std::sqrt(m_rmsLGMap[key]/m_counterMap[key]-m_meanLGMap[key]/m_counterMap[key]*m_meanLGMap[key]/m_counterMap[key]);
      hgMeanMap[ 100*iboard+it ]->Fill(iux/2 , iuy, hgMean );
      lgMeanMap[ 100*iboard+it ]->Fill(iux/2 , iuy, lgMean );
      hgRMSMap[ 100*iboard+it ]->Fill(iux/2 , iuy, hgRMS );
      lgRMSMap[ 100*iboard+it ]->Fill(iux/2 , iuy, lgRMS );
    }
  }
}

void RawHitPlotter::fillDescriptions(edm::ConfigurationDescriptions& descriptions)
{
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

DEFINE_FWK_MODULE(RawHitPlotter);
