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

const static size_t N_SKIROC_PER_HEXA = 4;
//const static size_t HGCAL_TB_GEOMETRY::N_CHANNELS_PER_SKIROC = 64;

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
  std::map<int,TH1F*> m_h_adcHigh;
  std::map<int,TH1F*> m_h_adcLow;
  std::map<int,TH2F*> m_h_pulseHigh;
  std::map<int,TH2F*> m_h_pulseLow;

  std::map<int,TH1F*> m_h_cmHigh;
  std::map<int,TH1F*> m_h_cmLow;

  TH2F* m_h_tot_vs_low[HGCAL_TB_GEOMETRY::NUMBER_OF_HEXABOARD][HGCAL_TB_GEOMETRY::N_SKIROC_PER_HEXA];
  TH2F* m_h_low_vs_high[HGCAL_TB_GEOMETRY::NUMBER_OF_HEXABOARD][HGCAL_TB_GEOMETRY::N_SKIROC_PER_HEXA];

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
  m_subtractCommonMode(iConfig.getUntrackedParameter<bool>("SubtractCommonMode",false)),
  m_commonModeThreshold(iConfig.getUntrackedParameter<double>("CommonModeThreshold",100))
{
  usesResource("TFileService");
  edm::Service<TFileService> fs;

  m_HGCalTBRawHitCollection = consumes<HGCalTBRawHitCollection>(iConfig.getParameter<edm::InputTag>("InputCollection"));

  m_evtID=0;
  
  std::ostringstream os( std::ostringstream::ate );
  TH2F* htmp2;
  TH1F* htmp1;
  for(size_t ib = 0; ib<HGCAL_TB_GEOMETRY::NUMBER_OF_HEXABOARD; ib++) {
    for( size_t iski=0; iski<HGCAL_TB_GEOMETRY::N_SKIROC_PER_HEXA; iski++ ){
      os.str("");os<<"HexaBoard"<<ib<<"_Skiroc"<<iski;
      TFileDirectory dir = fs->mkdir( os.str().c_str() );
      for( size_t it=0; it<NUMBER_OF_TIME_SAMPLES; it++ ){
	os.str("");
	os << "CommonModeHigh_HexaBoard" << ib << "_Chip" << iski << "_Sample" << it ;
	htmp1=dir.make<TH1F>(os.str().c_str(),os.str().c_str(),4000,-500,3500);
	m_h_cmHigh.insert( std::pair<int,TH1F*>(ib*100000+iski*10000+it, htmp1) );
	os.str("");
	os << "CommonModeLow_HexaBoard" << ib << "_Chip" << iski << "_Sample" << it ;
	htmp1=dir.make<TH1F>(os.str().c_str(),os.str().c_str(),4000,-500,3500);
	m_h_cmLow.insert( std::pair<int,TH1F*>(ib*100000+iski*10000+it, htmp1) );
      }
      os.str("");
      os << "ToTVsLowGain_Hexa" << ib << "_Chip" << iski ;
      htmp2=dir.make<TH2F>(os.str().c_str(),os.str().c_str(),2600,-100,2500,1000,0,1000);
      m_h_tot_vs_low[ib][iski]=htmp2;
      os.str("");
      os << "LowGainVsHighGain_Hexa" << ib << "_Chip" << iski ;
      htmp2=dir.make<TH2F>(os.str().c_str(),os.str().c_str(),2600,-100,2500,3500,-500,3000);
      m_h_low_vs_high[ib][iski]=htmp2;
      for( size_t ichan=0; ichan<HGCAL_TB_GEOMETRY::N_CHANNELS_PER_SKIROC; ichan++ ){
	for( size_t it=0; it<NUMBER_OF_TIME_SAMPLES; it++ ){
	  os.str("");
	  os << "HighGain_HexaBoard" << ib << "_Chip" << iski << "_Channel" << ichan << "_Sample" << it ;
	  htmp1=dir.make<TH1F>(os.str().c_str(),os.str().c_str(),4000,-500,3500);
	  m_h_adcHigh.insert( std::pair<int,TH1F*>(ib*100000+iski*10000+ichan*100+it, htmp1) );
	  os.str("");
	  os << "LowGain_HexaBoard" << ib << "_Chip" << iski << "_Channel" << ichan << "_Sample" << it ;
	  htmp1=dir.make<TH1F>(os.str().c_str(),os.str().c_str(),4000,-500,3500);
	  m_h_adcLow.insert( std::pair<int,TH1F*>(ib*100000+iski*10000+ichan*100+it, htmp1) );
	}
	os.str("");
	os << "PulseHighGain_Hexa" << ib << "_Chip" << iski << "_Channel" << ichan;
	htmp2=dir.make<TH2F>(os.str().c_str(),os.str().c_str(),NUMBER_OF_TIME_SAMPLES,0,(NUMBER_OF_TIME_SAMPLES)*25,4000,-500,3500);
	m_h_pulseHigh.insert( std::pair<int,TH2F*>(ib*1000+iski*100+ichan, htmp2) );
	os.str("");
	os << "PulseLowGain_Hexa" << ib << "_Chip" << iski << "_Channel" << ichan;
	htmp2=dir.make<TH2F>(os.str().c_str(),os.str().c_str(),NUMBER_OF_TIME_SAMPLES,0,(NUMBER_OF_TIME_SAMPLES)*25,4000,-500,3500);
	m_h_pulseLow.insert( std::pair<int,TH2F*>(ib*1000+iski*100+ichan, htmp2) );
	os.str("");
      }
    }
  }
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
  for( std::map<int,commonModeNoise>::iterator it=cmMap.begin(); it!=cmMap.end(); ++it ){
    uint16_t iboard=(it->first-1)/HGCAL_TB_GEOMETRY::N_SKIROC_PER_HEXA;
    uint16_t iski;
    if(iboard%2==0)
      iski=(HGCAL_TB_GEOMETRY::N_SKIROC_PER_HEXA-(it->first-1))%HGCAL_TB_GEOMETRY::N_SKIROC_PER_HEXA;
    else
      iski=(HGCAL_TB_GEOMETRY::N_SKIROC_PER_HEXA-(it->first%HGCAL_TB_GEOMETRY::N_SKIROC_PER_HEXA))%HGCAL_TB_GEOMETRY::N_SKIROC_PER_HEXA;
    for( size_t ts=0; ts<NUMBER_OF_TIME_SAMPLES; ts++ ){
      m_h_cmHigh[iboard*100000+iski*10000+ts]->Fill( it->second.fullHG[ts] );
      m_h_cmLow[iboard*100000+iski*10000+ts]->Fill( it->second.fullLG[ts] );
    }
  }

  for( auto hit : *hits ){
    for( size_t it=0; it<NUMBER_OF_TIME_SAMPLES; it++ ){
      float highGain,lowGain;
      if( m_subtractCommonMode ){
	HGCalTBElectronicsId eid( essource_.emap_.detId2eid(hit.detid().rawId()) );
	if( !essource_.emap_.existsEId(eid) ) continue;
  	int iski = eid.iskiroc();
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
      int iboard=hit.skiroc()/N_SKIROC_PER_HEXA;
      int iski=hit.skiroc()%N_SKIROC_PER_HEXA;
      int ichan=hit.channel();
      m_h_adcHigh[iboard*100000+iski*10000+ichan*100+it]->Fill(highGain);
      m_h_adcLow[iboard*100000+iski*10000+ichan*100+it]->Fill(lowGain);
      m_h_pulseHigh[iboard*1000+iski*100+ichan]->Fill(it*25,highGain);
      m_h_pulseLow[iboard*1000+iski*100+ichan]->Fill(it*25,lowGain);
      HGCalTBElectronicsId eid( essource_.emap_.detId2eid(hit.detid().rawId()) );
      if( !essource_.emap_.existsEId(eid) ) continue;
      std::pair<int,HGCalTBDetId> p( iboard*1000+iski*100+ichan,hit.detid() );
      setOfConnectedDetId.insert(p);
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
  TFileDirectory dir = fs->mkdir( "HexagonalPlotter" );
  std::map<int,TH2Poly*>  pedPolyMap;
  std::map<int,TH2Poly*>  pedPolyMapLG;
  std::map<int,TH2Poly*>  noisePolyMap;
  std::map<int,TH2Poly*>  noisePolyMapLG;
  std::map<int,TH2Poly*>  chanMap;
  std::ostringstream os( std::ostringstream::ate );
  TH2Poly *h;
  for(size_t ib = 0; ib<HGCAL_TB_GEOMETRY::NUMBER_OF_HEXABOARD; ib++) {
    for( size_t it=0; it<NUMBER_OF_TIME_SAMPLES; it++ ){
      h=dir.make<TH2Poly>();
      os.str("");
      os<<"HighGain_HexaBoard"<<ib<<"_TimeSample"<<it;
      h->SetName(os.str().c_str());
      h->SetTitle(os.str().c_str());
      InitTH2Poly(*h, (int)ib, 0, 0);
      pedPolyMap.insert( std::pair<int,TH2Poly*>(100*ib+it,h) );
      h=dir.make<TH2Poly>();
      os.str("");
      os<<"LowGain_HexaBoard"<<ib<<"_TimeSample"<<it;
      h->SetName(os.str().c_str());
      h->SetTitle(os.str().c_str());
      InitTH2Poly(*h, (int)ib, 0, 0);
      pedPolyMapLG.insert( std::pair<int,TH2Poly*>(100*ib+it,h) );
      h=dir.make<TH2Poly>();
      os.str("");
      os<<"Noise_HighGain_HexaBoard"<<ib<<"_TimeSample"<<it;
      h->SetName(os.str().c_str());
      h->SetTitle(os.str().c_str());
      InitTH2Poly(*h, (int)ib, 0, 0);
      noisePolyMap.insert( std::pair<int,TH2Poly*>(100*ib+it,h) );
      h=dir.make<TH2Poly>();
      os.str("");
      os<<"Noise_LowGain_HexaBoard"<<ib<<"_TimeSample"<<it;
      h->SetName(os.str().c_str());
      h->SetTitle(os.str().c_str());
      InitTH2Poly(*h, (int)ib, 0, 0);
      noisePolyMapLG.insert( std::pair<int,TH2Poly*>(100*ib+it,h) );
    }
  }

  for( std::set< std::pair<int,HGCalTBDetId> >::iterator it=setOfConnectedDetId.begin(); it!=setOfConnectedDetId.end(); ++it ){
    int iboard=(*it).first/1000;
    int iski=((*it).first%1000)/100;
    int ichan=(*it).first%100;
    HGCalTBDetId detid=(*it).second;
    CellCentreXY = TheCell.GetCellCentreCoordinatesForPlots( detid.layer(), detid.sensorIU(), detid.sensorIV(), detid.iu(), detid.iv(), m_sensorsize );
    double iux = (CellCentreXY.first < 0 ) ? (CellCentreXY.first + delta) : (CellCentreXY.first - delta) ;
    double iuy = (CellCentreXY.second < 0 ) ? (CellCentreXY.second + delta) : (CellCentreXY.second - delta);
    for( size_t it=0; it<NUMBER_OF_TIME_SAMPLES; it++ ){
      pedPolyMap[ 100*iboard+it ]->Fill(iux/2 , iuy, m_h_adcHigh[iboard*100000+iski*10000+ichan*100+it]->GetMean() );
      pedPolyMapLG[ 100*iboard+it ]->Fill(iux/2 , iuy, m_h_adcLow[iboard*100000+iski*10000+ichan*100+it]->GetMean() );
      noisePolyMap[ 100*iboard+it ]->Fill(iux/2 , iuy, m_h_adcHigh[iboard*100000+iski*10000+ichan*100+it]->GetRMS() );
      noisePolyMapLG[ 100*iboard+it ]->Fill(iux/2 , iuy, m_h_adcLow[iboard*100000+iski*10000+ichan*100+it]->GetRMS() );
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
