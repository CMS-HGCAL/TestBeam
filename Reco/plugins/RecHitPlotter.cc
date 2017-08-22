#include <iostream>
#include "TH1F.h"
#include "TH2F.h"
#include "TH2Poly.h"
#include <fstream>
#include <sstream>

#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/one/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ServiceRegistry/interface/Service.h"

#include "CommonTools/UtilAlgos/interface/TFileService.h"

#include "HGCal/DataFormats/interface/HGCalTBRecHitCollections.h"
#include "HGCal/DataFormats/interface/HGCalTBDetId.h"
#include "HGCal/CondObjects/interface/HGCalElectronicsMap.h"
#include "HGCal/CondObjects/interface/HGCalCondObjectTextIO.h"
#include "HGCal/DataFormats/interface/HGCalTBElectronicsId.h"
#include "HGCal/Geometry/interface/HGCalTBCellVertices.h"
#include "HGCal/Geometry/interface/HGCalTBTopology.h"
#include "HGCal/Geometry/interface/HGCalTBGeometryParameters.h"

#include <iomanip>
#include <set>

struct channelInfo{
  channelInfo(){;}
  int key;
  TH1F* h_adcHigh;
  TH1F* h_adcLow;
  TH2F* h_high_vs_low;
};

class RecHitPlotter : public edm::one::EDAnalyzer<edm::one::SharedResources>
{
public:
  explicit RecHitPlotter(const edm::ParameterSet&);
  ~RecHitPlotter();
  static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);
private:
  virtual void beginJob() override;
  virtual void endJob() override;
  void analyze(const edm::Event& , const edm::EventSetup&) override;
  void InitTH2Poly(TH2Poly& poly, int layerID, int sensorIU, int sensorIV);

  std::string m_electronicMap;

  struct {
    HGCalElectronicsMap emap_;
  } essource_;

  int m_sensorsize;
  bool m_eventPlotter;
  int m_evtID;
  double m_mipThreshold;
  double m_noiseThreshold;

  TH1F* m_h_mip[4];
  TH1F* m_h_hgSum;
  TH1F* m_h_lgSum;
  TH1F* m_h_enSum;
  std::map<uint32_t,channelInfo*> m_channelInfoMap;

  edm::EDGetTokenT<HGCalTBRecHitCollection> m_HGCalTBRecHitCollection;

  HGCalTBTopology IsCellValid;
  HGCalTBCellVertices TheCell;
  std::vector<std::pair<double, double>> CellXY;
  std::pair<double, double> CellCentreXY;
};

RecHitPlotter::RecHitPlotter(const edm::ParameterSet& iConfig) :
  m_electronicMap(iConfig.getUntrackedParameter<std::string>("ElectronicMap","HGCal/CondObjects/data/map_CERN_Hexaboard_28Layers_AllFlipped.txt")),
  m_sensorsize(iConfig.getUntrackedParameter<int>("SensorSize",128)),
  m_eventPlotter(iConfig.getUntrackedParameter<bool>("EventPlotter",false)),
  m_mipThreshold(iConfig.getUntrackedParameter<double>("MipThreshold",200.)),
  m_noiseThreshold(iConfig.getUntrackedParameter<double>("NoiseThreshold",10.))
{
  m_HGCalTBRecHitCollection = consumes<HGCalTBRecHitCollection>(iConfig.getParameter<edm::InputTag>("InputCollection"));

  m_evtID=0;
  
  std::cout << iConfig.dump() << std::endl;
}


RecHitPlotter::~RecHitPlotter()
{

}

void RecHitPlotter::beginJob()
{
  HGCalCondObjectTextIO io(0);
  edm::FileInPath fip(m_electronicMap);
  if (!io.load(fip.fullPath(), essource_.emap_)) {
    throw cms::Exception("Unable to load electronics map");
  };
  
  usesResource("TFileService");
  edm::Service<TFileService> fs;
  
  m_h_mip[0]=fs->make<TH1F>("Mip_Ski0","Mip_Ski0",180,20,200);
  m_h_mip[1]=fs->make<TH1F>("Mip_Ski1","Mip_Ski1",180,20,200);
  m_h_mip[2]=fs->make<TH1F>("Mip_Ski2","Mip_Ski2",180,20,200);
  m_h_mip[3]=fs->make<TH1F>("Mip_Ski3","Mip_Ski3",180,20,200);
  m_h_hgSum=fs->make<TH1F>("HighGainSum","HighGainSum",5000,0,50000);
  m_h_lgSum=fs->make<TH1F>("LowGainSum","LowGainSum",5000,0,50000);
  m_h_enSum=fs->make<TH1F>("EnergySum","EnergySum",5000,0,50000);

  std::ostringstream os( std::ostringstream::ate );
  for(size_t ib = 0; ib<HGCAL_TB_GEOMETRY::NUMBER_OF_HEXABOARD; ib++) {
    for( size_t iski=0; iski<HGCAL_TB_GEOMETRY::N_SKIROC_PER_HEXA; iski++ ){
      os.str("");os<<"HexaBoard"<<ib<<"_Skiroc"<<iski;
      TFileDirectory dir = fs->mkdir( os.str().c_str() );
      for( size_t ichan=0; ichan<HGCAL_TB_GEOMETRY::N_CHANNELS_PER_SKIROC; ichan++ ){
	int skiId=ib*HGCAL_TB_GEOMETRY::N_SKIROC_PER_HEXA+(HGCAL_TB_GEOMETRY::N_SKIROC_PER_HEXA-iski)%HGCAL_TB_GEOMETRY::N_SKIROC_PER_HEXA+1;
	HGCalTBElectronicsId eid(skiId,ichan);      
	if( !essource_.emap_.existsEId(eid.rawId()) ) continue;
	channelInfo *chan=new channelInfo();
	chan->key=ib*1000+iski*100+ichan;
	os.str("");
	os << "HighGain_Channel" << ichan  ;
	chan->h_adcHigh=dir.make<TH1F>(os.str().c_str(),os.str().c_str(),3600,-100,3500);
	os.str("");
	os << "LowGain_Channel" << ichan ;
	chan->h_adcLow=dir.make<TH1F>(os.str().c_str(),os.str().c_str(),3600,-100,3500);
	os.str("");
	os << "HighGainVsLowGain_HexaBoard" << ib << "_Chip" << iski << "_Channel" << ichan;
	chan->h_high_vs_low=dir.make<TH2F>(os.str().c_str(),os.str().c_str(),360,-100,3500,400,-500,3500);
	m_channelInfoMap.insert( std::pair<uint32_t,channelInfo*>(chan->key,chan) );
      }
    }
  }
}

void RecHitPlotter::analyze(const edm::Event& event, const edm::EventSetup& setup)
{
  usesResource("TFileService");
  edm::Service<TFileService> fs;

  edm::Handle<HGCalTBRecHitCollection> hits;
  event.getByToken(m_HGCalTBRecHitCollection, hits);

  std::map<int,TH2Poly*>  polyMapHG;
  std::map<int,TH2Poly*>  polyMapLG;
  if( m_eventPlotter ){
    std::ostringstream os( std::ostringstream::ate );
    os << "Event" << event.id().event();
    TFileDirectory dir = fs->mkdir( os.str().c_str() );
    for(size_t ib = 0; ib<HGCAL_TB_GEOMETRY::NUMBER_OF_HEXABOARD; ib++) {
      TH2Poly *h=dir.make<TH2Poly>();
      os.str("");
      os<<"HexaBoard"<<ib<<"_HighGain";
      h->SetName(os.str().c_str());
      h->SetTitle(os.str().c_str());
      InitTH2Poly(*h, (int)ib, 0, 0);
      polyMapHG.insert( std::pair<int,TH2Poly*>(ib,h) );
      h=dir.make<TH2Poly>();
      os.str("");
      os<<"HexaBoard"<<ib<<"_LowGain";
      h->SetName(os.str().c_str());
      h->SetTitle(os.str().c_str());
      InitTH2Poly(*h, (int)ib, 0, 0);
      polyMapLG.insert( std::pair<int,TH2Poly*>(ib,h) );
    }
  }
  
  float energyHighSum(0),energyLowSum(0),energySum(0);
  for( auto hit : *hits ){
    HGCalTBElectronicsId eid( essource_.emap_.detId2eid( hit.id().rawId() ) );
    int iski=(HGCAL_TB_GEOMETRY::N_SKIROC_PER_HEXA-(eid.iskiroc()-1)%HGCAL_TB_GEOMETRY::N_SKIROC_PER_HEXA)%HGCAL_TB_GEOMETRY::N_SKIROC_PER_HEXA;
    int key=(hit.id().layer()-1)*1000+iski*100+eid.ichan();
    if( m_channelInfoMap.find(key)==m_channelInfoMap.end() )
      std::cout << key << "\t" << hit.id().layer()<< " " << iski << " " << eid.ichan() << " "
		<< (eid.iskiroc()-1)%HGCAL_TB_GEOMETRY::N_SKIROC_PER_HEXA 
		<< std::endl;
    if( !hit.isUnderSaturationForHighGain() ) m_channelInfoMap[ key ]->h_adcHigh->Fill( hit.energyHigh() );
    if( !hit.isUnderSaturationForLowGain() ) m_channelInfoMap[ key ]->h_adcLow->Fill( hit.energyLow() );
    if( hit.energyHigh()>m_noiseThreshold && !hit.isUnderSaturationForLowGain() && !hit.isUnderSaturationForHighGain() ){
      m_channelInfoMap[ key ]->h_high_vs_low->Fill( hit.energyLow(), hit.energyHigh() );
      energyHighSum+=hit.energyHigh();
      energyLowSum+=hit.energyLow();
      energySum+=hit.energy();
      //if( hit.energyHigh()<m_mipThreshold )
      //	m_h_mip[(4-(eid.iskiroc()-1))%4]->Fill(hit.energyHigh());
    }
    if(m_eventPlotter){
      if(!IsCellValid.iu_iv_valid(hit.id().layer(),hit.id().sensorIU(),hit.id().sensorIV(),hit.id().iu(),hit.id().iv(),m_sensorsize))  continue;
      CellCentreXY=TheCell.GetCellCentreCoordinatesForPlots(hit.id().layer(),hit.id().sensorIU(),hit.id().sensorIV(),hit.id().iu(),hit.id().iv(),m_sensorsize);
      double iux = (CellCentreXY.first < 0 ) ? (CellCentreXY.first + HGCAL_TB_GEOMETRY::DELTA) : (CellCentreXY.first - HGCAL_TB_GEOMETRY::DELTA) ;
      double iuy = (CellCentreXY.second < 0 ) ? (CellCentreXY.second + HGCAL_TB_GEOMETRY::DELTA) : (CellCentreXY.second - HGCAL_TB_GEOMETRY::DELTA);
      polyMapHG[ (eid.iskiroc()-1)/HGCAL_TB_GEOMETRY::N_SKIROC_PER_HEXA ]->Fill(iux , iuy, hit.energyHigh());
      polyMapLG[ (eid.iskiroc()-1)/HGCAL_TB_GEOMETRY::N_SKIROC_PER_HEXA ]->Fill(iux , iuy, hit.energyLow());
    }
  }
  m_h_hgSum->Fill( energyHighSum );
  m_h_lgSum->Fill( energyLowSum );
  m_h_enSum->Fill( energySum );
}

void RecHitPlotter::InitTH2Poly(TH2Poly& poly, int layerID, int sensorIU, int sensorIV)
{
  double HexX[HGCAL_TB_GEOMETRY::MAXVERTICES] = {0.};
  double HexY[HGCAL_TB_GEOMETRY::MAXVERTICES] = {0.};
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

void RecHitPlotter::endJob()
{
}

void RecHitPlotter::fillDescriptions(edm::ConfigurationDescriptions& descriptions)
{
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

DEFINE_FWK_MODULE(RecHitPlotter);
