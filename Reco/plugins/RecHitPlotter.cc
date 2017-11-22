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
#include "HGCal/CondObjects/interface/HGCalTBDetectorLayout.h"
#include "HGCal/DataFormats/interface/HGCalTBElectronicsId.h"
#include "HGCal/Geometry/interface/HGCalTBCellVertices.h"
#include "HGCal/Geometry/interface/HGCalTBTopology.h"
#include "HGCal/Geometry/interface/HGCalTBGeometryParameters.h"

#include <iomanip>
#include <set>

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
  void InitTH2Poly(TH2Poly& poly, int det, int layerID);

  std::string m_electronicMap;
  std::string m_detectorLayoutFile;

  struct {
    HGCalElectronicsMap emap_;
    HGCalTBDetectorLayout layout_;
  } essource_;

  int m_sensorsize;
  bool m_eventPlotter;
  int m_evtID;
  double m_noiseThreshold;

  TH1F* m_h_hgSum;
  TH1F* m_h_lgSum;
  TH1F* m_h_enSum;

  edm::EDGetTokenT<HGCalTBRecHitCollection> m_HGCalTBRecHitCollection;

  HGCalTBTopology IsCellValid;
  HGCalTBCellVertices TheCell;
  std::vector<std::pair<double, double>> CellXY;
  std::pair<double, double> CellCentreXY;
};

RecHitPlotter::RecHitPlotter(const edm::ParameterSet& iConfig) :
  m_electronicMap(iConfig.getUntrackedParameter<std::string>("ElectronicMap","HGCal/CondObjects/data/map_CERN_Hexaboard_28Layers_AllFlipped.txt")),
  m_detectorLayoutFile(iConfig.getUntrackedParameter<std::string>("DetectorLayout","HGCal/CondObjects/data/layerGeom_oct2017_h2_17layers.txt")),
  m_sensorsize(iConfig.getUntrackedParameter<int>("SensorSize",128)),
  m_eventPlotter(iConfig.getUntrackedParameter<bool>("EventPlotter",false)),
  m_noiseThreshold(iConfig.getUntrackedParameter<double>("NoiseThreshold",0.5))
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
  fip=edm::FileInPath(m_detectorLayoutFile);
  if (!io.load(fip.fullPath(), essource_.layout_)) {
    throw cms::Exception("Unable to load detector layout file");
  };
 
  usesResource("TFileService");
  edm::Service<TFileService> fs;
  
  m_h_hgSum=fs->make<TH1F>("HighGainSum","HighGainSum",5000,0,10000);
  m_h_lgSum=fs->make<TH1F>("LowGainSum","LowGainSum",5000,0,10000);
  m_h_enSum=fs->make<TH1F>("EnergySum","EnergySum",5000,0,10000);

}

void RecHitPlotter::analyze(const edm::Event& event, const edm::EventSetup& setup)
{
  usesResource("TFileService");
  edm::Service<TFileService> fs;

  edm::Handle<HGCalTBRecHitCollection> hits;
  event.getByToken(m_HGCalTBRecHitCollection, hits);

  std::map<int,TH2Poly*>  polyMap;
  std::map<int,TH2Poly*>  polyMapHG;
  std::map<int,TH2Poly*>  polyMapLG;
  if( m_eventPlotter ){
    std::ostringstream os( std::ostringstream::ate );
    os << "Event" << event.id().event();
    TFileDirectory dir = fs->mkdir( os.str().c_str() );
    for(int il = 0; il<essource_.layout_.nlayers(); il++) {
      int subdetId = essource_.layout_.at(il).subdet();
      TH2Poly *h=dir.make<TH2Poly>();
      os.str("");
      os<<"Energy_Layer"<<il;
      h->SetName(os.str().c_str());
      h->SetTitle(os.str().c_str());
      InitTH2Poly(*h, subdetId, il);
      polyMap.insert( std::pair<int,TH2Poly*>(il,h) );
      h=dir.make<TH2Poly>();
      os.str("");
      os<<"HighGain_Layer"<<il;
      h->SetName(os.str().c_str());
      h->SetTitle(os.str().c_str());
      InitTH2Poly(*h, subdetId, il);
      polyMapHG.insert( std::pair<int,TH2Poly*>(il,h) );
      h=dir.make<TH2Poly>();
      os.str("");
      os<<"LowGain_Layer"<<il;
      h->SetName(os.str().c_str());
      h->SetTitle(os.str().c_str());
      InitTH2Poly(*h, subdetId, il);
      polyMapLG.insert( std::pair<int,TH2Poly*>(il,h) );
    }
  }
  float energyHighSum(0),energyLowSum(0),energySum(0);
  for( auto hit : *hits ){
    HGCalTBElectronicsId eid( essource_.emap_.detId2eid( hit.id().rawId() ) );
    if( hit.energy()>m_noiseThreshold && !hit.isUnderSaturationForLowGain() && !hit.isUnderSaturationForHighGain() && hit.id().cellType()!=5 ){
      energyHighSum+=hit.energyHigh();
      energyLowSum+=hit.energyLow();
      energySum+=hit.energy();
    }
    if(m_eventPlotter){
      if(!IsCellValid.iu_iv_valid(hit.id().layer(),hit.id().sensorIU(),hit.id().sensorIV(),hit.id().iu(),hit.id().iv(),m_sensorsize))  continue;
      CellCentreXY=TheCell.GetCellCentreCoordinatesForPlots(hit.id().layer(),hit.id().sensorIU(),hit.id().sensorIV(),hit.id().iu(),hit.id().iv(),m_sensorsize);
      polyMap[ hit.id().layer()-1 ]->Fill(CellCentreXY.first , CellCentreXY.second, hit.energy());
      polyMapHG[ hit.id().layer()-1 ]->Fill(CellCentreXY.first , CellCentreXY.second, hit.energyHigh());
      polyMapLG[ hit.id().layer()-1 ]->Fill(CellCentreXY.first , CellCentreXY.second, hit.energyLow());
    }
  }
  m_h_hgSum->Fill( energyHighSum );
  m_h_lgSum->Fill( energyLowSum );
  m_h_enSum->Fill( energySum );
}

void RecHitPlotter::InitTH2Poly(TH2Poly& poly, int det, int layerID)
{
  double HexX[HGCAL_TB_GEOMETRY::MAXVERTICES] = {0.};
  double HexY[HGCAL_TB_GEOMETRY::MAXVERTICES] = {0.};
  if( det==0 ){
    for(int iv = -7; iv < 8; iv++) {
      for(int iu = -7; iu < 8; iu++) {
	if(!IsCellValid.iu_iv_valid(layerID, 0, 0, iu, iv, m_sensorsize)) 
	  continue;
	CellXY = TheCell.GetCellCoordinatesForPlots(layerID, 0, 0, iu, iv, m_sensorsize);
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
  else if( det==1 ){
    for(int sensorIV = -1; sensorIV <= 1; sensorIV++){
      for(int sensorIU = -1; sensorIU <= 1; sensorIU++){
	for(int iv = -7; iv < 8; iv++) {
	  for(int iu = -7; iu < 8; iu++) {
	    if(!IsCellValid.iu_iv_valid(layerID, sensorIU, sensorIV, iu, iv, m_sensorsize)) 
	      continue;
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
    }// loop over Sensor_Iu ends here
  }// loop over Sensor_Iv ends here
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
