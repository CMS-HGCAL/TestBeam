#include <iostream>
#include "TTree.h"
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
#include "HGCal/DataFormats/interface/HGCalTBClusterCollection.h"
#include "HGCal/DataFormats/interface/HGCalTBDetId.h"
#include "HGCal/CondObjects/interface/HGCalCondObjectTextIO.h"
#include "HGCal/CondObjects/interface/HGCalTBDetectorLayout.h"
#include "HGCal/Geometry/interface/HGCalTBCellVertices.h"
#include "HGCal/Geometry/interface/HGCalTBTopology.h"
#include "HGCal/Geometry/interface/HGCalTBGeometryParameters.h"

#include <iomanip>
#include <set>

class ShowerAnalyzer : public edm::one::EDAnalyzer<edm::one::SharedResources>
{
public:
  explicit ShowerAnalyzer(const edm::ParameterSet&);
  ~ShowerAnalyzer();
  static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);
private:
  virtual void beginJob() override;
  virtual void endJob() override;
  void analyze(const edm::Event& , const edm::EventSetup&) override;

  std::string m_detectorLayoutFile;

  struct {
    HGCalTBDetectorLayout layout_;
  } essource_;

  int m_sensorsize;
  int m_evtID;
  double m_noiseThreshold;

  edm::EDGetTokenT<HGCalTBRecHitCollection> m_HGCalTBRecHitCollection;
  edm::EDGetTokenT<reco::HGCalTBClusterCollection> m_HGCalTBClusterCollection;
  edm::EDGetTokenT<reco::HGCalTBClusterCollection> m_HGCalTBClusterCollection7;
  edm::EDGetTokenT<reco::HGCalTBClusterCollection> m_HGCalTBClusterCollection19;

  HGCalTBTopology IsCellValid;
  HGCalTBCellVertices TheCell;
  std::vector<std::pair<double, double>> CellXY;
  std::pair<double, double> CellCentreXY;


  TTree* m_tree;
  std::vector<double> m_energyLayer;
  std::vector<double> m_energyLayerHG;
  std::vector<double> m_e1Layer;
  std::vector<double> m_e7Layer;
  std::vector<double> m_e19Layer;
  std::vector<int> m_n7Layer;
  std::vector<int> m_n19Layer;
  std::vector<double> m_x19Layer;
  std::vector<double> m_y19Layer;
  double m_energyRing[20];
  double m_energyEE;
  double m_energyFH;
  double m_energyTot;
  std::vector<int> m_detIdEmax;
  std::vector<int> m_nCellsAboveThr;
  std::vector<int> m_nCellsHGUnderSat;
  std::vector<int> m_nCellsLGUnderSat;
};

ShowerAnalyzer::ShowerAnalyzer(const edm::ParameterSet& iConfig) :
  m_detectorLayoutFile(iConfig.getUntrackedParameter<std::string>("DetectorLayout","HGCal/CondObjects/data/layerGeom_oct2017_h2_17layers.txt")),
  m_sensorsize(iConfig.getUntrackedParameter<int>("SensorSize",128)),
  m_noiseThreshold(iConfig.getUntrackedParameter<double>("NoiseThreshold",0.5))
{
  m_HGCalTBRecHitCollection = consumes<HGCalTBRecHitCollection>(iConfig.getParameter<edm::InputTag>("InputCollection"));
  m_HGCalTBClusterCollection = consumes<reco::HGCalTBClusterCollection>(iConfig.getParameter<edm::InputTag>("InputClusterCollection"));
  m_HGCalTBClusterCollection7 = consumes<reco::HGCalTBClusterCollection>(iConfig.getParameter<edm::InputTag>("InputCluster7Collection"));
  m_HGCalTBClusterCollection19 = consumes<reco::HGCalTBClusterCollection>(iConfig.getParameter<edm::InputTag>("InputCluster19Collection"));

  m_evtID=0;

  usesResource("TFileService");
  edm::Service<TFileService> fs;
  
  m_tree = fs->make<TTree>("tree","");
  m_tree->Branch("energyTot",&m_energyTot);
  m_tree->Branch("energyEE",&m_energyEE);
  m_tree->Branch("energyFH",&m_energyFH);
  m_tree->Branch("energyLayer","std::vector<double>",&m_energyLayer);
  m_tree->Branch("e1Layer","std::vector<double>",&m_e1Layer);
  m_tree->Branch("e7Layer","std::vector<double>",&m_e7Layer);
  m_tree->Branch("e19Layer","std::vector<double>",&m_e19Layer);
  m_tree->Branch("x19Layer","std::vector<double>",&m_x19Layer);
  m_tree->Branch("y19Layer","std::vector<double>",&m_y19Layer);
  m_tree->Branch("n7Layer","std::vector<int>",&m_n7Layer);
  m_tree->Branch("n19Layer","std::vector<int>",&m_n19Layer);
  m_tree->Branch("energyRing",&m_energyRing,"energyRing[20]/D");
  m_tree->Branch("detIdEmax","std::vector<int>",&m_detIdEmax);
  m_tree->Branch("nCellsAboveThr","std::vector<int>",&m_nCellsAboveThr);
  m_tree->Branch("nCellsHGUnderSat","std::vector<int>",&m_nCellsHGUnderSat);
  m_tree->Branch("nCellsLGUnderSat","std::vector<int>",&m_nCellsLGUnderSat);
  std::cout << iConfig.dump() << std::endl;
}


ShowerAnalyzer::~ShowerAnalyzer()
{

}

void ShowerAnalyzer::beginJob()
{
  HGCalCondObjectTextIO io(0);
  edm::FileInPath fip(m_detectorLayoutFile);
  if (!io.load(fip.fullPath(), essource_.layout_)) {
    throw cms::Exception("Unable to load detector layout file");
  };
}

void ShowerAnalyzer::analyze(const edm::Event& event, const edm::EventSetup& setup)
{
  m_energyTot = m_energyEE = m_energyFH = 0.0;
  m_energyLayer.clear();
  m_energyLayerHG.clear();
  m_e1Layer.clear();
  m_e7Layer.clear();
  m_e19Layer.clear();
  m_x19Layer.clear();
  m_y19Layer.clear();
  m_n7Layer.clear();
  m_n19Layer.clear();
  m_detIdEmax.clear();
  m_nCellsAboveThr.clear();
  m_nCellsHGUnderSat.clear();
  m_nCellsLGUnderSat.clear();
  std::vector<HGCalTBRecHit> maxEnergyHits;
  for( int ilayer=0; ilayer<essource_.layout_.nlayers(); ilayer++ ){
    m_energyLayer.push_back(0.0);
    m_energyLayerHG.push_back(0.0);
    m_e1Layer.push_back(0.0);
    m_e7Layer.push_back(0.0);
    m_e19Layer.push_back(0.0);
    m_x19Layer.push_back(0.0);
    m_y19Layer.push_back(0.0);
    m_n7Layer.push_back(0);
    m_n19Layer.push_back(0);
    m_detIdEmax.push_back(0);
    m_nCellsAboveThr.push_back(0);
    m_nCellsHGUnderSat.push_back(0);
    m_nCellsLGUnderSat.push_back(0);
    HGCalTBRecHit emptyHit;
    maxEnergyHits.push_back(emptyHit);
  }
  for( int i=0; i<20; i++ )
    m_energyRing[i]=0.0;

  edm::Handle<HGCalTBRecHitCollection> hits;
  event.getByToken(m_HGCalTBRecHitCollection, hits);
  edm::Handle<reco::HGCalTBClusterCollection> clusters;
  event.getByToken(m_HGCalTBClusterCollection, clusters);
  edm::Handle<reco::HGCalTBClusterCollection> clusters7;
  event.getByToken(m_HGCalTBClusterCollection7, clusters7);
  edm::Handle<reco::HGCalTBClusterCollection> clusters19;
  event.getByToken(m_HGCalTBClusterCollection19, clusters19);


  for( auto cluster : *clusters7 ){
    m_e7Layer.at( cluster.layer()-1 ) = cluster.energy();
    m_n7Layer.at( cluster.layer()-1 ) = cluster.size();
  }
  for( auto cluster : *clusters19 ){
    m_e19Layer.at( cluster.layer()-1 ) = cluster.energy();
    m_n19Layer.at( cluster.layer()-1 ) = cluster.size();
    m_x19Layer[cluster.layer()-1] = cluster.x();
    m_y19Layer[cluster.layer()-1] = cluster.y();
  }
  HGCalTBCellVertices cellVertice;
  std::pair<double, double> CellCentreXY;
  for( auto hit : *hits ){
    int layerId=hit.id().layer()-1;
    if( hit.isUnderSaturationForLowGain() )
      m_nCellsLGUnderSat.at(layerId)++;
    if( hit.isUnderSaturationForHighGain() )
      m_nCellsHGUnderSat.at(layerId)++;
    if( hit.energy()>m_noiseThreshold && hit.id().cellType()!=5 ){
      m_nCellsAboveThr.at(layerId)++;
      if( hit.energy()>maxEnergyHits.at(layerId).energy() )
	maxEnergyHits.at(layerId)=hit;

      m_energyTot+=hit.energy();
      if( essource_.layout_.at( layerId ).subdet()==0 )
	m_energyEE+=hit.energy();
      else
	m_energyFH+=hit.energy();
      
      m_energyLayer.at( layerId )+=hit.energy();
      m_energyLayerHG.at( layerId )+=hit.energyHigh();

      CellCentreXY=cellVertice.GetCellCentreCoordinatesForPlots( hit.id().layer(), hit.id().sensorIU(), hit.id().sensorIV(), hit.id().iu(), hit.id().iv(), m_sensorsize);
      math::XYZPoint hitpos = math::XYZPoint( CellCentreXY.first,CellCentreXY.second,0 );
      math::XYZPoint cog = math::XYZPoint( m_x19Layer[hit.id().layer()-1],m_y19Layer[hit.id().layer()-1],0 );
      math::XYZPoint diff=math::XYZPoint( hitpos.x()-cog.x(),hitpos.y()-cog.y(),hitpos.z()-cog.z() );
      int ring=(int)(std::sqrt(diff.mag2())/HGCAL_TB_CELL::FULL_CELL_SIDE/std::sqrt(3));
      if(ring<20)
	m_energyRing[ring]+=hit.energy();
    }
  }
  for( int ilayer=0; ilayer<essource_.layout_.nlayers(); ilayer++ ){
    m_e1Layer.at(ilayer)=maxEnergyHits.at(ilayer).energy();
    m_detIdEmax.at(ilayer)=maxEnergyHits.at(ilayer).id();
  }

  m_tree->Fill();
}

void ShowerAnalyzer::endJob()
{
}

void ShowerAnalyzer::fillDescriptions(edm::ConfigurationDescriptions& descriptions)
{
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

DEFINE_FWK_MODULE(ShowerAnalyzer);
