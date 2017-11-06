#include <memory>
#include <iostream>
#include <fstream>
#include <sstream>
#include <iterator>
#include "TH2Poly.h"
#include "TProfile.h"
#include "TH1F.h"
#include "TH2F.h"
// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/one/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ServiceRegistry/interface/Service.h"

#include "DataFormats/Math/interface/Error.h"
#include "DataFormats/Math/interface/Vector3D.h"
#include "DataFormats/Math/interface/Point3D.h"

#include "HGCal/DataFormats/interface/HGCalTBDetId.h"
#include "HGCal/DataFormats/interface/HGCalTBClusterCollection.h"
#include "HGCal/DataFormats/interface/HGCalTBElectronicsId.h"

#include "HGCal/Geometry/interface/HGCalTBCellVertices.h"
#include "HGCal/Geometry/interface/HGCalTBTopology.h"
#include "HGCal/Geometry/interface/HGCalTBGeometryParameters.h"
#include "HGCal/Geometry/interface/HGCalTBSpillParameters.h"

#include "CommonTools/UtilAlgos/interface/TFileService.h"

#include "HGCal/CondObjects/interface/HGCalElectronicsMap.h"
#include "HGCal/CondObjects/interface/HGCalCondObjectTextIO.h"

#define MAXVERTICES 6
static const double delta = 0.00001;//Add/subtract delta = 0.00001 to x,y of a cell centre so the TH2Poly::Fill doesnt have a problem at the edges where the centre of a half-hex cell passes through the sennsor boundary line.

using namespace std;

class EventDisplay : public edm::one::EDAnalyzer<edm::one::SharedResources>
{

public:
  explicit EventDisplay(const edm::ParameterSet&);
  ~EventDisplay();
  static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);

private:
  virtual void beginJob() override;
  void analyze(const edm::Event& , const edm::EventSetup&) override;
  virtual void endJob() override;
  void InitTH2Poly(TH2Poly& poly, int layerID, int sensorIU, int sensorIV);

  // ----------member data ---------------------------
  //  edm::EDGetToken HGCalTBRecHitCollection_;
  edm::EDGetToken HGCalTBClusterCollection_;
  edm::EDGetToken HGCalTBClusterCollection7_;
  edm::EDGetToken HGCalTBClusterCollection19_;
  std::string mapfile_ = "HGCal/CondObjects/data/map_CERN_8Layers_Sept2016.txt";
  struct {
    HGCalElectronicsMap emap_;
  } essource_;
  int nlayers;
  int sensorsize;
  int _evtID;
  HGCalTBTopology IsCellValid;
  HGCalTBCellVertices TheCell;
  std::vector<std::pair<double, double>> CellXY;
  std::pair<double, double> CellCentreXY;
  edm::Service<TFileService> fs;
};

EventDisplay::EventDisplay(const edm::ParameterSet& iConfig) :
  nlayers( iConfig.getUntrackedParameter<int>("Nlayers",8) ),
  sensorsize( iConfig.getUntrackedParameter<int>("SensorSize",128) )
{
  usesResource("TFileService");
  HGCalTBClusterCollection_ = consumes<reco::HGCalTBClusterCollection>(iConfig.getParameter<edm::InputTag>("HGCALTBCLUSTERS"));
  HGCalTBClusterCollection7_ = consumes<reco::HGCalTBClusterCollection>(iConfig.getParameter<edm::InputTag>("HGCALTBCLUSTERS7"));
  HGCalTBClusterCollection19_ = consumes<reco::HGCalTBClusterCollection>(iConfig.getParameter<edm::InputTag>("HGCALTBCLUSTERS19"));
  _evtID = 0;
}


EventDisplay::~EventDisplay()
{
}

void EventDisplay::InitTH2Poly(TH2Poly& poly, int layerID, int sensorIU, int sensorIV)
{
  double HexX[MAXVERTICES] = {0.};
  double HexY[MAXVERTICES] = {0.};

  for(int iv = -7; iv < 8; iv++) {
    for(int iu = -7; iu < 8; iu++) {
      if(!IsCellValid.iu_iv_valid(layerID, sensorIU, sensorIV, iu, iv, sensorsize)) continue;
      CellXY = TheCell.GetCellCoordinatesForPlots(layerID, sensorIU, sensorIV, iu, iv, sensorsize);
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
//

// ------------ method called for each event  ------------
void
EventDisplay::analyze(const edm::Event& event, const edm::EventSetup& setup)
{

  edm::Handle<reco::HGCalTBClusterCollection> clusters;
  event.getByToken(HGCalTBClusterCollection_, clusters);
  edm::Handle<reco::HGCalTBClusterCollection> clusters7;
  event.getByToken(HGCalTBClusterCollection7_, clusters7);
  edm::Handle<reco::HGCalTBClusterCollection> clusters19;
  event.getByToken(HGCalTBClusterCollection19_, clusters19);

  std::ostringstream os( std::ostringstream::ate );
  TH2Poly *h_RecHit_layer[nlayers];
  TH2Poly *h_Cluster_layer[nlayers];
  TH2Poly *h_Cluster7_layer[nlayers];
  TH2Poly *h_Cluster19_layer[nlayers];
  for( int ilayer=0; ilayer<nlayers; ilayer++ ){
    os.str("");
    os << "Layer" << ilayer << "_Event" << _evtID;
    h_RecHit_layer[ilayer] = fs->make<TH2Poly>();
    h_RecHit_layer[ilayer]->SetName( os.str().c_str() );
    h_RecHit_layer[ilayer]->SetTitle( os.str().c_str() );
    InitTH2Poly(*h_RecHit_layer[ilayer],ilayer,0,0);
    os.str("");
    os << "Layer" << ilayer << "_Event" << _evtID << "_cluster";
    h_Cluster_layer[ilayer] = fs->make<TH2Poly>();
    h_Cluster_layer[ilayer]->SetName( os.str().c_str() );
    h_Cluster_layer[ilayer]->SetTitle( os.str().c_str() );
    InitTH2Poly(*h_Cluster_layer[ilayer],ilayer,0,0);
    os.str("");
    os << "Layer" << ilayer << "_Event" << _evtID << "_cluster7";
    h_Cluster7_layer[ilayer] = fs->make<TH2Poly>();
    h_Cluster7_layer[ilayer]->SetName( os.str().c_str() );
    h_Cluster7_layer[ilayer]->SetTitle( os.str().c_str() );
    InitTH2Poly(*h_Cluster7_layer[ilayer],ilayer,0,0);
    os.str("");
    os << "Layer" << ilayer << "_Event" << _evtID << "_cluster19";
    h_Cluster19_layer[ilayer] = fs->make<TH2Poly>();
    h_Cluster19_layer[ilayer]->SetName( os.str().c_str() );
    h_Cluster19_layer[ilayer]->SetTitle( os.str().c_str() );
    InitTH2Poly(*h_Cluster19_layer[ilayer],ilayer,0,0);

  }
  _evtID++;
  
  float clusterID[nlayers];
  for(int i=0; i<nlayers; i++) 
    clusterID[i]=0.;

  for( auto cluster : *clusters ){
    clusterID[ cluster.layer()-1 ]+=1.0;
    for( std::vector< std::pair<DetId,float> >::const_iterator it=cluster.hitsAndFractions().begin(); it!=cluster.hitsAndFractions().end(); ++it ){
      HGCalTBDetId detID=(*it).first;
      if(!IsCellValid.iu_iv_valid( detID.layer(),
				   detID.sensorIU(), detID.sensorIV(), 
				   detID.iu(), detID.iv(), sensorsize ) 
	 )  continue;
      CellCentreXY = TheCell.GetCellCentreCoordinatesForPlots( detID.layer(), 
							       detID.sensorIU(), detID.sensorIV(), 
							       detID.iu(), detID.iv(), sensorsize );
      double iux = (CellCentreXY.first < 0 ) ? (CellCentreXY.first + delta) : (CellCentreXY.first - delta) ;
      double iuy = (CellCentreXY.second < 0 ) ? (CellCentreXY.second + delta) : (CellCentreXY.second - delta);
      if( detID.cellType()!=1 )
	h_Cluster_layer[ detID.layer()-1 ]->Fill(iux , iuy, clusterID[ detID.layer()-1 ]);
      h_RecHit_layer[ detID.layer()-1 ]->Fill(iux , iuy, (*it).second*cluster.energy());
    }
  }

  for(int i=0; i<nlayers; i++) 
    clusterID[i]=0.;
  for( auto cluster : *clusters7 ){
    clusterID[ cluster.layer()-1 ]+=1.0;
    for( std::vector< std::pair<DetId,float> >::const_iterator it=cluster.hitsAndFractions().begin(); it!=cluster.hitsAndFractions().end(); ++it ){
      HGCalTBDetId detID=(*it).first;
      if(!IsCellValid.iu_iv_valid( detID.layer(),
  				   detID.sensorIU(), detID.sensorIV(), 
  				   detID.iu(), detID.iv(), sensorsize ) 
  	 )  continue;
      CellCentreXY = TheCell.GetCellCentreCoordinatesForPlots( detID.layer(), 
  							       detID.sensorIU(), detID.sensorIV(), 
  							       detID.iu(), detID.iv(), sensorsize );
      double iux = (CellCentreXY.first < 0 ) ? (CellCentreXY.first + delta) : (CellCentreXY.first - delta) ;
      double iuy = (CellCentreXY.second < 0 ) ? (CellCentreXY.second + delta) : (CellCentreXY.second - delta);
      if( detID.cellType()!=1 )
	h_Cluster7_layer[ detID.layer()-1 ]->Fill(iux , iuy, clusterID[ detID.layer()-1 ]);
    }
  }
  
  for(int i=0; i<nlayers; i++) 
    clusterID[i]=0.;
  for( auto cluster : *clusters19 ){
    clusterID[ cluster.layer()-1 ]+=1.0;
    for( std::vector< std::pair<DetId,float> >::const_iterator it=cluster.hitsAndFractions().begin(); it!=cluster.hitsAndFractions().end(); ++it ){
      HGCalTBDetId detID=(*it).first;
      if(!IsCellValid.iu_iv_valid( detID.layer(),
  				   detID.sensorIU(), detID.sensorIV(), 
  				   detID.iu(), detID.iv(), sensorsize ) 
  	 )  continue;
      CellCentreXY = TheCell.GetCellCentreCoordinatesForPlots( detID.layer(), 
  							       detID.sensorIU(), detID.sensorIV(), 
  							       detID.iu(), detID.iv(), sensorsize );
      double iux = (CellCentreXY.first < 0 ) ? (CellCentreXY.first + delta) : (CellCentreXY.first - delta) ;
      double iuy = (CellCentreXY.second < 0 ) ? (CellCentreXY.second + delta) : (CellCentreXY.second - delta);
      if( detID.cellType()!=1 )
	h_Cluster19_layer[ detID.layer()-1 ]->Fill(iux , iuy, clusterID[ detID.layer()-1 ]);
    }
  }

}

void
EventDisplay::beginJob()
{
  HGCalCondObjectTextIO io(0);
  edm::FileInPath fip(mapfile_);
  if (!io.load(fip.fullPath(), essource_.emap_)) {
    throw cms::Exception("Unable to load electronics map");
  }
}

void
EventDisplay::endJob()
{

}

void
EventDisplay::fillDescriptions(edm::ConfigurationDescriptions& descriptions)
{
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(EventDisplay);
