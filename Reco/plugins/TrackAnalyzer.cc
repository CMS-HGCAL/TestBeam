#include <memory>
#include <iostream>
#include "TTree.h"
#include "TH1F.h"
#include <sstream>
#include <cmath>
#include <limits>

#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/one/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ServiceRegistry/interface/Service.h"

#include "CommonTools/UtilAlgos/interface/TFileService.h"

#include "DataFormats/Math/interface/Point3D.h"

#include "HGCal/DataFormats/interface/HGCalTBRecHitCollections.h"
#include "HGCal/DataFormats/interface/HGCalTBDetId.h"
#include "HGCal/DataFormats/interface/HGCalTBCaloTrackCollection.h"
#include "HGCal/DataFormats/interface/HGCalTBElectronicsId.h"

#include "HGCal/CondObjects/interface/HGCalElectronicsMap.h"
#include "HGCal/CondObjects/interface/HGCalCondObjectTextIO.h"

#include "HGCal/Geometry/interface/HGCalTBCellVertices.h"
#include "HGCal/Reco/interface/HGCalTBCaloTrackingUtil.h"

using namespace std;

class TrackAnalyzer : public edm::one::EDAnalyzer<edm::one::SharedResources>
{
public:
  explicit TrackAnalyzer(const edm::ParameterSet&);
  ~TrackAnalyzer();
  static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);

private:
  virtual void beginJob() override;
  void analyze(const edm::Event& , const edm::EventSetup&) override;
  virtual void endJob() override;

  // ----------member data ---------------------------
  edm::EDGetToken HGCalTBRecHitCollection_;
  edm::EDGetToken HGCalTBCaloTrackCollection_;
  int _nlayers;
  int _nskirocsperlayer;
  int _nchannelsperskiroc;
  int _sensorSize;
  double _maxChi2;
  double _noiseEnergyThreshold;
  std::vector<double> _layerZPositions;

  TTree* tree;
  int _evtID;
  float _chi2;
  float _theta;
  float _phi;
  int _nhit;
  int _tracknhit;
  float _vx0;
  float _vy0;
  float _px;
  float _py;
  bool _trackSuccess;

  std::vector<double> _deltas;

  std::vector<int> _nhitlayer;
  std::vector<double> _energylayer;

  std::map<int,TH1F*> h_delta_layer;
  std::map<int, TH1F*> h_mip_map;
  std::map<int, TH1F*> h_noise_map;

  const float PI = 3.1415927;

  std::string mapfile_ = "HGCal/CondObjects/data/map_CERN_8Layers_Sept2016.txt";
  struct {
    HGCalElectronicsMap emap_;
  } essource_;

};


TrackAnalyzer::TrackAnalyzer(const edm::ParameterSet& iConfig) :
  _nlayers( iConfig.getUntrackedParameter<int>("NLayers",8) ),
  _nskirocsperlayer( iConfig.getUntrackedParameter<int>("NSkirocsPerLayer",2) ),
  _nchannelsperskiroc( iConfig.getUntrackedParameter<int>("NChannelsPerSkiroc",64) ),
  _sensorSize( iConfig.getUntrackedParameter<int>("SensorSize",128) ),
  _maxChi2( iConfig.getUntrackedParameter<double>("maxChi2",9.48) ),
  _noiseEnergyThreshold( iConfig.getUntrackedParameter<double>("noiseEnergyThreshold",9.0) )
{
  HGCalTBRecHitCollection_ = consumes<HGCalTBRecHitCollection>(iConfig.getParameter<edm::InputTag>("HGCALTBRECHITS"));
  HGCalTBCaloTrackCollection_ = consumes<reco::HGCalTBCaloTrackCollection>(iConfig.getParameter<edm::InputTag>("HGCALTBCALOTRACKS"));
  float sum=0.;
  std::vector<double> vec;
  sum+=0.0 ; vec.push_back( sum );
  sum+=5.35; vec.push_back( sum );
  sum+=5.17; vec.push_back( sum );
  sum+=3.92; vec.push_back( sum );
  sum+=4.08; vec.push_back( sum );
  sum+=1.15; vec.push_back( sum );
  sum+=4.11; vec.push_back( sum );
  sum+=2.14; vec.push_back( sum );
  //cern config 1 (5X0->15X0) is default
  _layerZPositions = iConfig.getUntrackedParameter< std::vector<double> >("LayerZPositions",vec);
  std::cout << iConfig.dump() << std::endl;

  usesResource("TFileService");
  edm::Service<TFileService> fs;
  tree = fs->make<TTree>("tree", "HGCAL TB variables tree");
  tree->Branch( "evtID",&_evtID ); 
  tree->Branch( "chi2",&_chi2 ); 
  tree->Branch( "theta",&_theta );
  tree->Branch( "phi",&_phi );
  tree->Branch( "nhit",&_nhit );
  tree->Branch( "tracknhit",&_tracknhit );
  tree->Branch( "vx0",&_vx0 );
  tree->Branch( "vy0",&_vy0 );
  tree->Branch( "px",&_px );
  tree->Branch( "py",&_py );
  tree->Branch( "trackSuccess",&_trackSuccess );
  tree->Branch( "nhitlayer","std::vector<int>",&_nhitlayer);
  tree->Branch( "energylayer","std::vector<double>",&_energylayer);
  tree->Branch( "deltas","std::vector<double>",&_deltas);

  _evtID=0;

  std::ostringstream os( std::ostringstream::ate );
  for( int ilayer=0; ilayer<_nlayers; ilayer++ ){
    os.str("");
    os << "Delta_Layer" << ilayer;
    TH1F* h = fs->make<TH1F>(os.str().c_str(), os.str().c_str(), 1000, 0, 10 );
    h_delta_layer[ilayer]=h;
    for( int iskiroc=0; iskiroc<_nskirocsperlayer; iskiroc++ ){
      for( int ichannel=0; ichannel< _nchannelsperskiroc; ichannel++ ){
	os.str("");
	os << "Layer" << ilayer << "_Skiroc" << iskiroc << "_Channel" << ichannel << "_MIP";
	int key=ilayer*1000+iskiroc*100+ichannel;
	h = fs->make<TH1F>(os.str().c_str(), os.str().c_str(), 8000, -20, 60 );
	h_mip_map[key]=h;
	os.str("");
	os << "Layer" << ilayer << "_Skiroc" << iskiroc << "_Channel" << ichannel << "_Noise";
	h = fs->make<TH1F>(os.str().c_str(), os.str().c_str(), 10000, -50, 50 );
	h_noise_map[key]=h;
      }
    }
  }

  _nhitlayer = std::vector<int>(_nlayers,0);
  _energylayer = std::vector<double>(_nlayers,0);

}

TrackAnalyzer::~TrackAnalyzer()
{
}

void
TrackAnalyzer::analyze(const edm::Event& event, const edm::EventSetup& setup)
{
  _deltas.clear();
  _tracknhit = _chi2 = _theta = _phi = _nhit = _vx0 = _vy0 = _px = _py = 0;
  _trackSuccess=false;

  for( int i=0; i<_nlayers; i++ ){
    _nhitlayer.at(i)=0;
    _energylayer.at(i)=0;
  }

  edm::Handle<HGCalTBRecHitCollection> Rechits;
  event.getByToken(HGCalTBRecHitCollection_, Rechits);
  edm::Handle<reco::HGCalTBCaloTrackCollection> TracksColl;
  event.getByToken(HGCalTBCaloTrackCollection_, TracksColl);
  
  reco::HGCalTBCaloTrackCollection tracks=(*TracksColl);
  if( tracks.size()>0 ){
    reco::HGCalTBCaloTrack track = (*tracks.begin());
    track.isNull() ? _trackSuccess=false : _trackSuccess=true ;
    _tracknhit = track.getDetIds().size();
    _chi2 = track.normalisedChi2();
    if( track.isNull()==false && _chi2<_maxChi2 ){
      _theta = track.momentum().theta()*180/PI;
      _phi = track.momentum().phi()*180/PI;
      _vx0 = track.vertex().x();
      _vy0 = track.vertex().y();
      _px = track.momentum().x();
      _py = track.momentum().y();

      HGCalTBCellVertices cellVertice;
      DistanceBetweenTrackAndPoint dist;
      for( std::vector<HGCalTBDetId>::iterator it=track.getDetIds().begin(); it!=track.getDetIds().end(); ++it ){
	uint32_t EID = essource_.emap_.detId2eid( *it );
	HGCalTBElectronicsId eid(EID);
	int key=( (*it).layer()-1 )*1000 + (eid.iskiroc() - 1)%2*100 + eid.ichan();
	HGCalTBRecHit hit=(*(*Rechits).find(*it));
	h_mip_map[key]->Fill( hit.energy() );
	_nhitlayer.at( (*it).layer()-1 )+=1;
	_energylayer.at( (*it).layer()-1 )+=hit.energy();	
	std::pair<double,double> xy=cellVertice.GetCellCentreCoordinatesForPlots( (*it).layer(), (*it).sensorIU(), (*it).sensorIV(), (*it).iu(), (*it).iv(), _sensorSize);
	math::XYZPoint xyz(xy.first,xy.second,_layerZPositions.at( (*it).layer()-1 ));
	h_delta_layer[ (*it).layer()-1 ]->Fill( dist.distance(track,xyz) );
      }
    }
  }
  _nhit = (int)(*Rechits).size();
  for( auto hit : *Rechits ){
    if( hit.energy() < _noiseEnergyThreshold ){
      uint32_t EID = essource_.emap_.detId2eid( hit.id() );
      HGCalTBElectronicsId eid(EID);
      int key=( hit.id().layer()-1 )*1000 + (eid.iskiroc() - 1)%2*100 + eid.ichan();
      h_noise_map[key]->Fill( hit.energy() );
    }
  }

  tree->Fill();
  _evtID++;
}

void
TrackAnalyzer::beginJob()
{
  HGCalCondObjectTextIO io(0);
  edm::FileInPath fip(mapfile_);
  if (!io.load(fip.fullPath(), essource_.emap_)) {
    throw cms::Exception("Unable to load electronics map");
  }
}

void
TrackAnalyzer::endJob()
{

}

void
TrackAnalyzer::fillDescriptions(edm::ConfigurationDescriptions& descriptions)
{
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(TrackAnalyzer);
