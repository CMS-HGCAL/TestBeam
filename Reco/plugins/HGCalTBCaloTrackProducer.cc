#include "HGCal/Reco/plugins/HGCalTBCaloTrackProducer.h"
#include "HGCal/Reco/interface/HGCalTBCaloTrackingUtil.h"
#include "HGCal/DataFormats/interface/HGCalTBDetId.h"

HGCalTBCaloTrackProducer::HGCalTBCaloTrackProducer(const edm::ParameterSet& cfg) : 
  _outputCollectionName(cfg.getParameter<std::string>("OutputCollectionName")),
  _doTrackCleaning(cfg.getUntrackedParameter<bool>("doTrackCleaning",false)),
  _maxDistanceToRecoTrack(cfg.getUntrackedParameter<double>("maxDistanceToRecoTrack",1.3)),
  _minTouchedLayers(cfg.getUntrackedParameter<int>("minTouchedLayers",4))
{
  HGCalTBClusterCollection_ = consumes<reco::HGCalTBClusterCollection>(cfg.getParameter<edm::InputTag>("HGCALTBCLUSTERS"));
  std::cout << cfg.dump() << std::endl;

  produces <reco::HGCalTBCaloTrackCollection>(_outputCollectionName);
}

void HGCalTBCaloTrackProducer::produce(edm::Event& event, const edm::EventSetup& iSetup)
{
  std::auto_ptr<reco::HGCalTBCaloTrackCollection> tracks(new reco::HGCalTBCaloTrackCollection);
  edm::Handle<reco::HGCalTBClusterCollection> clusters;
  event.getByToken(HGCalTBClusterCollection_, clusters);

  std::set<int> touchedLayers;
  std::vector<HGCalTBDetId> detIds;
  for( auto cluster : *clusters ){
    touchedLayers.insert( cluster.layer() );
    for( std::vector< std::pair<DetId,float> >::const_iterator it=cluster.hitsAndFractions().begin(); it!=cluster.hitsAndFractions().end(); ++it ){
      detIds.push_back( HGCalTBDetId((*it).first) );
    }
  }
  if( (int)touchedLayers.size()>_minTouchedLayers ){
    reco::WeightedLeastSquare<reco::HGCalTBClusterCollection> wls;
    std::vector<float> trackPar;
    std::vector<float> trackParError;
    wls.run( (*clusters), trackPar, trackParError);
    float chi2 = wls.chi2( (*clusters), trackPar);
    int ndof = (*clusters).size();
    Vector momentum = Vector(-1., 0., trackPar[1]).Cross( Vector(0., -1., trackPar[3]) );
    Point vertex = Point( trackPar[0], trackPar[2], 0.0 );
    reco::HGCalTBCaloTrack track( chi2, ndof, vertex, momentum,/* cov,*/ detIds);
    bool success=true;
    if( _doTrackCleaning ){
      reco::HGCalTBClusterCollection cleancol;
      reco::TrackCleaner cleaner;
      cleaner.clean( (*clusters), cleancol, track, _maxDistanceToRecoTrack );
      touchedLayers.clear();
      detIds.clear();
      for( auto cluster : cleancol ){
	touchedLayers.insert( cluster.layer() );
	for( std::vector< std::pair<DetId,float> >::const_iterator it=cluster.hitsAndFractions().begin(); it!=cluster.hitsAndFractions().end(); ++it )
	  detIds.push_back( HGCalTBDetId((*it).first) );
      }
      if( (int)touchedLayers.size()>_minTouchedLayers ){
	wls.run( cleancol, trackPar, trackParError);
	chi2 = wls.chi2( cleancol, trackPar);
	ndof = cleancol.size();
	momentum = Vector(-1., 0., trackPar[1]).Cross( Vector(0., -1., trackPar[3]) );
	vertex = Point( trackPar[0], trackPar[2], 0.0 );
	track = reco::HGCalTBCaloTrack( chi2, ndof, vertex, momentum,/* cov,*/ detIds);
      }
      else success=false;
    }
    if( success==false )
      track=reco::HGCalTBCaloTrack();
    tracks->push_back(track);
  }
  event.put(tracks,_outputCollectionName);
}


DEFINE_FWK_MODULE(HGCalTBCaloTrackProducer);
