#include "HGCal/Reco/plugins/HGCalTBClusterProducer.h"
#include "HGCal/Geometry/interface/HGCalTBTopology.h"
#include "HGCal/Geometry/interface/HGCalTBCellVertices.h"
#include "HGCal/CondObjects/interface/HGCalCondObjectTextIO.h"
#include <map>
#include <algorithm>
#include <sstream>

HGCalTBClusterProducer::HGCalTBClusterProducer(const edm::ParameterSet& cfg) : 
  _elecMapFile(cfg.getUntrackedParameter<std::string>("ElectronicMapFile",std::string("HGCal/CondObjects/data/map_CERN_8Layers_Sept2016.txt"))),
  _outputCollectionName(cfg.getParameter<std::string>("OutputCollectionName")),
  _rechitToken(consumes<HGCalTBRecHitCollection>(cfg.getParameter<edm::InputTag>("rechitCollection"))),
  _runDynamicCluster(cfg.getUntrackedParameter<bool>("runDynamicCluster",true)),
  _runCluster7(cfg.getUntrackedParameter<bool>("runCluster7",true)),
  _runCluster19(cfg.getUntrackedParameter<bool>("runCluster19",true)),
  _sensorSize(cfg.getUntrackedParameter<int>("sensorSize",128)),
  _minEnergy(cfg.getUntrackedParameter<double>("minEnergy",0.0)),
  _rmSpecialCells(cfg.getUntrackedParameter<bool>("RemoveSpecialCells",false)),
  _positionWeights(cfg.getUntrackedParameter<std::string>("PositionWeightsOption",std::string("linear")))
{
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
  _layerZPositions = cfg.getUntrackedParameter< std::vector<double> >("LayerZPositions",vec);

  vec.clear(); vec.push_back( 5 ); vec.push_back( 1 );
  _logParams = cfg.getUntrackedParameter< std::vector<double> >("LogWeightParams",vec);
  if( _logParams.size()<2 )
    throw cms::Exception("Wrong vector size for _logParams in HGCalTBCluster");

  std::cout << cfg.dump() << std::endl;
  produces <reco::HGCalTBClusterCollection>(_outputCollectionName);
  if( _runCluster7 ){
    std::ostringstream os( std::ostringstream::ate );
    os.str("");os << _outputCollectionName << 7;
    _outputCollectionName7=os.str();
    produces <reco::HGCalTBClusterCollection>(_outputCollectionName7);
  }
  if( _runCluster19 ){
    std::ostringstream os( std::ostringstream::ate );
    os.str("");os << _outputCollectionName << 19;
    _outputCollectionName19=os.str();
    produces <reco::HGCalTBClusterCollection>(_outputCollectionName19);
  }
  
  HGCalCondObjectTextIO io(0);
  edm::FileInPath fip(_elecMapFile);
  if (!io.load(fip.fullPath(), _elecMap)) {
    throw cms::Exception("Unable to load electronics map");
  }
  
}

void HGCalTBClusterProducer::produce(edm::Event& event, const edm::EventSetup& iSetup)
{

  std::auto_ptr<reco::HGCalTBClusterCollection> clusters(new reco::HGCalTBClusterCollection);
  std::auto_ptr<reco::HGCalTBClusterCollection> clusters7(new reco::HGCalTBClusterCollection);
  std::auto_ptr<reco::HGCalTBClusterCollection> clusters19(new reco::HGCalTBClusterCollection);

  edm::Handle<HGCalTBRecHitCollection> rechits;
  event.getByToken(_rechitToken, rechits);
  std::map<int, HGCalTBRecHitCollection> hitmap;

  for(auto hit : *rechits ){
    if( hit.energy()<_minEnergy )
       continue;
    if( _rmSpecialCells && hit.id().cellType()!=0 && hit.id().cellType()!=4 )
      continue;
    if( hitmap.find( hit.id().layer() )!=hitmap.end() )
      hitmap[ hit.id().layer() ].push_back(hit);
    else{
      HGCalTBRecHitCollection hitcol;
      hitcol.push_back(hit);
      std::pair< int, HGCalTBRecHitCollection > p( hit.id().layer(),hitcol );
      hitmap.insert(p);
    }
  }
  for( std::map<int,HGCalTBRecHitCollection>::iterator it=hitmap.begin(); it!=hitmap.end(); ++it ){
    if( _runDynamicCluster ){
      std::vector<reco::HGCalTBCluster> vec;
      createDynamicClusters(it->second, vec);
      for( std::vector<reco::HGCalTBCluster>::iterator jt=vec.begin(); jt!=vec.end(); ++jt )
	clusters->push_back(*jt);
    }
    if( _runCluster7 ){
      reco::HGCalTBCluster cluster;
      createSeededClusters(it->second,cluster,1);
      clusters7->push_back(cluster);
    }
    if( _runCluster19 ){
      reco::HGCalTBCluster cluster;
      createSeededClusters(it->second,cluster,2);
      clusters19->push_back(cluster);
    }
  }
  if( _runDynamicCluster )
    event.put(clusters, _outputCollectionName);
  if( _runCluster7 )
    event.put(clusters7, _outputCollectionName7);
  if( _runCluster19 )
    event.put(clusters19, _outputCollectionName19);
}

void HGCalTBClusterProducer::createDynamicClusters(HGCalTBRecHitCollection rechits, std::vector<reco::HGCalTBCluster> &clusterCol)
{
  _maxTransverse=1;
  std::vector<HGCalTBDetId> temp;
  for( std::vector<HGCalTBRecHit>::iterator it=rechits.begin(); it!=rechits.end(); ++it){
    if( std::find( temp.begin(),temp.end(),(*it).id() )!=temp.end() ) continue;
    temp.push_back( (*it).id() );
    std::vector<HGCalTBDetId> clusterDetIDs;
    clusterDetIDs.push_back( (*it).id() );
    buildCluster(rechits, temp, clusterDetIDs);            
    reco::HGCalTBCluster cluster;
    cluster.setLayer( (*it).id().layer() );
    float energyHigh=0.;
    float energyLow=0.;
    float energy=0.;
    HGCalTBRecHit seed=(*rechits.find(*clusterDetIDs.begin()));
    for( std::vector<HGCalTBDetId>::iterator jt=clusterDetIDs.begin(); jt!=clusterDetIDs.end(); ++jt){
      HGCalTBRecHit hit=(*rechits.find(*jt));
      energyHigh+=hit.energyHigh();
      energyLow+=hit.energyLow();
      energy+=hit.energy();
      if( hit.energy()<seed.energy() ) seed=hit;
    } 
    cluster.setSeed(seed.id());
    cluster.setEnergyLow(energyLow);
    cluster.setEnergyHigh(energyHigh);
    cluster.setEnergy(energy);
    for( std::vector<HGCalTBDetId>::iterator jt=clusterDetIDs.begin(); jt!=clusterDetIDs.end(); ++jt)
      cluster.addHitAndFraction( (*jt), (*rechits.find(*jt)).energy()/energy );
    
    cluster.setPosition( clusterPosition(cluster) );
    clusterCol.push_back(cluster);
  }
}

void HGCalTBClusterProducer::buildCluster(  HGCalTBRecHitCollection rechits,
					    std::vector<HGCalTBDetId> &temp,
					    std::vector<HGCalTBDetId> &clusterDetIDs
					    )
{
  HGCalTBTopology top;
  HGCalTBDetId detID=clusterDetIDs.back();
  std::set<HGCalTBDetId> neighbors=top.getNeighboringCellsDetID( detID, _sensorSize , _maxTransverse, _elecMap );
  for( std::set<HGCalTBDetId>::const_iterator it=neighbors.begin(); it!=neighbors.end(); ++it){
    if( std::find(temp.begin(), temp.end(), (*it))!=temp.end() || rechits.find(*it)==rechits.end() )
      continue;
    temp.push_back( (*it) );
    clusterDetIDs.push_back( (*it) );
    buildCluster(rechits, temp, clusterDetIDs);
  }
}

void HGCalTBClusterProducer::createSeededClusters(HGCalTBRecHitCollection rechits, reco::HGCalTBCluster &cluster, int maxDist)
{
  if( rechits.size()==0 ) return;
  _maxTransverse=maxDist;
  HGCalTBTopology top;
  HGCalTBRecHit seed=(*rechits.begin());
  for( std::vector<HGCalTBRecHit>::iterator it=rechits.begin(); it!=rechits.end(); ++it){
    if( (*it).energy() > seed.energy() )
      seed=(*it);
  }
  cluster.setLayer( seed.id().layer() );
  cluster.setSeed( seed.id() );
  float energyHigh=seed.energyHigh();
  float energyLow=seed.energyLow();
  float energy=seed.energy();
  
  std::set<HGCalTBDetId> neighbors=top.getNeighboringCellsDetID( seed.id(), _sensorSize , _maxTransverse, _elecMap );
  for( std::set<HGCalTBDetId>::iterator jt=neighbors.begin(); jt!=neighbors.end(); ++jt){
    if( rechits.find(*jt) != rechits.end() ){
      HGCalTBRecHit hit=(*rechits.find(*jt));
      energyHigh+=hit.energyHigh();
      energyLow+=hit.energyLow();
      energy+=hit.energy();
    }
  }
  cluster.setEnergyLow(energyLow);
  cluster.setEnergyHigh(energyHigh);
  cluster.setEnergy(energy);

  cluster.addHitAndFraction( seed.id(), seed.energy()/energy );
  for( std::set<HGCalTBDetId>::iterator jt=neighbors.begin(); jt!=neighbors.end(); ++jt)
    if( rechits.find(*jt) != rechits.end() )
      cluster.addHitAndFraction( (*jt), (*rechits.find(*jt)).energy()/energy );
  
  cluster.setPosition( clusterPosition(cluster) );
}

math::XYZPoint HGCalTBClusterProducer::clusterPosition(reco::HGCalTBCluster cluster)
{
  HGCalTBCellVertices cellVertice;
  std::pair<double, double> xy;
  double sumweight=0.;
  double x=0.;
  double y=0.;
  for( std::vector< std::pair<DetId,float> >::const_iterator it=cluster.hitsAndFractions().begin(); it!=cluster.hitsAndFractions().end(); ++it){
    HGCalTBDetId detId=(*it).first;
    float weight;
    if( _positionWeights==std::string("logarithmic") )
      weight=_logParams[0]+_logParams[1]*std::log( (*it).second );
    else 
      weight=(*it).second;
    if( weight<0 ) continue;
    xy=cellVertice.GetCellCentreCoordinatesForPlots( detId.layer(), detId.sensorIU(), detId.sensorIV(), detId.iu(), detId.iv(), _sensorSize);
    x+=xy.first*weight;
    y+=xy.second*weight;
    sumweight+=weight;
  }
  x/=sumweight;
  y/=sumweight;
  return math::XYZPoint( x,y,_layerZPositions.at( cluster.layer()-1 ) );
}

DEFINE_FWK_MODULE(HGCalTBClusterProducer);
