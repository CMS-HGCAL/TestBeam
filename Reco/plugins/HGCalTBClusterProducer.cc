#include "HGCal/Reco/plugins/HGCalTBClusterProducer.h"
#include "HGCal/Geometry/interface/HGCalTBTopology.h"
#include "HGCal/Geometry/interface/HGCalTBCellVertices.h"
#include "HGCal/CondObjects/interface/HGCalCondObjectTextIO.h"
#include <map>
#include <algorithm>
#include <sstream>

HGCalTBClusterProducer::HGCalTBClusterProducer(const edm::ParameterSet& cfg) : 
  m_electronicMap(cfg.getUntrackedParameter<std::string>("ElectronicMapFile",std::string("HGCal/CondObjects/data/map_CERN_Hexaboard_October_20Sensors_5EELayers_7FHLayers_V0.txt"))),
  m_detectorLayoutFile(cfg.getUntrackedParameter<std::string>("DetectorLayout","HGCal/CondObjects/data/layerGeom_oct2017_h2_17layers.txt")),
  m_sensorSize(cfg.getUntrackedParameter<int>("SensorSize",128)),
  m_outputCollectionName(cfg.getParameter<std::string>("OutputCollectionName")),
  m_runDynamicCluster(cfg.getUntrackedParameter<bool>("runDynamicCluster",true)),
  m_runCluster7(cfg.getUntrackedParameter<bool>("runCluster7",true)),
  m_runCluster19(cfg.getUntrackedParameter<bool>("runCluster19",true)),
  m_minEnergy(cfg.getUntrackedParameter<double>("minEnergy",0.0))
{
  m_HGCalTBRecHitCollection = consumes<HGCalTBRecHitCollection>(cfg.getParameter<edm::InputTag>("InputCollection"));

  produces <reco::HGCalTBClusterCollection>(m_outputCollectionName);
  if( m_runCluster7 ){
    std::ostringstream os( std::ostringstream::ate );
    os.str("");os << m_outputCollectionName << 7;
    m_outputCollectionName7=os.str();
    produces <reco::HGCalTBClusterCollection>(m_outputCollectionName7);
  }
  if( m_runCluster19 ){
    std::ostringstream os( std::ostringstream::ate );
    os.str("");os << m_outputCollectionName << 19;
    m_outputCollectionName19=os.str();
    produces <reco::HGCalTBClusterCollection>(m_outputCollectionName19);
  }
  
  std::cout << cfg.dump() << std::endl;
}

void HGCalTBClusterProducer::beginJob()
{
  HGCalCondObjectTextIO io(0);
  edm::FileInPath fip(m_electronicMap);
  if (!io.load(fip.fullPath(), m_essource.emap)) {
    throw cms::Exception("Unable to load electronics map");
  }
  fip=edm::FileInPath(m_detectorLayoutFile);
  if (!io.load(fip.fullPath(), m_essource.layout)) {
    throw cms::Exception("Unable to load detector layout file");
  };
  for( auto layer : m_essource.layout.layers() )
    layer.print();
}

void HGCalTBClusterProducer::produce(edm::Event& event, const edm::EventSetup& iSetup)
{

  std::auto_ptr<reco::HGCalTBClusterCollection> clusters(new reco::HGCalTBClusterCollection);
  std::auto_ptr<reco::HGCalTBClusterCollection> clusters7(new reco::HGCalTBClusterCollection);
  std::auto_ptr<reco::HGCalTBClusterCollection> clusters19(new reco::HGCalTBClusterCollection);

  edm::Handle<HGCalTBRecHitCollection> rechits;
  event.getByToken(m_HGCalTBRecHitCollection, rechits);
  std::map<int, HGCalTBRecHitCollection> hitmap;

  for(auto hit : *rechits ){
    if( hit.energy()<m_minEnergy )
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
    if( m_runDynamicCluster ){
      std::vector<reco::HGCalTBCluster> vec;
      createDynamicClusters(it->second, vec);
      for( std::vector<reco::HGCalTBCluster>::iterator jt=vec.begin(); jt!=vec.end(); ++jt ){
	clusters->push_back(*jt);
      }
    }
    if( m_runCluster7 ){
      reco::HGCalTBCluster cluster;
      createSeededClusters(it->second,cluster,1);
      clusters7->push_back(cluster);
    }
    if( m_runCluster19 ){
      reco::HGCalTBCluster cluster;
      createSeededClusters(it->second,cluster,2);
      clusters19->push_back(cluster);
    }
  }
  if( m_runDynamicCluster )
    event.put(clusters, m_outputCollectionName);
  if( m_runCluster7 )
    event.put(clusters7, m_outputCollectionName7);
  if( m_runCluster19 )
    event.put(clusters19, m_outputCollectionName19);
}

void HGCalTBClusterProducer::createDynamicClusters(HGCalTBRecHitCollection rechits, std::vector<reco::HGCalTBCluster> &clusterCol)
{
  int maxDistance=1;
  HGCalTBCellVertices cellVertice;
  std::pair<double, double> CellCentreXY;
  std::vector<HGCalTBDetId> temp;
  for( std::vector<HGCalTBRecHit>::iterator it=rechits.begin(); it!=rechits.end(); ++it){
    if( std::find( temp.begin(),temp.end(),(*it).id() )!=temp.end() ) continue;
    temp.push_back( (*it).id() );
    std::vector<HGCalTBDetId> clusterDetIDs;
    clusterDetIDs.push_back( (*it).id() );
    buildCluster(rechits, temp, clusterDetIDs, maxDistance);
    reco::HGCalTBCluster cluster;
    cluster.setLayer( (*it).id().layer() );
    float energyHigh=0.;
    float energyLow=0.;
    float energy=0.;
    float x,y,z;
    x = y = z = 0.0;
    for( std::vector<HGCalTBDetId>::iterator jt=clusterDetIDs.begin(); jt!=clusterDetIDs.end(); ++jt){
      HGCalTBRecHit hit=(*rechits.find(*jt));
      energyHigh+=hit.energyHigh();
      energyLow+=hit.energyLow();
      energy+=hit.energy();
      CellCentreXY=cellVertice.GetCellCentreCoordinatesForPlots( (*jt).layer(), (*jt).sensorIU(), (*jt).sensorIV(), (*jt).iu(), (*jt).iv(), m_sensorSize);
      x += CellCentreXY.first*hit.energy();
      y += CellCentreXY.second*hit.energy();
      z += m_essource.layout.at( hit.id().layer()-1 ).z()*hit.energy();
    } 
    cluster.setPosition( math::XYZPoint(x,y,z)/energy );
    cluster.setEnergyLow(energyLow);
    cluster.setEnergyHigh(energyHigh);
    cluster.setEnergy(energy);
    for( std::vector<HGCalTBDetId>::iterator jt=clusterDetIDs.begin(); jt!=clusterDetIDs.end(); ++jt)
      cluster.addHitAndFraction( (*jt), (*rechits.find(*jt)).energy()/energy );
  
    clusterCol.push_back(cluster);
  }
}

void HGCalTBClusterProducer::buildCluster(  HGCalTBRecHitCollection rechits,
					    std::vector<HGCalTBDetId> &temp,
					    std::vector<HGCalTBDetId> &clusterDetIDs,
					    int maxDistance
					    )
{
  HGCalTBTopology top;
  HGCalTBDetId detID=clusterDetIDs.back();
  std::set<HGCalTBDetId> neighbors=top.getNeighboringCellsDetID( detID, m_sensorSize , maxDistance, m_essource.emap );
  for( std::set<HGCalTBDetId>::const_iterator it=neighbors.begin(); it!=neighbors.end(); ++it){
    if( std::find(temp.begin(), temp.end(), (*it))!=temp.end() || rechits.find(*it)==rechits.end() )
      continue;
    temp.push_back( (*it) );
    clusterDetIDs.push_back( (*it) );
    buildCluster(rechits, temp, clusterDetIDs,maxDistance);
  }
}

void HGCalTBClusterProducer::createSeededClusters(HGCalTBRecHitCollection rechits, reco::HGCalTBCluster &cluster, int maxDistance)
{
  if( rechits.size()==0 ) return;
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
  HGCalTBCellVertices cellVertice;
  std::pair<double, double> CellCentreXY;
  CellCentreXY=cellVertice.GetCellCentreCoordinatesForPlots( seed.id().layer(), seed.id().sensorIU(), seed.id().sensorIV(), seed.id().iu(), seed.id().iv(), m_sensorSize);
  float x = CellCentreXY.first*seed.energy();
  float y = CellCentreXY.second*seed.energy();
  float z = m_essource.layout.at( seed.id().layer()-1 ).z();
  
  std::set<HGCalTBDetId> neighbors=top.getNeighboringCellsDetID( seed.id(), m_sensorSize , maxDistance, m_essource.emap );
  for( std::set<HGCalTBDetId>::iterator jt=neighbors.begin(); jt!=neighbors.end(); ++jt){
    if( rechits.find(*jt) != rechits.end() ){
      HGCalTBRecHit hit=(*rechits.find(*jt));
      energyHigh+=hit.energyHigh();
      energyLow+=hit.energyLow();
      energy+=hit.energy();
      CellCentreXY=cellVertice.GetCellCentreCoordinatesForPlots( (*jt).layer(), (*jt).sensorIU(), (*jt).sensorIV(), (*jt).iu(), (*jt).iv(), m_sensorSize);
      x += CellCentreXY.first*hit.energy();
      y += CellCentreXY.second*hit.energy();
    }
  }
  cluster.setPosition( math::XYZPoint(x/energy,y/energy,z) );
  cluster.setEnergyLow(energyLow);
  cluster.setEnergyHigh(energyHigh);
  cluster.setEnergy(energy);

  cluster.addHitAndFraction( seed.id(), seed.energy()/energy );
  for( std::set<HGCalTBDetId>::iterator jt=neighbors.begin(); jt!=neighbors.end(); ++jt)
    if( rechits.find(*jt) != rechits.end() )
      cluster.addHitAndFraction( (*jt), (*rechits.find(*jt)).energy()/energy );
  
}

DEFINE_FWK_MODULE(HGCalTBClusterProducer);
