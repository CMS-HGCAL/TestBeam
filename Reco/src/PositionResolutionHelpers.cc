#include "HGCal/Reco/interface/PositionResolutionHelpers.h"



//Particle_Track implementations
Particle_Track::Particle_Track(edm::Handle<HGCalTBRecHitCollection> Rechits, HGCalElectronicsMap& emap) {
}

Particle_Track::~Particle_Track() {
}




//****   Sensor Hit Maps    ****//

//public functions
SensorHitMap::SensorHitMap(){
  centralHitPoint = std::make_pair(0., 0.);
  layerZ = 0;
  sensorSize = 128;
}

void SensorHitMap::setSensorSize(int s){
  sensorSize = s;
}

void SensorHitMap::setZ(double z){
  this->layerZ = z;
}

//reduces the information from the Rechit towards what is necessary for the impact point calculation
void SensorHitMap::addHit(HGCalTBRecHit Rechit) {
  CellCenterXY = TheCell.GetCellCentreCoordinatesForPlots((Rechit.id()).layer(), (Rechit.id()).sensorIU(), (Rechit.id()).sensorIV(), (Rechit.id()).iu(), (Rechit.id()).iv(), sensorSize);
  double iux = CellCenterXY.first;
  double ivy = CellCenterXY.second;
  //Rechit.setCartesianCoordinates(iux, ivy, Layer_Z_Positions[layer]);   //in cm,   is not really necessary

  HitTriple* hit = new HitTriple;
  hit->x = iux;
  hit->y = ivy;
  hit->I = Rechit.energy();
  Hits.push_back(hit);
}

std::pair<double, double> SensorHitMap::calculateCenterPosition(WeightingMethod method) {
  switch(method){
    case SQUAREDWEIGHTING:
      SensorHitMap::squaredWeighting();
    case LINEARWEIGHTING:
    default:
      SensorHitMap::linearWeighting();
  }
  return centralHitPoint;
}

std::pair<double, double> SensorHitMap::getCenterPosition() {
  return centralHitPoint;
}

//private functions
void SensorHitMap::squaredWeighting() {
  //todo
  centralHitPoint.first = 1.0;
  centralHitPoint.second = 1.0;
  return;
}

void SensorHitMap::linearWeighting() {
  //todo
  centralHitPoint.first = 2.0;
  centralHitPoint.second = 2.0;
  return;
}

//debug function
void SensorHitMap::printHits() {
  
  for(std::vector<HitTriple*>::iterator it=Hits.begin(); it!=Hits.end(); it++){
    std::cout<<(*it)->x<<"  "<<(*it)->y<<"  "<<(*it)->I<<std::endl;
  }
  
}
