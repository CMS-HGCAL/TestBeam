#include "HGCal/Reco/interface/PositionResolutionHelpers.h"



//Particle_Track implementations
Particle_Track::Particle_Track(edm::Handle<HGCalTBRecHitCollection> Rechits){;//, HGCalElectronicsMap& emap) {
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
  int ID = (Rechit.id()).cellType();
  double energy = Rechit.energy(); //input to centrum calculation ?
  double energyHigh = Rechit.energyHigh(); //used for pedestal subtraction analogous to the event display

  //Rechit.setCartesianCoordinates(iux, ivy, Layer_Z_Positions[layer]);   //in cm,   is not really necessary

  HitTriple* hit = new HitTriple;
  hit->ID = ID;
  hit->x = iux;
  hit->y = ivy;
  hit->I = energy;
  Hits.push_back(hit);

  if (energyHigh <= 30.0) {
    if (cellTypeCount.find(ID) == cellTypeCount.end()) {
      cellTypeCount[ID] = 0;
      pedestalCount[ID] = 0;
    }
    cellTypeCount[ID] += 1;
    pedestalCount[ID] += energyHigh;
  }
}

void SensorHitMap::subtractPedestals() {
  for(std::vector<HitTriple*>::iterator hit=Hits.begin(); hit!=Hits.end(); hit++){
    //analogous treatment of pedestals as in the RecHitPlotter_HighGain_New plugin
    switch((*hit)->ID) {
      case 0:
      case 4:
        (*hit)->I -= (pedestalCount[0]+pedestalCount[4])/(cellTypeCount[0]+cellTypeCount[4]);
      case 2:
        (*hit)->I -= pedestalCount[2]/cellTypeCount[2];
      case 3:
        (*hit)->I -= pedestalCount[3]/cellTypeCount[3];
      case 5:
        (*hit)->I -= pedestalCount[5]/cellTypeCount[5];
      case 1:
      default:
        continue;
    }
  }
}


std::pair<double, double> SensorHitMap::calculateCenterPosition(WeightingMethod method) {
  switch(method){
    case SQUAREDWEIGHTING:
      SensorHitMap::poweredWeighting(2);
    case LINEARWEIGHTING:
    default:
      SensorHitMap::poweredWeighting(1);
  }
  return centralHitPoint;
}

std::pair<double, double> SensorHitMap::getCenterPosition() {
  return centralHitPoint;
}

//private functions
void SensorHitMap::poweredWeighting(int exponent) {
  double numerator_x, numerator_y, denominator;
  numerator_x = numerator_y = denominator = 0; 
  double w;
  for(std::vector<HitTriple*>::iterator hit=Hits.begin(); hit!=Hits.end(); hit++){
    w = (*hit)->I >= 0. ? pow((*hit)->I, exponent) : 0.0;
    denominator += w;
    numerator_x += w*(*hit)->x;
    numerator_y += w*(*hit)->y;
  }
  centralHitPoint.first = numerator_x/denominator;
  centralHitPoint.second = numerator_y/denominator;
}


//debug function
void SensorHitMap::printHits() {
  
  for(std::vector<HitTriple*>::iterator it=Hits.begin(); it!=Hits.end(); it++){
    std::cout<<(*it)->x<<"  "<<(*it)->y<<"  "<<(*it)->I<<std::endl;
  }
  
}
