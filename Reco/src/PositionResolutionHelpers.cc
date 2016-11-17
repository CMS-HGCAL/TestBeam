#include "HGCal/Reco/interface/PositionResolutionHelpers.h"


//****   Sensor Hit Maps    ****//

//public functions
SensorHitMap::SensorHitMap(){
  centralHitPoint = std::make_pair(0., 0.);
  layerZ = 0;
  sensorSize = 128;
}

void SensorHitMap::setSensorSize(int s) {
  sensorSize = s;
}

void SensorHitMap::setZ(double z) {
  this->layerZ = z;
}

double SensorHitMap::getZ() {
  return this->layerZ;
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

  //analogous to RecHitPlotter_HighGain_New, only add to pedestals if energyHigh exceeds 30
  if (energyHigh <= 30.0) { //TODO: make a threshold
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


void SensorHitMap::calculateCenterPosition(WeightingMethod method) {
  switch(method){
    case SQUAREDWEIGHTING:
      SensorHitMap::poweredWeighting(2);
      break;
    case LINEARWEIGHTING:
      SensorHitMap::poweredWeighting(1);
      break;
    //case NEWMETHOD:
      //SensorHitMap::newWeightingFunction()
      //break;
    default:
      SensorHitMap::poweredWeighting(2);
  }
}

std::pair<double, double> SensorHitMap::getCenterPosition() {
  return centralHitPoint;
}
std::pair<double, double> SensorHitMap::getCenterPositionError() {
  return centralHitPointError;
}


//private functions
void SensorHitMap::poweredWeighting(int exponent) {
  double threshold = 0.0;     //TODO: maybe as parameter

  double numerator_x, numerator_y, denominator;
  numerator_x = numerator_y = denominator = 0; 
  double w;
  for(std::vector<HitTriple*>::iterator hit=Hits.begin(); hit!=Hits.end(); hit++){
    w = (*hit)->I >= threshold ? pow((*hit)->I, exponent) : 0.0;    //0.0 --> not included in the sum
    denominator += w;
    numerator_x += w*(*hit)->x;
    numerator_y += w*(*hit)->y;
  }
  centralHitPoint.first = numerator_x/denominator;
  centralHitPoint.second = numerator_y/denominator;

  //calculate the RMs
  numerator_x = numerator_y = 0.0;
  for(std::vector<HitTriple*>::iterator hit=Hits.begin(); hit!=Hits.end(); hit++){ 
    w = (*hit)->I >= threshold ? pow((*hit)->I, exponent) : 0.0; 
    numerator_x += w*pow((*hit)->x - centralHitPoint.first, 2);
    numerator_y += w*pow((*hit)->y - centralHitPoint.second, 2);
  }
  centralHitPointError.first = sqrt(numerator_x/denominator);
  centralHitPointError.second = sqrt(numerator_y/denominator);
}




//****   Particle Tracks    ****//

//public functions
ParticleTrack::ParticleTrack(){
  lastAppliedMethod = DEFAULTFITTING;
  ROOTpol_x = ROOTpol_y = 0;

};
ParticleTrack::~ParticleTrack(){

};

void ParticleTrack::addFitPoint(SensorHitMap* sensor){
  x.push_back(sensor->getCenterPosition().first);  
  x_err.push_back(sensor->getCenterPositionError().first);  
  y.push_back(sensor->getCenterPosition().second);  
  y_err.push_back(sensor->getCenterPositionError().second);  
  z.push_back(sensor->getZ());  
  z_err.push_back(0.0);
};


void ParticleTrack::fitTrack(TrackFittingMethod method){
  try {
    switch(method) {
      case LINEFITTGRAPHERRORS:
        lineFitTGraphErrors();
        break;
      default:
        lastAppliedMethod = DEFAULTFITTING;
        break;
    }
    lastAppliedMethod = method;      
  } catch(cms::Exception& iException) {
    std::cout<<"Fitting method has failed with error code "<< iException <<". Attempting the default fitting."<<std::endl<<std::endl;
    fitTrack(DEFAULTFITTING);
  }
}

std::pair<double, double> ParticleTrack::calculatePositionXY(double z) {
  switch(lastAppliedMethod) {
    case LINEFITTGRAPHERRORS:
      return positionFromLineFitTGraphErrors(z);
    default:
      return std::make_pair(0.,0.); 
  }
}

//private functions
void ParticleTrack::lineFitTGraphErrors(){  
  if (ROOTpol_x == 0) {
    ROOTpol_x = new TF1("ROOTpol_x", "pol1", *min_element(z.begin(), z.end())-1.0, *max_element(z.begin(), z.end())+1.0);
  }
  if (ROOTpol_y == 0) {
    ROOTpol_y = new TF1("ROOTpol_y", "pol1", *min_element(z.begin(), z.end())-1.0, *max_element(z.begin(), z.end())+1.0);
  }
  TGraph* tmp_graph_x = new TGraphErrors(z.size(), &(z[0]), &(x[0]), &(z_err[0]), &(x_err[0])); //z_err should be filled with zeros
  TGraph* tmp_graph_y = new TGraphErrors(z.size(), &(z[0]), &(y[0]), &(z_err[0]), &(y_err[0]));
  
  tmp_graph_x->Fit(ROOTpol_x, "Q");
  tmp_graph_y->Fit(ROOTpol_y, "Q");

}; 

std::pair<double, double> ParticleTrack::positionFromLineFitTGraphErrors(double z) {
  return std::make_pair(ROOTpol_x->Eval(z), ROOTpol_y->Eval(z));
}



//debug function
void SensorHitMap::printHits() {
  
  for(std::vector<HitTriple*>::iterator it=Hits.begin(); it!=Hits.end(); it++){
    std::cout<<(*it)->x<<"  "<<(*it)->y<<"  "<<(*it)->I<<std::endl;
  }
  
}
