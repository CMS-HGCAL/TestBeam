#include "HGCal/Reco/interface/Sensors.h"


//****   Sensor Hit Maps    ****//

//public functions
SensorHitMap::SensorHitMap(int l){
  _label = l;

  mostSignificantHit = NULL;  //points to the most significant hit

  centralHitPoint = std::make_pair(0., 0.);
  centralHitZ = 0;
  centralHitPointError = std::make_pair(sqrt(12.), sqrt(12.));
  CM_threshold = 2.;
  layerLabZ = 0;
  layerX0 = 0;
  particleEnergy = 0;
  ADC_per_MIP = 1.;
  sensorSize = 128;

  d_alpha = 0., d_beta = 0., d_gamma = 0., d_x0 = 0., d_y0 = 0., d_z0 = 0.;
  residualResolution = -1;

  CM_cells_count = 0;
  CM_sum = 0.;

  totalWeight = 1.; 
  totalEnergy = 0.;
}

SensorHitMap::~SensorHitMap(){
  for(std::map<int, HitData*>::iterator hit=Hits.begin(); hit!=Hits.end(); hit++){
    delete (*hit).second;
  }
  Hits.clear();
}

void SensorHitMap::setSensorSize(int s) {
  sensorSize = s;
}

void SensorHitMap::setCenterHitPosition(double x, double y, double x_err, double y_err) {
  centralHitPoint.first = x;
  centralHitPoint.second = y;
  centralHitPointError.first = x_err;
  centralHitPointError.second = y_err;

  centralHitZ = -d_gamma*centralHitPoint.first - d_beta*centralHitPoint.second;
}

void SensorHitMap::setLabZ(double z_cm, double X0) {
  this->layerLabZ = z_cm;
  this->layerX0 = X0;
}

void SensorHitMap::setParticleEnergy(double e) {
  this->particleEnergy = e;
}

void SensorHitMap::setAlignmentParameters(double d_alpha, double d_beta, double d_gamma, double d_x0, double d_y0, double d_z0) {
  this->d_alpha = d_alpha;
  this->d_beta = d_beta;
  this->d_gamma = d_gamma;
  this->d_x0 = d_x0;
  this->d_y0 = d_y0;
  this->d_z0 = d_z0;
}

void SensorHitMap::setResidualResolution(double r) {
  residualResolution = r;
}

void SensorHitMap::setADCPerMIP(double ADC_per_MIP) {
  this->ADC_per_MIP = ADC_per_MIP;
}

void SensorHitMap::setPedestalThreshold(double t) {
  this->CM_threshold = t; 
}

double SensorHitMap::getLabZ() {
  return this->layerLabZ + this->d_z0;
}

double SensorHitMap::getParticleEnergy() {
  return this->particleEnergy;
}

double SensorHitMap::getIntrinsicHitZPosition() { //this value is always zero if rotational angles are all zero
  return this->centralHitZ;
}

double SensorHitMap::getX0() {
  return this->layerX0;
}
//reduces the information from the Rechit towards what is necessary for the impact point calculation
//here all hits independent from the cellType are included
void SensorHitMap::addHit(HGCalTBRecHit Rechit) {
  int uniqueID = (Rechit.id()).rawId();

  //CellCenterXY = TheCell.GetCellCentreCoordinatesForPlots((Rechit.id()).layer(), (Rechit.id()).sensorIU(), (Rechit.id()).sensorIV(), (Rechit.id()).iu(), (Rechit.id()).iv(), sensorSize);
  double iux = Rechit.getCellCenterCartesianCoordinate(0);  //CellCenterXY.first;
  double ivy = Rechit.getCellCenterCartesianCoordinate(1); //CellCenterXY.second;
  int ID = (Rechit.id()).cellType();

  if (filterByCellType(ID)) return; //returns false so far

  double energy = Rechit.energy() / ADC_per_MIP;  //the LayerSumAnalyzer also energy deposits in MIP units

  Hits[uniqueID] = new HitData;
  Hits[uniqueID]->ID = ID;
  Hits[uniqueID]->x = iux;
  Hits[uniqueID]->y = ivy;
  Hits[uniqueID]->I = energy;
  Hits[uniqueID]->E = energy;
  totalEnergy += energy;
  
  if (mostSignificantHit==NULL || energy > mostSignificantHit->E) {
    mostSignificantHit = Hits[uniqueID];
  }

  //analogous to RecHitPlotter_HighGain_New, only add to pedestals if energyHigh exceeds a threshold (default is 2. if not set in the setPedestalThreshold)
  if (energy <= CM_threshold) { //also analogous to the implementation in the LayerSumAnalyzer
    CM_cells_count++;
    CM_sum += energy;
  }
}

void SensorHitMap::registerClusterHit(HGCalTBDetId hit, int N_considered) {  //requires that all hits have been added to the layer
  int uniqueID = hit.rawId();

  if(Hits.find(uniqueID) == Hits.end()) return; //if cluster hit does not represent a valid cell type, do not add it

  if (totalClusterEnergy.find(N_considered) == totalClusterEnergy.end()) {
    totalClusterEnergy[N_considered] = 0;
  }
  totalClusterEnergy[N_considered] += Hits[uniqueID]->E;

  clusterIndexes[N_considered].push_back(uniqueID);
}

std::pair<int, double> SensorHitMap::subtractCM() {
  totalEnergy = 0.;
  double cm_subtraction = CM_cells_count > 0. ? CM_sum/CM_cells_count : 0.;

  for(std::map<int, HitData*>::iterator hit=Hits.begin(); hit!=Hits.end(); hit++){
    if (filterByCellType((*hit).second->ID)) continue;
    // we want: - all cells that were input to the cm (common mode) calculation get weight 0
    //          - all the others are corrected by the cm
    if ((*hit).second->E > CM_threshold) {    
      (*hit).second->I = std::max((*hit).second->I - cm_subtraction, 0.);
    } else {
      (*hit).second->I = 0.;//std::max((*hit).second->I - cm_subtraction, 0.);  //everything below the given threshold is set to zero
    }
    totalEnergy += (*hit).second->I;
  }

  return std::make_pair(CM_cells_count, CM_sum);
}

void SensorHitMap::calculateCenterPosition(ConsiderationMethod considerationMethod, WeightingMethod weightingMethod) {
  switch(considerationMethod){
    case CONSIDERALL:
      SensorHitMap::considerNClosest(-1);
      break;
    case CONSIDERSEVEN:
      //analogous to the LayerSumAnalyzer (25.11.16)
      SensorHitMap::considerNClosest(7); 
      break;
    case CONSIDERNINETEEN:
      //analogous to the LayerSumAnalyzer (25.11.16)
      SensorHitMap::considerNClosest(19); 
      break;
    case CONSIDERCLUSTERSALL:
      SensorHitMap::considerClusters(-1);
      break;
    case CONSIDERCLUSTERSSEVEN:
      SensorHitMap::considerClusters(7);
      break;
    case CONSIDERCLUSTERSNINETEEN:
      SensorHitMap::considerClusters(19);
      break;
    default:
      SensorHitMap::considerNClosest(-1.);
      break;
  }

  switch(weightingMethod){
    case SQUAREDWEIGHTING:
      SensorHitMap::poweredWeighting(2);
      break;
    case LINEARWEIGHTING:
      SensorHitMap::poweredWeighting(1);
      break;
    case LOGWEIGHTING_30_10:
      SensorHitMap::logWeighting(3.0, 1.0);
      break;
    case LOGWEIGHTING_30_15:
      SensorHitMap::logWeighting(3.0, 1.5);
      break;
    case LOGWEIGHTING_40_10:
      SensorHitMap::logWeighting(4.0, 1.0);
      break;
    case LOGWEIGHTING_40_15:
      SensorHitMap::logWeighting(4.0, 1.5);
      break;
    case LOGWEIGHTING_50_10:
      SensorHitMap::logWeighting(5.0, 1.0);
      break;
    case LOGWEIGHTING_50_15:
      SensorHitMap::logWeighting(5.0, 1.5);
      break;
    case LOGWEIGHTING_60_10:
      SensorHitMap::logWeighting(6.0, 1.0);
      break;
    case LOGWEIGHTING_60_15:
      SensorHitMap::logWeighting(6.0, 1.5);
      break;
    case LOGWEIGHTING_70_10:
      SensorHitMap::logWeighting(7.0, 1.0);
      break;
    case LOGWEIGHTING_70_15:
      SensorHitMap::logWeighting(7.0, 1.5);
      break;
    //case NEWMETHOD:
      //SensorHitMap::newWeightingFunction()
      //break;
    default:
      SensorHitMap::poweredWeighting(2);
  }

  //set the displacement in z w.r.t. the rotational angles only
  centralHitZ = -d_gamma*centralHitPoint.first - d_beta*centralHitPoint.second;
}

double SensorHitMap::getResidualResolution() {
  return residualResolution;
};

std::pair<double, double> SensorHitMap::getHitPosition() {
  return centralHitPoint;
}

std::pair<double, double> SensorHitMap::getLabHitPosition() {
  double x_lab =   centralHitPoint.first - d_alpha * centralHitPoint.second - d_x0; 
  double y_lab = + d_alpha * centralHitPoint.first + centralHitPoint.second - d_y0;

  return std::make_pair(x_lab, y_lab);
}

std::pair<double, double> SensorHitMap::getHitPositionError() {
  return centralHitPointError;
}

std::pair<double, double> SensorHitMap::getCenterOfClosestCell(std::pair<double, double> X_ref) {
  std::vector<std::pair<double, HitData*>> to_sort;
  for(std::map<int, HitData*>::iterator hit=Hits.begin(); hit!=Hits.end(); hit++){
    if (filterByCellType((*hit).second->ID)) continue;
    double current_radius = sqrt(pow((*hit).second->x - X_ref.first,2) + pow((*hit).second->y - X_ref.second,2));
    to_sort.push_back(std::make_pair(current_radius, (*hit).second));
  }
  //sort the hits by their radial distance
  std::sort(to_sort.begin(), to_sort.end(), 
    [](const std::pair<double, HitData*>& a, const std::pair<double, HitData*> b){
      return a.first < b.first;
    }
  );  
  return std::make_pair(to_sort[0].second->x, to_sort[0].second->y);
}

double SensorHitMap::getTotalEnergy() {
  return this->totalEnergy;
}

double SensorHitMap::getTotalClusterEnergy(int N_considered) {
  return this->totalClusterEnergy[N_considered];
}

double SensorHitMap::getTotalWeight() {
  return this->totalWeight;
};

//private functions 

//only consider certain cellTypes for the N closest approach (analogous to the LayerSumAnalyzer)
//TODO
bool SensorHitMap::filterByCellType(int ID) {  
  /*
  ID = 0 : full cell,
  ID = 1 : calibration pad,
  ID = 2 : half cell,
  ID = 3 : mouse bite cell,
  ID = 4 : outer calibration pad cell,
  ID = 5 : merged cell
  */
  
  if (ID!=0 )  //we only want to consider the main cells in the middle for first estimation
    return true;      
  
  return false;
}

void SensorHitMap::considerNClosest(int N_considered) {     
  //better: ranking by radial distance and then take the closest ones
  HitsForPositioning.clear();
  if (N_considered < 0) {
    for(std::map<int, HitData*>::iterator hit=Hits.begin(); hit!=Hits.end(); hit++){
      if (filterByCellType((*hit).second->ID)) continue;
      HitsForPositioning.push_back((*hit).second);
    }
    return;
  }

  //calculate radial distance for all pairs to the most significant hit
  std::vector<std::pair<double, HitData*>> to_sort;
  for(std::map<int, HitData*>::iterator hit=Hits.begin(); hit!=Hits.end(); hit++){
    if (filterByCellType((*hit).second->ID)) continue;
    double current_radius = sqrt(pow((*hit).second->x - mostSignificantHit->x,2) + pow((*hit).second->y - mostSignificantHit->y,2));
    to_sort.push_back(std::make_pair(current_radius, (*hit).second));
  }

  //sort the hits by their radial distance
  std::sort(to_sort.begin(), to_sort.end(), 
    [](const std::pair<double, HitData*>& a, const std::pair<double, HitData*> b){
      return a.first < b.first;
    }
  );

  double maxDist = 0;
  if (N_considered == 7) maxDist = (0.5 + sqrt(3)) * HGCAL_TB_CELL::FULL_CELL_SIDE;
  else if (N_considered == 7) maxDist = (0.5 + sqrt(3) * 2.) * HGCAL_TB_CELL::FULL_CELL_SIDE;

  int considerCounter = 0;
  for (std::vector<std::pair<double, HitData*>>::iterator sorted_hit = to_sort.begin(); 
    sorted_hit != to_sort.end(); sorted_hit++) {
    if ((*sorted_hit).first > maxDist && considerCounter >= N_considered){ //fill at least Nconsider but also mind radial distance
      break;
    }
    considerCounter++;
    HitsForPositioning.push_back((*sorted_hit).second);
  }

  //additional cleanup to be safe
  to_sort.clear();
}

void SensorHitMap::considerClusters(int N_considered) {
  for (std::vector<int>::iterator clusterIndex = clusterIndexes[N_considered].begin();
  clusterIndex != clusterIndexes[N_considered].end(); clusterIndex++){
    HitsForPositioning.push_back(Hits[*clusterIndex]);
  }
}

void SensorHitMap::poweredWeighting(int exponent) {
  double numerator_x, numerator_y;
  double w;
  
  numerator_x = numerator_y = totalWeight = 0; 
  for(std::vector<HitData*>::iterator hit=HitsForPositioning.begin(); hit!=HitsForPositioning.end(); hit++){
    w = (*hit)->I >= 0.0 ? pow((*hit)->I, exponent) : 0.0;    //0.0 --> not included in the sum
    totalWeight += w;
    numerator_x += w*(*hit)->x;
    numerator_y += w*(*hit)->y;
  }

  //prevent divisions through zero
  if (totalWeight == 0.0) {   //equivalent to the statement that no sensor has fired 
    return; //defaults are preserved
  }

  centralHitPoint.first = numerator_x/totalWeight;
  centralHitPoint.second = numerator_y/totalWeight;

  //calculate the RMs
  numerator_x = numerator_y = 0.0;
  for(std::vector<HitData*>::iterator hit=HitsForPositioning.begin(); hit!=HitsForPositioning.end(); hit++){ 
    w = (*hit)->I >= 0.0 ? pow((*hit)->I, exponent) : 0.0; 
    numerator_x += w*pow((*hit)->x - centralHitPoint.first, 2);
    numerator_y += w*pow((*hit)->y - centralHitPoint.second, 2);
  }
  centralHitPointError.first = sqrt(numerator_x/totalWeight) > 1./sqrt(12.) ? sqrt(numerator_x/totalWeight): 1./sqrt(12.);
  centralHitPointError.second = sqrt(numerator_y/totalWeight) > 1./sqrt(12.) ? sqrt(numerator_y/totalWeight): 1./sqrt(12.);
};

void SensorHitMap::logWeighting(double log_a, double log_b) {
  double I_max = 0;   //determine the 'intensity' maximum first
  double I_i = 0;
  for(std::vector<HitData*>::iterator hit=HitsForPositioning.begin(); hit!=HitsForPositioning.end(); hit++){
    I_i = (*hit)->I >= 0.0 ? (*hit)->I : 0.0;    
    I_max += I_i;
  }
  

  double numerator_x, numerator_y;
  double w;
  numerator_x = numerator_y = totalWeight = 0; 
  
  int counter = 0;
  for(std::vector<HitData*>::iterator hit=HitsForPositioning.begin(); hit!=HitsForPositioning.end(); hit++){
    counter++;
    I_i = (*hit)->I >= 0.0 ? (*hit)->I : 0.0;    
    if (I_i == 0.) continue;  
    w = std::max(log_a + log_b*log(I_i/I_max), 0.0);
    totalWeight += w;
    numerator_x += w*(*hit)->x;
    numerator_y += w*(*hit)->y;    
  }

  //prevent divisions through zero
  if (totalWeight == 0.0) {   //equivalent to the statement that no sensor has fired 
    return; //defaults are preserved
  }

  centralHitPoint.first = numerator_x/totalWeight;
  centralHitPoint.second = numerator_y/totalWeight;

  //calculate the RMs
  numerator_x = numerator_y = 0.0;
  for(std::vector<HitData*>::iterator hit=HitsForPositioning.begin(); hit!=HitsForPositioning.end(); hit++){ 
    I_i = (*hit)->I >= 0.0 ? (*hit)->I : 0.0;    
    if (I_i == 0.) continue;  
    w = std::max(log_a + log_b*log(I_i/I_max), 0.0);
    numerator_x += w*pow((*hit)->x - centralHitPoint.first, 2);
    numerator_y += w*pow((*hit)->y - centralHitPoint.second, 2);
  }
  centralHitPointError.first = sqrt(numerator_x/totalWeight) > 1./sqrt(12.) ? sqrt(numerator_x/totalWeight): 1./sqrt(12.);
  centralHitPointError.second = sqrt(numerator_y/totalWeight) > 1./sqrt(12.) ? sqrt(numerator_y/totalWeight): 1./sqrt(12.);
};