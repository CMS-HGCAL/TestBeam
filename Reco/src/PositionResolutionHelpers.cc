#include "HGCal/Reco/interface/PositionResolutionHelpers.h"


//****   Sensor Hit Maps    ****//

//public functions
SensorHitMap::SensorHitMap(){
  mostSignificantHit = NULL;  //will point to the most significant hit
  
  centralHitPoint = std::make_pair(0., 0.);
  centralHitPointError = std::make_pair(sqrt(12.), sqrt(12.));
  CM_threshold = 2.;
  layerZ_cm = 0;
  layerZ_X0 = 0;
  ADC_per_MIP = 1.;
  sensorSize = 128;

  CM_cells_count = 0;
  CM_sum = 0;

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

void SensorHitMap::setZ(double z_cm, double z_X0) {
  this->layerZ_cm = z_cm;
  this->layerZ_X0 = z_X0;
}

void SensorHitMap::setADCPerMIP(double ADC_per_MIP) {
  this->ADC_per_MIP = ADC_per_MIP;
}

void SensorHitMap::setPedestalThreshold(double t) {
  this->CM_threshold = t; 
}

double SensorHitMap::getZ_cm() {
  return this->layerZ_cm;
}

double SensorHitMap::getZ_X0() {
  return this->layerZ_X0;
}
//reduces the information from the Rechit towards what is necessary for the impact point calculation
//here all hits independent from the cellType are included
void SensorHitMap::addHit(HGCalTBRecHit Rechit) {
  int uniqueID = (Rechit.id()).rawId();

  CellCenterXY = TheCell.GetCellCentreCoordinatesForPlots((Rechit.id()).layer(), (Rechit.id()).sensorIU(), (Rechit.id()).sensorIV(), (Rechit.id()).iu(), (Rechit.id()).iv(), sensorSize);
  double iux = CellCenterXY.first;
  double ivy = CellCenterXY.second;
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

  if (totalClusterEnergy.find(N_considered) == totalClusterEnergy.end()) {
    totalClusterEnergy[N_considered] = 0;
  }
  totalClusterEnergy[N_considered] += Hits[uniqueID]->E;

  clusterIndexes[N_considered].push_back(uniqueID);
}

void SensorHitMap::subtractCM() {
  totalEnergy = 0.;
  double cm_subtraction = CM_cells_count > 0. ? CM_sum/CM_cells_count : 0.;

  for(std::map<int, HitData*>::iterator hit=Hits.begin(); hit!=Hits.end(); hit++){
    if (filterByCellType((*hit).second->ID)) continue;
    // we want: - all cells that were input to the cm (common mode) calculation get weight 0
    //          - all the others are corrected by the cm
    if ((*hit).second->E > CM_threshold) {    
      (*hit).second->I = (*hit).second->I - cm_subtraction;
    } else {
      (*hit).second->I = std::max((*hit).second->I - cm_subtraction, 0.);
    }
    totalEnergy += (*hit).second->I;
  }
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
}

std::pair<double, double> SensorHitMap::getCenterPosition() {
  return centralHitPoint;
}
std::pair<double, double> SensorHitMap::getCenterPositionError() {
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
  if (ID!=0 )  //we only want to consider the main cells in the middle for first estimation
    return true;      //TODO
  */
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
  centralHitPointError.first = sqrt(numerator_x/totalWeight);
  centralHitPointError.second = sqrt(numerator_y/totalWeight);
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
  centralHitPointError.first = sqrt(numerator_x/totalWeight);
  centralHitPointError.second = sqrt(numerator_y/totalWeight);
};


//****   Particle Tracks    ****//

//public functions
ParticleTrack::ParticleTrack(){
  lastAppliedMethod = DEFAULTFITTING;
  ROOTpol_x = ROOTpol_y = 0;

  N_points = 0;
};

ParticleTrack::~ParticleTrack(){
  delete ROOTpol_x;
  delete ROOTpol_y;
  x.clear(); x_err.clear(); y.clear(); y_err.clear(); z.clear(); z_err.clear();
};

void ParticleTrack::addFitPoint(SensorHitMap* sensor){
  N_points++;
  x.push_back(sensor->getCenterPosition().first);  
  x_err.push_back(sensor->getCenterPositionError().first);  
  y.push_back(sensor->getCenterPosition().second);  
  y_err.push_back(sensor->getCenterPositionError().second);  
  z.push_back(sensor->getZ_cm());  
  z_err.push_back(0.0);
  Energies.push_back(sensor->getTotalEnergy());
};

void ParticleTrack::weightFitPoints(FitPointWeightingMethod method) {
  //calculate the sum of all Energies first 
  double sum_Energies = 0.0;
  for (int i=0; i<N_points; i++) {
    sum_Energies += Energies[i];
  }

  //we want to assure that the sum of Energies is conserved to not majorly influence the goodness of the fit
  double nom_x=0.; double nom_y=0.; double denom_x=0.; double denom_y=0.; double new_weight = 1.0;
  for (int i=0; i<N_points; i++) {
    new_weight = weightToFitPointWeight(Energies[i], sum_Energies, method);
    nom_x += x_err[i];
    nom_y += y_err[i];
    denom_x  += x_err[i]*new_weight; //original Energies are untouched
    denom_y  += y_err[i]*new_weight;
  }
  for (int i=0; i<N_points; i++) {
    new_weight = weightToFitPointWeight(Energies[i], sum_Energies, method);
    x_err[i] = (nom_x/denom_x)*x_err[i]*new_weight;
    y_err[i] = (nom_y/denom_y)*y_err[i]*new_weight;
  }
}

void ParticleTrack::fitTrack(TrackFittingMethod method){
  try {
    switch(method) {
      case LINEFITTGRAPHERRORS:
        polFitTGraphErrors(1);
        break;
      case POL2TGRAPHERRORS:
        polFitTGraphErrors(2);
        break;
      case POL3TGRAPHERRORS:
        polFitTGraphErrors(3);
        break;
      default:
        lastAppliedMethod = DEFAULTFITTING;
        break;
    }
    lastAppliedMethod = method;      
  } catch(cms::Exception& iException) {
    //std::cout<<"Fitting method has failed with error code "<< iException <<". Attempting the default fitting."<<std::endl<<std::endl;
    fitTrack(DEFAULTFITTING);
  }
}

std::pair<double, double> ParticleTrack::calculatePositionXY(double z) {
  switch(lastAppliedMethod) {
    case LINEFITTGRAPHERRORS:
      return positionFromPolFitTGraphErrors(1, z);
    case POL2TGRAPHERRORS:
      return positionFromPolFitTGraphErrors(2, z);
    case POL3TGRAPHERRORS:
      return positionFromPolFitTGraphErrors(3, z);
    default:
      return std::make_pair(0.,0.); 
  }
}
std::pair<double, double> ParticleTrack::calculatePositionErrorXY(double z) {
  switch(lastAppliedMethod) {
    case LINEFITTGRAPHERRORS:
      return positionErrorFromPolFitTGraphErrors(1, z);
    case POL2TGRAPHERRORS:
      return positionErrorFromPolFitTGraphErrors(2, z);
    case POL3TGRAPHERRORS:
      return positionErrorFromPolFitTGraphErrors(3, z);
    default:
      return std::make_pair(0.,0.); 
  }
}

double ParticleTrack::getSumOfEnergies() {  //returns the original Energies (i.e. sum(layers for fit) sum(hits in layer) of WeightingMethod(intensity))
  double sum_Energies = 0.0;
  for (int i=0; i<N_points; i++) {
    sum_Energies += Energies[i];
  }
  return sum_Energies;
}

//private functions
void ParticleTrack::polFitTGraphErrors(int degree){  
  //todo: clear existing Polynomial pointers if existing
  if (ROOTpol_x != 0)
    delete ROOTpol_x;
  if (ROOTpol_y != 0)
    delete ROOTpol_y;

  ROOTpol_x = new TF1("ROOTpol_x", ("pol"+std::to_string(degree)).c_str(), *min_element(z.begin(), z.end())-1.0, *max_element(z.begin(), z.end())+1.0);
  ROOTpol_y = new TF1("ROOTpol_y", ("pol"+std::to_string(degree)).c_str(), *min_element(z.begin(), z.end())-1.0, *max_element(z.begin(), z.end())+1.0);
  
  TGraphErrors* tmp_graph_x = new TGraphErrors(z.size(), &(z[0]), &(x[0]), &(z_err[0]), &(x_err[0])); //z_err should be filled with zeros
  TGraphErrors* tmp_graph_y = new TGraphErrors(z.size(), &(z[0]), &(y[0]), &(z_err[0]), &(y_err[0]));
  
  fit_result_x = tmp_graph_x->Fit(ROOTpol_x, "QFS");    //for some ROOT versions, the option 'S' is essential to retrieve the fit results object
  fit_result_y = tmp_graph_y->Fit(ROOTpol_y, "QFS");

  delete tmp_graph_x;
  delete tmp_graph_y;
  
}; 

std::pair<double, double> ParticleTrack::positionFromPolFitTGraphErrors(int degree, double z) {
  return std::make_pair(ROOTpol_x->Eval(z), ROOTpol_y->Eval(z));
}

std::pair<double, double> ParticleTrack::positionErrorFromPolFitTGraphErrors(int degree, double z) {
  double err_x = 0.; double err_y = 0.;
  for (int i=0; i<=degree; i++){
    err_x += pow( ROOTpol_x->GetParError(i) * pow(z, i) ,2);
    err_y += pow( ROOTpol_y->GetParError(i) * pow(z, i) ,2);
    
    for (int j=i+1; j<=degree; j++) {
      err_x += 2*fit_result_x->Correlation(i, j) * ROOTpol_x->GetParError(i) * ROOTpol_x->GetParError(j) * fabs(pow(z, i+j));
      err_y += 2*fit_result_y->Correlation(i, j) * ROOTpol_y->GetParError(i) * ROOTpol_y->GetParError(j) * fabs(pow(z, i+j));
    }
  }
  return std::make_pair(err_x, err_y);
}


double weightToFitPointWeight(double e, double sum_e, FitPointWeightingMethod m) {
  switch(m) {
    case NONE:
      return 1.0;
    case LINEAR:
      return 1.0/(1.0+e/sum_e);
    case SQUARED:
      return pow(1.0/(1.0+e/sum_e) ,2);
    case LOGARITHMIC:
      return 1.0/(1.0+log(1.0+e/sum_e));
    case EXPONENTIAL:
      return exp(-e/sum_e);
    default:
      return 1.0;
  }  
}