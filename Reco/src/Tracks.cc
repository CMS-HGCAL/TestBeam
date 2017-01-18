#include "HGCal/Reco/interface/Tracks.h"

//****   Particle Tracks    ****//

//public functions
ParticleTrack::ParticleTrack(){
  lastAppliedMethod = DEFAULTFITTING;
  ROOTpol_x = ROOTpol_y = NULL;
  linefit_x = linefit_y = NULL;

  N_points = 0;
};

ParticleTrack::~ParticleTrack(){
  delete ROOTpol_x;
  delete ROOTpol_y;
  delete linefit_x;
  delete linefit_y;
  x.clear(); x_err.clear(); y.clear(); y_err.clear(); z.clear(); z_err.clear();
};

void ParticleTrack::addFitPoint(SensorHitMap* sensor){
  N_points++;
  x.push_back(sensor->getLabHitPosition().first);  
  x_err.push_back(sensor->getHitPositionError().first);  
  y.push_back(sensor->getLabHitPosition().second);  
  y_err.push_back(sensor->getHitPositionError().second);  
  z.push_back(sensor->getLabZ()+sensor->getIntrinsicHitZPosition());  
  z_err.push_back(0.0);
  particleEnergy.push_back(sensor->getParticleEnergy());
  preAbsorber.push_back(sensor->getX0()); 

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
      case LINEFITANALYTICAL: 
        analyticalStraightLineFit();
        break;
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
    case LINEFITANALYTICAL:
      return positionFromAnalyticalStraightLine(z);
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
    case LINEFITANALYTICAL:
      return positionFromAnalyticalStraightLineErrors(z);
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
void ParticleTrack::gblTrackFit() {};
std::pair<double, double> ParticleTrack::positionFromGblTrackFit(double z) {
  return std::make_pair(0.0, 0.0);
};

std::pair<double, double> ParticleTrack::positionFromGblTrackFitErrors(double z) {
  return std::make_pair(0.0, 0.0);
};





void ParticleTrack::analyticalStraightLineFit() {
  if (linefit_x != NULL)
    delete linefit_x;
  if (linefit_y != NULL)
    delete linefit_y;

  linefit_x = new LineFitter(z, x, x_err);
  linefit_y = new LineFitter(z, y, y_err);
  linefit_x->fit();
  linefit_y->fit();

  polFitTGraphErrors(1);
};

std::pair<double, double> ParticleTrack::positionFromAnalyticalStraightLine(double z) {
  return std::make_pair(linefit_x->eval(z), linefit_y->eval(z));
};
std::pair<double, double> ParticleTrack::positionFromAnalyticalStraightLineErrors(double z) {
  return std::make_pair(linefit_x->evalError(z), linefit_y->evalError(z));
};



void ParticleTrack::polFitTGraphErrors(int degree){  
  if (ROOTpol_x != NULL)
    delete ROOTpol_x;
  if (ROOTpol_y != NULL)
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
  return std::make_pair(sqrt(err_x), sqrt(err_y));
}


double weightToFitPointWeight(double e, double sum_e, FitPointWeightingMethod m) {
  switch(m) {
    case NONE:
      return 1.0;
    case LINEAR:
      return 1.0/(0.01+e/sum_e);
    case SQUARED:
      return 1.0/(0.01+pow(e/sum_e,2));
    case LOGARITHMIC:
      return 1.0/(0.01+log(1.0+e/sum_e));
    case EXPONENTIAL:
      return exp(-e/sum_e);
    default:
      return 1.0;
  }  
}