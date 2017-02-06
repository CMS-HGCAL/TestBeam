#include "HGCal/Reco/interface/Tracks.h"

//****   gblhelpers    ****//
double gblhelpers::showerProfile(double *t, double *par) {
  return par[0]*par[2]*(pow(par[2]*t[0], par[1]-1)*exp(-par[2]*t[0]))/TMath::Gamma(par[1]);
}

double gblhelpers::computeEnergyLoss(double t, double E0) {
  double tmax = 2.6648 * pow(10., -7) * pow(E0, 3) + -1.72 * pow(10., -4) * pow(E0, 2) + 0.05003425 * E0 + 4.09805383;

  TF1 *profile = new TF1("profile", showerProfile, 0., 50., 3);
  profile->SetParameter(0, E0);
  profile->SetParameter(1, 0.5*tmax+1); //valid approximation for b ~0.5 yields this a
  profile->SetParameter(2, 0.5);

  double integral = profile->Integral(0, t);
  delete profile;
  return integral;
}


//****   Particle Tracks    ****//

//public functions
ParticleTrack::ParticleTrack(){
  lastAppliedMethod = DEFAULTFITTING;
  ROOTpol_x = ROOTpol_y = NULL;
  linefit_x = linefit_y = NULL;
  gblTrajectory = NULL;

  N_points = 0;
};

ParticleTrack::~ParticleTrack(){
  delete ROOTpol_x;
  delete ROOTpol_y;
  delete linefit_x;
  delete linefit_y;
  delete gblTrajectory;
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
  preAbsorber.push_back(sensor->getX0()); 
  particleEnergy.push_back(sensor->getParticleEnergy());

  Energies.push_back(sensor->getTotalEnergy());

  //This is for gbl: Check why this has to be set to exactly those values for it to work
  std::pair<double, double> measurementError = sensor->getHitPositionError();
  measurementError.first = sensor->getResidualResolution() ==-1 ? measurementError.first : sensor->getResidualResolution();
  measurementError.second = sensor->getResidualResolution() ==-1 ? measurementError.second : sensor->getResidualResolution();
  gblLayers.push_back(gblhelpers::layer(sensor->label(), sensor->getLabZ()+sensor->getIntrinsicHitZPosition(), sensor->getX0(), sensor->getParticleEnergy(), false, sensor->getLabHitPosition(), measurementError));
};

void ParticleTrack::addReferenceSensor(SensorHitMap* sensor) {
  referenceSensor = sensor;
  gblLayers.push_back(gblhelpers::layer(sensor->label(), sensor->getLabZ()+sensor->getIntrinsicHitZPosition(), sensor->getX0(), sensor->getParticleEnergy(), true));
}

void ParticleTrack::addDummySensor(SensorHitMap* sensor) {
  gblLayers.push_back(gblhelpers::layer(sensor->label(), sensor->getLabZ()+sensor->getIntrinsicHitZPosition(), sensor->getX0(), sensor->getParticleEnergy(), true));
}

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
      case GBLTRACK:
        gblTrackFit();
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

std::pair<double, double> ParticleTrack::calculateReferenceXY() {
  return ParticleTrack::calculatePositionXY(referenceSensor->getLabZ() + referenceSensor->getIntrinsicHitZPosition(), referenceSensor->label());
};

std::pair<double, double> ParticleTrack::calculateReferenceErrorXY() {
  return ParticleTrack::calculatePositionErrorXY(referenceSensor->getLabZ() + referenceSensor->getIntrinsicHitZPosition(), referenceSensor->label());
};

std::pair<double, double> ParticleTrack::calculatePositionXY(double z, int layerLabel) {
  switch(lastAppliedMethod) {
    case LINEFITANALYTICAL:
      return positionFromAnalyticalStraightLine(z);
    case LINEFITTGRAPHERRORS:
      return positionFromPolFitTGraphErrors(1, z);
    case POL2TGRAPHERRORS:
      return positionFromPolFitTGraphErrors(2, z);
    case POL3TGRAPHERRORS:
      return positionFromPolFitTGraphErrors(3, z);
    case GBLTRACK:
      return positionFromGblTrackFit(layerLabel);
    default:
      return std::make_pair(0.,0.); 
  }
}
std::pair<double, double> ParticleTrack::calculatePositionErrorXY(double z, int layerLabel) {
  switch(lastAppliedMethod) {
    case LINEFITANALYTICAL:
      return positionFromAnalyticalStraightLineErrors(z);
    case LINEFITTGRAPHERRORS:
      return positionErrorFromPolFitTGraphErrors(1, z);
    case POL2TGRAPHERRORS:
      return positionErrorFromPolFitTGraphErrors(2, z);
    case POL3TGRAPHERRORS:
      return positionErrorFromPolFitTGraphErrors(3, z);
    case GBLTRACK:
      return positionFromGblTrackFitErrors(layerLabel);
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

void ParticleTrack::gblTrackToMilleBinary(gbl::MilleBinary *milleBinary) {
  if (gblTrajectory != NULL) {
    //std::cout<<"Writing to MilleBinary"<<std::endl;
    gblTrajectory->milleOut(*milleBinary);
  }
}

//private functions
void ParticleTrack::gblTrackFit() {
  std::vector<gbl::GblPoint> listOfGblPoints;
  std::sort(gblLayers.begin(), gblLayers.end());
  
  double dz, particleEnergy, layerX0, thetaRMS_abs;//, thetaRMS_sensor;
  int label;

  double z_prev = gblLayers[0].z() - 1.;    //1cm to incorporate possible absorbers

  for (size_t i=0; i<gblLayers.size(); i++) {
    dz = gblLayers[i].z() - z_prev;
    z_prev = gblLayers[i].z();
    particleEnergy = gblLayers[i].particleEnergy();
    layerX0 = gblLayers[i].absorberX0();
    label = gblLayers[i].label();

    thetaRMS_abs = (0.0136*sqrt(layerX0)/particleEnergy*(1+0.038*log(layerX0)));     //according to PDG formula

    TMatrixD jac_abs1(5, 5);
    jac_abs1.UnitMatrix();
    jac_abs1[3][1] = (0.5-sqrt(1/12.))*dz;
    jac_abs1[4][2] = (0.5-sqrt(1/12.))*dz;
    gbl::GblPoint point_abs1(jac_abs1);
    TVectorD scat_mean_abs1(2); //mean of kink angle variation
    scat_mean_abs1.Zero();
    TVectorD scat_invResRMS_abs1(2);
    scat_invResRMS_abs1[0] = 1./(2.*thetaRMS_abs*thetaRMS_abs);
    scat_invResRMS_abs1[1] = 1./(2.*thetaRMS_abs*thetaRMS_abs);
    point_abs1.addScatterer(scat_mean_abs1, scat_invResRMS_abs1);
    listOfGblPoints.push_back(point_abs1);

    
    TMatrixD jac_abs2(5, 5);
    jac_abs2.UnitMatrix();
    jac_abs2[3][1] = 2*sqrt(1./12.)*dz;
    jac_abs2[4][2] = 2*sqrt(1./12.)*dz;
    gbl::GblPoint point_abs2(jac_abs2);
    TVectorD scat_mean_abs2(2); //mean of kink angle variation
    scat_mean_abs2.Zero();
    TVectorD scat_invResRMS_abs2(2);
    scat_invResRMS_abs2[0] = 1./(2.*thetaRMS_abs*thetaRMS_abs);
    scat_invResRMS_abs2[1] = 1./(2.*thetaRMS_abs*thetaRMS_abs);
    point_abs2.addScatterer(scat_mean_abs2, scat_invResRMS_abs2);
    listOfGblPoints.push_back(point_abs2);
    
   
    //1. Make Jacobian (5x5) from plane i-1 to i  
    TMatrixD jac(5, 5);
    jac.UnitMatrix();
    jac[3][1] = (0.5-sqrt(1/12.))*dz;
    jac[4][2] = (0.5-sqrt(1/12.))*dz;
    //2. initiate the point
    gbl::GblPoint point(jac);

    if (!gblLayers[i].isReference()) {
      //3. add measurement
      TVectorD meas(2);
      meas[0] = (gblLayers[i].measurement()).first;
      meas[1] = (gblLayers[i].measurement()).second;
      //4. add measurement precision
      TVectorD meas_prec(2);
      meas_prec[0] = 1./pow((gblLayers[i].measurementError()).first, 2);
      meas_prec[1] = 1./pow((gblLayers[i].measurementError()).second,2);
      //5. add the measurement to the point
      TMatrixD proL2m(2,2);
      proL2m.UnitMatrix();
      point.addMeasurement(proL2m, meas, meas_prec);
      //6. add global alignment parameters
      std::vector<int> labels; labels.push_back(100*label + 11); labels.push_back(100*label + 12); labels.push_back(100*label + 21);
      TMatrixD glder(2, 3); glder[0][0] = glder[1][1] = 1.;  glder[0][1] = glder[1][0] = 0.; glder[0][2] = meas[1]; glder[1][2] = -meas[0];
      point.addGlobals(labels, glder);
    }

    listOfGblPoints.push_back(point);
    gblPointLabels[label] = listOfGblPoints.size();
  }


  gblTrajectory = new gbl::GblTrajectory(listOfGblPoints, false);
  

  double chi2, lw; int ndf;
  gblTrajectory->fit(chi2, ndf, lw);
  //std::cout<<"chi2: "<<chi2<<"   lw: "<<lw<<"    ndf: "<<ndf<<"  listOfGblPoints: "<<listOfGblPoints.size()<<std::endl;
};

std::pair<double, double> ParticleTrack::positionFromGblTrackFit(int layerLabel) {
  TVectorD aCorr(5);
  TMatrixDSym aCov(5);
  gblTrajectory->getResults(gblPointLabels[layerLabel], aCorr, aCov);
  return std::make_pair(aCorr(3), aCorr(4));
};

std::pair<double, double> ParticleTrack::positionFromGblTrackFitErrors(int layerLabel) {
  TVectorD aCorr(5);
  TMatrixDSym aCov(5);
  gblTrajectory->getResults(gblPointLabels[layerLabel], aCorr, aCov);
  return std::make_pair(sqrt(aCov(3,3)), sqrt(aCov(4,4)));
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