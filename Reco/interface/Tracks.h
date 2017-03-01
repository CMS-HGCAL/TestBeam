#ifndef HGCAL_RECO_TRACK_H
#define HGCAL_RECO_TRACK_H

#include <iostream>
#include <utility>
#include <algorithm>
#include <cmath>
#include <fstream>

#include "TF1.h"
#include "TGraphErrors.h"
#include "TFitResult.h"
#include "TProfile.h"
#include "TMath.h"

#include "HGCal/Reco/interface/PositionResolutionHelpers.h"
#include "HGCal/Reco/interface/Sensors.h"
#include "FWCore/Utilities/interface/Exception.h"

#include "Alignment/ReferenceTrajectories/interface/GblTrajectory.h"
#include "Alignment/ReferenceTrajectories/interface/GblPoint.h"
#include "Alignment/ReferenceTrajectories/interface/MilleBinary.h"

enum TrackFittingMethod {
  DEFAULTFITTING,
  LINEFITANALYTICAL,
  LINEFITTGRAPHERRORS,
  POL2TGRAPHERRORS,
  POL3TGRAPHERRORS,
  GBLTRACK
};

enum FitPointWeightingMethod {
  NONE, 
  LINEAR,
  SQUARED,
  LOGARITHMIC,
  EXPONENTIAL
};

namespace gblhelpers {
  class layer {
    public: 
      layer(int label, double z0, double absorberX0, double particleEnergy, bool isReference, 
      std::pair<double, double> measurement = std::make_pair(-1, -1), std::pair<double, double> measurementError = std::make_pair(-1, -1)) : 
        _label(label), _z0(z0), _absorberX0(absorberX0), _particleEnergy(particleEnergy), _isReference(isReference), _measurement(measurement), _measurementError(measurementError) {}; 
          
      int label() {return _label;}
      double z() {return _z0;}
      double absorberX0() {return _absorberX0;}
      double particleEnergy() {return _particleEnergy;}
      bool isReference() {return _isReference;}
      std::pair<double, double> measurement() {return _measurement;};
      std::pair<double, double> measurementError() {return _measurementError;};

      bool operator < (const layer& lay) const {
        return (_z0 < lay._z0);
      }

    private:
      int _label;
      double _z0;
      double _absorberX0;
      double _particleEnergy;
      bool _isReference;
      std::pair<double, double> _measurement;
      std::pair<double, double> _measurementError;
 
  };
  double showerProfile(double *t, double *par);
  double computeEnergyLoss(double t, double E0);
}

//
// class declarations
//

class ParticleTrack{
  public:
    ParticleTrack();
    ~ParticleTrack();
    void addFitPoint(SensorHitMap* sensor);
    void addReferenceSensor(SensorHitMap* reference);
    void addDummySensor(SensorHitMap* dummy);
    void weightFitPoints(FitPointWeightingMethod method);
    void fitTrack(TrackFittingMethod method);
    std::pair<double, double> calculateReferenceXY();
    std::pair<double, double> calculateReferenceErrorXY();
    std::pair<double, double> calculatePositionXY(double z, int layerLabel);
    std::pair<double, double> calculatePositionErrorXY(double z, int layerLabel);
    double getSumOfEnergies();
  
    void gblTrackToMilleBinary(gbl::MilleBinary *milleBinary);
  private:
    int N_points; //cross check, should be identical to x.size() etc.
    SensorHitMap* referenceSensor;
    
    //general information, the fit points
    std::vector<double> x;  
    std::vector<double> x_err;  
    std::vector<double> y;  
    std::vector<double> y_err;  
    std::vector<double> z;  
    std::vector<double> z_err;  
    std::vector<double> particleEnergy;   //approximated electron energy left over at the individual layers
    std::vector<double> preAbsorber;      //absorber thickness in radiation lengths deposited directly in front of the sensor
    std::vector<double> Energies; //sum of energies deposited in the layers that are added to the track
    
    std::vector<gblhelpers::layer> gblLayers;
    std::map<int, int> gblPointLabels;

    TrackFittingMethod lastAppliedMethod;

    //different fit functions
    void gblTrackFit();
    std::pair<double, double> positionFromGblTrackFit(int layerLabel);
    std::pair<double, double> positionFromGblTrackFitErrors(int layerLabel);
    gbl::GblTrajectory* gblTrajectory;

    void analyticalStraightLineFit();
    std::pair<double, double> positionFromAnalyticalStraightLine(double z);
    std::pair<double, double> positionFromAnalyticalStraightLineErrors(double z);
    LineFitter* linefit_x;
    LineFitter* linefit_y;

    void polFitTGraphErrors(int degree);  
    std::pair<double, double> positionFromPolFitTGraphErrors(int degree, double z);
    std::pair<double, double> positionErrorFromPolFitTGraphErrors(int degree, double z);
    TF1* ROOTpol_x;
    TF1* ROOTpol_y;
    TFitResultPtr fit_result_x;
    TFitResultPtr fit_result_y;
};

double weightToFitPointWeight(double w, double sum_w, FitPointWeightingMethod m);




#endif