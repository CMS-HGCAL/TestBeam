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
  POL3TGRAPHERRORS
};

enum FitPointWeightingMethod {
  NONE, 
  LINEAR,
  SQUARED,
  LOGARITHMIC,
  EXPONENTIAL
};

//
// class declarations
//

class ParticleTrack{
  public:
    ParticleTrack();
    ~ParticleTrack();
    void addFitPoint(SensorHitMap* sensor);
    void weightFitPoints(FitPointWeightingMethod method);
    void fitTrack(TrackFittingMethod method);
    std::pair<double, double> calculatePositionXY(double z);
    std::pair<double, double> calculatePositionErrorXY(double z);
    double getSumOfEnergies();
  private:
    int N_points; //cross check, should be identical to x.size() etc.
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
    TrackFittingMethod lastAppliedMethod;

    //different fit functions
    void gblTrackFit();
    std::pair<double, double> positionFromGblTrackFit(double z);
    std::pair<double, double> positionFromGblTrackFitErrors(double z);
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