#ifndef HGCAL_RECO_PARTICLETRACK_H
#define HGCAL_RECO_PARTICLETRACK_H

#include <iostream>
#include <utility>
#include <algorithm>
#include <cmath>

#include "TF1.h"
#include "TGraphErrors.h"
#include "TFitResult.h"

#include "HGCal/DataFormats/interface/HGCalTBRecHitCollections.h"
#include "HGCal/DataFormats/interface/HGCalTBRecHit.h"
#include "HGCal/DataFormats/interface/HGCalTBClusterCollection.h"
#include "HGCal/DataFormats/interface/HGCalTBDetId.h"
#include "HGCal/Geometry/interface/HGCalTBCellVertices.h"
#include "HGCal/Geometry/interface/HGCalTBCellParameters.h" //e.g. to get the cell's dimensions

#include "FWCore/Utilities/interface/Exception.h"


enum ConsiderationMethod {
  CONSIDERALL,
  CONSIDERSEVEN,
  CONSIDERNINETEEN,
  CONSIDERCLUSTERSALL,
  CONSIDERCLUSTERSSEVEN,
  CONSIDERCLUSTERSNINETEEN
};

enum WeightingMethod {
  DEFAULTWEIGHTING,
  SQUAREDWEIGHTING,
  LINEARWEIGHTING,
  LOGWEIGHTING_50_10,
  LOGWEIGHTING_50_05,
  LOGWEIGHTING_70_10
};

enum TrackFittingMethod {
  DEFAULTFITTING,
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

struct HitData {
  double x; double y; 
  double I; //intensity that is input to the weight calculation
  double E; //actual energy of the hit
  int ID;   //the ID corresponds to the cell ID, it is necessary for the pedestal subtraction
};

//
// class declarations
//
class SensorHitMap {
  private:
    std::pair<double, double> centralHitPoint;
    std::pair<double, double> centralHitPointError;
    double layerZ_cm;
    double layerZ_X0;
    int sensorSize;
    double CM_threshold;
    double ADC_per_MIP;
    std::map<int, HitData*> Hits; //those are all the hits in the layer
    std::map<int, std::vector<int>> clusterIndexes;
    std::vector<HitData*> HitsForPositioning;
    HitData* mostSignificantHit;

    //helpers to obtain the x-y coordinate
    HGCalTBCellVertices TheCell;
    std::pair<double, double> CellCenterXY;
    int CM_cells_count;
    double CM_sum;

    double totalWeight; //equivalent to the denominator in the according weighting method

    bool filterByCellType(HitData* hit);
    void considerNClosest(int N_considered);
    void considerClusters(int N_considered);
    void poweredWeighting(int exponent);
    void logWeighting(double log_a, double log_b);

  public:
    SensorHitMap();
    ~SensorHitMap();
    void setZ(double z_cm, double z_X0);
    void setADCPerMIP(double ADC_per_MIP);
    void setSensorSize(int s);
    void setPedestalThreshold(double t);
    //reduces the information from the Rechit towards what is necessary for the impact point calculation
    void addHit(HGCalTBRecHit Rechit);
    void addClusterHit(HGCalTBDetId hit, int N_considered);
    void subtractCM();
    void calculateCenterPosition(ConsiderationMethod considerationMethod, WeightingMethod weightingMethod);
    double getTotalWeight();
    double getZ_cm();
    double getZ_X0();
    std::pair<double, double> getCenterPosition();
    std::pair<double, double> getCenterPositionError(); //calculated via RMS

    //debug
    void printHits();
};


class ParticleTrack{
  public:
    ParticleTrack();
    ~ParticleTrack();
    void addFitPoint(SensorHitMap* sensor);
    void weightFitPoints(FitPointWeightingMethod method);
    void fitTrack(TrackFittingMethod method);
    std::pair<double, double> calculatePositionXY(double z);
    std::pair<double, double> calculatePositionErrorXY(double z);
    double getSumOfWeights();
  private:
    int N_points; //cross check, should be identical to x.size() etc.
    //general information, the fit points
    std::vector<double> x;  
    std::vector<double> x_err;  
    std::vector<double> y;  
    std::vector<double> y_err;  
    std::vector<double> z;  
    std::vector<double> z_err;  
    std::vector<double> weights;
    TrackFittingMethod lastAppliedMethod;

    //different fit functions
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
