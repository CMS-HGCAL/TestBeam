#ifndef HGCAL_RECO_PARTICLETRACK_H
#define HGCAL_RECO_PARTICLETRACK_H

#include <iostream>
#include <utility>
#include <algorithm>
#include <cmath>

#include "TF1.h"
#include "TGraphErrors.h"

#include "HGCal/DataFormats/interface/HGCalTBRecHitCollections.h"
#include "HGCal/DataFormats/interface/HGCalTBRecHit.h"
#include "HGCal/Geometry/interface/HGCalTBCellVertices.h"

#include "FWCore/Utilities/interface/Exception.h"


enum ConsiderationMethod {
  CONSIDERALL,
  CONSIDERSEVEN,
  CONSIDERNINETEEN,
  CONSIDERCLUSTERS
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

struct HitData {
  double x; double y; 
  double I; //intensity that is input to the weight calculation
  double E; //actual energy of the hit
  int ID;   //the ID corresponds to the cell ID, it is necessary for the pedestal subtraction
};

struct DeviationTriple {
  double deviation, predicted_x, predicted_y;
};

//
// class declarations
//
class SensorHitMap {
  private:
    std::pair<double, double> centralHitPoint;
    std::pair<double, double> centralHitPointError;
    double layerZ;
    int sensorSize;
    double CM_threshold;
    double ADC_per_MIP;
    std::vector<HitData*> Hits; //those are all the hits in the layer
    std::vector<HitData*> HitsForPositioning;
    HitData* mostSignificantHit;

    //helpers to obtain the x-y coordinate
    HGCalTBCellVertices TheCell;
    std::pair<double, double> CellCenterXY;
    int CM_cells_count;
    double CM_sum;

    void fillHitsForPositioningByRadius(double R);
    void poweredWeighting(int exponent);
    void logWeighting(double log_a, double log_b);

  public:
    SensorHitMap();
    ~SensorHitMap();
    void setZ(double z);
    void setADCPerMIP(double ADC_per_MIP);
    void setSensorSize(int s);
    void setPedestalThreshold(double t);
    //reduces the information from the Rechit towards what is necessary for the impact point calculation
    void addHit(HGCalTBRecHit Rechit);
    void subtractCM();
    void calculateCenterPosition(ConsiderationMethod considerationMethod, WeightingMethod weightingMethod);
    double getZ();
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
    void fitTrack(TrackFittingMethod method);
    std::pair<double, double> calculatePositionXY(double z);
  private:
    //general information, the fit points
    std::vector<double> x;  
    std::vector<double> x_err;  
    std::vector<double> y;  
    std::vector<double> y_err;  
    std::vector<double> z;  
    std::vector<double> z_err;  
    TrackFittingMethod lastAppliedMethod;

    //different fit functions
    void polFitTGraphErrors(int degree);  
    std::pair<double, double> positionFromPolFitTGraphErrors(double z);
    TF1* ROOTpol_x;
    TF1* ROOTpol_y;
    TGraph* tmp_graph_x;
    TGraph* tmp_graph_y;
};


#endif
