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


//
// class declaration
//

enum WeightingMethod {
  DEFAULTWEIGHTING,
  SQUAREDWEIGHTING,
  LINEARWEIGHTING 
};

enum TrackFittingMethod {
  DEFAULTFITTING,
  LINEFITTGRAPHERRORS
};

struct HitTriple {
  double x; double y; 
  double I;
  int ID;   //the ID corresponds to the cell ID, it is necessary for the pedestal subtraction
};

class SensorHitMap {
  private:
    std::pair<double, double> centralHitPoint;
    std::pair<double, double> centralHitPointError;
    double layerZ;
    int sensorSize;
    std::vector<HitTriple*> Hits;
    //helpers to obtain the x-y coordinate
    HGCalTBCellVertices TheCell;
    std::pair<double, double> CellCenterXY;
    std::map<int, int> cellTypeCount;
    std::map<int, double> pedestalCount;

    void poweredWeighting(int exponent);

  public:
    SensorHitMap();
    ~SensorHitMap();
    void setZ(double z);
    void setSensorSize(int s);
    //reduces the information from the Rechit towards what is necessary for the impact point calculation
    void addHit(HGCalTBRecHit Rechit);
    void subtractPedestals();
    void calculateCenterPosition(WeightingMethod method);
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
    void lineFitTGraphErrors();  
    std::pair<double, double> positionFromLineFitTGraphErrors(double z);
    TF1* ROOTpol_x;
    TF1* ROOTpol_y;
    TGraph* tmp_graph_x;
    TGraph* tmp_graph_y;
};


#endif
