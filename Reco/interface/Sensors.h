#ifndef HGCAL_RECO_SENSORS_H
#define HGCAL_RECO_SENSORS_H

#include <iostream>
#include <utility>
#include <algorithm>
#include <cmath>
#include <fstream>

#include "HGCal/DataFormats/interface/HGCalTBRecHitCollections.h"
#include "HGCal/DataFormats/interface/HGCalTBRecHit.h"
#include "HGCal/DataFormats/interface/HGCalTBClusterCollection.h"
#include "HGCal/DataFormats/interface/HGCalTBDetId.h"
#include "HGCal/Geometry/interface/HGCalTBCellParameters.h" //e.g. to get the cell's dimensions


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
  LOGWEIGHTING_20_10,
  LOGWEIGHTING_21_10,
  LOGWEIGHTING_22_10,
  LOGWEIGHTING_23_10,
  LOGWEIGHTING_24_10,
  LOGWEIGHTING_25_10,
  LOGWEIGHTING_26_10,
  LOGWEIGHTING_27_10,
  LOGWEIGHTING_28_10,
  LOGWEIGHTING_29_10,
  LOGWEIGHTING_30_10,
  LOGWEIGHTING_31_10,
  LOGWEIGHTING_32_10,
  LOGWEIGHTING_33_10,
  LOGWEIGHTING_34_10,
  LOGWEIGHTING_35_10,
  LOGWEIGHTING_36_10,
  LOGWEIGHTING_37_10,
  LOGWEIGHTING_38_10,
  LOGWEIGHTING_39_10,
  LOGWEIGHTING_40_10,
  LOGWEIGHTING_41_10,
  LOGWEIGHTING_42_10,
  LOGWEIGHTING_43_10,
  LOGWEIGHTING_44_10,
  LOGWEIGHTING_45_10,
  LOGWEIGHTING_46_10,
  LOGWEIGHTING_47_10,
  LOGWEIGHTING_48_10,
  LOGWEIGHTING_49_10,
  LOGWEIGHTING_50_10,
  LOGWEIGHTING_60_10,
  LOGWEIGHTING_70_10,
  LOGWEIGHTING_205_10,
  LOGWEIGHTING_215_10,
  LOGWEIGHTING_225_10,
  LOGWEIGHTING_235_10,
  LOGWEIGHTING_245_10,
  LOGWEIGHTING_255_10,
  LOGWEIGHTING_265_10,
  LOGWEIGHTING_275_10,
  LOGWEIGHTING_285_10,
  LOGWEIGHTING_295_10,
  LOGWEIGHTING_305_10,
  LOGWEIGHTING_315_10,
  LOGWEIGHTING_325_10,
  LOGWEIGHTING_335_10,
  LOGWEIGHTING_345_10,
  LOGWEIGHTING_355_10,
  LOGWEIGHTING_365_10,
  LOGWEIGHTING_375_10,
  LOGWEIGHTING_385_10,
  LOGWEIGHTING_395_10,
  LOGWEIGHTING_405_10,
  LOGWEIGHTING_415_10,
  LOGWEIGHTING_425_10,
  LOGWEIGHTING_435_10,
  LOGWEIGHTING_445_10,
  LOGWEIGHTING_455_10,
  LOGWEIGHTING_465_10,
  LOGWEIGHTING_475_10,
  LOGWEIGHTING_485_10,
  LOGWEIGHTING_495_10
};

struct HitData {
  double x; double y; 
  double I; //intensity that is input to the weight calculation
  double E; //actual energy of the hit
  int ID;   //the ID corresponds to the cell ID, it is necessary for the pedestal subtraction
};


class SensorHitMap {
  private:
    int _label;
    std::pair<double, double> centralHitPoint;
    double centralHitZ;
    std::pair<double, double> centralHitPointError;
    double layerLabZ;
    double layerX0;
    double particleEnergy;  //approximated energy of the particle at that layer  
    double d_alpha, d_beta, d_gamma, d_x0, d_y0, d_z0; //alignment parameters
    double residualResolution;    //externally set, from residual histograms
    int sensorSize;
    double CM_threshold;
    std::map<int, HitData*> Hits; //those are all the hits in the layer
    std::map<int, std::vector<int>> clusterIndexes;
    std::vector<HitData*> HitsForPositioning;
    HitData* mostSignificantHit;

    //helpers to obtain the x-y coordinate
    //HGCalTBCellVertices TheCell;
    //std::pair<double, double> CellCenterXY;
    int CM_cells_count;
    double CM_sum;

    double totalWeight; //equivalent to the denominator in the according weighting method
    double totalEnergy;
    std::map<int, double> totalClusterEnergy;

    bool filterByCellType(int ID);
    void considerNClosest(int N_considered);
    void considerClusters(int N_considered);
    void poweredWeighting(int exponent);
    void logWeighting(double log_a, double log_b);

  public:
    SensorHitMap(int _label);
    ~SensorHitMap();
    int label() {return _label;};
    void setLabZ(double z_cm, double X0);
    void setCenterHitPosition(double x, double y, double x_err, double y_err);
    void setParticleEnergy(double e);
    void setAlignmentParameters(double d_alpha, double d_beta, double d_gamma, double d_x0, double d_y0, double d_z0);
    void setResidualResolution(double r);
    void setSensorSize(int s);
    void setPedestalThreshold(double t);
    //reduces the information from the Rechit towards what is necessary for the impact point calculation
    void addHit(HGCalTBRecHit Rechit, double ADC_per_MIP);
    void registerClusterHit(HGCalTBDetId hit, int N_considered);
    std::pair<int, double> subtractCM();  //returns the sum of Common mode noise and the number of cells that enter the calculation
    void calculateCenterPosition(ConsiderationMethod considerationMethod, WeightingMethod weightingMethod);
    double getTotalEnergy();
    double getTotalClusterEnergy(int N_considered);
    double getTotalWeight();
    double getLabZ();
    double getX0();
    double getParticleEnergy();
    double getIntrinsicHitZPosition();
    double getResidualResolution();
    std::pair<double, double> getHitPosition(); //returns central hit in layer's own frame
    std::pair<double, double> getLabHitPosition();  //returns central hit in lab frame
    std::pair<double, double> getHitPositionError(); //calculated via RMS
    std::pair<double, double> getCenterOfClosestCell(std::pair<double, double> X_ref);
    
    //debug
    void printHits();
};

#endif
