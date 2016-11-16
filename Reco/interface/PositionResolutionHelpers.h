#ifndef HGCAL_RECO_PARTICLETRACK_H
#define HGCAL_RECO_PARTICLETRACK_H

#include <iostream>
#include <utility>
#include <cmath>
//#include "FWCore/Framework/interface/Frameworkfwd.h"
//#include "FWCore/Framework/interface/Event.h"
//#include "FWCore/ServiceRegistry/interface/Service.h"
#include "HGCal/DataFormats/interface/HGCalTBRecHitCollections.h"
//#include "HGCal/DataFormats/interface/HGCalTBDetId.h"
#include "HGCal/DataFormats/interface/HGCalTBRecHit.h"
#include "HGCal/Geometry/interface/HGCalTBCellVertices.h"
//#include "HGCal/Geometry/interface/HGCalTBTopology.h"
//#include "HGCal/CondObjects/interface/HGCalElectronicsMap.h"
//#include "HGCal/DataFormats/interface/HGCalTBElectronicsId.h"
//#include "HGCal/Geometry/interface/HGCalTBGeometryParameters.h"


//
// class declaration
//

enum WeightingMethod {
  SQUAREDWEIGHTING,
  LINEARWEIGHTING, 
  DEFAULT
};

class Particle_Track{
  public:
    Particle_Track(edm::Handle<HGCalTBRecHitCollection> Rechits);//, HGCalElectronicsMap& emap);
    ~Particle_Track();
  private:
    edm::Handle<HGCalTBRecHitCollection> Rechits_;
};


struct HitTriple {
  double x; double y; 
  double I;
  int ID;   //the ID corresponds to the cell ID, it is necessary for the pedestal subtraction
};

class SensorHitMap {
  private:
    std::pair<double, double> centralHitPoint;
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
    void setZ(double z);
    void setSensorSize(int s);
    //reduces the information from the Rechit towards what is necessary for the impact point calculation
    void addHit(HGCalTBRecHit Rechit);
    void subtractPedestals();
    std::pair<double, double> calculateCenterPosition(WeightingMethod method);
    std::pair<double, double> getCenterPosition();

    //debug
    void printHits();
};



#endif
