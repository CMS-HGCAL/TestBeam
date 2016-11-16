#ifndef HGCAL_RECO_PARTICLETRACK_H
#define HGCAL_RECO_PARTICLETRACK_H


//#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/Event.h"
//#include "FWCore/ServiceRegistry/interface/Service.h"
#include "HGCal/DataFormats/interface/HGCalTBRecHitCollections.h"
//#include "HGCal/DataFormats/interface/HGCalTBDetId.h"
#include "HGCal/DataFormats/interface/HGCalTBRecHit.h"
//#include "HGCal/Geometry/interface/HGCalTBCellVertices.h"
//#include "HGCal/Geometry/interface/HGCalTBTopology.h"
#include "HGCal/CondObjects/interface/HGCalElectronicsMap.h"
#include "HGCal/DataFormats/interface/HGCalTBElectronicsId.h"
//#include "HGCal/Geometry/interface/HGCalTBGeometryParameters.h"

using namespace std;
//
// class declaration
//

class Particle_Track{

  public:
    Particle_Track(edm::Handle<HGCalTBRecHitCollection> Rechits, HGCalElectronicsMap& emap);
    ~Particle_Track();

  private:
    edm::Handle<HGCalTBRecHitCollection> Rechits_;
    struct {
      HGCalElectronicsMap emap_;
    } essource_;

};

#endif
