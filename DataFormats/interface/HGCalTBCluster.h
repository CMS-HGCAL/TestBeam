#ifndef DATAFORMATS_HGCALTBCLUSTER_H
#define DATAFORMATS_HGCALTBCLUSTER_H 1

#include "DataFormats/CaloRecHit/interface/CaloCluster.h"

/** \class HGCalTBCluster
 *
 * \author Arnaud Steen
 * 2D hgcal cluster class
 */

namespace reco{

  class HGCalTBCluster : public reco::CaloCluster
    {
    public:
      HGCalTBCluster();
      HGCalTBCluster(int layer, float energy, float energyLow, float energyHigh );
	
      int   layer() const { return _layer; }
      float energyLow() const {	return _energyLow; }
      float energyHigh() const { return _energyHigh; }
      float recHitEnergyHigh( int i) const{ return hitsAndFractions_.at(i).second*_energyHigh; }
      float recHitEnergyLow( int i) const{ return hitsAndFractions_.at(i).second*_energyLow; }

      void setLayer(int val){ _layer=val; }
      void setEnergyLow(float val){ _energyLow=val; }
      void setEnergyHigh(float val){ _energyHigh=val; }
      
      
    protected : 
      int _layer;
      float _energyLow;
      float _energyHigh;
    };

  std::ostream& operator<<(std::ostream& s, const HGCalTBCluster& cluster);
}

#endif
