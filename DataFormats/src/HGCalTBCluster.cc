#include "HGCal/DataFormats/interface/HGCalTBCluster.h"
//#include "DataFormats/Math/interface/Point3D.h"
#include <iostream>
namespace reco
{

  HGCalTBCluster::HGCalTBCluster() : reco::CaloCluster()
  {

  }

  HGCalTBCluster::HGCalTBCluster(int layer, float energy, float energyLow, float energyHigh) :
    reco::CaloCluster(),
    _layer(layer),
    _energyLow(energyLow),
    _energyHigh(energyHigh)
  {
		this->setEnergy(energy);
  }

  std::ostream& operator<<(std::ostream& out, const HGCalTBCluster& cluster)
  {
 
    if(!out) return out;

    out<<"CaloCluster , algoID="<<cluster.algoID()
       <<", Layer="<<cluster.layer()    
       <<", Ehigh="<<cluster.energyHigh()
       <<", Elow="<<cluster.energyLow();
      //<<", Position="<<cluster.position().x()<<","<<cluster.position().y()<<","<<cluster.position().z();
    if( cluster.correctedEnergy() != -1.0 ) {
      out << ", E_corr="<<cluster.correctedEnergy();
    }
    out<<", nhits="<<cluster.hitsAndFractions().size()<<std::endl;
    for(unsigned i=0; i<cluster.hitsAndFractions().size(); i++ ) {
      out<<""<<cluster.printHitAndFraction(i)<<", ";
    }
    return out;
  }

}
