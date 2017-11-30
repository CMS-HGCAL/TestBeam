#include "HGCal/DataFormats/interface/HGCalTBCluster.h"
#include "HGCal/DataFormats/interface/HGCalTBDetId.h"
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
       <<", E="<<cluster.energy();
      //       <<", Elow="<<cluster.energyLow();
      //<<", Position="<<cluster.position().x()<<","<<cluster.position().y()<<","<<cluster.position().z();
    // if( cluster.correctedEnergy() != -1.0 ) {
    //   out << ", E_corr="<<cluster.correctedEnergy();
    // }
    out<<", nhits="<<cluster.hitsAndFractions().size()<<std::endl;
    for( std::vector< std::pair<DetId,float> >::const_iterator it=cluster.hitsAndFractions().begin(); it!=cluster.hitsAndFractions().end(); ++it ) {
      out<<"("<<HGCalTBDetId((*it).first)<<", "<<(*it).second*cluster.energy()<<")\t";
    }
    return out;
  }

}
