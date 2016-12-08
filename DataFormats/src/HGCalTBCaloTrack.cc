#include "HGCal/DataFormats/interface/HGCalTBCaloTrack.h"
#include <iostream>
namespace reco
{

  HGCalTBCaloTrack::HGCalTBCaloTrack(double _chi2, double _ndof, const Point &_vertex,
				     const Vector &_momentum,/* const CovarianceMatrix &_cov, */
				     std::vector<HGCalTBDetId> &_decIds) 
  {
    if( _chi2!=_chi2 )
      return;

    isNull_ = false;
    chi2_ = _chi2;
    ndof_ = _ndof; 
    
    vertex_ = _vertex;
    momentum_ = _momentum; 
    //cov_ = _cov; 

    detIds_ = _decIds;
  }

  std::ostream& operator<<(std::ostream& out, const HGCalTBCaloTrack& caloTrack)
  {
 
    if(!out) return out;

    out<<"CaloCaloTrack , chi2="<<caloTrack.chi2()
       <<",\t ndof="<<caloTrack.ndof()    
       <<",\t first point="<<caloTrack.vertex().x()<<";"<<caloTrack.vertex().y()<<";"<<caloTrack.vertex().z()
       <<",\t momentum="<<caloTrack.momentum().x()<<";"<<caloTrack.momentum().y()<<";"<<caloTrack.momentum().z()
      /*<<",\t nhits="<<caloTrack.getDetIds().size() << */ << std::endl;
    return out;
  }

}
