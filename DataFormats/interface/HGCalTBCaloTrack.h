#ifndef DATAFORMATS_HGCALTBCALOTRACK_H
#define DATAFORMATS_HGCALTBCALOTRACK_H 1

#include "HGCal/DataFormats/interface/HGCalTBDetId.h"
#include "DataFormats/Math/interface/Vector3D.h"
#include "DataFormats/Math/interface/Point3D.h"
#include "DataFormats/Math/interface/Error.h"

/** \class HGCalTBCaloTrack
 *
 * \author Arnaud Steen
 * hgcal caloTrack class
 */

typedef math::XYZVector Vector;
typedef math::XYZPoint Point;

namespace reco{

  class HGCalTBCaloTrack 
    {
    public:
      //typedef math::Error<dimension>::type CovarianceMatrix;
      HGCalTBCaloTrack( ){ isNull_=true;}
      HGCalTBCaloTrack( double _chi2, double _ndof, const Point &_vertex,
			const Vector &_momentum,/* const CovarianceMatrix &_cov, */
			std::vector<HGCalTBDetId> &_detIds);
      
      bool isNull() const { return isNull_; }
      double chi2() const { return chi2_; }
      double ndof() const { return ndof_; }
      double normalisedChi2() const { return chi2_/ndof_; }

      HGCalTBDetId getDetId(int i) const { return detIds_.at(i); }
      std::vector<HGCalTBDetId> &getDetIds() { return detIds_; }

      Point vertex() const { return vertex_; }
      Vector momentum() const { return momentum_; }
      //CovarianceMatrix covarianceMatrix() const { return cov_; }; 

      Point expectedTrackProjection(float z) const { return Point( vertex_.x()+momentum_.x()*z,vertex_.y()+momentum_.y()*z,z); }

    private:
      bool isNull_;
      double chi2_;
      double ndof_; 
      
      Point vertex_;
      Vector momentum_; 
      //CovarianceMatrix cov_; 
      
      std::vector<HGCalTBDetId> detIds_;
    };

  std::ostream& operator<<(std::ostream& s, const HGCalTBCaloTrack& caloTrack);
}

#endif
