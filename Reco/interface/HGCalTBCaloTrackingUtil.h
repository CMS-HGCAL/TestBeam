#ifndef HGCALTBCALOTRACKINGUTIL_HH
#define HGCALTBCALOTRACKINGUTIL_HH

#include <iostream>
#include <vector>
#include <cmath>
#include <Eigen/Dense>

#include "DataFormats/Math/interface/Vector3D.h"
#include "DataFormats/Math/interface/Point3D.h"
#include "DataFormats/Math/interface/Error.h"

#include "HGCal/DataFormats/interface/HGCalTBClusterCollection.h"
#include "HGCal/DataFormats/interface/HGCalTBRecHitCollections.h"
#include "HGCal/DataFormats/interface/HGCalTBCaloTrackCollection.h"

#include "HGCal/Reco/interface/Distance.h"

namespace reco{
  template <typename T>
    class LeastSquare
    {
    public : 
      LeastSquare(){;}
      ~LeastSquare(){;}
    
      void run(T t, std::vector<float> &par, std::vector<float> &parerror)
      { 
	par.clear();
	parerror.clear();
	float xsum = 0.0;
	float ysum = 0.0;
	float zsum = 0.0;
	float zzsum = 0.0;
	float xzsum = 0.0;
	float yzsum = 0.0;
	for( auto p : t ){  
	  xsum = xsum + p.x();
	  ysum = ysum + p.y();
	  zsum = zsum + p.z();
	  xzsum = xzsum + p.x()*p.z();
	  yzsum = yzsum + p.y()*p.z();
	  zzsum = zzsum + p.z()*p.z();
	}

	par.push_back( (zzsum*xsum-xzsum*zsum)/(t.size()*zzsum-zsum*zsum) );
	par.push_back( (xzsum*t.size()-xsum*zsum)/(t.size()*zzsum-zsum*zsum) );
	par.push_back( (zzsum*ysum-yzsum*zsum)/(t.size()*zzsum-zsum*zsum) );
	par.push_back( (yzsum*t.size()-ysum*zsum)/(t.size()*zzsum-zsum*zsum) );
  
	parerror.push_back( std::sqrt( zzsum/(t.size()*zzsum-zsum*zsum) ) );
	parerror.push_back( std::sqrt( t.size()/(t.size()*zzsum-zsum*zsum) ) );
	parerror.push_back( std::sqrt( zzsum/(t.size()*zzsum-zsum*zsum) ) );
	parerror.push_back( std::sqrt( t.size()/(t.size()*zzsum-zsum*zsum) ) );
      }

      float chi2(T t, std::vector<float> &par)
      {
	float _chi2=0.0;
	for( auto p : t ){
	  Point p_track(  par[0]+par[1]*p.z(), par[2]+par[3]*p.z(), p.z() );
	  _chi2 += (p_track-p.position()).Mag2()*12;
	}
	return _chi2/t.size();
      }
    };

  template <typename T>
    class WeightedLeastSquare
    {
    public : 
      WeightedLeastSquare(){;}
      ~WeightedLeastSquare(){;}
    
      void run(T t, std::vector<float> &par, std::vector<float> &parerror)
      { 
	par.clear();
	parerror.clear();
	float xsum = 0.0;
	float ysum = 0.0;
	float zsum = 0.0;
	float zzsum = 0.0;
	float xzsum = 0.0;
	float yzsum = 0.0;
	float esum = 0.0;
	for( auto p : t ){  
	  xsum = xsum + p.x()*p.energy();
	  ysum = ysum + p.y()*p.energy();
	  zsum = zsum + p.z()*p.energy();
	  xzsum = xzsum + p.x()*p.z()*p.energy();
	  yzsum = yzsum + p.y()*p.z()*p.energy();
	  zzsum = zzsum + p.z()*p.z()*p.energy();
	  esum = esum + p.energy();
	}

	par.push_back( (zzsum*xsum-xzsum*zsum)/(esum*zzsum-zsum*zsum) );
	par.push_back( (xzsum*esum-xsum*zsum)/(esum*zzsum-zsum*zsum) );
	par.push_back( (zzsum*ysum-yzsum*zsum)/(esum*zzsum-zsum*zsum) );
	par.push_back( (yzsum*esum-ysum*zsum)/(esum*zzsum-zsum*zsum) );
  
	parerror.push_back( std::sqrt( zzsum/(esum*zzsum-zsum*zsum) ) );
	parerror.push_back( std::sqrt( esum/(esum*zzsum-zsum*zsum) ) );
	parerror.push_back( std::sqrt( zzsum/(esum*zzsum-zsum*zsum) ) );
	parerror.push_back( std::sqrt( esum/(esum*zzsum-zsum*zsum) ) );
      }

      float chi2(T t, std::vector<float> &par)
      {
	float _chi2=0.0;
	for( auto p : t ){
	  Point p_track(  par[0]+par[1]*p.z(), par[2]+par[3]*p.z(), p.z() );
	  _chi2 += (p_track-p.position()).Mag2()*12;
	}
	return _chi2;
      }
    };

  template <typename T>
    class PrincipalComponentAnalysis
    {
    public : 
      PrincipalComponentAnalysis(){;}
      ~PrincipalComponentAnalysis(){;}
      
      Eigen::MatrixXd covarianceMatrix() const {return _covarianceMatrix;};
      Eigen::MatrixXd eigenVectors() const {return _eigenVectors;}
      std::vector<double> &eigenValues() const {return _eigenValues; }
    private :
      Eigen::MatrixXd _covarianceMatrix;
      Eigen::MatrixXd _eigenVectors;
      std::vector<double> _eigenValues;
    public:
      void runPCA(T t)
      {
	covarianceMatrix = Eigen::MatrixXd(t.size(),3);
	float x,y,z; x=y=z=0.0;
	for( auto p : t ){
	  x+=p.x();
	  y+=p.y();
	  z+=p.z();
	}
	math::XYZPoint mean=math::XYZPoint(x,y,z)/t.size();
	int count=0;
	for( auto p : t ){
	  covarianceMatrix(count,0) = p.x()-mean.x();
	  covarianceMatrix(count,1) = p.y()-mean.y();
	  covarianceMatrix(count,2) = p.z()-mean.z();
	  count++;
	}    
	Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> eigensolver( _covarianceMatrix.transpose()*_covarianceMatrix );
	_eigenVectors=eigensolver.eigenvectors();
	for(unsigned int i=0; i<eigensolver.eigenvalues().size(); i++){
	  _eigenValues.push_back(eigensolver.eigenvalues()[i]);
	}
      }
    };

  class TrackCleaner
  {
  public: 
    void clean( HGCalTBRecHitCollection col, HGCalTBRecHitCollection &cleancol, reco::HGCalTBCaloTrack &track, double maxDistance )
    {
      Distance<HGCalTBRecHit,reco::HGCalTBCaloTrack> dist;
      for( auto p : col )
	if( dist.distance( p,track )<maxDistance )
	  cleancol.push_back(p);
    }
    void clean( HGCalTBClusterCollection col, HGCalTBClusterCollection &cleancol, reco::HGCalTBCaloTrack &track, double maxDistance )
    {
      Distance<HGCalTBCluster,reco::HGCalTBCaloTrack> dist;
      for( auto p : col )
	if( dist.distance( p,track )<maxDistance )
	  cleancol.push_back(p);
    }
  };
};
#endif
