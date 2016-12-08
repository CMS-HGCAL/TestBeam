#ifndef DISTANCE_HH
#define DISTANCE_HH

#include <iostream>
#include <cmath>

#include "HGCal/DataFormats/interface/HGCalTBCluster.h"
#include "HGCal/DataFormats/interface/HGCalTBRecHit.h"
#include "HGCal/DataFormats/interface/HGCalTBCaloTrack.h"

#include "DataFormats/Math/interface/Point3D.h"
#include "DataFormats/Math/interface/Vector3D.h"

template <typename T, typename S>
  class Distance
{
 public : 
  Distance(){;}
  ~Distance(){;}
    
  float distance(T t, S s){ return std::sqrt( (t.position()-s.position()).mag2() ); }
};

template <typename T>
class Distance<T, math::XYZPoint >
{
 public : 
  Distance(){;}
  ~Distance(){;}
    
  float distance(T t, math::XYZPoint p){ return std::sqrt( (t.position()-p).mag2() ); }
};

template <typename S> 
class Distance<S, reco::HGCalTBCaloTrack>
{ 
 public :  
  Distance(){;} 
  ~Distance(){;} 
      
  float distance(S s, reco::HGCalTBCaloTrack t)
  {
    /*
      point H(x,y,z) (S s)
      track T, orientation vector u(tx,ty,tz)
      d(H,T) = || vec(BH) * u || / || u || where B is a point from the track
    */

    math::XYZPoint H=s.position();

    //B : a track point 
    math::XYZPoint B=t.vertex();

    math::XYZVector u=t.momentum();
    math::XYZVector v=B-H;
	
    return std::sqrt( (v.Cross(u)).mag2()/u.mag2() );
  }

  float distanceInLayer(S s,reco::HGCalTBCaloTrack t)
  {
    math::XYZPoint impact=t.expectedTrackProjection(s.z());
    return std::sqrt( (impact-s.position()).mag2() );
  }
}; 

#endif

