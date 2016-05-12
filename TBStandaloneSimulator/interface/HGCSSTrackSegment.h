#ifndef HGCSSTRACKSEGMENT_H
#define HGCSSTRACKSEGMENT_H

#include "Rtypes.h"
#include "G4ThreeVector.hh"

class HGCSSTrackSegment
{
public:
	typedef int key_type;

	HGCSSTrackSegment()
		: track(0),
		  pdgid(0),
		  time(0),
		  x1(0),
		  y1(0),
		  z1(0),
		  x2(0),
		  y2(0),
		  z2(0),
		  px(0),
		  py(0),
		  pz(0)
	{
	};

	HGCSSTrackSegment(int track_, int pdgid_, double time_,
	                  const G4ThreeVector& p1,
	                  const G4ThreeVector& p2,
	                  const G4ThreeVector& p)
		: track(track_),
		  pdgid(pdgid_),
		  time(time_),
		  x1(p1[0]),
		  y1(p1[1]),
		  z1(p1[2]),
		  x2(p2[0]),
		  y2(p2[1]),
		  z2(p2[2]),
		  px(p[0]),
		  py(p[1]),
		  pz(p[2])
	{
	}

	int id() const
	{
		return track;
	}

	virtual ~HGCSSTrackSegment() {}

private:
	int track, pdgid;
	double time;
	double x1, y1, z1;
	double x2, y2, z2;
	double px, py, pz;

	ClassDef(HGCSSTrackSegment, 1);
};


typedef std::vector<HGCSSTrackSegment> HGCSSTrackSegmentVec;

#endif
