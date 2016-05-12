#include <vector>
#include "HGCal/TBStandaloneSimulator/interface/HGCSSInfo.hh"
#include "HGCal/TBStandaloneSimulator/interface/HGCSSEvent.hh"
#include "HGCal/TBStandaloneSimulator/interface/HGCSSSamplingSection.hh"
#include "HGCal/TBStandaloneSimulator/interface/HGCSSSimHit.hh"
#include "HGCal/TBStandaloneSimulator/interface/HGCSSTrackSegment.h"
#include "HGCal/TBStandaloneSimulator/interface/HGCSSGenParticle.hh"
#include "HGCal/TBStandaloneSimulator/interface/HGCSSRecoHit.hh"
#include "HGCal/TBStandaloneSimulator/interface/HGCSSGeometryConversion.hh"
#include "HGCal/TBStandaloneSimulator/interface/HGCCellMap.h"
#include "DataFormats/Common/interface/Wrapper.h"

namespace HGCSS
{
struct dictionary {
	HGCSSSamplingSection _a1;
	std::vector<HGCSSSamplingSection> _v1;
	edm::Wrapper<std::vector<HGCSSSamplingSection> > _w1;

	HGCSSSimHit _a2;
	std::vector<HGCSSSimHit> _v2;
	edm::Wrapper<std::vector<HGCSSSimHit> > _w2;

	HGCSSGenParticle _a3;
	std::vector<HGCSSGenParticle> _v3;
	edm::Wrapper<std::vector<HGCSSGenParticle> > _w3;

	HGCSSRecoHit _a4;
	std::vector<HGCSSRecoHit> _v4;
	edm::Wrapper<std::vector<HGCSSRecoHit> > _w4;

	HGCSSTrackSegment _a5;
	std::vector<HGCSSTrackSegment> _v5;
	edm::Wrapper<std::vector<HGCSSTrackSegment> > _w5;
};
}

