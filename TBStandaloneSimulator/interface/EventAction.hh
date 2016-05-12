#ifndef EventAction_h
#define EventAction_h 1



#include "G4ThreeVector.hh"
#include "G4UserEventAction.hh"
#include "globals.hh"

#include "TFile.h"
#include "TTree.h"

#include "HGCal/TBStandaloneSimulator/interface/SamplingSection.hh"
#include "HGCal/TBStandaloneSimulator/interface/G4SiHit.hh"
#include "HGCal/TBStandaloneSimulator/interface/HGCSSEvent.hh"
#include "HGCal/TBStandaloneSimulator/interface/HGCSSSamplingSection.hh"
#include "HGCal/TBStandaloneSimulator/interface/HGCSSSimHit.hh"
#include "HGCal/TBStandaloneSimulator/interface/HGCSSGenParticle.hh"
#include "HGCal/TBStandaloneSimulator/interface/HGCSSGeometryConversion.hh"
#include "HGCal/TBStandaloneSimulator/interface/HGCSSTrackSegment.h"

#include <vector>
#include <map>
#include "fstream"

class RunAction;
class EventActionMessenger;

class EventAction : public G4UserEventAction
{
public:
	EventAction();
	virtual ~EventAction();
	void BeginOfEventAction(const G4Event*);
	void EndOfEventAction(const G4Event*);

	void Store(G4int trackId, G4int pdgId, G4double globalTime,
	           const G4ThreeVector& p1, const G4ThreeVector& p2,
	           const G4ThreeVector& p);

	void Detect(G4double edep, G4double stepl, G4double globalTime, G4int pdgId,
	            G4VPhysicalVolume *volume, const G4ThreeVector & position,
	            G4int trackID, G4int parentID,
	            const HGCSSGenParticle & genPart);

	void SetPrintModulo(G4int    val)
	{
		printModulo = val;
	};
	void Add( std::vector<SamplingSection> *newDetector )
	{
		detector_ = newDetector;
	}

	bool isFirstVolume(const std::string volname) const;

private:
	RunAction*  runAct;
	std::vector<SamplingSection> *detector_;
	G4int     evtNb_, printModulo;

	bool saveTracks;

	HGCSSGeometryConversion* geomConv_;

	TFile *outF_;
	TTree *tree_;
	HGCSSEvent event_;
	HGCSSSamplingSectionVec ssvec_;
	HGCSSSimHitVec hitvec_;
	HGCSSGenParticleVec genvec_;
	HGCSSTrackSegmentVec trkvec_;
	EventActionMessenger*  eventMessenger;
};

#endif


