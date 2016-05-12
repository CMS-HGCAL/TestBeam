#include "HGCal/TBStandaloneSimulator/interface/SteppingAction.hh"
#include "HGCal/TBStandaloneSimulator/interface/DetectorConstruction.hh"
#include "HGCal/TBStandaloneSimulator/interface/EventAction.hh"
#include "HGCal/TBStandaloneSimulator/interface/HGCSSGenParticle.hh"
#include "G4Step.hh"
#include "G4RunManager.hh"

//
SteppingAction::SteppingAction()
{
	eventAction_ = (EventAction*)G4RunManager::GetRunManager()->
	               GetUserEventAction();
	eventAction_->Add(  ((DetectorConstruction*)G4RunManager::GetRunManager()->
	                     GetUserDetectorConstruction())->getStructure() );
	saturationEngine = new G4EmSaturation();
	// A hack to avoid compiler warning
	int hack = CLHEP::HepRandomGenActive;
	timeLimit_ = hack;
	timeLimit_ = 100;//ns
}

//
SteppingAction::~SteppingAction()
{ }

//
void SteppingAction::UserSteppingAction(const G4Step* aStep)
{
	// get PreStepPoint
	const G4StepPoint *thePreStepPoint = aStep->GetPreStepPoint();
	const G4StepPoint *thePostStepPoint = aStep->GetPostStepPoint();
	// get TouchableHandle
	const G4Track* lTrack = aStep->GetTrack();
	G4int trackID = lTrack->GetTrackID();
	G4int parentID = lTrack->GetParentID();

	G4VPhysicalVolume* volume = thePreStepPoint->GetPhysicalVolume();
	std::string thePrePVname("null");
	if(volume == 0) {
	} else {
		thePrePVname = volume->GetName();
	}
	G4VPhysicalVolume* postvolume = thePostStepPoint->GetPhysicalVolume();
	std::string thePostPVname("null");
	if(postvolume == 0) {
	} else {
		thePostPVname = postvolume->GetName();
	}

	G4double edep = aStep->GetTotalEnergyDeposit();

	//correct with Birk's law for scintillator material
	if (volume->GetName().find("Scint") != volume->GetName().npos) {
		G4double attEdep
		    = saturationEngine->
		      VisibleEnergyDeposition(lTrack->GetDefinition(),
		                              lTrack->GetMaterialCutsCouple(),
		                              aStep->GetStepLength(), edep, 0.);
		edep = attEdep;
	}

	G4double stepl = 0.;
	if (lTrack->GetDefinition()->GetPDGCharge() != 0.)
		stepl = aStep->GetStepLength();

	G4int pdgId = lTrack->GetDefinition()->GetPDGEncoding();
	G4double globalTime = lTrack->GetGlobalTime();

	// store track segment
	const G4ThreeVector& position = thePreStepPoint->GetPosition();
	const G4ThreeVector& postposition = thePostStepPoint->GetPosition();
	const G4ThreeVector& p = lTrack->GetMomentum();
	eventAction_->Store(trackID, pdgId, globalTime, p, position, postposition);

	HGCSSGenParticle genPart;
	if (globalTime < timeLimit_ &&
	        thePrePVname == "Wphys"
	        && eventAction_->isFirstVolume(thePostPVname)) {
		const G4ThreeVector&  postposition = thePostStepPoint->GetPosition();
		const G4ThreeVector&  p  = lTrack->GetMomentum();
		G4ParticleDefinition* pd = lTrack->GetDefinition();
		genPart.setPosition(postposition[0],
		                    postposition[1],
		                    postposition[2]);
		genPart.setMomentum(p[0], p[1], p[2]);
		genPart.mass(pd->GetPDGMass());
		genPart.time(globalTime);
		genPart.pdgid(pdgId);
		genPart.charge(pd->GetPDGCharge());
		genPart.trackID(trackID);
	}

	eventAction_->Detect(edep,
	                     stepl,
	                     globalTime,
	                     pdgId,
	                     volume,
	                     position,
	                     trackID,
	                     parentID,
	                     genPart);
}
