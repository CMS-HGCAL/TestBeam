#include "HGCal/TBStandaloneSimulator/interface/EventAction.hh"
#include "HGCal/TBStandaloneSimulator/interface/RunAction.hh"
#include "HGCal/TBStandaloneSimulator/interface/EventActionMessenger.hh"
#include "HGCal/TBStandaloneSimulator/interface/DetectorConstruction.hh"
#include "HGCal/TBStandaloneSimulator/interface/HGCSSInfo.hh"
#include "HGCal/TBStandaloneSimulator/interface/HGCSSSimHit.hh"
#include "HGCal/TBStandaloneSimulator/interface/HGCSSTrackSegment.h"

#include "G4RunManager.hh"
#include "G4Event.hh"
#include "G4UnitsTable.hh"

#include <iomanip>

//
EventAction::EventAction()
{
	runAct = (RunAction*)G4RunManager::GetRunManager()->GetUserRunAction();
	eventMessenger = new EventActionMessenger(this);
	printModulo = 100;
	outF_ = TFile::Open("PFcal.root", "RECREATE");
	outF_->cd();

	DetectorConstruction* detector
	    = (DetectorConstruction*)(G4RunManager::GetRunManager()->
	                              GetUserDetectorConstruction());

	// cache save tracks flag

	saveTracks = detector->saveTracks();

	// A hack to avoid compiler warning
	int hack = CLHEP::HepRandomGenActive;
	int model = hack;

	model          = detector->getModel();
	int version    = detector->getVersion();
	double xysize  = detector->GetCalorSizeXY();
	double cellsize = detector->GetCellSize();

	//save some info
	HGCSSInfo* info = new HGCSSInfo();
	info->calorSizeXY(xysize);
	info->cellSize(cellsize);
	info->model(model);
	info->version(version);

	G4cout << " -- check Info: version = " << info->version()
	       << " model = " << info->model() << G4endl;
	outF_->WriteObjectAny(info, "HGCSSInfo", "Info");

	//honeycomb
	geomConv_ = new HGCSSGeometryConversion(info->model(),
	                                        xysize,
	                                        cellsize);
	geomConv_->initialiseHoneyComb(xysize, cellsize);

	//square map
	geomConv_->initialiseSquareMap(xysize, cellsize);

	tree_ = new TTree("HGCSSTree", "HGC Standalone simulation");
	tree_->Branch("HGCSSEvent", "HGCSSEvent", &event_);
	tree_->Branch("HGCSSSamplingSectionVec",
	              "std::vector<HGCSSSamplingSection>", &ssvec_);
	tree_->Branch("HGCSSSimHitVec", "std::vector<HGCSSSimHit>", &hitvec_);
	tree_->Branch("HGCSSGenParticleVec",
	              "std::vector<HGCSSGenParticle>", &genvec_);
	if ( saveTracks )
		tree_->Branch("HGCSSTrackSegmentVec",
		              "std::vector<HGCSSTrackSegment>", &trkvec_);
}

//
EventAction::~EventAction()
{
	outF_->cd();
	tree_->Write();
	outF_->Close();
	delete eventMessenger;
}

//
void EventAction::BeginOfEventAction(const G4Event* evt)
{
	evtNb_ = evt->GetEventID();
	if (evtNb_ % printModulo == 0) {
		G4cout << "\n---> Start event: " << evtNb_ << G4endl;
		CLHEP::HepRandom::showEngineStatus();
	}
}

//
void EventAction::Store(G4int trackId, G4int pdgId, G4double globalTime,
                        const G4ThreeVector& p1, const G4ThreeVector& p2,
                        const G4ThreeVector& p)
{
	if ( saveTracks )
		trkvec_.push_back(HGCSSTrackSegment(trackId, pdgId, globalTime,
		                                    p1, p2, p));
}

//
void EventAction::Detect(G4double edep,
                         G4double stepl,
                         G4double globalTime,
                         G4int pdgId,
                         G4VPhysicalVolume* volume,
                         const G4ThreeVector& position,
                         G4int trackID,
                         G4int parentID,
                         const HGCSSGenParticle& genPart)
{
	for(size_t i = 0; i < detector_->size(); i++)
		(*detector_)[i].add(edep,
		                    stepl,
		                    globalTime,
		                    pdgId,
		                    volume,
		                    position,
		                    trackID,
		                    parentID,
		                    i);
	if (genPart.isIncoming()) genvec_.push_back(genPart);
}

bool EventAction::isFirstVolume(const std::string volname) const
{
	if ( ! detector_ ) return  false;
	if ( detector_->size() == 0 ) return  false;

	// get the first sampling section.
	// loop over all of its volumes
	// and check if one of them matches the
	// given volume name
	SamplingSection& section = (*detector_)[0];
	if ( section.n_elements == 0 ) return false;

	bool found = false;
	for(size_t c = 0; c < section.ele_vol.size(); c++) {
		std::string name = std::string(section.ele_vol[c]->GetName());
		if ( name == volname ) {
			found = true;
			break;
		}
	}
	return found;
}

//
void EventAction::EndOfEventAction(const G4Event* g4evt)
{
	assert(g4evt != 0);
	bool debug(evtNb_ % printModulo == 0);

	hitvec_.clear();
	trkvec_.clear();
	genvec_.clear();
	ssvec_.clear();

	event_.eventNumber(evtNb_);

	//G4cout << " EndOfEventAction - Number of primary vertices: "
	// << g4evt->GetNumberOfPrimaryVertex() << G4endl;
	assert(g4evt->GetNumberOfPrimaryVertex() > 0);

	if ( debug )
		G4cout << " -- vtx pos x=" << g4evt->GetPrimaryVertex(0)->GetX0()
		       << " y=" << g4evt->GetPrimaryVertex(0)->GetY0()
		       << " z=" << g4evt->GetPrimaryVertex(0)->GetZ0()
		       << " t=" << g4evt->GetPrimaryVertex(0)->GetT0()
		       << G4endl;

	event_.vtx_x(g4evt->GetPrimaryVertex(0)->GetX0());
	event_.vtx_y(g4evt->GetPrimaryVertex(0)->GetY0());
	event_.vtx_z(g4evt->GetPrimaryVertex(0)->GetZ0());

	ssvec_.clear();
	ssvec_.reserve(detector_->size());

	for(size_t i = 0; i < detector_->size(); i++) {
		HGCSSSamplingSection lSec;

		// NOTE: get a reference, NOT a copy!
		SamplingSection& section = (*detector_)[i];

		lSec.volNb(i);
		lSec.volX0trans(section.getAbsorberX0());
		lSec.voldEdx(section.getAbsorberdEdx());
		lSec.volLambdatrans(section.getAbsorberLambda());
		lSec.absorberE(section.getAbsorbedEnergy());
		lSec.measuredE(section.getMeasuredEnergy(false));
		lSec.totalE(section.getTotalEnergy());
		lSec.gFrac(section.getPhotonFraction());
		lSec.eFrac(section.getElectronFraction());
		lSec.muFrac(section.getMuonFraction());
		lSec.neutronFrac(section.getNeutronFraction());
		lSec.hadFrac(section.getHadronicFraction());
		lSec.avgTime(section.getAverageTime());
		lSec.nSiHits(section.getTotalSensHits());
		ssvec_.push_back(lSec);

		if (evtNb_ == 1) std::cout << "if (layer==" << i << ") return "
			                           <<  lSec.voldEdx() << ";"
			                           << std::endl;
		bool is_scint = section.hasScintillator;
		for (unsigned idx(0); idx < section.n_sens_elements; ++idx) {
			std::map<unsigned, HGCSSSimHit> lHitMap;
			std::pair<std::map<unsigned, HGCSSSimHit>::iterator, bool> isInserted;

			if (i > 0) {
				// look at previous section, but check that it has
				// sensitive elements.
				SamplingSection& prevsection = (*detector_)[i - 1];
				if ( prevsection.SiHitVecSize() > 0 ) {
					const G4SiHitVec vhit = prevsection.getSiHitVec(idx);
					section.trackParticleHistory(idx, vhit);
				}
			}

			for (unsigned iSiHit(0);
			        iSiHit < section.getSiHitVec(idx).size();
			        ++iSiHit) {
				G4SiHit lSiHit = section.getSiHitVec(idx)[iSiHit];
				HGCSSSimHit lHit(lSiHit, section.sens_layer[idx], is_scint
				                 ? geomConv_->squareMap()
				                 : geomConv_->hexagonMap());

				isInserted = lHitMap.
				             insert(std::pair<unsigned, HGCSSSimHit>(lHit.cellid(), lHit));
				if (!isInserted.second) isInserted.first->second.Add(lSiHit);
			}
			std::map<unsigned, HGCSSSimHit>::iterator lIter = lHitMap.begin();
			hitvec_.reserve(hitvec_.size() + lHitMap.size());
			for (; lIter != lHitMap.end(); ++lIter) {
				(lIter->second).calculateTime();
				hitvec_.push_back(lIter->second);
			}

		}//loop on sensitive layers

		if(debug) {
			section.report( (i == 0) );
		}
		//if (i==0) G4cout << " ** evt " << evt->GetEventID() << G4endl;
		section.resetCounters();
	}
	if(debug) {
		G4cout << " -- Number of track segments    = " << trkvec_.size() << G4endl
		       << " -- Number of truth particles   = " << genvec_.size() << G4endl
		       << " -- Number of simhits           = " << hitvec_.size() << G4endl
		       << " -- Number of sampling sections = " << ssvec_.size()  << G4endl;
	}

	tree_->Fill();

	//G4cout << "END(EventAction::EndOfEventAction" << G4endl;
}
