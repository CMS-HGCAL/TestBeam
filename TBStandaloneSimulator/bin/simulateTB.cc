// ----------------------------------------------------------------------------
// read geometry from a config file 04-28-2016 HBP
// ----------------------------------------------------------------------------
#include "G4RunManager.hh"
#include "G4UImanager.hh"
#include "G4VisExecutive.hh"
#include "G4UIExecutive.hh"

#include "HGCal/TBStandaloneSimulator/interface/DetectorConstruction.hh"
#include "HGCal/TBStandaloneSimulator/interface/PhysicsList.hh"
#include "HGCal/TBStandaloneSimulator/interface/PrimaryGeneratorAction.hh"
#include "HGCal/TBStandaloneSimulator/interface/RunAction.hh"
#include "HGCal/TBStandaloneSimulator/interface/EventAction.hh"
#include "HGCal/TBStandaloneSimulator/interface/SteppingAction.hh"
#include "HGCal/TBStandaloneSimulator/interface/SteppingVerbose.hh"
#include "HGCal/TBStandaloneSimulator/interface/TBGeometry.h"

#include "TSystem.h"
#include "TString.h"

int main(int argc, char** argv)
{
	// A hack to avoid compiler warning
	int hack = CLHEP::HepRandomGenActive;
	hack = argc;
	argc = hack;

	if ( argc < 3 ) {
		std::cout << "Usage: " << std::endl
		          << "  simulateTB \x1b[1;31;48m<geometry file>\x1b[0m"
		          << " \x1b[1;32;48m<macro file>\x1b[0m"
		          << std::endl;
		exit(0);
	}

	G4String geomFile  = gSystem->ExpandPathName(argv[1]);
	G4String macroFile = gSystem->ExpandPathName(argv[2]);
	TBConfig config;
	config.macro = macroFile;
	config.geometry = geomFile;
	if ( argc > 3 )
		config.savetracks = atoi(argv[3]);
	else
		config.savetracks = 0;
	G4cout << config << G4endl;

	TBGeometry geometry(geomFile);

	// Choose the Random engine
	CLHEP::HepRandom::setTheEngine(new CLHEP::RanecuEngine);

	// User Verbose output class
	G4VSteppingVerbose::SetInstance(new SteppingVerbose);

	// Construct the default run manager
	G4RunManager * runManager = new G4RunManager;

	// Set mandatory initialization classes

	runManager->SetUserInitialization(new DetectorConstruction(geometry,
	                                  config));
	runManager->SetUserInitialization(new PhysicsList);

	// Set user action classes
	int model = geometry.model();
	double eta = 0;
	runManager->SetUserAction(new PrimaryGeneratorAction(model, eta));
	runManager->SetUserAction(new RunAction);
	runManager->SetUserAction(new EventAction);
	runManager->SetUserAction(new SteppingAction);

	// Initialize G4 kernel
	runManager->Initialize();

	// Initialize visualization
	G4VisManager* visManager = new G4VisExecutive;
	visManager->Initialize();

	// Start User Interface manager
	G4UImanager* UImanager = G4UImanager::GetUIpointer();
	G4String command = "/control/execute " + macroFile;
	UImanager->ApplyCommand(command);

	delete visManager;
	delete runManager;
	return 0;
}

