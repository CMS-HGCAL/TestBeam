//
// ********************************************************************
// * License and Disclaimer                                           *
// *                                                                  *
// * The  Geant4 software  is  copyright of the Copyright Holders  of *
// * the Geant4 Collaboration.  It is provided  under  the terms  and *
// * conditions of the Geant4 Software License,  included in the file *
// * LICENSE and available at  http://cern.ch/geant4/license .  These *
// * include a list of copyright holders.                             *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.  Please see the license in the file  LICENSE  and URL above *
// * for the full disclaimer and the limitation of liability.         *
// *                                                                  *
// * This  code  implementation is the result of  the  scientific and *
// * technical work of the GEANT4 collaboration.                      *
// * By using,  copying,  modifying or  distributing the software (or *
// * any work based  on the software)  you  agree  to acknowledge its *
// * use  in  resulting  scientific  publications,  and indicate your *
// * acceptance of all terms of the Geant4 Software license.          *
// ********************************************************************
//
//
// $Id$
//
//

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include "HGCal/TBStandaloneSimulator/interface/PrimaryGeneratorMessenger.hh"
#include "HGCal/TBStandaloneSimulator/interface/PrimaryGeneratorAction.hh"

#include "G4UIdirectory.hh"
#include "G4UIcmdWithAString.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

PrimaryGeneratorMessenger::PrimaryGeneratorMessenger(
    PrimaryGeneratorAction* Gun)
	: Action(Gun)
{
	dir = new G4UIdirectory("/generator/");
	dir-> SetGuidance("Control commands for primary generator");

	select = new G4UIcmdWithAString("/generator/select", this);
	select-> SetGuidance("select generator type");
	select-> SetParameterName("generator_type", false, false);
	select-> SetCandidates("particleGun pythia hepmcAscii");
	select-> SetDefaultValue("particleGun");

	RndmCmd = new G4UIcmdWithAString("/generator/particleGun/rndm", this);
	RndmCmd->SetGuidance("Shoot randomly the incident particle.");
	RndmCmd->SetGuidance("  Choice : on(default), off");
	RndmCmd->SetParameterName("choice", true);
	RndmCmd->SetDefaultValue("on");
	RndmCmd->SetCandidates("on off");
	RndmCmd->AvailableForStates(G4State_PreInit, G4State_Idle);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

PrimaryGeneratorMessenger::~PrimaryGeneratorMessenger()
{
	delete select;
	delete RndmCmd;
	delete dir;
	//delete gunDir;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void PrimaryGeneratorMessenger::SetNewValue(
    G4UIcommand* command, G4String newValue)
{
	if ( command == select) {
		Action->SetGenerator(newValue);
		G4cout << "current generator type: "
		       << Action-> GetGeneratorName() << G4endl;
	} else if( command == RndmCmd ) {
		Action->SetRndmFlag(newValue);
	} else {}
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......


//////////////////////////////////////////////////////////////////////////////
G4String PrimaryGeneratorMessenger::GetCurrentValue(G4UIcommand* command)
//////////////////////////////////////////////////////////////////////////////
{
	G4String cv, st;
	if (command == select) {
		cv = Action-> GetGeneratorName();
	}

	return cv;
}
