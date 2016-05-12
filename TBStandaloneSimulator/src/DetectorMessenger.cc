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

#include "HGCal/TBStandaloneSimulator/interface/DetectorMessenger.hh"
#include "HGCal/TBStandaloneSimulator/interface/DetectorConstruction.hh"

#include "G4UIdirectory.hh"
#include "G4UIcmdWithAString.hh"
#include "G4UIcmdWithAnInteger.hh"
#include "G4UIcmdWithADoubleAndUnit.hh"
#include "G4UIcmdWithoutParameter.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

DetectorMessenger::DetectorMessenger(
    DetectorConstruction* Det)
	: Detector(Det)
{
	N03Dir = new G4UIdirectory("/N03/");
	N03Dir->SetGuidance("UI commands of this example");

	detDir = new G4UIdirectory("/N03/det/");
	detDir->SetGuidance("detector control");

	MagFieldCmd = new G4UIcmdWithADoubleAndUnit("/N03/det/setField", this);
	MagFieldCmd->SetGuidance("Define magnetic field.");
	MagFieldCmd->SetGuidance("Magnetic field will be in Z direction.");
	MagFieldCmd->SetParameterName("Bz", false);
	MagFieldCmd->SetUnitCategory("Magnetic flux density");
	MagFieldCmd->AvailableForStates(G4State_PreInit, G4State_Idle);

	SetModelCmd = new G4UIcmdWithAnInteger("/N03/det/setModel", this);
	SetModelCmd->SetGuidance("Define detector model");
	SetModelCmd->SetParameterName("Model", true);
	SetModelCmd->SetDefaultValue(0);
	SetModelCmd->AvailableForStates(G4State_PreInit, G4State_Idle);

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

DetectorMessenger::~DetectorMessenger()
{
	delete MagFieldCmd;
	delete SetModelCmd;
	delete detDir;
	delete N03Dir;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void DetectorMessenger::SetNewValue(G4UIcommand* command, G4String newValue)
{
	if( command == MagFieldCmd ) {
		Detector->SetMagField(MagFieldCmd->GetNewDoubleValue(newValue));
	}
	if (command == SetModelCmd ) {
		Detector->SetDetModel(SetModelCmd->GetNewIntValue(newValue));
	}

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
