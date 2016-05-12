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

#ifndef PrimaryGeneratorAction_h
#define PrimaryGeneratorAction_h 1

#include "G4VUserPrimaryGeneratorAction.hh"
#include "globals.hh"
#include <map>

class G4Event;
class DetectorConstruction;
class PrimaryGeneratorMessenger;
class G4VPrimaryGenerator;
class G4ParticleGun;
class HepMCG4AsciiReader;
class HepMCG4PythiaInterface;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

class PrimaryGeneratorAction : public G4VUserPrimaryGeneratorAction
{
public:
	PrimaryGeneratorAction(G4int mod = 0, double eta = 0);
	virtual ~PrimaryGeneratorAction();

	void GeneratePrimaries(G4Event*);
	void SetRndmFlag(G4String val)
	{
		rndmFlag = val;
	}

	void SetGenerator(G4VPrimaryGenerator* gen);
	void SetGenerator(G4String genname);

	G4VPrimaryGenerator* GetGenerator() const;
	G4String GetGeneratorName() const;

private:
	int model_;
	double eta_;
	G4ParticleGun* particleGun;
	HepMCG4AsciiReader* hepmcAscii;
	HepMCG4PythiaInterface* pythiaGen;

	G4VPrimaryGenerator* currentGenerator;
	G4String currentGeneratorName;
	std::map<G4String, G4VPrimaryGenerator*> gentypeMap;

	DetectorConstruction*    Detector;     //pointer to the geometry
	PrimaryGeneratorMessenger* gunMessenger; //messenger of this class
	G4String                   rndmFlag;     //flag for a rndm impact point
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
// ====================================================================
// inline functions
// ====================================================================

inline void PrimaryGeneratorAction::SetGenerator
(G4VPrimaryGenerator* gen)
{
	currentGenerator = gen;
}

inline void PrimaryGeneratorAction::SetGenerator(G4String genname)
{
	std::map<G4String, G4VPrimaryGenerator*>::iterator pos = gentypeMap.find(genname);
	if(pos != gentypeMap.end()) {
		currentGenerator = pos->second;
		currentGeneratorName = genname;
	}
}

inline G4VPrimaryGenerator* PrimaryGeneratorAction::GetGenerator() const
{
	return currentGenerator;
}

inline G4String PrimaryGeneratorAction::GetGeneratorName() const
{
	return currentGeneratorName;
}
#endif


