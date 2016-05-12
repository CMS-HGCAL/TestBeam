#ifndef RunAction_h
#define RunAction_h 1

#include "G4UserRunAction.hh"
#include "globals.hh"

class G4Run;

class RunAction : public G4UserRunAction
{
public:
	RunAction();
	virtual ~RunAction();

	void BeginOfRunAction(const G4Run*);
	void   EndOfRunAction(const G4Run*);

	void fillPerEvent(G4double, G4double, G4double, G4double);

private:
	G4double sumEAbs, sum2EAbs;
	G4double sumEGap, sum2EGap;

	G4double sumLAbs, sum2LAbs;
	G4double sumLGap, sum2LGap;
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif

