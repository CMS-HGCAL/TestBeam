#ifndef _samplingsection_hh_
#define _samplingsection_hh_

#include "G4VUserDetectorConstruction.hh"
#include "globals.hh"
#include "G4ThreeVector.hh"
#include "G4Colour.hh"

#include <iomanip>
#include <vector>

#include "G4SiHit.hh"
#include "HGCal/TBStandaloneSimulator/interface/TBGeometry.h"

class SamplingSection
{
public:
	SamplingSection(std::vector<TBGeometry::Element>& elements_);

	//DTOR
	~SamplingSection() { };

	TBGeometry::Element getElement(int aEle)
	{
		if ( aEle > (int)elements.size() - 1 ) return TBGeometry::Element();
		if ( aEle < -(int)elements.size() )  return TBGeometry::Element();
		if ( aEle >= 0 )
			return elements[aEle];
		else
			return elements[elements.size() + aEle];
	}

	//
	void add(G4double den, G4double dl,
	         G4double globalTime, G4int pdgId,
	         G4VPhysicalVolume* vol,
	         const G4ThreeVector & position,
	         G4int trackID, G4int parentID,
	         G4int layerId);

	inline bool isSensitiveElement(size_t aEle)
	{
		if ( aEle < n_elements )
			return elements[aEle].imap["sensitive"];
		else
			return false;
	};

	inline unsigned getSensitiveLayerIndex(std::string astr)
	{
		if (astr.find("_") == astr.npos) return 0;
		size_t pos = astr.find("phys");
		if (pos != astr.npos && pos > 1) {
			unsigned idx = 0;
			std::istringstream(astr.substr(pos - 1, 1)) >> idx;
			return idx;
		}
		return 0;
	};

	inline size_t SiHitVecSize()
	{
		return sens_HitVec.size();
	}

	inline G4Colour g4Colour(size_t aEle)
	{
		if (isSensitiveElement(aEle))return G4Colour::Red();
		if (ele_name[aEle] == "Cu")  return G4Colour::Brown();
		if (isAbsorberElement(aEle)) return G4Colour::Gray();
		if (ele_name[aEle] == "PCB") return G4Colour::Blue();
		if (ele_name[aEle] == "Air") return G4Colour::Cyan();
		return G4Colour::Yellow();
	};

	inline bool isAbsorberElement(size_t aEle)
	{
		if (aEle < n_elements &&
		        (
		            ele_name[aEle] == "Pb"  || ele_name[aEle] == "Cu" ||
		            ele_name[aEle] == "W"   || ele_name[aEle] == "Brass" ||
		            ele_name[aEle] == "Fe"  || ele_name[aEle] == "Steel" ||
		            ele_name[aEle] == "SSteel" || ele_name[aEle] == "Al" ||
		            ele_name[aEle] == "WCu" || ele_name[aEle] == "NeutMod"
		        )
		   ) return true;
		return false;
	};

	//reset
	inline void resetCounters()
	{
		ele_den.resize(n_elements, 0);
		ele_dl.resize(n_elements, 0);
		sens_time.resize(n_sens_elements, 0);
		sens_gFlux.resize(n_sens_elements, 0);
		sens_eFlux.resize(n_sens_elements, 0);
		sens_muFlux.resize(n_sens_elements, 0);
		sens_neutronFlux.resize(n_sens_elements, 0);
		sens_hadFlux.resize(n_sens_elements, 0);
		//reserve some space based on first event....
		for (unsigned idx(0); idx < n_sens_elements; ++idx) {
			if (sens_HitVec[idx].size() > sens_HitVec_size_max) {
				sens_HitVec_size_max = 2 * sens_HitVec[idx].size();
				G4cout << "-- SamplingSection::resetCounters(), space reserved for HitVec vector increased to " << sens_HitVec_size_max << G4endl;
			}
			sens_HitVec[idx].clear();
			sens_HitVec[idx].reserve(sens_HitVec_size_max);
		}
	}

	//
	G4double getMeasuredEnergy(bool weighted = true);
	G4double getAbsorbedEnergy();
	G4double getTotalEnergy();
	G4double getAbsorberX0();
	G4double getAbsorberdEdx();
	G4double getAbsorberLambda();
	G4double getHadronicFraction();
	G4double getNeutronFraction();
	G4double getMuonFraction();
	G4double getPhotonFraction();
	G4double getElectronFraction();
	G4double getAverageTime();
	G4int    getTotalSensHits();
	G4double getTotalSensE();

	const G4SiHitVec & getSiHitVec(const unsigned & idx) const;
	void trackParticleHistory(const unsigned & idx, const G4SiHitVec & incoming);

	//
	void report(bool header = false);

	//members
	std::vector<TBGeometry::Element> elements;

	size_t n_elements;
	size_t n_sens_elements;
	std::vector<std::string>        ele_name;
	std::vector<G4double>           ele_thick;
	std::vector<G4double>           ele_X0;
	std::vector<G4double>           ele_dEdx;
	std::vector<G4double>           ele_L0;
	std::vector<G4double>           ele_den;
	std::vector<G4double>           ele_dl;
	std::vector<G4VPhysicalVolume*> ele_vol;
	std::vector<G4double>
	sens_gFlux, sens_eFlux, sens_muFlux,
	            sens_neutronFlux, sens_hadFlux, sens_time;
	G4double Total_thick;
	std::vector<G4SiHitVec> sens_HitVec;
	unsigned sens_HitVec_size_max;
	bool hasScintillator;
	std::vector<G4int> sens_layer;

	static G4int sens_layer_count;
};

#endif
