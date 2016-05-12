#include <cassert>
#include "G4VPhysicalVolume.hh"
#include "G4SystemOfUnits.hh"
#include "HGCal/TBStandaloneSimulator/interface/SamplingSection.hh"

using namespace std;

int SamplingSection::sens_layer_count = 0;

SamplingSection::SamplingSection(std::vector<TBGeometry::Element>&
                                 elements_)
	: elements(elements_)
{
	G4cout << " -- BEGIN(sampling section)" << G4endl;

	Total_thick = 0;
	n_sens_elements = 0;
	n_elements = 0;
	ele_name.clear();
	ele_X0.clear();
	ele_L0.clear();
	ele_vol.clear();
	sens_layer.clear();
	hasScintillator = false;
	n_elements = elements.size();

	for(size_t c = 0; c < elements.size(); c++) {
		TBGeometry::Element& e = elements[c];
		G4double units = mm;
		string theunits = e.smap["units"];
		if      ( theunits == "m" )
			units = m;
		else if ( theunits == "cm" )
			units = cm;

		bool isSensitive = static_cast<bool>(e.imap["sensitive"]);
		string material  = e.smap["material"];
		double thickness = e.dmap["thickness"] * units;
		double zpos      = e.dmap["z"] * units;

		if ( material == "Scintillator" ) hasScintillator = true;

		ele_name.push_back(material);
		ele_thick.push_back(thickness);
		ele_X0.push_back(0);
		ele_dEdx.push_back(0);
		ele_L0.push_back(0);
		ele_vol.push_back(0);

		Total_thick += thickness;

		if ( isSensitive ) {
			G4SiHitVec lVec;
			sens_HitVec.push_back(lVec);
			sens_layer.push_back(sens_layer_count);
			sens_layer_count++;
		}

		cout << "\tmaterial: " << material
		     << "\tthickness: " << thickness << theunits
		     << "\tz: " << zpos;
		if ( isSensitive ) cout << " <=== sensitive";
		cout << endl;
	}
	n_sens_elements = sens_HitVec.size();

	sens_HitVec_size_max = 0;
	resetCounters();

	cout << "        number of elements:           " << n_elements  << endl
	     << "        number of sensitive elements: " << n_sens_elements << endl
	     << " -- END(sampling section)" << endl
	     << endl;
}

//
void SamplingSection::add(G4double den, G4double dl,
                          G4double globalTime, G4int pdgId,
                          G4VPhysicalVolume* vol,
                          const G4ThreeVector & position,
                          G4int trackID, G4int parentID,
                          G4int layerId)
{
	std::string lstr = vol->GetName();
	for (size_t ie(0); ie < n_elements; ++ie) {
		if(ele_vol[ie] && lstr == ele_vol[ie]->GetName()) {
			size_t eleidx = ie;
			ele_den[eleidx] += den;
			ele_dl[eleidx] += dl;
			if (isSensitiveElement(eleidx)) {
				size_t idx = getSensitiveLayerIndex(lstr);
				sens_time[idx] += den * globalTime;

				//discriminate further by particle type
				if(abs(pdgId) == 22)      sens_gFlux[idx] += den;
				else if(abs(pdgId) == 11) sens_eFlux[idx] += den;
				else if(abs(pdgId) == 13) sens_muFlux[idx] += den;
				else if (abs(pdgId) == 2112) sens_neutronFlux[idx] += den;
				else {
					sens_hadFlux[idx] += den;
				}

				//add hit
				G4SiHit lHit;
				lHit.energy = den;
				lHit.time  = globalTime;
				lHit.pdgId = pdgId;
				lHit.layer = layerId;
				lHit.hit_x = position.x();
				lHit.hit_y = position.y();
				lHit.hit_z = position.z();
				lHit.trackId = trackID;
				lHit.parentId = parentID;
				sens_HitVec[idx].push_back(lHit);
			}//if Si
		}//if in right material

	}//loop on available materials
}

//
void SamplingSection::report(bool header)
{
	if(header) G4cout << "E/[MeV]\t  Si\tAbsorber\tTotal\tSi g frac\tSi e frac\tSi mu frac\tSi had frac\tSi <t> \t nG4SiHits" << G4endl;
	G4cout << std::setprecision(3) << "\t  " << getMeasuredEnergy(false) << "\t" << getAbsorbedEnergy() << "\t\t" << getTotalEnergy() << "\t"
	       << getPhotonFraction() << "\t" << getElectronFraction() << "\t" << getMuonFraction() << "\t" << getHadronicFraction() << "\t"
	       << getAverageTime() << "\t"
	       << getTotalSensHits() << "\t"
	       << G4endl;
}

G4int SamplingSection::getTotalSensHits()
{
	G4int tot = 0;
	for (unsigned ie(0); ie < n_sens_elements; ++ie) {
		tot += sens_HitVec[ie].size();
	}
	return tot;
}

G4double SamplingSection::getTotalSensE()
{
	double etot = 0;
	for (unsigned ie(0); ie < n_elements; ++ie) {
		if (isSensitiveElement(ie)) etot += ele_den[ie];
	}
	return etot;
}

G4double SamplingSection::getAverageTime()
{
	double etot = getTotalSensE();
	double time = 0;
	for (unsigned ie(0); ie < n_sens_elements; ++ie) {
		time += sens_time[ie];
	}
	return etot > 0 ? time / etot : 0 ;
}

//
G4double SamplingSection::getPhotonFraction()
{
	double etot = getTotalSensE();
	double val = 0;
	for (unsigned ie(0); ie < n_sens_elements; ++ie) {
		val += sens_gFlux[ie];
	}
	return etot > 0 ? val / etot : 0 ;
}

//
G4double SamplingSection::getElectronFraction()
{
	double etot = getTotalSensE();
	double val = 0;
	for (unsigned ie(0); ie < n_sens_elements; ++ie) {
		val += sens_eFlux[ie];
	}
	return etot > 0 ? val / etot : 0 ;
}

//
G4double SamplingSection::getMuonFraction()
{
	double etot = getTotalSensE();
	double val = 0;
	for (unsigned ie(0); ie < n_sens_elements; ++ie) {
		val += sens_muFlux[ie];
	}
	return etot > 0 ? val / etot : 0 ;
}

//
G4double SamplingSection::getNeutronFraction()
{
	double etot = getTotalSensE();
	double val = 0;
	for (unsigned ie(0); ie < n_sens_elements; ++ie) {
		val += sens_neutronFlux[ie];
	}
	return etot > 0 ? val / etot : 0 ;
}

//
G4double SamplingSection::getHadronicFraction()
{
	double etot = getTotalSensE();
	double val = 0;
	for (unsigned ie(0); ie < n_sens_elements; ++ie) {
		val += sens_hadFlux[ie];
	}
	return etot > 0 ? val / etot : 0 ;
}

//
G4double SamplingSection::getMeasuredEnergy(bool weighted)
{
	G4double weight = (weighted ? getAbsorberX0() : 1.0);
	return weight * getTotalSensE();
}

//
G4double SamplingSection::getAbsorberX0()
{
	double val = 0;
	for (unsigned ie(0); ie < n_elements; ++ie) {
		if (isAbsorberElement(ie))
			if (ele_X0[ie] > 0) val += ele_thick[ie] / ele_X0[ie];
	}
	return val;
}
//
G4double SamplingSection::getAbsorberdEdx()
{
	double val = 0;
	for (unsigned ie(0); ie < n_elements; ++ie) {
		if (isAbsorberElement(ie))
			val += ele_thick[ie] * ele_dEdx[ie];
	}
	return val;
}

//
G4double SamplingSection::getAbsorberLambda()
{
	double val = 0;
	for (unsigned ie(0); ie < n_elements; ++ie) {
		if (isAbsorberElement(ie) && ele_L0[ie] > 0)
			val += ele_thick[ie] / ele_L0[ie];
	}
	return val;
}

//
G4double SamplingSection::getAbsorbedEnergy()
{
	double val = 0;
	for (unsigned ie(0); ie < n_elements; ++ie) {
		if (isAbsorberElement(ie)) val += ele_den[ie];
	}
	return val;
}

//
G4double SamplingSection::getTotalEnergy()
{
	double val = 0;
	for (unsigned ie(0); ie < n_elements; ++ie) {
		val += ele_den[ie];
	}
	return val;
}

const G4SiHitVec & SamplingSection::getSiHitVec(const unsigned & idx) const
{
	assert(sens_HitVec.size() > 0);
	return sens_HitVec[idx];
}

void SamplingSection::trackParticleHistory(const unsigned & idx,
        const G4SiHitVec & incoming)
{
	for (unsigned iP(0); iP < sens_HitVec[idx].size(); ++iP) { //loop on g4hits
		G4int parId = sens_HitVec[idx][iP].parentId;
		for (unsigned iI(0); iI < incoming.size(); ++iI) { //loop on previous layer
			G4int trId = incoming[iI].trackId;
			if (trId == parId)
				sens_HitVec[idx][iP].parentId = incoming[iI].parentId;
		}//loop on previous layer
	}//loop on g4hits
}
