#include <iomanip>
#include <cmath>
#include <stdlib.h>

#include "HGCal/TBStandaloneSimulator/interface/HGCSSSimHit.hh"
#include "HGCal/TBStandaloneSimulator/interface/HGCSSGeometryConversion.hh"
#include "HGCal/TBStandaloneSimulator/interface/G4SiHit.hh"
#include "TH2Poly.h"


HGCSSSimHit::HGCSSSimHit(const G4SiHit & aSiHit,
                         const unsigned & asilayer,
                         TH2Poly* map,
                         const float )
{
	energy_ = aSiHit.energy;
	//energy weighted time
	//PS: need to call calculateTime() after all hits
	//have been added to have divided by totalE!!
	time_ = aSiHit.time * aSiHit.energy;
	zpos_ = aSiHit.hit_z;
	setLayer(aSiHit.layer, asilayer);

	//coordinates in mm
	double x = aSiHit.hit_x;
	double y = aSiHit.hit_y;
	//cellid encoding:
	cellid_ = map->FindBin(x, y);

	//std::cout << "aSiHit.layer = " << aSiHit.layer
	//	    << " asilayer = " << asilayer
	//	    << std::endl;

	nGammas_ = 0;
	nElectrons_ = 0;
	nMuons_ = 0;
	nNeutrons_ = 0;
	nProtons_ = 0;
	nHadrons_ = 0;
	if(abs(aSiHit.pdgId) == 22) nGammas_++;
	else if(abs(aSiHit.pdgId) == 11) nElectrons_++;
	else if(abs(aSiHit.pdgId) == 13) nMuons_++;
	else if(abs(aSiHit.pdgId) == 2112) nNeutrons_++;
	else if(abs(aSiHit.pdgId) == 2212) nProtons_++;
	else nHadrons_++;

	trackIDMainParent_ = aSiHit.parentId;
	energyMainParent_ = aSiHit.energy;

}

void HGCSSSimHit::Add(const G4SiHit & aSiHit)
{

	time_ = time_ + aSiHit.time * aSiHit.energy;
	//PS: need to call calculateTime() after all hits
	//have been added to have divided by totalE!!

	if(abs(aSiHit.pdgId) == 22) nGammas_++;
	else if(abs(aSiHit.pdgId) == 11) nElectrons_++;
	else if(abs(aSiHit.pdgId) == 13) nMuons_++;
	else if(abs(aSiHit.pdgId) == 2112) nNeutrons_++;
	else if(abs(aSiHit.pdgId) == 2212) nProtons_++;
	else nHadrons_++;

	energy_ += aSiHit.energy;
	if (aSiHit.energy > energyMainParent_) {
		trackIDMainParent_ = aSiHit.parentId;
		energyMainParent_ = aSiHit.energy;
	}

}


std::pair<double, double> HGCSSSimHit::get_xy(const bool isScintillator,
        const HGCSSGeometryConversion & aGeom) const
{
	if (isScintillator) return aGeom.squareGeom.find(cellid_)->second;
	else return aGeom.hexaGeom.find(cellid_)->second;

}

ROOT::Math::XYZPoint HGCSSSimHit::position(const bool isScintillator,
        const HGCSSGeometryConversion & aGeom) const
{
	std::pair<double, double> xy = get_xy(isScintillator, aGeom);
	return ROOT::Math::XYZPoint(xy.first / 10., xy.second / 10., zpos_ / 10.);
}

double HGCSSSimHit::theta(const bool isScintillator,
                          const HGCSSGeometryConversion & aGeom) const
{
	return 2 * atan(exp(-1.*eta(isScintillator, aGeom)));
}

double HGCSSSimHit::eta(const bool isScintillator,
                        const HGCSSGeometryConversion & aGeom) const
{
	return position(isScintillator, aGeom).eta();
}

double HGCSSSimHit::phi(const bool isScintillator,
                        const HGCSSGeometryConversion & aGeom) const
{
	return position(isScintillator, aGeom).phi();
}

void HGCSSSimHit::Print(std::ostream & aOs) const
{
	aOs << "====================================" << std::endl
	    << " = Layer " << layer() << " siLayer " << silayer() << " cellid " << cellid_ << std::endl
	    << " = Energy " << energy_ << " time " << time_ << std::endl
	    << " = g " << nGammas_
	    << " e " << nElectrons_
	    << " mu " << nMuons_
	    << " neutron " << nNeutrons_
	    << " proton " << nProtons_
	    << " had " << nHadrons_
	    << std::endl
	    << " = main parent: trackID " << trackIDMainParent_ << " efrac " << mainParentEfrac()
	    << std::endl
	    << "====================================" << std::endl;

}

