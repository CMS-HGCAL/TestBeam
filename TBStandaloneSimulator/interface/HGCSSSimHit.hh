#ifndef _hgcsssimhit_hh_
#define _hgcsssimhit_hh_

#include <iomanip>
#include <vector>
#include <sstream>
#include <iostream>
#include <map>

#include "Rtypes.h"
#include "Math/Point3D.h"
#include "Math/Point3Dfwd.h"

class TH2Poly;
class G4SiHit;
class HGCSSGeometryConversion;

//for hexagons: side size.
static const float CELL_SIZE_X = 6.496345; //2.5;//mm
static const float CELL_SIZE_Y = CELL_SIZE_X;

class HGCSSSimHit
{

public:

	typedef unsigned key_type;

	HGCSSSimHit():
		energy_(0),
		time_(0),
		zpos_(0),
		layer_(0),
		cellid_(0),
		nGammas_(0),
		nElectrons_(0),
		nMuons_(0),
		nNeutrons_(0),
		nProtons_(0),
		nHadrons_(0),
		trackIDMainParent_(0),
		energyMainParent_(0)
	{

	};

	HGCSSSimHit(const G4SiHit & aSiHit,
	            const unsigned & asilayer,
	            TH2Poly* map,
	            const float cellSize = CELL_SIZE_X);

	virtual ~HGCSSSimHit() {}

	unsigned id() const
	{
		return cellid_;
	}

	inline double energy() const
	{
		return energy_;
	};

	inline double time() const
	{
		return time_;
	};

	inline void calculateTime()
	{
		if (energy_ > 0) time_ = time_ / energy_;
	};

	inline unsigned layer() const
	{
		return layer_;
	};

	inline unsigned silayer() const
	{
		return silayer_;
	};

	inline void setLayer(const unsigned & layer, const unsigned & silayer)
	{
		layer_ = layer;
		silayer_ = silayer;
	};

	inline unsigned cellid() const
	{
		return cellid_;
	};

	inline unsigned nGammas() const
	{
		return nGammas_;
	};

	inline unsigned nElectrons() const
	{
		return nElectrons_;
	};

	inline unsigned nMuons() const
	{
		return nMuons_;
	};

	inline unsigned nNeutrons() const
	{
		return nNeutrons_;
	};

	inline unsigned nProtons() const
	{
		return nProtons_;
	};
	inline unsigned nHadrons() const
	{
		return nHadrons_;
	};
	inline unsigned numberOfParticles() const
	{
		return nGammas_ + nElectrons_ + nMuons_ + nNeutrons_ + nProtons_ + nHadrons_;
	};

	inline double gFrac() const
	{
		return nGammas_ / numberOfParticles();
	};

	inline double eFrac() const
	{
		return nElectrons_ / numberOfParticles();
	};

	inline  double muFrac() const
	{
		return nMuons_ / numberOfParticles();
	};

	inline double neutronFrac() const
	{
		return nNeutrons_ / numberOfParticles();
	};

	inline double protonFrac() const
	{
		return nProtons_ / numberOfParticles();
	};

	inline double hadFrac() const
	{
		return nHadrons_ / numberOfParticles();
	};

	void Add(const G4SiHit & aSiHit);

	std::pair<double, double> get_xy(const bool isScintillator,
	                                 const HGCSSGeometryConversion & aGeom) const;

	ROOT::Math::XYZPoint position(const bool isScintillator,
	                              const HGCSSGeometryConversion & aGeom) const;


	inline double get_z() const
	{
		return zpos_;
	};

	double eta(const bool isScintillator,
	           const HGCSSGeometryConversion & aGeom) const;
	double theta(const bool isScintillator,
	             const HGCSSGeometryConversion & aGeom) const;
	double phi(const bool isScintillator,
	           const HGCSSGeometryConversion & aGeom) const;

	inline unsigned getGranularity() const
	{
		return 1;
	};

	inline int mainParentTrackID() const
	{
		return trackIDMainParent_;
	};

	inline double mainParentEfrac() const
	{
		return energyMainParent_ / energy_;
	};

	void Print(std::ostream & aOs) const ;

private:

	double energy_;
	double time_;
	double zpos_;
	unsigned layer_;
	unsigned silayer_;
	unsigned cellid_;
	unsigned nGammas_;
	unsigned nElectrons_;
	unsigned nMuons_;
	unsigned nNeutrons_;
	unsigned nProtons_;
	unsigned nHadrons_;
	int trackIDMainParent_;
	double energyMainParent_;

	ClassDef(HGCSSSimHit, 1);



};


typedef std::vector<HGCSSSimHit> HGCSSSimHitVec;



#endif
