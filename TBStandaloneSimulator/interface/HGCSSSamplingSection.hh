#ifndef _hgcsssamplingsection_hh_
#define _hgcsssamplingsection_hh_

#include <iomanip>
#include <vector>
#include <sstream>
#include <map>

#include "Rtypes.h"

class HGCSSSamplingSection
{


public:

	typedef unsigned key_type;

	HGCSSSamplingSection():
		volNb_(0),
		volX0trans_(0),
		voldEdx_(0),
		volLambdatrans_(0),
		measuredE_(0),
		absorberE_(0),
		totalE_(0),
		gFrac_(0),
		eFrac_(0),
		muFrac_(0),
		neutronFrac_(0),
		hadFrac_(0),
		avgTime_(0),
		nSiHits_(0)
	{

	};

	virtual ~HGCSSSamplingSection() {};

	unsigned id() const
	{
		return volNb_;
	}

	//getters
	inline unsigned volNb() const
	{
		return volNb_;
	};

	inline double volX0trans() const
	{
		return volX0trans_;
	};

	inline double voldEdx() const
	{
		return voldEdx_;
	};

	inline double volLambdatrans() const
	{
		return volLambdatrans_;
	};
	inline double measuredE() const
	{
		return measuredE_;
	};
	inline double absorberE() const
	{
		return absorberE_;
	};
	inline double totalE() const
	{
		return totalE_;
	};
	inline double gFrac() const
	{
		return gFrac_;
	};
	inline double eFrac() const
	{
		return eFrac_;
	};
	inline double muFrac() const
	{
		return muFrac_;
	};
	inline double neutronFrac() const
	{
		return neutronFrac_;
	};
	inline double hadFrac() const
	{
		return hadFrac_;
	};
	inline double avgTime() const
	{
		return avgTime_;
	};
	inline unsigned nSiHits() const
	{
		return nSiHits_;
	};

	//setters
	inline void volNb(const unsigned & aVal)
	{
		volNb_ = aVal;
	};
	inline void volX0trans(const double & aVal)
	{
		volX0trans_ = aVal;
	};
	inline void voldEdx(const double & aVal)
	{
		voldEdx_ = aVal;
	};
	inline void volLambdatrans(const double & aVal)
	{
		volLambdatrans_ = aVal;
	};
	inline void measuredE(const double & aVal)
	{
		measuredE_ = aVal;
	};
	inline void absorberE(const double & aVal)
	{
		absorberE_ = aVal;
	};
	inline void totalE(const double & aVal)
	{
		totalE_ = aVal;
	};
	inline void gFrac(const double & aVal)
	{
		gFrac_ = aVal;
	};
	inline void eFrac(const double & aVal)
	{
		eFrac_ = aVal;
	};
	inline void muFrac(const double & aVal)
	{
		muFrac_ = aVal;
	};
	inline void neutronFrac(const double & aVal)
	{
		neutronFrac_ = aVal;
	};
	inline void hadFrac(const double & aVal)
	{
		hadFrac_ = aVal;
	};
	inline void avgTime(const double & aVal)
	{
		avgTime_ = aVal;
	};
	inline void nSiHits(const unsigned & aVal)
	{
		nSiHits_ = aVal;
	};

private:
	unsigned volNb_;
	double volX0trans_;
	double voldEdx_;
	double volLambdatrans_;
	double measuredE_;
	double absorberE_;
	double totalE_;
	double gFrac_;
	double eFrac_;
	double muFrac_;
	double neutronFrac_;
	double hadFrac_;
	double avgTime_;
	unsigned nSiHits_;

	ClassDef(HGCSSSamplingSection, 1);



};

typedef std::vector<HGCSSSamplingSection> HGCSSSamplingSectionVec;

#endif













