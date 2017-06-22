#ifndef DATAFORMATS_HGCALTBRECHIT_H
#define DATAFORMATS_HGCALTBRECHIT_H 1

#include "DataFormats/CaloRecHit/interface/CaloRecHit.h"
#include "HGCal/DataFormats/interface/HGCalTBDetId.h"
#include <vector>

/** \class HGCalTBRecHit
 *
 * \author Jeremy Mans
 *
 * \todo fix the energy threshold for low gain saturation in a different way: now it's hardcoded
 */
#define _lowGainSaturationThreshold 2000

class HGCalTBRecHit : public CaloRecHit
{

public:
	typedef DetId key_type;

	enum Flags {
		kGood = 0,
		kHighGainSaturated,
		kLowGainSaturated,
		kTotGainSaturated,
		kThirdSample
	};


	HGCalTBRecHit();
	// by default a recHit is greated with no flag
	//	HGCalTBRecHit(const DetId& id, float energyLow, float energyHigh, float time, uint32_t flags = 0); // when constructing from digis using 2 gains for the ADC
	HGCalTBRecHit(const DetId& id, float energy, float energyLow, float energyHigh, float energyToT, float time, uint32_t flags = 0); // when constructing from digis using 2 gains for the ADC
	
	/// get the id
	HGCalTBDetId id() const
	{
		return HGCalTBDetId(detid());
	};
	/////  bool isRecovered() const;
	float _energyLow, _energyHigh, _energyTot;
	float cellCenter_x;
	float cellCenter_y;
	float _time;


	float energyLow() const
	{
		return _energyLow;
	};

	float energyHigh() const
	{
		return _energyHigh;
	};

	float energyTot() const
	{
		return _energyTot;
	};

	void setTime(float time);

	float time(){ 
	  return _time; 
	} 

	// set the flags
	void setFlag(int flag)
	{
		setFlagField(1, flag, 1);
	}; // flagBits_|= (0x1 << flag);}
	void unsetFlag(int flag)
	{
		setFlagField(0, flag, 1);
	}; //_ &= ~(0x1 << flag);}

	// check if the flag is true
	bool checkFlag(int flag) const
	{
		return flagField(flag, 1);
	}; //flagBits_ & ( 0x1<<flag);}

	void setCellCenterCoordinate(float x, float y);

	float getCellCenterCartesianCoordinate(int index);	//index of the access
};

std::ostream& operator<<(std::ostream& s, const HGCalTBRecHit& hit);

#endif
