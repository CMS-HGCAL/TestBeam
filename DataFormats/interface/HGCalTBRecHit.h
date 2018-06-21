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
	float _energyTSLow, _energyTSHigh;
	float cellCenter_x;
	float cellCenter_y;
	float _time;
  	float _timeMaxHG;
  	float _timeMaxLG;

	int _toaRise, _toaFall;

	float energyLow() const
	{
		return _energyLow;
	};

	float energyHigh() const
	{
		return _energyHigh;
	};

	float energyTSLow() const
	{
		return _energyTSLow;
	};

	float energyTSHigh() const
	{
		return _energyTSHigh;
	};

	float energyTot() const
	{
		return _energyTot;
	};

	void setTime(float time){_time = time;return;}
	void setTimeMaxHG(float time){_timeMaxHG = time;return;}
	void setTimeMaxLG(float time){_timeMaxLG = time;return;}
	void setToaRise(float toaRise) { _toaRise = toaRise; } ;
	void setToaFall(float toaFall) { _toaFall = toaFall; } ;

	void setEnergyTOT(float _energy) {_energyTot=_energy;};
	void setEnergyLow(float _energy) {_energyLow=_energy;};
	void setEnergyHigh(float _energy) {_energyHigh=_energy;};
	void setEnergyTSLow(float _energy) {_energyTSLow=_energy;};
	void setEnergyTSHigh(float _energy) {_energyTSHigh=_energy;};

	float time(){ 
	  return _time; 
	} 
	float timeMaxHG(){ 
	  return _timeMaxHG; 
	} 	
	float timeMaxLG(){ 
	  return _timeMaxLG; 
	} 	
	float toaRise() const { return _toaRise; };
	float toaFall() const { return _toaFall; };

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

	void setUnderSaturationForHighGain(){ m_underSaturationHG=true; }
	void setUnderSaturationForLowGain(){ m_underSaturationLG=true; }
	bool isUnderSaturationForHighGain(){ return m_underSaturationHG; }
	bool isUnderSaturationForLowGain(){ return m_underSaturationLG; }
 private:
	bool m_underSaturationHG;
	bool m_underSaturationLG;
	
};

std::ostream& operator<<(std::ostream& s, const HGCalTBRecHit& hit);

#endif
