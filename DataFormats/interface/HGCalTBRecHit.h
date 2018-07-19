#ifndef DATAFORMATS_HGCALTBRECHIT_H
#define DATAFORMATS_HGCALTBRECHIT_H 1

#include "DataFormats/CaloRecHit/interface/CaloRecHit.h"
#include "HGCal/DataFormats/interface/HGCalTBDetId.h"
#include <cstdlib>
#include <vector>
#include "FWCore/Utilities/interface/TypeWithDict.h"

/** \class HGCalTBRecHit
 *
 * \author Jeremy Mans
 *
 * \todo fix the energy threshold for low gain saturation in a different way: now it's hardcoded
 */

class HGCalTBRecHit : public CaloRecHit
{

public:
	typedef DetId key_type;

	enum Flags {
		kGood = 0,
		kHGFitFailed,
		kLGFitFailed,
		kHighGainSaturated,
		kLowGainSaturated,
		kTotGainSaturated,
		kFullyCalibrated
	};


	HGCalTBRecHit();
	// by default a recHit is greated with no flag
	HGCalTBRecHit(const DetId& id, Float16_t energy, Float16_t energyLow, Float16_t energyHigh, unsigned int short energyToT, Float16_t time, uint32_t flags = 0); // when constructing from digis using 2 gains for the ADC
	
	/// get the id
	HGCalTBDetId id() const
	{
		return HGCalTBDetId(detid());
	};
	/////  bool isRecovered() const;
	Float16_t _energyLow, _energyHigh;
	unsigned int short _energyTot;
	std::pair<Float16_t, Float16_t> _energyTSLow, _energyTSHigh;
	Float16_t _energy_HGExcl;
	Float16_t cellCenter_x;
	Float16_t cellCenter_y;
	Float16_t _time;
  	Float16_t _timeMaxHG;
  	Float16_t _timeMaxLG;

	unsigned int short _toaRise, _toaFall;

	Float16_t energyLow() const
	{
		return _energyLow;
	};

	Float16_t energyHigh() const
	{
		return _energyHigh;
	};

	std::pair<Float16_t, Float16_t> energyTSLow() const
	{
		return _energyTSLow;
	};

	std::pair<Float16_t, Float16_t> energyTSHigh() const
	{
		return _energyTSHigh;
	};

	unsigned int short energyTot() const
	{
		return _energyTot;
	};

	Float16_t energy_HGExcl() const
	{
		return _energy_HGExcl;
	};

	void setTime(Float16_t time){_time = time;return;}
	void setTimeMaxHG(Float16_t time){_timeMaxHG = time;return;}
	void setTimeMaxLG(Float16_t time){_timeMaxLG = time;return;}
	void setToaRise(unsigned int short toaRise) { _toaRise = toaRise; } ;
	void setToaFall(unsigned int short toaFall) { _toaFall = toaFall; } ;

	//all ADC
	void setEnergyTOT(unsigned int short _energy) {_energyTot=_energy;};
	void setEnergyLow(Float16_t _energy) {_energyLow=_energy;};
	void setEnergyHigh(Float16_t _energy) {_energyHigh=_energy;};
	void setEnergyTSLow(Float16_t _energy1, Float16_t _energy2) {_energyTSLow=std::make_pair(_energy1, _energy2);};
	void setEnergyTSHigh(Float16_t _energy1, Float16_t _energy2) {_energyTSHigh=std::make_pair(_energy1, _energy2);};
	void setEnergy_HGExcl(Float16_t _energy) {_energy_HGExcl=_energy;};

	Float16_t time(){ 
	  return _time; 
	} 
	Float16_t timeMaxHG(){ 
	  return _timeMaxHG; 
	} 	
	Float16_t timeMaxLG(){ 
	  return _timeMaxLG; 
	} 	
	unsigned int short toaRise() const { return _toaRise; };
	unsigned int short toaFall() const { return _toaFall; };

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

	void setCellCenterCoordinate(Float16_t x, Float16_t y);

	Float16_t getCellCenterCartesianCoordinate(int index);	//index of the access

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
