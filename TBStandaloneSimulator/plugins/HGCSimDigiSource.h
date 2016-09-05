#ifndef HGCSIMDIGISOURCE_H
#define HGCSIMDIGISOURCE_H
/** \class

	\author Harrison B. Prosper
 */
#include <stdio.h>
#include <iostream>
#include <vector>

#include "TChain.h"
#include "TTree.h"
#include "TRandom3.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Sources/interface/ProducerSourceFromFiles.h"
#include "FWCore/ParameterSet/interface/ConfigurationDescriptions.h"
#include "FWCore/ParameterSet/interface/ParameterSetDescription.h"
#include "FWCore/Utilities/interface/Exception.h"

#include "DataFormats/FEDRawData/interface/FEDRawDataCollection.h"
#include "HGCal/CondObjects/interface/HGCalElectronicsMap.h"
#include "HGCal/TBStandaloneSimulator/interface/HGCSSSimHit.hh"
#include "HGCal/TBStandaloneSimulator/interface/HGCCellMap.h"
#include "HGCal/DataFormats/interface/HGCalTBDataFrameContainers.h"

class HGCSimDigiSource : public edm::ProducerSourceFromFiles
{
public:
	explicit HGCSimDigiSource(const edm::ParameterSet& pset,
	                          edm::InputSourceDescription
	                          const& desc);

	virtual ~HGCSimDigiSource();

	static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);

	struct Cell {
		int skiroc;      // SKIROC number
		int channel;     // channel number
		DetId detid;     // detector id
		uint16_t ADClow;
		uint16_t ADChigh;
		uint16_t TDC;
		int layer;
		int sensor_u;
		int sensor_v;
		int u;
		int v;
		int celltype;
		double x;
		double y;
		double energy;

		bool operator<(HGCSimDigiSource::Cell& o)
		{
			// SKIROC 2 occurs before SKIROC 1
			int ski   = skiroc % 2;
			int o_ski = o.skiroc % 2;
			int lhs   = 64 * ski + channel;
			int rhs   = 64 * o_ski + o.channel;
			return lhs < rhs;
		}
	};

private:
	virtual bool setRunAndEventInfo(edm::EventID& id,
	                                edm::TimeValue_t& time,
	                                edm::EventAuxiliary::ExperimentType&);

	virtual void produce(edm::Event& e);

	virtual void digitize(std::vector<HGCSimDigiSource::Cell>& channels);

	int _run;
	int _maxevents;
	int _minadccount;
	double _adcpermev;
	std::vector<std::string> _filenames;       ///<name of input sim files
	std::vector<std::string> _noisefilenames;  ///<name of input pedestal files

	/// Sim objects
	TChain* _chain;
	TTree*  _tree;
	size_t  _entries;
	size_t  _entry;
	HGCCellMap  _cellmap;
	HGCalElectronicsMap _emap;
	HGCSSSimHitVec*  _simhits;

	/// noise objects
	TChain* _noisechain;
	TTree*  _noisetree;
	size_t  _noiseentries;
	size_t  _noiseentry;
	TRandom3 _random;
	std::vector<std::map<uint32_t, uint16_t> > _noise;
};


#endif
