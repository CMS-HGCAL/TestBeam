#include "FWCore/Framework/interface/Event.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Sources/interface/ProducerSourceFromFiles.h"
#include "FWCore/ParameterSet/interface/ConfigurationDescriptions.h"
#include "FWCore/ParameterSet/interface/ParameterSetDescription.h"

#include "HGCal/DataFormats/interface/HGCalTBRunData.h"
#include "HGCal/DataFormats/interface/HGCalTBWireChamberData.h"

#include <fstream>
#include <iostream>
#include <stdio.h>
#include <string>
#include <sstream>
#include "stdlib.h"
#include <vector>
#include <cmath>

#include "TFile.h"
#include "TBranch.h"
#include "TTree.h"



enum DWC_to_TDC_MAP {
	DWC1_LEFT = 0,
	DWC1_RIGHT = 1,
	DWC1_DOWN = 2,
	DWC1_UP = 3,
	DWC2_LEFT = 4,
	DWC2_RIGHT = 5,
	DWC2_DOWN = 6,
	DWC2_UP = 7,
	DWC3_LEFT = 14,		//default: 15 --> but mind the flip for these chambers
	DWC3_RIGHT = 15,	//default: 14 --> but mind the flip for these chambers
	DWC3_DOWN = 13,
	DWC3_UP = 12,
	DWC4_LEFT = 10,		//defalt: 10 --> but mind the flip for these chambers
	DWC4_RIGHT = 11,	//default: 11 --> but mind the flip for these chambers
	DWC4_DOWN = 9,
	DWC4_UP = 8
};

//numbers to be validated through more precise measurements
double dwc_z1 = -103;	//z=0 is the HGCal table, unit is cm
double dwc_z2 = -231.;
double dwc_z3 = -1479.;
double dwc_z4 = -1784.;

//to the EDM::Event via auxiliary information
class HGCalTBWireChamberSource : public edm::ProducerSourceFromFiles {
	private:
		bool setRunAndEventInfo(edm::EventID& id, edm::TimeValue_t& time, edm::EventAuxiliary::ExperimentType& type);
		virtual void produce(edm::Event & e);
	  	virtual void beginJob() override;
	  	void ReadAlignmentParameters();
	  	void ReadTimingFile(std::string);
		
		std::string outputCollectionName;
		std::vector<double> slope_x;
		std::vector<double> slope_y;

		bool performAlignment;
	  	std::string alignmentParamaterFile;
		std::map<int, double> alignmentParameters;

		std::vector<std::string> timingFileNames;
		std::vector<int> skipFirstEventInDWCProducer;
		std::vector<int> runType;
		std::map<int, int> trigger_to_event_table;

		int rootTreeIndex;
		int fileCounter;
		int nextFileIndex;
		
		TFile *rootFile;
	  	TTree *tree;   //!pointer to the analyzed TTree or TChain

	  	unsigned int n_run; unsigned int n_trigger; 
	  	std::vector<int> *channels;
	  	std::vector<int> *dwc_timestamps;
	 

	  	TBranch                   *b_run;   
	  	TBranch                   *b_trigger;   
	  	TBranch                   *b_channels;   
	  	TBranch                   *b_dwc_timestamps;   


	public:
		explicit HGCalTBWireChamberSource(const edm::ParameterSet & pset, edm::InputSourceDescription const& desc);
		virtual ~HGCalTBWireChamberSource() {
		}
	
};
