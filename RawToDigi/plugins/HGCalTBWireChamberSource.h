#include "FWCore/Framework/interface/Event.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Sources/interface/ProducerSourceFromFiles.h"
#include "FWCore/ParameterSet/interface/ConfigurationDescriptions.h"
#include "FWCore/ParameterSet/interface/ParameterSetDescription.h"

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
	DWC3_LEFT = 15,
	DWC3_RIGHT = 14,
	DWC3_DOWN = 13,
	DWC3_UP = 12,
	DWC4_LEFT = 11,
	DWC4_RIGHT = 10,
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
	  	virtual void endJob() override;
	  	void ReadAlignmentParameters();
		
		std::string outputCollectionName;
		std::vector<double> slope_x;
		std::vector<double> slope_y;

		bool performAlignment;
	  	std::string alignmentParamaterFile;
		std::map<int, double> alignmentParameters;

		int eventCounter;
		int fileCounter;
		int nextFileIndex;
		
		TFile *rootFile;
	  	TTree *tree;   //!pointer to the analyzed TTree or TChain

	  	unsigned int n_run; unsigned int n_trigger; 
	  	std::vector<int> *channels;
	  	std::vector<int> *dwc_timestamps;
	 
	  	bool makeTree;
	  	TTree* outTree;
	  	double time_DWC1;
	  	double time_DWC2;
	  	double time_DWC3;
	  	double time_DWC4;

	  	TBranch                   *b_run;   
	  	TBranch                   *b_trigger;   
	  	TBranch                   *b_channels;   
	  	TBranch                   *b_dwc_timestamps;   


	public:
		explicit HGCalTBWireChamberSource(const edm::ParameterSet & pset, edm::InputSourceDescription const& desc);
		virtual ~HGCalTBWireChamberSource() {
		}
	
};
