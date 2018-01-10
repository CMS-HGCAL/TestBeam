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


double dwc_z1_H2 = -109.;	//z=0 is the HGCal table, unit is cm
double dwc_z2_H2 = -235.;
double dwc_z3_H2 = -1509.;
double dwc_z4_H2 = -1769.;

double dwc_z1_H6A = -500.;	//z=0 is the HGCal table, unit is cm


//to the EDM::Event via auxiliary information
class HGCalTBWireChamberSource : public edm::ProducerSourceFromFiles {
	private:
		bool setRunAndEventInfo(edm::EventID& id, edm::TimeValue_t& time, edm::EventAuxiliary::ExperimentType& type);
		virtual void produce(edm::Event & e);
	  	void ReadAlignmentParameters(int);
	  	void ReadTimingFile(std::string, bool);
		
		std::string outputCollectionName;
		std::vector<double> slope_x;
		std::vector<double> slope_y;

		std::vector<double> wc_resolutions;

		double triggerTimeDifferenceTolerance;	//indicated in ms
		double TDCTriggerTimeStampConversionToMs; //conversion from TDC trigger time stamp to ms

		bool performAlignment;
	  	std::vector<std::string> alignmentParamaterFiles;
		
	  	std::map<std::pair<int, int> ,std::map<int, double> >loadedAlignmentParameters;
		std::map<int, double> currentAlignmentParameters;

		std::vector<std::string> timingFileNames;
		std::vector<int> sumTriggerTimes;
		std::vector<int> skipFirstNEvents;
		std::vector<int> triggerCountOffsets;
		std::vector<int> allowForTDCEventSkipping;
		std::vector<int> setupIDs;
		std::vector<int> pdgIDs;
		std::vector<double> beamEnergies;
		std::vector<int> triggerTimingFormat; 	//default: ms, 1: micro seconds
		std::vector<int> hitsPerChannelStored; 	//default: 0, 1: hits per Channel are Stored and in principle a more sophisticated analysis can be run

		std::string areaSpecification;

		std::map<int, int> trigger_to_event_table;
		std::map<int, double> event_trigger_time;

		double dwc_z1, dwc_z2, dwc_z3, dwc_z4;

		int rootTreeIndex;
		int fileCounter;
		int nextFileIndex;
		int eventCounter, goodEventCounter;
		std::vector<int> syncCounter;

		TFile *rootFile;
	  	TTree *tree;   //!pointer to the analyzed TTree or TChain

	  	unsigned int n_run; unsigned int n_trigger_tdc; unsigned int n_trigger_orm; unsigned int skippedTDCTriggers;
	  	std::vector<int> *channels;
	  	std::vector<int> *dwc_timestamps;
	  	std::map<int, std::vector<int>* > hits;
	  	int timeSinceStart;
	  	Long64_t timeSinceStart_long;
	 
	  	double ref_time_sync, ref_time_dwc, delta_T_priorDWCTrigger;


	  	TBranch                   *b_run;   
	  	TBranch                   *b_trigger;   
	  	TBranch                   *b_channels;   
	  	TBranch                   *b_dwc_timestamps;   
	  	TBranch                   *b_timeSinceStart;   
	  	std::map<int, TBranch*>	  b_hits;


	public:
		explicit HGCalTBWireChamberSource(const edm::ParameterSet & pset, edm::InputSourceDescription const& desc);
		virtual ~HGCalTBWireChamberSource() {
		}
	
};
