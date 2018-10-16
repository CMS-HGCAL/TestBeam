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


//DWC_to_TDC_MAP
int DWC1_LEFT;
int DWC1_RIGHT;
int DWC1_DOWN;
int DWC1_UP;
int DWC2_LEFT;
int DWC2_RIGHT;
int DWC2_DOWN;
int DWC2_UP;
int DWC3_LEFT;
int DWC3_RIGHT;
int DWC3_DOWN;
int DWC3_UP;
int DWC4_LEFT;
int DWC4_RIGHT;
int DWC4_DOWN;
int DWC4_UP;


//indication in cm
double dwc_z1_H2_Summer2017 = -109.;	//z=0 is the HGCal table, unit is cm
double dwc_z2_H2_Summer2017 = -235.;
double dwc_z3_H2_Summer2017 = -1509.;
double dwc_z4_H2_Summer2017 = -1769.;

double dwc_z1_H6A_October2017 = -500.;	//z=0 is the HGCal table, unit is cm

double dwc_z1_H2_June2018 = -120.;	//z=0 is the HGCal table, unit is cm
double dwc_z2_H2_June2018 = -246.;
double dwc_z3_H2_June2018 = -1520.;
double dwc_z4_H2_June2018 = -1780.;


double dwc_z1_H2_October2018 = -3300.;	//z=0 is the HGCal table, unit is cm
double dwc_z2_H2_October2018 = -2900.;
double dwc_z3_H2_October2018 = -900.;
double dwc_z4_H2_October2018 = -200.;


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


	  	int N_TDC_channels;
	  	int fileFormat;

	public:
		explicit HGCalTBWireChamberSource(const edm::ParameterSet & pset, edm::InputSourceDescription const& desc);
		virtual ~HGCalTBWireChamberSource() {
		}
	
};
