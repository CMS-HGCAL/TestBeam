#include "HGCal/RawToDigi/plugins/HGCalTBWireChamberSource.h"

#define DEBUG


bool validTimestamp(int ts) {
	return (ts>=0);
}

HGCalTBWireChamberSource::HGCalTBWireChamberSource(const edm::ParameterSet & pset, edm::InputSourceDescription const& desc) :  edm::ProducerSourceFromFiles(pset, desc, true),
	syncCounter(3, 0), rootFile(NULL)
{
	//find and fill the configured runs
	outputCollectionName = pset.getParameter<std::string>("OutputCollectionName");

	triggerTimeDifferenceTolerance = pset.getUntrackedParameter<double>("triggerTimeDifferenceTolerance", 0.2);		//given in ms
	TDCTriggerTimeStampConversionToMs = pset.getUntrackedParameter<double>("TDCTriggerTimeStampConversionToMs", 1.0/1000);		//given in ms

  	std::vector<double> v0(4, 0.2/40);		//0.2mm/ns and TDC binning is 25ps
	slope_x = pset.getUntrackedParameter<std::vector<double> >("slope_x", v0);
	slope_y = pset.getUntrackedParameter<std::vector<double> >("slope_y", v0);
  	std::vector<double> v1(4, 1.0);
	wc_resolutions = pset.getUntrackedParameter<std::vector<double> >("wc_resolutions", v1);

	performAlignment = pset.getUntrackedParameter<bool>("performAlignment", false);
	alignmentParamaterFiles = pset.getParameter<std::vector<std::string> >("alignmentParamaterFiles");

	timingFileNames = pset.getParameter<std::vector<std::string> >("timingFileNames");
	sumTriggerTimes = pset.getParameter<std::vector<int> >("sumTriggerTimes");
	triggerCountOffsets = pset.getParameter<std::vector<int> >("triggerCountOffsets");
	skipFirstNEvents = pset.getParameter<std::vector<int> >("skipFirstNEvents");
	allowForTDCEventSkipping = pset.getParameter<std::vector<int> >("allowForTDCEventSkipping");

	setupIDs = pset.getParameter<std::vector<int> >("setupIDs");
	pdgIDs = pset.getParameter<std::vector<int> >("pdgIDs");
	beamEnergies = pset.getParameter<std::vector<double> >("beamEnergies");
	triggerTimingFormat = pset.getParameter<std::vector<int> >("triggerTimingFormat");
	hitsPerChannelStored = pset.getParameter<std::vector<int> >("hitsPerChannelStored"); 

	areaSpecification = pset.getUntrackedParameter<std::string>("areaSpecification", "H2");

	if (areaSpecification=="H2_Summer2017") {
		dwc_z1 = dwc_z1_H2_Summer2017;
		dwc_z2 = dwc_z2_H2_Summer2017;
		dwc_z3 = dwc_z3_H2_Summer2017;
		dwc_z4 = dwc_z4_H2_Summer2017;
		DWC1_LEFT = 0;
		DWC1_RIGHT = 1;
		DWC1_DOWN = 2;
		DWC1_UP = 3;
		DWC2_LEFT = 4;
		DWC2_RIGHT = 5;
		DWC2_DOWN = 6;
		DWC2_UP = 7;
		DWC3_LEFT = 14;
		DWC3_RIGHT = 15;
		DWC3_DOWN = 13;
		DWC3_UP = 12;
		DWC4_LEFT = 10;
		DWC4_RIGHT = 11;
		DWC4_DOWN = 9;
		DWC4_UP = 8;
		N_TDC_channels = 16;
		fileFormat=1;
	} else if (areaSpecification=="H6A_October2017") {
		dwc_z1 = dwc_z1_H6A_October2017;
		dwc_z2 = dwc_z3 = dwc_z4 = -1.;
		DWC1_LEFT = 0;
		DWC1_RIGHT = 1;
		DWC1_DOWN = 2;
		DWC1_UP = 3;
		DWC2_LEFT = 4;
		DWC2_RIGHT = 5;
		DWC2_DOWN = 6;
		DWC2_UP = 7;
		DWC3_LEFT = 14;
		DWC3_RIGHT = 15;
		DWC3_DOWN = 13;
		DWC3_UP = 12;
		DWC4_LEFT = 10;
		DWC4_RIGHT = 11;
		DWC4_DOWN = 9;
		DWC4_UP = 8;	
		N_TDC_channels = 16;
		fileFormat=1;	
	} else if (areaSpecification=="H2_June2018") {
		dwc_z1 = dwc_z1_H2_June2018;
		dwc_z2 = dwc_z2_H2_June2018;
		dwc_z3 = dwc_z3_H2_June2018;
		dwc_z4 = dwc_z4_H2_June2018;
		DWC1_LEFT = 0;
		DWC1_RIGHT = 1;
		DWC1_DOWN = 2;
		DWC1_UP = 3;
		DWC2_LEFT = 4;
		DWC2_RIGHT = 5;
		DWC2_DOWN = 6;
		DWC2_UP = 7;
		DWC3_LEFT = 9;	
		DWC3_RIGHT = 8;	
		DWC3_DOWN = 10;
		DWC3_UP = 11;
		DWC4_LEFT = 13;	
		DWC4_RIGHT = 12;
		DWC4_DOWN = 14;
		DWC4_UP = 15;
		N_TDC_channels = 16;
		fileFormat=1;
	} else if (areaSpecification=="H2_October2018") {
		dwc_z1 = dwc_z1_H2_October2018;
		dwc_z2 = dwc_z2_H2_October2018;
		dwc_z3 = dwc_z3_H2_October2018;
		dwc_z4 = dwc_z4_H2_October2018;
		DWC1_LEFT = 4;
		DWC1_RIGHT = 5;
		DWC1_DOWN = 6;
		DWC1_UP = 7;
		DWC2_LEFT = 9;
		DWC2_RIGHT = 8;
		DWC2_DOWN = 10;
		DWC2_UP = 11;
		DWC3_LEFT = 13;	
		DWC3_RIGHT = 12;	
		DWC3_DOWN = 14;
		DWC3_UP = 15;
		DWC4_LEFT = 1;	
		DWC4_RIGHT = 0;
		DWC4_DOWN = 2;
		DWC4_UP = 3;
		N_TDC_channels=32;
		fileFormat=2;
	}

	produces<std::map<int, WireChamberData> >(outputCollectionName);
	produces<RunData>("RunData");

	tree = NULL;

	//values from the 
	n_run=0;
	n_trigger_tdc=n_trigger_orm=0;
	channels=0;
	dwc_timestamps=0;

	for (int ch=0; ch<N_TDC_channels; ch++) {
		if (hitsPerChannelStored[ch]==1) {
			hits[ch] = 0;
		}
	}
	
	fileCounter = -1;
	rootTreeIndex = 0;
	nextFileIndex = 0;
	rootFile = NULL;
	tree = NULL;

	if (performAlignment) ReadAlignmentParameters(0);
}


bool HGCalTBWireChamberSource::setRunAndEventInfo(edm::EventID& id, edm::TimeValue_t& time, edm::EventAuxiliary::ExperimentType& evType) {	

	if (fileCounter != nextFileIndex || rootTreeIndex == tree->GetEntries()) {		//initial loading of a file

		if (tree!=NULL && rootTreeIndex == tree->GetEntries()) {
			nextFileIndex++;
			fileCounter = -1;
			rootFile->Close();
			std::cout<<"Number of good (DWC) events: "<<goodEventCounter<<" / "<<eventCounter<<std::endl;
			std::cout<<"Number of synchronised events: "<<syncCounter[0]<< " while "<<syncCounter[1]<<" are rejected of which "<<syncCounter[2]<<" with major discrepancy"<<std::endl<<std::endl;
			ref_time_sync = ref_time_dwc = 0;
		}

		if (nextFileIndex == (int)fileNames().size()) {
			return false; 		//end of files is reached
		}

		//set the root file
		fileCounter = nextFileIndex;
		rootTreeIndex = 0;
  		eventCounter = goodEventCounter = 0;
		syncCounter[0] = syncCounter[1] =  syncCounter[2] = 0;

		std::cout<<"Opening "<<fileNames()[fileCounter].c_str()<<std::endl;
		rootFile = new TFile(fileNames()[fileCounter].c_str());	
		
		if (fileFormat==1) {
			tree = (TTree*)rootFile->Get("DelayWireChambers");
			tree->SetBranchAddress("run", &n_run, &b_run);
			tree->SetBranchAddress("event", &n_trigger_tdc, &b_trigger);
			tree->SetBranchAddress("channels", &channels, &b_channels);
			tree->SetBranchAddress("dwc_timestamps", &dwc_timestamps, &b_dwc_timestamps);
			
			if (triggerTimingFormat[fileCounter]==0) {
				tree->SetBranchAddress("timeSinceStart", &timeSinceStart, &b_timeSinceStart);
			} else {
				tree->SetBranchAddress("timeSinceStart", &timeSinceStart_long, &b_timeSinceStart);
			}


			for (int ch=0; ch<N_TDC_channels; ch++) {
				if (hitsPerChannelStored[ch]==1) {
					b_hits[ch] = 0;
					tree->SetBranchAddress(("dwc_hits_ch"+std::to_string(ch)).c_str(), &(hits.at(ch)), &(b_hits.at(ch)));
				}
			}
		} else if (fileFormat==2) {
			tree = (TTree*)rootFile->Get("CAENData");
			tree->SetBranchAddress("run", &n_run, &b_run);
			tree->SetBranchAddress("event", &n_trigger_tdc, &b_trigger);
			tree->SetBranchAddress("tdc_first_hits", &dwc_timestamps, &b_dwc_timestamps);
			tree->SetBranchAddress("trigger_timestamp_cpu_TDC", &timeSinceStart_long, &b_timeSinceStart);
		
			for (int ch=0; ch<N_TDC_channels; ch++) {
				if (hitsPerChannelStored[ch]==1) {
					b_hits[ch] = 0;
					tree->SetBranchAddress(("tdc_hits_ch"+std::to_string(ch)).c_str(), &(hits.at(ch)), &(b_hits.at(ch)));
				}
			}			
		}

		skippedTDCTriggers = 0;
		ReadTimingFile(timingFileNames[fileCounter], sumTriggerTimes[fileCounter]);
	}

	tree->GetEntry(rootTreeIndex);
	rootTreeIndex++;
	

	std::map<std::pair<int, int> ,std::map<int, double> >::iterator alignmentParameterIterator;
	for (alignmentParameterIterator=loadedAlignmentParameters.begin(); alignmentParameterIterator!=loadedAlignmentParameters.end(); alignmentParameterIterator++) {
		int run_min = alignmentParameterIterator->first.first;
		int run_max = alignmentParameterIterator->first.second;
		int this_run = n_run;

		if (this_run>=run_min && (this_run<=run_max || run_max==-1) ) {
			currentAlignmentParameters = alignmentParameterIterator->second;
			break;
		}

	}

	return true;
}


void HGCalTBWireChamberSource::produce(edm::Event & event) {	
	if (fileCounter==-1) return;
	
	std::unique_ptr<std::map<int, WireChamberData> > dwcs(new std::map<int, WireChamberData>);

	//make the wire chambers
	int N_DWC_points = 0;

	//DWC 1
	WireChamberData dwc1;
	dwc1.ID = 1;
	dwc1.recordedTimeStamps=0;
	dwc1.averagedTimeStamp=0;
	dwc1.recordedTimeStamps+=validTimestamp(dwc_timestamps->at(DWC1_LEFT)) ? 1 : 0;
	dwc1.averagedTimeStamp+=validTimestamp(dwc_timestamps->at(DWC1_LEFT)) ? dwc_timestamps->at(DWC1_LEFT) : 0;
	dwc1.recordedTimeStamps+=validTimestamp(dwc_timestamps->at(DWC1_RIGHT)) ? 1 : 0;
	dwc1.averagedTimeStamp+=validTimestamp(dwc_timestamps->at(DWC1_RIGHT)) ? dwc_timestamps->at(DWC1_RIGHT) : 0;
	dwc1.recordedTimeStamps+=validTimestamp(dwc_timestamps->at(DWC1_DOWN)) ? 1 : 0;
	dwc1.averagedTimeStamp+=validTimestamp(dwc_timestamps->at(DWC1_DOWN)) ? dwc_timestamps->at(DWC1_DOWN) : 0;
	dwc1.recordedTimeStamps+=validTimestamp(dwc_timestamps->at(DWC1_UP)) ? 1 : 0;
	dwc1.averagedTimeStamp+=validTimestamp(dwc_timestamps->at(DWC1_UP)) ? dwc_timestamps->at(DWC1_UP) : 0;
	if (dwc1.recordedTimeStamps>0) dwc1.averagedTimeStamp /= dwc1.recordedTimeStamps;

	dwc1.goodMeasurement_X = validTimestamp(dwc_timestamps->at(DWC1_LEFT)) && validTimestamp(dwc_timestamps->at(DWC1_RIGHT));
	dwc1.goodMeasurement_Y = validTimestamp(dwc_timestamps->at(DWC1_DOWN)) && validTimestamp(dwc_timestamps->at(DWC1_UP));
	dwc1.x = dwc1.goodMeasurement_X ? slope_x.at(0) * (dwc_timestamps->at(DWC1_LEFT)-dwc_timestamps->at(DWC1_RIGHT)): -999;
	dwc1.res_x = wc_resolutions[0];
	dwc1.y = dwc1.goodMeasurement_Y ? slope_y.at(0) * (dwc_timestamps->at(DWC1_DOWN)-dwc_timestamps->at(DWC1_UP)): -999;
	dwc1.res_y = wc_resolutions[0];
	dwc1.z = dwc_z1;


	//DWC 2
	WireChamberData dwc2;
	dwc2.ID = 2;
	dwc2.recordedTimeStamps=0;
	dwc2.averagedTimeStamp=0;
	dwc2.recordedTimeStamps+=validTimestamp(dwc_timestamps->at(DWC2_LEFT)) ? 1 : 0;
	dwc2.averagedTimeStamp+=validTimestamp(dwc_timestamps->at(DWC2_LEFT)) ? dwc_timestamps->at(DWC2_LEFT) : 0;
	dwc2.recordedTimeStamps+=validTimestamp(dwc_timestamps->at(DWC2_RIGHT)) ? 1 : 0;
	dwc2.averagedTimeStamp+=validTimestamp(dwc_timestamps->at(DWC2_RIGHT)) ? dwc_timestamps->at(DWC2_RIGHT) : 0;
	dwc2.recordedTimeStamps+=validTimestamp(dwc_timestamps->at(DWC2_DOWN)) ? 1 : 0;
	dwc2.averagedTimeStamp+=validTimestamp(dwc_timestamps->at(DWC2_DOWN)) ? dwc_timestamps->at(DWC2_DOWN) : 0;
	dwc2.recordedTimeStamps+=validTimestamp(dwc_timestamps->at(DWC2_UP)) ? 1 : 0;
	dwc2.averagedTimeStamp+=validTimestamp(dwc_timestamps->at(DWC2_UP)) ? dwc_timestamps->at(DWC2_UP) : 0;
	if (dwc2.recordedTimeStamps>0) dwc2.averagedTimeStamp /= dwc2.recordedTimeStamps;

	dwc2.goodMeasurement_X = validTimestamp(dwc_timestamps->at(DWC2_LEFT)) && validTimestamp(dwc_timestamps->at(DWC2_RIGHT));
	dwc2.goodMeasurement_Y = validTimestamp(dwc_timestamps->at(DWC2_DOWN)) && validTimestamp(dwc_timestamps->at(DWC2_UP));
	dwc2.x = dwc2.goodMeasurement_X ? slope_x.at(1) * (dwc_timestamps->at(DWC2_LEFT)-dwc_timestamps->at(DWC2_RIGHT)): -999;
	dwc2.res_x = wc_resolutions[1];
	dwc2.y = dwc2.goodMeasurement_Y ? slope_y.at(1) * (dwc_timestamps->at(DWC2_DOWN)-dwc_timestamps->at(DWC2_UP)): -999;
	dwc2.res_y = wc_resolutions[1];
	dwc2.z = dwc_z2;


	//DWC 3
	WireChamberData dwc3;
	dwc3.ID = 3;
	dwc3.recordedTimeStamps=0;
	dwc3.averagedTimeStamp=0;
	dwc3.recordedTimeStamps+=validTimestamp(dwc_timestamps->at(DWC3_LEFT)) ? 1 : 0;
	dwc3.averagedTimeStamp+=validTimestamp(dwc_timestamps->at(DWC3_LEFT)) ? dwc_timestamps->at(DWC3_LEFT) : 0;
	dwc3.recordedTimeStamps+=validTimestamp(dwc_timestamps->at(DWC3_RIGHT)) ? 1 : 0;
	dwc3.averagedTimeStamp+=validTimestamp(dwc_timestamps->at(DWC3_RIGHT)) ? dwc_timestamps->at(DWC3_RIGHT) : 0;
	dwc3.recordedTimeStamps+=validTimestamp(dwc_timestamps->at(DWC3_DOWN)) ? 1 : 0;
	dwc3.averagedTimeStamp+=validTimestamp(dwc_timestamps->at(DWC3_DOWN)) ? dwc_timestamps->at(DWC3_DOWN) : 0;
	dwc3.recordedTimeStamps+=validTimestamp(dwc_timestamps->at(DWC3_UP)) ? 1 : 0;
	dwc3.averagedTimeStamp+=validTimestamp(dwc_timestamps->at(DWC3_UP)) ? dwc_timestamps->at(DWC3_UP) : 0;
	if (dwc3.recordedTimeStamps>0) dwc3.averagedTimeStamp /= dwc3.recordedTimeStamps;

	dwc3.goodMeasurement_X = validTimestamp(dwc_timestamps->at(DWC3_LEFT)) && validTimestamp(dwc_timestamps->at(DWC3_RIGHT));
	dwc3.goodMeasurement_Y = validTimestamp(dwc_timestamps->at(DWC3_DOWN)) && validTimestamp(dwc_timestamps->at(DWC3_UP));
	dwc3.x = dwc3.goodMeasurement_X ? slope_x.at(2) * (dwc_timestamps->at(DWC3_LEFT)-dwc_timestamps->at(DWC3_RIGHT)): -999;
	dwc3.res_x = wc_resolutions[2];
	dwc3.y = dwc3.goodMeasurement_Y ? slope_y.at(2) * (dwc_timestamps->at(DWC3_DOWN)-dwc_timestamps->at(DWC3_UP)): -999;
	dwc3.res_y = wc_resolutions[2];
	dwc3.z = dwc_z3;
	
	if ((n_run>=1195) && (n_run<=1333)) {//from run 1195, x-coordinate of DWC A was not connected anymore. TDC channels 14 and 15 were input by trigger signals.
						//also the channels for the y-coordinate of DWCA must have been flipped
		dwc3.goodMeasurement_X = dwc3.goodMeasurement = false;	
		dwc3.x = -999;
		dwc3.y = -dwc3.y;
	} 



	//DWC 4
	WireChamberData dwc4;
	dwc4.ID = 4;
	dwc4.recordedTimeStamps=0;
	dwc4.averagedTimeStamp=0;
	dwc4.recordedTimeStamps+=validTimestamp(dwc_timestamps->at(DWC4_LEFT)) ? 1 : 0;
	dwc4.averagedTimeStamp+=validTimestamp(dwc_timestamps->at(DWC4_LEFT)) ? dwc_timestamps->at(DWC4_LEFT) : 0;
	dwc4.recordedTimeStamps+=validTimestamp(dwc_timestamps->at(DWC4_RIGHT)) ? 1 : 0;
	dwc4.averagedTimeStamp+=validTimestamp(dwc_timestamps->at(DWC4_RIGHT)) ? dwc_timestamps->at(DWC4_RIGHT) : 0;
	dwc4.recordedTimeStamps+=validTimestamp(dwc_timestamps->at(DWC4_DOWN)) ? 1 : 0;
	dwc4.averagedTimeStamp+=validTimestamp(dwc_timestamps->at(DWC4_DOWN)) ? dwc_timestamps->at(DWC4_DOWN) : 0;
	dwc4.recordedTimeStamps+=validTimestamp(dwc_timestamps->at(DWC4_UP)) ? 1 : 0;
	dwc4.averagedTimeStamp+=validTimestamp(dwc_timestamps->at(DWC4_UP)) ? dwc_timestamps->at(DWC4_UP) : 0;
	if (dwc4.recordedTimeStamps>0) dwc4.averagedTimeStamp /= dwc4.recordedTimeStamps;
	
	dwc4.goodMeasurement_X = validTimestamp(dwc_timestamps->at(DWC4_LEFT)) && validTimestamp(dwc_timestamps->at(DWC4_RIGHT));
	dwc4.goodMeasurement_Y = validTimestamp(dwc_timestamps->at(DWC4_DOWN)) && validTimestamp(dwc_timestamps->at(DWC4_UP));
	dwc4.x = dwc4.goodMeasurement_X ? slope_x.at(2) * (dwc_timestamps->at(DWC4_LEFT)-dwc_timestamps->at(DWC4_RIGHT)): -999;
	dwc4.res_x = wc_resolutions[3];
	dwc4.y = dwc4.goodMeasurement_Y ? slope_y.at(2) * (dwc_timestamps->at(DWC4_DOWN)-dwc_timestamps->at(DWC4_UP)): -999;
	dwc4.res_y = wc_resolutions[3];
	dwc4.z = dwc_z4;


	//store multiplicities:
	int N(0), sumHits(0);
	if (hitsPerChannelStored[DWC1_LEFT]==1) {N++;sumHits+=hits.at(DWC1_LEFT)->size();}
	if (hitsPerChannelStored[DWC1_RIGHT]==1) {N++;sumHits+=hits.at(DWC1_RIGHT)->size();}
	if (hitsPerChannelStored[DWC1_DOWN]==1) {N++;sumHits+=hits.at(DWC1_DOWN)->size();}
	if (hitsPerChannelStored[DWC1_UP]==1) {N++;sumHits+=hits.at(DWC1_UP)->size();}
	dwc1.averageHitMultiplicty = (N > 0) ? sumHits*1.*1./N : -999; 

	N = 0; sumHits = 0;
	if (hitsPerChannelStored[DWC2_LEFT]==1) {N++;sumHits+=hits.at(DWC2_LEFT)->size();}
	if (hitsPerChannelStored[DWC2_RIGHT]==1) {N++;sumHits+=hits.at(DWC2_RIGHT)->size();}
	if (hitsPerChannelStored[DWC2_DOWN]==1) {N++;sumHits+=hits.at(DWC2_DOWN)->size();}
	if (hitsPerChannelStored[DWC2_UP]==1) {N++;sumHits+=hits.at(DWC2_UP)->size();}
	dwc2.averageHitMultiplicty = (N > 0) ? sumHits*1./N : -999; 

	N = 0; sumHits = 0;
	if (hitsPerChannelStored[DWC3_LEFT]==1) {N++;sumHits+=hits.at(DWC3_LEFT)->size();}
	if (hitsPerChannelStored[DWC3_RIGHT]==1) {N++;sumHits+=hits.at(DWC3_RIGHT)->size();}
	if (hitsPerChannelStored[DWC3_DOWN]==1) {N++;sumHits+=hits.at(DWC3_DOWN)->size();}
	if (hitsPerChannelStored[DWC3_UP]==1) {N++;sumHits+=hits.at(DWC3_UP)->size();}
	dwc3.averageHitMultiplicty = (N > 0) ? sumHits*1./N : -999; 

	N = 0; sumHits = 0;
	if (hitsPerChannelStored[DWC4_LEFT]==1) {N++;sumHits+=hits.at(DWC4_LEFT)->size();}
	if (hitsPerChannelStored[DWC4_RIGHT]==1) {N++;sumHits+=hits.at(DWC4_RIGHT)->size();}
	if (hitsPerChannelStored[DWC4_DOWN]==1) {N++;sumHits+=hits.at(DWC4_DOWN)->size();}
	if (hitsPerChannelStored[DWC4_UP]==1) {N++;sumHits+=hits.at(DWC4_UP)->size();}
	dwc4.averageHitMultiplicty = (N > 0) ? sumHits*1./N : -999; 

	//perform alignment
	if (fabs(dwc1.x) <= 50.) dwc1.x = dwc1.x - currentAlignmentParameters[11] - dwc1.y*currentAlignmentParameters[21]; else dwc1.goodMeasurement_X=false;
	if (fabs(dwc1.y) <= 50.) dwc1.y = dwc1.y - currentAlignmentParameters[12] +  dwc1.x*currentAlignmentParameters[21]; else dwc1.goodMeasurement_Y=false;
	if (fabs(dwc2.x) <= 50.) dwc2.x = dwc2.x - currentAlignmentParameters[111] - dwc2.y*currentAlignmentParameters[121]; else dwc2.goodMeasurement_X=false;
	if (fabs(dwc2.y) <= 50.) dwc2.y = dwc2.y - currentAlignmentParameters[112] +  dwc2.x*currentAlignmentParameters[121]; else dwc2.goodMeasurement_Y=false;
	if (fabs(dwc3.x) <= 50.) dwc3.x = dwc3.x - currentAlignmentParameters[211] - dwc3.y*currentAlignmentParameters[221]; else dwc3.goodMeasurement_X=false;
	if (fabs(dwc3.y) <= 50.) dwc3.y = dwc3.y - currentAlignmentParameters[212] +  dwc3.x*currentAlignmentParameters[221]; else dwc3.goodMeasurement_Y=false;
	if (fabs(dwc4.x) <= 50.) dwc4.x = dwc4.x - currentAlignmentParameters[311] - dwc4.y*currentAlignmentParameters[321]; else dwc4.goodMeasurement_X=false;
	if (fabs(dwc4.y) <= 50.) dwc4.y = dwc4.y - currentAlignmentParameters[312] +  dwc4.x*currentAlignmentParameters[321]; else dwc4.goodMeasurement_Y=false;

	//set the good measurement flags
	dwc1.goodMeasurement = (dwc1.goodMeasurement_X && dwc1.goodMeasurement_Y);
	N_DWC_points = dwc1.goodMeasurement ? N_DWC_points+1 : N_DWC_points;
	dwc2.goodMeasurement = (dwc2.goodMeasurement_X && dwc2.goodMeasurement_Y);
	N_DWC_points = dwc2.goodMeasurement ? N_DWC_points+1 : N_DWC_points;
	dwc3.goodMeasurement = (dwc3.goodMeasurement_X && dwc3.goodMeasurement_Y);
	N_DWC_points = dwc3.goodMeasurement ? N_DWC_points+1 : N_DWC_points;
	dwc4.goodMeasurement = (dwc4.goodMeasurement_X && dwc4.goodMeasurement_Y);
	N_DWC_points = dwc4.goodMeasurement ? N_DWC_points+1 : N_DWC_points;


	(*dwcs)[0] = dwc1;
	(*dwcs)[1] = dwc2;
	(*dwcs)[2] = dwc3;
	(*dwcs)[3] = dwc4;

	event.put(std::move(dwcs), outputCollectionName);		

	bool oneHit = false;
	for (int index=0; index<N_TDC_channels; index++)
		oneHit = oneHit || validTimestamp(dwc_timestamps->at(index));

	#ifdef DEBUG
		if (!oneHit) {
			std::cout<<"!!!!!!!!!!!!!!!"<<std::endl;
		}
	#endif

	//add the RunData
	std::unique_ptr<RunData> rd(new RunData);

	int n_trigger_orm = n_trigger_tdc-skipFirstNEvents[fileCounter]+skippedTDCTriggers;

	rd->configuration = setupIDs[fileCounter];
  	rd->energy = beamEnergies[fileCounter];
  	rd->pdgID = pdgIDs[fileCounter];
  	rd->runType = HGCAL_TB_BEAM;


	rd->run = n_run;
	rd->trigger = n_trigger_orm;
	rd->booleanUserRecords.add("hasValidDWCMeasurement", (N_DWC_points>=3));		//need two points for track extrapolation
	
	//do the matching to the HGCal events
	if (trigger_to_event_table.count(n_trigger_orm)==0) {
		rd->event=-1;
	} else {
		int event_candidate_index = trigger_to_event_table[n_trigger_orm];

		double timeSinceStart_ms = timeSinceStart;
		if (triggerTimingFormat[fileCounter]==1) timeSinceStart_ms = timeSinceStart_long * TDCTriggerTimeStampConversionToMs;

		double deltaTs = (event_trigger_time[event_candidate_index]-ref_time_sync) - (timeSinceStart_ms - ref_time_dwc);
		rd->doubleUserRecords.add("triggerDeltaT_to_TDC", deltaTs);

		if (allowForTDCEventSkipping[fileCounter]&&(deltaTs>-100.)&&(deltaTs<-15.)) {		
		//average time in between two events is around 20ms given by the sync board. So cutting on -15. is reasonable for this configuration (20 Oct 2017 in H6A)
			std::cout<<"Skipping one TDC trigger"<<std::endl;
			skippedTDCTriggers+=1;
		}

		#ifdef DEBUG
			std::cout<<"Event: "<<event_candidate_index<<"  tdc trigger: "<<n_trigger_tdc<<"  orm trigger: "<<n_trigger_orm<<": "<<deltaTs<<" = "<<(event_trigger_time[event_candidate_index]-ref_time_sync)<<" - "<<(timeSinceStart_ms - ref_time_dwc)<<std::endl;
		#endif

		//use absolute value for comparison
		deltaTs = fabs(deltaTs);

		if ((deltaTs <= triggerTimeDifferenceTolerance) ||  sumTriggerTimes[fileCounter]==-1){
			rd->event = event_candidate_index;

			if (rd->booleanUserRecords.has("hasValidDWCMeasurement")&&rd->booleanUserRecords.get("hasValidDWCMeasurement")) goodEventCounter++;
			syncCounter[0]++;
		} else {
			rd->event=-1;
			syncCounter[1]++;
		}

		if (deltaTs > 200. && deltaTs<100000 && eventCounter > 3 && sumTriggerTimes[fileCounter]!=-1) {			//the first line is not used to test synchronisation since a time offset is present for the two streams.
			//throw cms::Exception("EventAsynch") << "Trigger time interval differs by more than 200ms. The files are likely not synchronised. ";
			#ifdef DEBUG
				std::cout<<std::endl<<std::endl << "Trigger time interval differs by more than 200ms. The files are likely not synchronised. "<<std::endl<<std::endl<<std::endl;
			#endif
			syncCounter[2]++;
		} else {			
			ref_time_sync = event_trigger_time[event_candidate_index]; 
			ref_time_dwc = timeSinceStart_ms;	
		}
	
		delta_T_priorDWCTrigger = deltaTs;
		eventCounter++;
	}
	
	event.put(std::move(rd), "RunData");

}


void HGCalTBWireChamberSource::ReadTimingFile(std::string timingFilePath, bool sumTriggerTimes) {
	trigger_to_event_table.clear();
	event_trigger_time.clear();
		//Todo
	std::fstream file; 
	char fragment[100];
	int readCounter = -5, currentEvent = 0;

	file.open(timingFilePath.c_str(), std::fstream::in);

	double time_sinceStart = 0;

	std::cout<<"Reading file "<<timingFilePath<<" -open: "<<file.is_open()<<std::endl;
	while (file.is_open() && !file.eof()) {
		readCounter++;
		file >> fragment;
		if (readCounter==0) currentEvent++;
		if (readCounter==1) {
			trigger_to_event_table[atoi(fragment)-triggerCountOffsets[fileCounter]] = currentEvent;
		}
		if (readCounter==3) {
			std::istringstream reader(fragment);
			long time;					//must read this value as long since it is given in ns
			reader >> time;
			time_sinceStart =  sumTriggerTimes ? time_sinceStart + ((double) 25e-6 * time) : ((double) 25e-6 * time);
			event_trigger_time[currentEvent] = time_sinceStart; 		//one time stamp is 25ns --> conversion to ms for comparison.
			readCounter = -1;
		}
	}

	/*
	#ifdef DEBUG
		std::map<int, int>::iterator it;
		for (it=event_trigger_time.begin(); it!=event_trigger_time.end(); it++) {
			std::cout<<it->first<<"  -  "<<it->second<<std::endl;
		}
	#endif
	*/
	
	
}

void HGCalTBWireChamberSource::ReadAlignmentParameters(int fileIndex) {
  	if (fileIndex==(int)alignmentParamaterFiles.size()) return;

	std::fstream file; 
	char fragment[100];
	int readCounter = -2, currentParameter = 0;

	std::map<int, double> _parameters;

	if (readCounter==-2) {
		for (int i=0; i< 4; i++) {
		  _parameters[i*100+11] = 0.;
		  _parameters[i*100+12] = 0.;
		  _parameters[i*100+21] = 0.;
		}
	}	

	if (alignmentParamaterFiles[fileIndex]!=""){
		std::cout<<"Opening: "<<alignmentParamaterFiles[fileIndex]<<std::endl;
		file.open(alignmentParamaterFiles[fileIndex].c_str(), std::fstream::in);
	}

	int minRun, maxRun;
	if (file.is_open()) {
		file >> fragment;
		minRun = atoi(fragment);
		file >> fragment;
		maxRun = atoi(fragment);
	}

	while (file.is_open() && !file.eof()) {
		
		if (readCounter!=-2) readCounter++;
			file >> fragment;

		if (std::string(fragment)=="11") readCounter = 0;  //first parameter is read out

		if (readCounter==0) currentParameter = atoi(fragment);
		if (readCounter==1) currentParameter = _parameters[currentParameter] = atof(fragment); 
		if (readCounter==2) if (atof(fragment)==-1.) readCounter = -1;
		if (readCounter==4) readCounter = -1;
	}

	#ifdef DEBUG
		std::cout<<"Min run: "<<minRun<<"   Max run: "<<maxRun<<std::endl;
		for (int i=0; i<4; i++) {
		  std::cout<<"Alignment parameter: "<<i*100+11<<": "<<_parameters[i*100+11]<<std::endl;
		  std::cout<<"Alignment parameter: "<<i*100+12<<": "<<_parameters[i*100+12]<<std::endl;
		  std::cout<<"Alignment parameter: "<<i*100+21<<": "<<_parameters[i*100+21]<<std::endl;
		  std::cout<<"Alignment parameter: "<<i*100+22<<": "<<_parameters[i*100+22]<<std::endl;
		}
	#endif 
	loadedAlignmentParameters[std::make_pair(minRun, maxRun)] = _parameters;

	return ReadAlignmentParameters(fileIndex+1);
}

#include "FWCore/Framework/interface/InputSourceMacros.h"

DEFINE_FWK_INPUT_SOURCE(HGCalTBWireChamberSource);
