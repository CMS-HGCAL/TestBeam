#include "HGCal/RawToDigi/plugins/HGCalTBWireChamberSource.h"

//#define DEBUG


bool validTimestamp(int ts) {
	return (ts>=0);
}

HGCalTBWireChamberSource::HGCalTBWireChamberSource(const edm::ParameterSet & pset, edm::InputSourceDescription const& desc) :  edm::ProducerSourceFromFiles(pset, desc, true),
	rootFile(NULL)
{
	//find and fill the configured runs
	outputCollectionName = pset.getParameter<std::string>("OutputCollectionName");

  	std::vector<double> v0(4, 0.2/40);		//0.2mm/ns and TDC binning is 25ps
	triggerTimeDifferenceTolerance = pset.getUntrackedParameter<int>("triggerTimeDifferenceTolerance", 2);
	
	slope_x = pset.getUntrackedParameter<std::vector<double> >("slope_x", v0);
	slope_y = pset.getUntrackedParameter<std::vector<double> >("slope_y", v0);
	wc_resolution = pset.getUntrackedParameter<double>("wc_resolution", 0.5);

	performAlignment = pset.getUntrackedParameter<bool>("performAlignment", false);
	alignmentParamaterFile = pset.getUntrackedParameter<std::string>("alignmentParamaterFile", "");

	timingFileNames = pset.getParameter<std::vector<std::string> >("timingFileNames");
	sumTriggerTimes = pset.getParameter<std::vector<int> >("sumTriggerTimes");
	triggerCountOffsets = pset.getParameter<std::vector<int> >("triggerCountOffsets");
	skipFirstNEvents = pset.getParameter<std::vector<int> >("skipFirstNEvents");
	runType = pset.getParameter<std::vector<std::string> >("runType");


	produces<WireChambers>("WireChambers");
	produces<RunData>("RunData");

	tree = NULL;

	//values from the 
	n_run=0;
	n_trigger=0;
	channels=0;
	dwc_timestamps=0;

}

void HGCalTBWireChamberSource::beginJob() {
	fileCounter = -1;
	rootTreeIndex = 0;
	nextFileIndex = 0;
	rootFile = NULL;
	tree = NULL;

	if (performAlignment) ReadAlignmentParameters();
}


bool HGCalTBWireChamberSource::setRunAndEventInfo(edm::EventID& id, edm::TimeValue_t& time, edm::EventAuxiliary::ExperimentType& evType) {	


	if (fileCounter != nextFileIndex || rootTreeIndex == tree->GetEntries()) {		//initial loading of a file

		if (tree!=NULL && rootTreeIndex == tree->GetEntries()) {
			nextFileIndex++;
			fileCounter = -1;
			rootFile->Close();
			std::cout<<"Number of good (DWC) events: "<<goodEventCounter<<" / "<<eventCounter<<std::endl<<std::endl;
			ref_time_sync = ref_time_dwc = 0;
		}

		if (nextFileIndex == (int)fileNames().size()) {
			return false; 		//end of files is reached
		}

		//set the root file
		fileCounter = nextFileIndex;
		rootTreeIndex = 0;
  		eventCounter = goodEventCounter = 0;

		std::cout<<"Opening "<<fileNames()[fileCounter].c_str()<<std::endl;
		std::cout<<"Run type "<<runType[fileCounter]<<std::endl;
		rootFile = new TFile(fileNames()[fileCounter].c_str());	
		tree = (TTree*)rootFile->Get("DelayWireChambers");

		tree->SetBranchAddress("run", &n_run, &b_run);
		tree->SetBranchAddress("event", &n_trigger, &b_trigger);
		tree->SetBranchAddress("channels", &channels, &b_channels);
		tree->SetBranchAddress("dwc_timestamps", &dwc_timestamps, &b_dwc_timestamps);
		tree->SetBranchAddress("timeSinceStart", &timeSinceStart, &b_timeSinceStart);

		ReadTimingFile(timingFileNames[fileCounter], sumTriggerTimes[fileCounter]);
	}

	tree->GetEntry(rootTreeIndex);
	rootTreeIndex++;
	
	return true;
}


void HGCalTBWireChamberSource::produce(edm::Event & event) {	
	if (fileCounter==-1) return;
	
	std::auto_ptr<WireChambers> mwcs(new WireChambers);	

	//make the wire chambers
	int N_DWC_points = 0;


	//DWC 1
	WireChamberData* dwc1 = new WireChamberData();
	dwc1->ID = 1;
	dwc1->recordedTimeStamps=0;
	dwc1->averagedTimeStamp=0;
	dwc1->recordedTimeStamps+=validTimestamp(dwc_timestamps->at(DWC1_LEFT)) ? 1 : 0;
	dwc1->averagedTimeStamp+=validTimestamp(dwc_timestamps->at(DWC1_LEFT)) ? dwc_timestamps->at(DWC1_LEFT) : 0;
	dwc1->recordedTimeStamps+=validTimestamp(dwc_timestamps->at(DWC1_RIGHT)) ? 1 : 0;
	dwc1->averagedTimeStamp+=validTimestamp(dwc_timestamps->at(DWC1_RIGHT)) ? dwc_timestamps->at(DWC1_RIGHT) : 0;
	dwc1->recordedTimeStamps+=validTimestamp(dwc_timestamps->at(DWC1_DOWN)) ? 1 : 0;
	dwc1->averagedTimeStamp+=validTimestamp(dwc_timestamps->at(DWC1_DOWN)) ? dwc_timestamps->at(DWC1_DOWN) : 0;
	dwc1->recordedTimeStamps+=validTimestamp(dwc_timestamps->at(DWC1_UP)) ? 1 : 0;
	dwc1->averagedTimeStamp+=validTimestamp(dwc_timestamps->at(DWC1_UP)) ? dwc_timestamps->at(DWC1_UP) : 0;
	if (dwc1->recordedTimeStamps>0) dwc1->averagedTimeStamp /= dwc1->recordedTimeStamps;

	dwc1->goodMeasurement_X = validTimestamp(dwc_timestamps->at(DWC1_LEFT)) && validTimestamp(dwc_timestamps->at(DWC1_RIGHT));
	dwc1->goodMeasurement_Y = validTimestamp(dwc_timestamps->at(DWC1_DOWN)) && validTimestamp(dwc_timestamps->at(DWC1_UP));
	dwc1->goodMeasurement = (dwc1->goodMeasurement_X && dwc1->goodMeasurement_Y);
	dwc1->x = dwc1->goodMeasurement_X ? slope_x.at(0) * (dwc_timestamps->at(DWC1_LEFT)-dwc_timestamps->at(DWC1_RIGHT)): -999;
	dwc1->res_x = wc_resolution;
	dwc1->y = dwc1->goodMeasurement_Y ? slope_y.at(0) * (dwc_timestamps->at(DWC1_DOWN)-dwc_timestamps->at(DWC1_UP)): -999;
	dwc1->res_y = wc_resolution;
	dwc1->z = dwc_z1;
	N_DWC_points = dwc1->goodMeasurement ? N_DWC_points+1 : N_DWC_points;


	//DWC 2
	WireChamberData* dwc2 = new WireChamberData();
	dwc2->ID = 2;
	dwc2->recordedTimeStamps=0;
	dwc2->averagedTimeStamp=0;
	dwc2->recordedTimeStamps+=validTimestamp(dwc_timestamps->at(DWC2_LEFT)) ? 1 : 0;
	dwc2->averagedTimeStamp+=validTimestamp(dwc_timestamps->at(DWC2_LEFT)) ? dwc_timestamps->at(DWC2_LEFT) : 0;
	dwc2->recordedTimeStamps+=validTimestamp(dwc_timestamps->at(DWC2_RIGHT)) ? 1 : 0;
	dwc2->averagedTimeStamp+=validTimestamp(dwc_timestamps->at(DWC2_RIGHT)) ? dwc_timestamps->at(DWC2_RIGHT) : 0;
	dwc2->recordedTimeStamps+=validTimestamp(dwc_timestamps->at(DWC2_DOWN)) ? 1 : 0;
	dwc2->averagedTimeStamp+=validTimestamp(dwc_timestamps->at(DWC2_DOWN)) ? dwc_timestamps->at(DWC2_DOWN) : 0;
	dwc2->recordedTimeStamps+=validTimestamp(dwc_timestamps->at(DWC2_UP)) ? 1 : 0;
	dwc2->averagedTimeStamp+=validTimestamp(dwc_timestamps->at(DWC2_UP)) ? dwc_timestamps->at(DWC1_UP) : 0;
	if (dwc2->recordedTimeStamps>0) dwc2->averagedTimeStamp /= dwc2->recordedTimeStamps;

	dwc2->goodMeasurement_X = validTimestamp(dwc_timestamps->at(DWC2_LEFT)) && validTimestamp(dwc_timestamps->at(DWC2_RIGHT));
	dwc2->goodMeasurement_Y = validTimestamp(dwc_timestamps->at(DWC2_DOWN)) && validTimestamp(dwc_timestamps->at(DWC2_UP));
	dwc2->goodMeasurement = (dwc2->goodMeasurement_X && dwc2->goodMeasurement_Y);
	dwc2->x = dwc2->goodMeasurement_X ? slope_x.at(1) * (dwc_timestamps->at(DWC2_LEFT)-dwc_timestamps->at(DWC2_RIGHT)): -999;
	dwc2->res_x = wc_resolution;
	dwc2->y = dwc2->goodMeasurement_Y ? slope_y.at(1) * (dwc_timestamps->at(DWC2_DOWN)-dwc_timestamps->at(DWC2_UP)): -999;
	dwc2->res_y = wc_resolution;
	dwc2->z = dwc_z2;
	N_DWC_points = dwc2->goodMeasurement ? N_DWC_points+1 : N_DWC_points;


	//DWC 3
	WireChamberData* dwc3 = new WireChamberData();
	dwc3->ID = 3;
	dwc3->recordedTimeStamps=0;
	dwc3->averagedTimeStamp=0;
	dwc3->recordedTimeStamps+=validTimestamp(dwc_timestamps->at(DWC3_LEFT)) ? 1 : 0;
	dwc3->averagedTimeStamp+=validTimestamp(dwc_timestamps->at(DWC3_LEFT)) ? dwc_timestamps->at(DWC3_LEFT) : 0;
	dwc3->recordedTimeStamps+=validTimestamp(dwc_timestamps->at(DWC3_RIGHT)) ? 1 : 0;
	dwc3->averagedTimeStamp+=validTimestamp(dwc_timestamps->at(DWC3_RIGHT)) ? dwc_timestamps->at(DWC3_RIGHT) : 0;
	dwc3->recordedTimeStamps+=validTimestamp(dwc_timestamps->at(DWC3_DOWN)) ? 1 : 0;
	dwc3->averagedTimeStamp+=validTimestamp(dwc_timestamps->at(DWC3_DOWN)) ? dwc_timestamps->at(DWC3_DOWN) : 0;
	dwc3->recordedTimeStamps+=validTimestamp(dwc_timestamps->at(DWC3_UP)) ? 1 : 0;
	dwc3->averagedTimeStamp+=validTimestamp(dwc_timestamps->at(DWC3_UP)) ? dwc_timestamps->at(DWC1_UP) : 0;
	if (dwc3->recordedTimeStamps>0) dwc3->averagedTimeStamp /= dwc3->recordedTimeStamps;

	dwc3->goodMeasurement_X = validTimestamp(dwc_timestamps->at(DWC3_LEFT)) && validTimestamp(dwc_timestamps->at(DWC3_RIGHT));
	dwc3->goodMeasurement_Y = validTimestamp(dwc_timestamps->at(DWC3_DOWN)) && validTimestamp(dwc_timestamps->at(DWC3_UP));
	dwc3->goodMeasurement = (dwc3->goodMeasurement_X && dwc3->goodMeasurement_Y);
	dwc3->x = dwc3->goodMeasurement_X ? slope_x.at(2) * (dwc_timestamps->at(DWC3_LEFT)-dwc_timestamps->at(DWC3_RIGHT)): -999;
	dwc3->res_x = wc_resolution;
	dwc3->y = dwc3->goodMeasurement_Y ? slope_y.at(2) * (dwc_timestamps->at(DWC3_DOWN)-dwc_timestamps->at(DWC3_UP)): -999;
	dwc3->res_y = wc_resolution;
	dwc3->z = dwc_z3;
	N_DWC_points = dwc3->goodMeasurement ? N_DWC_points+1 : N_DWC_points;


	//DWC 4
	WireChamberData* dwc4 = new WireChamberData();
	dwc4->ID = 4;
	dwc4->recordedTimeStamps=0;
	dwc4->averagedTimeStamp=0;
	dwc4->recordedTimeStamps+=validTimestamp(dwc_timestamps->at(DWC4_LEFT)) ? 1 : 0;
	dwc4->averagedTimeStamp+=validTimestamp(dwc_timestamps->at(DWC4_LEFT)) ? dwc_timestamps->at(DWC4_LEFT) : 0;
	dwc4->recordedTimeStamps+=validTimestamp(dwc_timestamps->at(DWC4_RIGHT)) ? 1 : 0;
	dwc4->averagedTimeStamp+=validTimestamp(dwc_timestamps->at(DWC4_RIGHT)) ? dwc_timestamps->at(DWC4_RIGHT) : 0;
	dwc4->recordedTimeStamps+=validTimestamp(dwc_timestamps->at(DWC4_DOWN)) ? 1 : 0;
	dwc4->averagedTimeStamp+=validTimestamp(dwc_timestamps->at(DWC4_DOWN)) ? dwc_timestamps->at(DWC4_DOWN) : 0;
	dwc4->recordedTimeStamps+=validTimestamp(dwc_timestamps->at(DWC4_UP)) ? 1 : 0;
	dwc4->averagedTimeStamp+=validTimestamp(dwc_timestamps->at(DWC4_UP)) ? dwc_timestamps->at(DWC1_UP) : 0;
	if (dwc4->recordedTimeStamps>0) dwc4->averagedTimeStamp /= dwc4->recordedTimeStamps;
	
	dwc4->goodMeasurement_X = validTimestamp(dwc_timestamps->at(DWC4_LEFT)) && validTimestamp(dwc_timestamps->at(DWC4_RIGHT));
	dwc4->goodMeasurement_Y = validTimestamp(dwc_timestamps->at(DWC4_DOWN)) && validTimestamp(dwc_timestamps->at(DWC4_UP));
	dwc4->goodMeasurement = (dwc4->goodMeasurement_X && dwc4->goodMeasurement_Y);
	dwc4->x = dwc4->goodMeasurement_X ? slope_x.at(2) * (dwc_timestamps->at(DWC4_LEFT)-dwc_timestamps->at(DWC4_RIGHT)): -999;
	dwc4->res_x = wc_resolution;
	dwc4->y = dwc4->goodMeasurement_Y ? slope_y.at(2) * (dwc_timestamps->at(DWC4_DOWN)-dwc_timestamps->at(DWC4_UP)): -999;
	dwc4->res_y = wc_resolution;
	dwc4->z = dwc_z4;
	N_DWC_points = dwc4->goodMeasurement ? N_DWC_points+1 : N_DWC_points;


	if (performAlignment) {
		dwc1->x = dwc1->x - alignmentParameters[11];
		dwc1->y = dwc1->y - alignmentParameters[12];
		dwc2->x = dwc2->x - alignmentParameters[111];
		dwc2->y = dwc2->y - alignmentParameters[112];
		dwc3->x = dwc3->x - alignmentParameters[211];
		dwc3->y = dwc3->y - alignmentParameters[212];
		dwc4->x = dwc4->x - alignmentParameters[311];
		dwc4->y = dwc4->y - alignmentParameters[312];
	}

	mwcs->push_back(*dwc1);
	mwcs->push_back(*dwc2);
	mwcs->push_back(*dwc3);
	mwcs->push_back(*dwc4);

	event.put(std::move(mwcs), "WireChambers");		



	//add the RunData
	std::auto_ptr<RunData> rd(new RunData);

	int n_trigger_corrected = n_trigger-skipFirstNEvents[fileCounter];

	rd->configuration = -1;
	rd->runType = runType[fileCounter];
	switch(atoi(runType[fileCounter].c_str()) % 10) {
		case 1: rd->energy=100; break;
		case 2: rd->energy=150; break;
		case 3: rd->energy=200; break;
		case 4: rd->energy=250; break;
		case 5: rd->energy=300; break;
		default: rd->energy=-1;
	}


	rd->run = n_run;
	rd->trigger = n_trigger_corrected;
	rd->hasDanger = false;
	rd->hasValidMWCMeasurement = (N_DWC_points>=3);		//need two points for track extrapolation
	
	//do the matching to the HGCal events
	if (trigger_to_event_table.count(n_trigger_corrected)==0) {
		rd->event=-1;
	} else {
		int event_candidate_index = trigger_to_event_table[n_trigger_corrected];

		double deltaTs = fabs((event_trigger_time[event_candidate_index]-ref_time_sync) - (timeSinceStart - ref_time_dwc));

		#ifdef DEBUG
			std::cout<<"Event: "<<event_candidate_index<<"  trigger: "<<n_trigger<<": "<<deltaTs<<std::endl;
		#endif

		if (deltaTs < triggerTimeDifferenceTolerance) {
			rd->event = event_candidate_index;
			#ifdef DEBUG
				std::cout<<"good x1: "<<dwc1->goodMeasurement_X<<"   good y1: "<<dwc1->goodMeasurement_Y<<"     good x2: "<<dwc2->goodMeasurement_X<<"   good y2: "<<dwc2->goodMeasurement_Y;
				std::cout<<"   good x3: "<<dwc3->goodMeasurement_X<<"   good y3: "<<dwc3->goodMeasurement_Y<<"     good x4: "<<dwc4->goodMeasurement_X<<"   good y4: "<<dwc4->goodMeasurement_Y;
				std::cout<<std::endl;
			#endif
			if (rd->hasValidMWCMeasurement) goodEventCounter++;
		} else {
			rd->event=-1;
		}

		if (deltaTs > 200. && event_candidate_index != 1) {			//the first line is not used to test synchronisation since a time offset is present for the two streams.
			throw cms::Exception("EventAsynch") << "Trigger time interval differs by more than 200ms. The files are likely not synchronised. ";
		}

		ref_time_sync = event_trigger_time[event_candidate_index]; 
		ref_time_dwc = timeSinceStart;	
		


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

	int time_sinceStart = 0;

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

	#ifdef DEBUG
		std::map<int, int>::iterator it;
		for (it=event_trigger_time.begin(); it!=event_trigger_time.end(); it++) {
			std::cout<<it->first<<"  -  "<<it->second<<std::endl;
		}
	#endif
	
	
	
}

void HGCalTBWireChamberSource::ReadAlignmentParameters() {
  
	std::fstream file; 
	char fragment[100];
	int readCounter = -2, currentParameter = 0;

	if (readCounter==-2) {
		for (int i=0; i< 4; i++) {
		  alignmentParameters[i*100+11] = 0.;
		  alignmentParameters[i*100+12] = 0.;
		  alignmentParameters[i*100+21] = 1.;
		  alignmentParameters[i*100+22] = 1.;
		}
	}	

	if (alignmentParamaterFile!="")
		file.open(alignmentParamaterFile.c_str(), std::fstream::in);


	while (file.is_open() && !file.eof()) {
		if (readCounter!=-2) readCounter++;
			file >> fragment;
		if (std::string(fragment)=="11") readCounter = 0;  //first parameter is read out

		if (readCounter==0) currentParameter = atoi(fragment);
		if (readCounter==1) currentParameter = alignmentParameters[currentParameter] = atof(fragment); 
		if (readCounter==2) if (atof(fragment)==-1.) readCounter = -1;
		if (readCounter==4) readCounter = -1;
	}

	for (int i=0; i<4; i++) {
	  std::cout<<"Alignment parameter: "<<i*100+11<<": "<<alignmentParameters[i*100+11]<<std::endl;
	  std::cout<<"Alignment parameter: "<<i*100+12<<": "<<alignmentParameters[i*100+12]<<std::endl;
	  std::cout<<"Alignment parameter: "<<i*100+21<<": "<<alignmentParameters[i*100+21]<<std::endl;
	  std::cout<<"Alignment parameter: "<<i*100+22<<": "<<alignmentParameters[i*100+22]<<std::endl;
	}

}

#include "FWCore/Framework/interface/InputSourceMacros.h"

DEFINE_FWK_INPUT_SOURCE(HGCalTBWireChamberSource);
