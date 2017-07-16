#include "HGCal/RawToDigi/plugins/HGCalTBWireChamberSource.h"




HGCalTBWireChamberSource::HGCalTBWireChamberSource(const edm::ParameterSet & pset, edm::InputSourceDescription const& desc) :  edm::ProducerSourceFromFiles(pset, desc, true),
	rootFile(NULL)
{
	//find and fill the configured runs
	outputCollectionName = pset.getParameter<std::string>("OutputCollectionName");

  	std::vector<double> v0(4, 0.2/40);		//0.2mm/ns and TDC binning is 25ps
	slope_x = pset.getUntrackedParameter<std::vector<double> >("slope_x", v0);
	slope_y = pset.getUntrackedParameter<std::vector<double> >("slope_y", v0);
	performAlignment = pset.getUntrackedParameter<bool>("performAlignment", false);
	alignmentParamaterFile = pset.getUntrackedParameter<std::string>("alignmentParamaterFile", "");


	produces<WireChambers>("WireChambers");

	tree = NULL;

	//values from the 
	n_run=0;
	n_trigger=0;
	channels=0;
	dwc_timestamps=0;
}

void HGCalTBWireChamberSource::beginJob() {
	fileCounter = -1;
	eventCounter = 0;
	nextFileIndex = 0;
	rootFile = NULL;
	tree = NULL;

	if (performAlignment) ReadAlignmentParameters();

}


bool HGCalTBWireChamberSource::setRunAndEventInfo(edm::EventID& id, edm::TimeValue_t& time, edm::EventAuxiliary::ExperimentType& evType) {	
	if (nextFileIndex == (int)fileNames().size()) {

		return false; 		//end of files is reached
	}

	if (fileCounter != nextFileIndex) {		//initial loading of a file
		fileCounter = nextFileIndex;
		eventCounter = 0;
  
		std::cout<<"Opening "<<fileNames()[fileCounter].c_str()<<std::endl;
		rootFile = new TFile(fileNames()[fileCounter].c_str());	
		tree = (TTree*)rootFile->Get("DelayWireChambers");

		tree->SetBranchAddress("run", &n_run, &b_run);
		tree->SetBranchAddress("event", &n_trigger, &b_trigger);
		tree->SetBranchAddress("channels", &channels, &b_channels);
		tree->SetBranchAddress("dwc_timestamps", &dwc_timestamps, &b_dwc_timestamps);
	}


	if (eventCounter == tree->GetEntries()) {
		nextFileIndex++;
		fileCounter = -1;
		rootFile->Close();
		setRunAndEventInfo(id, time, evType);
	}
	
	tree->GetEntry(eventCounter);
	eventCounter++;

	return true;
}


void HGCalTBWireChamberSource::produce(edm::Event & event) {	
	std::auto_ptr<WireChambers> mwcs(new WireChambers);	

	//make the wire chambers

	//DWC 1
	WireChamberData* dwc1 = new WireChamberData();
	dwc1->ID = 1;
	dwc1->recordedTimeStamps=0;
	dwc1->averagedTimeStamp=0;
	dwc1->recordedTimeStamps+=(dwc_timestamps->at(DWC1_LEFT)>=0) ? 1 : 0;
	dwc1->averagedTimeStamp+=(dwc_timestamps->at(DWC1_LEFT)>=0) ? dwc_timestamps->at(DWC1_LEFT) : 0;
	dwc1->recordedTimeStamps+=(dwc_timestamps->at(DWC1_RIGHT)>=0) ? 1 : 0;
	dwc1->averagedTimeStamp+=(dwc_timestamps->at(DWC1_RIGHT)>=0) ? dwc_timestamps->at(DWC1_RIGHT) : 0;
	dwc1->recordedTimeStamps+=(dwc_timestamps->at(DWC1_DOWN)>=0) ? 1 : 0;
	dwc1->averagedTimeStamp+=(dwc_timestamps->at(DWC1_DOWN)>=0) ? dwc_timestamps->at(DWC1_DOWN) : 0;
	dwc1->recordedTimeStamps+=(dwc_timestamps->at(DWC1_UP)>=0) ? 1 : 0;
	dwc1->averagedTimeStamp+=(dwc_timestamps->at(DWC1_UP)>=0) ? dwc_timestamps->at(DWC1_UP) : 0;
	if (dwc1->recordedTimeStamps>0) dwc1->averagedTimeStamp /= dwc1->recordedTimeStamps;

	dwc1->x = slope_x.at(0) * (dwc_timestamps->at(DWC1_LEFT)-dwc_timestamps->at(DWC1_RIGHT));
	dwc1->y = slope_y.at(0) * (dwc_timestamps->at(DWC1_DOWN)-dwc_timestamps->at(DWC1_UP));
	dwc1->z = dwc_z1;
	dwc1->goodMeasurement_X = ((dwc_timestamps->at(DWC1_LEFT) >= 0) && (dwc_timestamps->at(DWC1_RIGHT) >=0));
	dwc1->goodMeasurement_Y = ((dwc_timestamps->at(DWC1_DOWN) >= 0) && (dwc_timestamps->at(DWC1_UP) >=0));
	dwc1->goodMeasurement = (dwc1->goodMeasurement_X && dwc1->goodMeasurement_Y);

	mwcs->push_back(*dwc1);


	//DWC 2
	WireChamberData* dwc2 = new WireChamberData();
	dwc2->ID = 2;
	dwc2->recordedTimeStamps=0;
	dwc2->averagedTimeStamp=0;
	dwc2->recordedTimeStamps+=(dwc_timestamps->at(DWC2_LEFT)>=0) ? 1 : 0;
	dwc2->averagedTimeStamp+=(dwc_timestamps->at(DWC2_LEFT)>=0) ? dwc_timestamps->at(DWC2_LEFT) : 0;
	dwc2->recordedTimeStamps+=(dwc_timestamps->at(DWC2_RIGHT)>=0) ? 1 : 0;
	dwc2->averagedTimeStamp+=(dwc_timestamps->at(DWC2_RIGHT)>=0) ? dwc_timestamps->at(DWC2_RIGHT) : 0;
	dwc2->recordedTimeStamps+=(dwc_timestamps->at(DWC2_DOWN)>=0) ? 1 : 0;
	dwc2->averagedTimeStamp+=(dwc_timestamps->at(DWC2_DOWN)>=0) ? dwc_timestamps->at(DWC2_DOWN) : 0;
	dwc2->recordedTimeStamps+=(dwc_timestamps->at(DWC2_UP)>=0) ? 1 : 0;
	dwc2->averagedTimeStamp+=(dwc_timestamps->at(DWC2_UP)>=0) ? dwc_timestamps->at(DWC1_UP) : 0;
	if (dwc2->recordedTimeStamps>0) dwc2->averagedTimeStamp /= dwc2->recordedTimeStamps;

	dwc2->x = slope_x.at(1) * (dwc_timestamps->at(DWC2_LEFT)-dwc_timestamps->at(DWC2_RIGHT));
	dwc2->y = slope_y.at(1) * (dwc_timestamps->at(DWC2_DOWN)-dwc_timestamps->at(DWC2_UP));
	dwc2->z = dwc_z2;
	dwc2->goodMeasurement_X = ((dwc_timestamps->at(DWC2_LEFT) >= 0) && (dwc_timestamps->at(DWC2_RIGHT) >=0));
	dwc2->goodMeasurement_Y = ((dwc_timestamps->at(DWC2_DOWN) >= 0) && (dwc_timestamps->at(DWC2_UP) >=0));
	dwc2->goodMeasurement = (dwc2->goodMeasurement_X && dwc2->goodMeasurement_Y);

	mwcs->push_back(*dwc2);


	//DWC 3
	WireChamberData* dwc3 = new WireChamberData();
	dwc3->ID = 3;
	dwc3->recordedTimeStamps=0;
	dwc3->averagedTimeStamp=0;
	dwc3->recordedTimeStamps+=(dwc_timestamps->at(DWC3_LEFT)>=0) ? 1 : 0;
	dwc3->averagedTimeStamp+=(dwc_timestamps->at(DWC3_LEFT)>=0) ? dwc_timestamps->at(DWC3_LEFT) : 0;
	dwc3->recordedTimeStamps+=(dwc_timestamps->at(DWC3_RIGHT)>=0) ? 1 : 0;
	dwc3->averagedTimeStamp+=(dwc_timestamps->at(DWC3_RIGHT)>=0) ? dwc_timestamps->at(DWC3_RIGHT) : 0;
	dwc3->recordedTimeStamps+=(dwc_timestamps->at(DWC3_DOWN)>=0) ? 1 : 0;
	dwc3->averagedTimeStamp+=(dwc_timestamps->at(DWC3_DOWN)>=0) ? dwc_timestamps->at(DWC3_DOWN) : 0;
	dwc3->recordedTimeStamps+=(dwc_timestamps->at(DWC3_UP)>=0) ? 1 : 0;
	dwc3->averagedTimeStamp+=(dwc_timestamps->at(DWC3_UP)>=0) ? dwc_timestamps->at(DWC1_UP) : 0;
	if (dwc3->recordedTimeStamps>0) dwc3->averagedTimeStamp /= dwc3->recordedTimeStamps;

	dwc3->x = slope_x.at(2) * (dwc_timestamps->at(DWC3_LEFT)-dwc_timestamps->at(DWC3_RIGHT));
	dwc3->y = slope_y.at(2) * (dwc_timestamps->at(DWC3_DOWN)-dwc_timestamps->at(DWC3_UP));
	dwc3->z = dwc_z3;
	dwc3->goodMeasurement_X = ((dwc_timestamps->at(DWC3_LEFT) >= 0) && (dwc_timestamps->at(DWC3_RIGHT) >=0));
	dwc3->goodMeasurement_Y = ((dwc_timestamps->at(DWC3_DOWN) >= 0) && (dwc_timestamps->at(DWC3_UP) >=0));
	dwc3->goodMeasurement = (dwc3->goodMeasurement_X && dwc3->goodMeasurement_Y);

	mwcs->push_back(*dwc3);


	//DWC 4
	WireChamberData* dwc4 = new WireChamberData();
	dwc4->ID = 4;
	dwc4->recordedTimeStamps=0;
	dwc4->averagedTimeStamp=0;
	dwc4->recordedTimeStamps+=(dwc_timestamps->at(DWC4_LEFT)>=0) ? 1 : 0;
	dwc4->averagedTimeStamp+=(dwc_timestamps->at(DWC4_LEFT)>=0) ? dwc_timestamps->at(DWC4_LEFT) : 0;
	dwc4->recordedTimeStamps+=(dwc_timestamps->at(DWC4_RIGHT)>=0) ? 1 : 0;
	dwc4->averagedTimeStamp+=(dwc_timestamps->at(DWC4_RIGHT)>=0) ? dwc_timestamps->at(DWC4_RIGHT) : 0;
	dwc4->recordedTimeStamps+=(dwc_timestamps->at(DWC4_DOWN)>=0) ? 1 : 0;
	dwc4->averagedTimeStamp+=(dwc_timestamps->at(DWC4_DOWN)>=0) ? dwc_timestamps->at(DWC4_DOWN) : 0;
	dwc4->recordedTimeStamps+=(dwc_timestamps->at(DWC4_UP)>=0) ? 1 : 0;
	dwc4->averagedTimeStamp+=(dwc_timestamps->at(DWC4_UP)>=0) ? dwc_timestamps->at(DWC1_UP) : 0;
	if (dwc4->recordedTimeStamps>0) dwc4->averagedTimeStamp /= dwc4->recordedTimeStamps;
	
	dwc4->x = slope_x.at(2) * (dwc_timestamps->at(DWC4_LEFT)-dwc_timestamps->at(DWC4_RIGHT));
	dwc4->y = slope_y.at(2) * (dwc_timestamps->at(DWC4_DOWN)-dwc_timestamps->at(DWC4_UP));
	dwc4->z = dwc_z3;
	dwc4->goodMeasurement_X = ((dwc_timestamps->at(DWC4_LEFT) >= 0) && (dwc_timestamps->at(DWC4_RIGHT) >=0));
	dwc4->goodMeasurement_Y = ((dwc_timestamps->at(DWC4_DOWN) >= 0) && (dwc_timestamps->at(DWC4_UP) >=0));
	dwc4->goodMeasurement = (dwc4->goodMeasurement_X && dwc4->goodMeasurement_Y);

	mwcs->push_back(*dwc4);


	if (performAlignment) {
		dwc1->x = dwc1->x * alignmentParameters[121] - alignmentParameters[111];
		dwc1->y = dwc1->y * alignmentParameters[122] - alignmentParameters[112];
		dwc2->x = dwc2->x * alignmentParameters[221] - alignmentParameters[211];
		dwc2->y = dwc2->y * alignmentParameters[222] - alignmentParameters[212];
		dwc3->x = dwc3->x * alignmentParameters[321] - alignmentParameters[311];
		dwc3->y = dwc3->y * alignmentParameters[322] - alignmentParameters[312];
		dwc4->x = dwc4->x * alignmentParameters[421] - alignmentParameters[411];
		dwc4->y = dwc4->y * alignmentParameters[422] - alignmentParameters[412];
	}




	event.put(std::move(mwcs), "WireChambers");		

}


void HGCalTBWireChamberSource::ReadAlignmentParameters() {
  
  std::fstream file; 
  char fragment[100];
  int readCounter = -2, currentParameter = 0;
  
  if (readCounter==-2) {
    for (int i=1; i<= 4; i++) {
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
    if (std::string(fragment)=="111") readCounter = 0;  //first parameter is read out

    if (readCounter==0) currentParameter = atoi(fragment);
    if (readCounter==1) currentParameter = alignmentParameters[currentParameter] = atof(fragment); 
    if (readCounter==2) if (atof(fragment)==-1.) readCounter = -1;
    if (readCounter==4) readCounter = -1;
  }

for (int i=1; i<= 4; i++) {
  std::cout<<"Alignment parameter: "<<i*100+11<<": "<<alignmentParameters[i*100+11]<<std::endl;
  std::cout<<"Alignment parameter: "<<i*100+12<<": "<<alignmentParameters[i*100+12]<<std::endl;
  std::cout<<"Alignment parameter: "<<i*100+21<<": "<<alignmentParameters[i*100+21]<<std::endl;
  std::cout<<"Alignment parameter: "<<i*100+22<<": "<<alignmentParameters[i*100+22]<<std::endl;
}

}

#include "FWCore/Framework/interface/InputSourceMacros.h"

DEFINE_FWK_INPUT_SOURCE(HGCalTBWireChamberSource);
