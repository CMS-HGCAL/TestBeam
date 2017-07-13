#include "HGCal/RawToDigi/plugins/HGCalTBWireChamberSource.h"




HGCalTBWireChamberSource::HGCalTBWireChamberSource(const edm::ParameterSet & pset, edm::InputSourceDescription const& desc) :  edm::ProducerSourceFromFiles(pset, desc, true),
	rootFile(NULL)
{
	//find and fill the configured runs
	outputCollectionName = pset.getParameter<std::string>("OutputCollectionName");

  	std::vector<double> v0(4, 0.2);
	slope_x = pset.getUntrackedParameter<std::vector<double> >("slope_x", v0);
	slope_y = pset.getUntrackedParameter<std::vector<double> >("slope_y", v0);

	produces<WireChambers>("WireChambers");

	tree = NULL;

	n_run=0;
	n_event=0;
	channels=0;
	dwc_timestamps=0;
}

void HGCalTBWireChamberSource::beginJob() {

	rootFile = new TFile(fileNames()[0].c_str());	
	tree = (TTree*)rootFile->Get("DelayWireChambers");

	tree->SetBranchAddress("run", &n_run, &b_run);
	tree->SetBranchAddress("event", &n_event, &b_event);
	tree->SetBranchAddress("channels", &channels, &b_channels);
	tree->SetBranchAddress("dwc_timestamps", &dwc_timestamps, &b_dwc_timestamps);
}


bool HGCalTBWireChamberSource::setRunAndEventInfo(edm::EventID& id, edm::TimeValue_t& time, edm::EventAuxiliary::ExperimentType& type) {
	
	if (n_event == tree->GetEntries()) {
		return false;
	}
	
	tree->GetEntry(n_event);
	
	return true;
}


void HGCalTBWireChamberSource::produce(edm::Event & event) {	
	std::auto_ptr<WireChambers> mwcs(new WireChambers);	

	//make the wire chambers
	WireChamberData* dwc1 = new WireChamberData();
	dwc1->ID = 1;
	dwc1->x = slope_x.at(0) * (dwc_timestamps->at(DWC1_LEFT)-dwc_timestamps->at(DWC1_RIGHT));
	dwc1->y = slope_y.at(0) * (dwc_timestamps->at(DWC1_DOWN)-dwc_timestamps->at(DWC1_UP));
	dwc1->z = dwc_z1;
	dwc1->goodMeasurement_X = ((dwc_timestamps->at(DWC1_LEFT) >= 0) && (dwc_timestamps->at(DWC1_RIGHT) >=0));
	dwc1->goodMeasurement_Y = ((dwc_timestamps->at(DWC1_DOWN) >= 0) && (dwc_timestamps->at(DWC1_UP) >=0));
	dwc1->goodMeasurement = (dwc1->goodMeasurement_X && dwc1->goodMeasurement_Y);
	mwcs->push_back(*dwc1);


	WireChamberData* dwc2 = new WireChamberData();
	dwc2->ID = 2;
	dwc2->x = slope_x.at(1) * (dwc_timestamps->at(DWC2_LEFT)-dwc_timestamps->at(DWC2_RIGHT));
	dwc2->y = slope_y.at(1) * (dwc_timestamps->at(DWC2_DOWN)-dwc_timestamps->at(DWC2_UP));
	dwc2->z = dwc_z2;
	dwc2->goodMeasurement_X = ((dwc_timestamps->at(DWC2_LEFT) >= 0) && (dwc_timestamps->at(DWC2_RIGHT) >=0));
	dwc2->goodMeasurement_Y = ((dwc_timestamps->at(DWC2_DOWN) >= 0) && (dwc_timestamps->at(DWC2_UP) >=0));
	dwc2->goodMeasurement = (dwc2->goodMeasurement_X && dwc2->goodMeasurement_Y);
	mwcs->push_back(*dwc2);


	WireChamberData* dwc3 = new WireChamberData();
	dwc3->ID = 3;
	dwc3->x = slope_x.at(2) * (dwc_timestamps->at(DWC3_LEFT)-dwc_timestamps->at(DWC3_RIGHT));
	dwc3->y = slope_y.at(2) * (dwc_timestamps->at(DWC3_DOWN)-dwc_timestamps->at(DWC3_UP));
	dwc3->z = dwc_z3;
	dwc3->goodMeasurement_X = ((dwc_timestamps->at(DWC3_LEFT) >= 0) && (dwc_timestamps->at(DWC3_RIGHT) >=0));
	dwc3->goodMeasurement_Y = ((dwc_timestamps->at(DWC3_DOWN) >= 0) && (dwc_timestamps->at(DWC3_UP) >=0));
	dwc3->goodMeasurement = (dwc3->goodMeasurement_X && dwc3->goodMeasurement_Y);
	mwcs->push_back(*dwc3);


	WireChamberData* dwc4 = new WireChamberData();
	dwc4->ID = 4;
	dwc4->x = slope_x.at(2) * (dwc_timestamps->at(DWC4_LEFT)-dwc_timestamps->at(DWC4_RIGHT));
	dwc4->y = slope_y.at(2) * (dwc_timestamps->at(DWC4_DOWN)-dwc_timestamps->at(DWC4_UP));
	dwc4->z = dwc_z3;
	dwc4->goodMeasurement_X = ((dwc_timestamps->at(DWC4_LEFT) >= 0) && (dwc_timestamps->at(DWC4_RIGHT) >=0));
	dwc4->goodMeasurement_Y = ((dwc_timestamps->at(DWC4_DOWN) >= 0) && (dwc_timestamps->at(DWC4_UP) >=0));
	dwc4->goodMeasurement = (dwc4->goodMeasurement_X && dwc4->goodMeasurement_Y);
	mwcs->push_back(*dwc4);

	/*
	if (dwc1->goodMeasurement && dwc2->goodMeasurement && dwc3->goodMeasurement && dwc4->goodMeasurement) {
		for (size_t i=0; i<16; i++) 
			std::cout<<dwc_timestamps->at(i)<<" ";
		std::cout<<std::endl;
	}
	*/

	event.put(std::move(mwcs), "WireChambers");		

}

void HGCalTBWireChamberSource::endJob() {
	delete tree;
	delete rootFile;
}


#include "FWCore/Framework/interface/InputSourceMacros.h"

DEFINE_FWK_INPUT_SOURCE(HGCalTBWireChamberSource);
