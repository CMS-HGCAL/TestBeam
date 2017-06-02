#include "HGCal/Reco/plugins/HGCalTBRecHitProducer_2017.h"
#include <iostream>

using namespace std;

#define DEBUG 0

HGCalTBRecHitProducer_2017::HGCalTBRecHitProducer_2017(const edm::ParameterSet& cfg) : outputCollectionName(cfg.getParameter<std::string>("OutputCollectionName")),
_rawHitsToken(consumes<HGCalTBRawHitCollection>(cfg.getParameter<edm::InputTag>("rawHitCollection"))),
_pedestalLow_filename(cfg.getParameter<std::string>("pedestalLow")),		//not used yet
_pedestalHigh_filename(cfg.getParameter<std::string>("pedestalHigh")),		//not used yet
_gainsLow_filename(cfg.getParameter<std::string>("gainLow")),				//not used yet	
_gainsHigh_filename(cfg.getParameter<std::string>("gainHigh")),				//not used yet
_adcSaturation(cfg.getParameter<int>("adcSaturation")),
_layers_config(cfg.getParameter<int>("layers_config")) ,
_N_timestamps(11) {
	
	produces <HGCalTBRecHitCollection>(outputCollectionName);

	if(_layers_config == 1){
		_LG2HG_value = cfg.getParameter<std::vector<double> >("LG2HG_May2017");
		_mapFile = cfg.getParameter<std::string>("mapFile_CERN_May2017");
	} else{
		throw cms::Exception("Not defined configuration");
	}

	HGCalCondObjectTextIO io(0);
	edm::FileInPath fip(_mapFile);
	if (!io.load(fip.fullPath(), essource_.emap_)) {
		throw cms::Exception("Unable to load electronics map");
	};
}

void HGCalTBRecHitProducer_2017::produce(edm::Event& event, const edm::EventSetup& iSetup){

	std::auto_ptr<HGCalTBRecHitCollection> rechits(new HGCalTBRecHitCollection); //auto_ptr are automatically deleted when out of scope

	edm::Handle<HGCalTBRawHitCollection> rawHitsHandle;
	event.getByToken(_rawHitsToken, rawHitsHandle);


	HGCalCondGains adcToGeV_low, adcToGeV_high;
	HGCalCondPedestals pedestals_low, pedestals_high;
	HGCalCondObjectTextIO condIO(HGCalTBNumberingScheme::scheme());

	if(_pedestalLow_filename != "") assert(condIO.load(_pedestalLow_filename, pedestals_low));
	if(_pedestalHigh_filename != "") 	assert(condIO.load(_pedestalHigh_filename, pedestals_high));
	if(_gainsLow_filename != "") 	assert(condIO.load(_gainsLow_filename, adcToGeV_low));
	if(_gainsHigh_filename != "") 	assert(condIO.load(_gainsHigh_filename, adcToGeV_high));



	//Core Assumptions (01. June 2017):
	//only one layer exsits
	//common mode is given as baseline from all other disconnected channels

	//Step 1: compute the baseline of high- and low-gain ADC counts
	double ave_lowADC_nonConnected = 0.;
	double ave_highADC_nonConnected = 0.;
	int sum_nonConnected = 0.;

	for( auto hit : *rawHitsHandle ) {
  		if (hit.detid().layer() != 127) continue;		//hits with layer index 127 are not connected
		if (DEBUG) std::cout << "[RECHIT PRODUCER (2017): raw hit] " << hit.detid()<<"   skiroc: "<< hit.skiroc() << "   channel: " << hit.channel() << std::endl;
  		for (int t_stamp=0; t_stamp<11; t_stamp++) {
  			if (DEBUG) std::cout<<"Time stamp: "<<t_stamp<<"   low: "<<hit.lowGainADC(t_stamp)<<"  high: "<<hit.highGainADC(t_stamp)<<std::endl;
  			ave_highADC_nonConnected += hit.highGainADC(t_stamp);
  			ave_lowADC_nonConnected += hit.lowGainADC(t_stamp);
			sum_nonConnected++;  			
  		}
	}
	ave_lowADC_nonConnected/=sum_nonConnected;
	ave_highADC_nonConnected/=sum_nonConnected;

	if (DEBUG) std::cout<<"Baseline for low gain: "<<ave_lowADC_nonConnected<<"    high gain: "<<ave_highADC_nonConnected<<std::endl;

	
	//Step 2: reconstruction algorithm
	for( auto hit : *rawHitsHandle ) {
  		int layer = hit.detid().layer();
  		int skiroc = hit.skiroc();
  		int channel = hit.channel();

  		if (layer == 127) continue;		//hits with layer index 127 are not connected

		double ADC_sum_low = 0.;
		double ADC_sum_high = 0.;
		double ADC_sum_low_baseline_subtracted = 0.;
		double ADC_sum_high_baseline_subtracted = 0.;

  		for (int t_stamp=0; t_stamp<_N_timestamps; t_stamp++) {
  			double low_ADC = hit.lowGainADC(t_stamp);
  			double high_ADC = hit.highGainADC(t_stamp);	
  			if (DEBUG) std::cout<<"Time stamp: "<<t_stamp<<"   low: "<<low_ADC<<"  high: "<<high_ADC<<std::endl;
  			ADC_sum_low += low_ADC;
  			ADC_sum_high += high_ADC;
  			ADC_sum_low_baseline_subtracted += (low_ADC > ave_lowADC_nonConnected) ? low_ADC - ave_lowADC_nonConnected : 0;
  			ADC_sum_high_baseline_subtracted += (high_ADC > ave_highADC_nonConnected) ? high_ADC - ave_highADC_nonConnected : 0;
  		}

	
  		HGCalTBRecHit recHit(hit.detid(), ADC_sum_high_baseline_subtracted, ADC_sum_low, ADC_sum_high, 0);
  		

  		CellCenterXY = TheCell.GetCellCentreCoordinatesForPlots((recHit.id()).layer(), (recHit.id()).sensorIU(), (recHit.id()).sensorIV(), (recHit.id()).iu(), (recHit.id()).iv(), 128); 
  		recHit.setCellCenterCoordinate(CellCenterXY.first, CellCenterXY.second);

  		if (DEBUG) std::cout<<"layer: "<<(recHit.id()).layer()<<"   skiroc: "<<skiroc<<"    channel: "<<channel<<"    x: "<<CellCenterXY.first<<"    y: "<<CellCenterXY.second<<"    ADC_sum_high_baseline_subtracted: "<<ADC_sum_high_baseline_subtracted<<std::endl;

  		rechits->push_back(recHit);
	}

	event.put(rechits, outputCollectionName);

}

DEFINE_FWK_MODULE(HGCalTBRecHitProducer_2017);
