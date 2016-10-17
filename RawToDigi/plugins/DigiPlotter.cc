// -*- C++ -*-
//
// Package:    HGCal/DigiPlotter
// Class:      DigiPlotter
//
/**\class DigiPlotter DigiPlotter.cc HGCal/DigiPlotter/plugins/DigiPlotter.cc

   Description: Plugin to make 2D and 1D histograms of digis
   
   Implementation:
   \author Rajdeep Mohan Chatterjee
   \author Shervin Nourbakhsh
*/
//
// Original Author:  Rajdeep Mohan Chatterjee
//         Created:  Mon, 15 Feb 2016 09:47:43 GMT
//         Modified by Shervin Nourbakhsh
//
//


// system include files
#include <memory>
#include <iostream>
#include "TH2Poly.h"
#include "TH1F.h"
#include "TProfile.h"
#include <fstream>
// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/one/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "HGCal/DataFormats/interface/HGCalTBDetId.h"
#include "HGCal/Geometry/interface/HGCalTBCellVertices.h"
#include "HGCal/Geometry/interface/HGCalTBCellParameters.h"
#include "HGCal/Geometry/interface/HGCalTBTopology.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "HGCal/CondObjects/interface/HGCalElectronicsMap.h"
#include "HGCal/CondObjects/interface/HGCalCondObjectTextIO.h"
#include "HGCal/DataFormats/interface/HGCalTBElectronicsId.h"
#include "HGCal/DataFormats/interface/HGCalTBDataFrameContainers.h"
#include "HGCal/Geometry/interface/HGCalTBGeometryParameters.h"


//
// class declaration
//

// If the analyzer does not use TFileService, please remove
// the template argument to the base class so the class inherits
// from  edm::one::EDAnalyzer<> and also remove the line from
// constructor "usesResource("TFileService");"
// This will improve performance in multithreaded jobs.

class DigiPlotter : public edm::one::EDAnalyzer<edm::one::SharedResources>
{

public:
	explicit DigiPlotter(const edm::ParameterSet&);
	~DigiPlotter();
	static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);
private:
	virtual void beginJob() override;
	void analyze(const edm::Event& , const edm::EventSetup&) override;
	virtual void endJob() override;
	// ----------member data ---------------------------
	HGCalTBTopology IsCellValid;
	HGCalTBCellVertices TheCell;
	std::string mapfile_ = "HGCal/CondObjects/data/map_CERN_8Layers_Sept2016.txt";
	struct {
		HGCalElectronicsMap emap_;
	} essource_;
	int sensorsize = 128;// The geometry for a 256 cell sensor hasnt been implemted yet. Need a picture to do this.

	std::vector<std::pair<double, double>>::const_iterator it;
	TH2Poly *h_digi_layer[SKIROC::MAXSAMPLES][MAXLAYERS];
	TH1F    *h_digi_layer_summed[SKIROC::MAXSAMPLES][MAXLAYERS];
	TProfile    *h_digi_layer_profile[SKIROC::MAXSAMPLES][MAXLAYERS];
	int Sensor_Iu = 0;
	int Sensor_Iv = 0;
	TH2F* Noise_2D_Profile[SKIROC::MAXSAMPLES][MAXLAYERS];
	TH1F  *h_digi_layer_channel[MAXSKIROCS][64][SKIROC::MAXSAMPLES];
//        TH1F  *h_digi_layer_cell_event[SKIROC::MAXSAMPLES][MAXLAYERS][cellx][celly][512];
	char name[50], title[50];
	double ADC_Sum_SKI_Layer[2][MAXLAYERS][2]; // 2 SKIROCs per layer, High gain and low gain ADC HARD CODED
	int Cell_Count_SKI_Layer[2][4]; // 2 SKIROCs per layer, High gain and low gain ADC HARD CODED
	std::string m_pedestalsHighGain;
	std::string m_pedestalsLowGain;
};

//
// constants, enums and typedefs
//

//
// static data member definitions
//

//
// constructors and destructor
//
DigiPlotter::DigiPlotter(const edm::ParameterSet& iConfig)
{
	//now do what ever initialization is needed
	usesResource("TFileService");
	edm::Service<TFileService> fs;
	consumesMany<SKIROC2DigiCollection>();

	double CellXs[HGCAL_TB_CELL::FullHexVertices] = {0.};
	double CellYs[HGCAL_TB_CELL::FullHexVertices] = {0.};

	for(int i_sample = 0; i_sample < SKIROC::MAXSAMPLES; i_sample++) {
		for(int nlayers = 0; nlayers < MAXLAYERS; nlayers++) {
			sprintf(name, "Noise_2D_Profile_ADC%i_Layer%i", i_sample, nlayers);
			sprintf(title, "Noise 2D Profile ADC%i Layer%i", i_sample, nlayers);
			Noise_2D_Profile[i_sample][nlayers] = fs->make<TH2F>(name, title, 128, 0, 127, 2000, -1000, 1000);
		}
	}

	for(int ISkiroc = 1; ISkiroc <= MAXSKIROCS; ISkiroc++) {
		for(int Channel = 0; Channel < SKIROC::NCHANNELS; Channel++) {
			for(int iii = 0; iii < SKIROC::MAXSAMPLES; iii++) {
				sprintf(name, "Ski_%i_Channel_%i_ADC%i", ISkiroc, Channel, iii);
				sprintf(title, "Ski %i Channel %i ADC%i", ISkiroc, Channel, iii);
				h_digi_layer_channel[ISkiroc - 1][Channel][iii] = fs->make<TH1F>(name, title, 4096, 0., 4095.);
			}
		}
	}

	for(int i_sample = 0; i_sample < SKIROC::MAXSAMPLES; i_sample++) {// one histogram for each sample (high and low gain)
		for(int nlayers = 0; nlayers < MAXLAYERS; nlayers++) {
//Booking a "hexagonal" histograms to display the sum of Digis for SKIROC::MAXSAMPLES, in 1 SKIROC in 1 layer. To include all layers soon. Also the 1D Digis per cell in a sensor is booked here for SKIROC::MAXSAMPLES.
			sprintf(name, "FullLayer_ADC%i_Layer%i", i_sample, nlayers + 1);
			sprintf(title, "Sum of adc counts per cell for ADC%i Layer%i", i_sample, nlayers + 1);
			h_digi_layer[i_sample][nlayers] = fs->make<TH2Poly>();
			h_digi_layer[i_sample][nlayers]->SetName(name);
			h_digi_layer[i_sample][nlayers]->SetTitle(title);
			sprintf(name, "FullLayer_ADC%i_Layer%i_summed", i_sample, nlayers + 1);
			sprintf(title, "Sum of adc counts for all cells in ADC%i Layer%i", i_sample, nlayers + 1);
			h_digi_layer_summed[i_sample][nlayers] = fs->make<TH1F>(name, title, 4096, 0., 4095.);
			h_digi_layer_summed[i_sample][nlayers]->GetXaxis()->SetTitle("Digis[adc counts]");
			sprintf(name, "FullLayer_ADC%i_Layer%i_profile", i_sample, nlayers + 1);
			sprintf(title, "profile of adc counts for all cells in ADC%i Layer%i", i_sample, nlayers + 1);
			h_digi_layer_profile[i_sample][nlayers] = fs->make<TProfile>(name, title, 128, 0, 127, 0., 4095.);
			h_digi_layer_profile[i_sample][nlayers]->GetXaxis()->SetTitle("Channel #");
			h_digi_layer_profile[i_sample][nlayers]->GetYaxis()->SetTitle("ADC counts");


			for(int iv = -7; iv < 8; iv++) {
				for(int iu = -7; iu < 8; iu++) {
					if(!IsCellValid.iu_iv_valid(nlayers, Sensor_Iu, Sensor_Iv, iu, iv, sensorsize)) continue;
					std::vector<std::pair<double, double>> CellXY = TheCell.GetCellCoordinatesForPlots(nlayers, Sensor_Iu, Sensor_Iv, iu, iv, sensorsize);
					size_t vertex_idx = 0;
					for(it = CellXY.begin(); it != CellXY.end(); it++) {
						CellXs[vertex_idx]   =  it->first;
						CellYs[vertex_idx++] =  it->second;
					}
//Somehow cloning of the TH2Poly was not working. Need to look at it. Currently physically booked another one.
					h_digi_layer[i_sample][nlayers]->AddBin(CellXY.size(), CellXs, CellYs);
					

				}//loop over iu
			}//loop over iv
		}//loop over nlayers
	}//loop over i_samples
	
	m_pedestalsHighGain = iConfig.getUntrackedParameter<std::string>("pedestalsHighGain", "");
	m_pedestalsLowGain = iConfig.getUntrackedParameter<std::string>("pedestalsLowGain", "");
}//contructor ends here


DigiPlotter::~DigiPlotter()
{

	// do anything here that needs to be done at desctruction time
	// (e.g. close files, deallocate resources etc.)

}


//
// member functions
//

// ------------ method called for each event  ------------
void
DigiPlotter::analyze(const edm::Event& event, const edm::EventSetup& setup)
{
	using namespace edm;
	std::vector<edm::Handle<SKIROC2DigiCollection> > ski;
	event.getManyByType(ski);
//        int Event = event.id().event();
	/*
	        for(int ski=0;ski<2;ski++){
	            for(int layers =0; layers<MAXLAYERS; layers++){
	                Cell_Count_SKI_Layer[ski][layers] = 0;
	                for(int samples =0; samples<2; samples++){
	                    ADC_Sum_SKI_Layer[ski][layers][samples] = 0.; // 2 SKIROCs per layer, High gain and low gain ADC HARD CODED
	                    }
	               }
	           }
	*/
	if(!ski.empty()) {

		std::vector<edm::Handle<SKIROC2DigiCollection> >::iterator i;
		int counter1 = 0, counter2 = 0;
		for(i = ski.begin(); i != ski.end(); i++) {
			const SKIROC2DigiCollection& Coll = *(*i);

//////////////////////////////////Evaluate average pedestal per event to subtract out//////////////////////////////////
			for(SKIROC2DigiCollection::const_iterator k = Coll.begin(); k != Coll.end(); k++) {
				const SKIROC2DataFrame& SKI = *k ;
				const HGCalTBDetId& detId = SKI.detid();
				uint32_t EID = essource_.emap_.detId2eid(SKI.detid());
				HGCalTBElectronicsId eid(EID);
				if(!IsCellValid.iu_iv_valid(detId, sensorsize))  continue;
;
//                                ADC_Sum_SKI_Layer[eid.iskiroc() - 2*(n_layer - 1) - 1][n_layer - 1][1] += SKI_1[0].adcHigh();
//                                ADC_Sum_SKI_Layer[eid.iskiroc() - 2*(n_layer - 1) - 1][n_layer - 1][0] += SKI_1[0].adcLow();
//                                Cell_Count_SKI_Layer[eid.iskiroc() - 2*(n_layer - 1) - 1][n_layer - 1] += 1;
			}



			for(SKIROC2DigiCollection::const_iterator j = Coll.begin(); j != Coll.end(); j++) {
				const SKIROC2DataFrame& SKI = *j ;
				const HGCalTBDetId& detId = SKI.detid();
				uint32_t EID = essource_.emap_.detId2eid(SKI.detid());
				HGCalTBElectronicsId eid(EID);
#ifdef DEBUG
				std::cout << std::endl << " Layer = " << detId.layer() << " Sensor IU = " << detId.sensorIU() << " Sensor IV = " << detId.sensorIV() << " Cell iu = " << detId.iu() << " Cell iu = " << detId.iv() << std::endl;
#endif
				if(!IsCellValid.iu_iv_valid(detId, sensorsize))  continue;
				std::pair<double, double> CellCentreXY = TheCell.GetCellCentreCoordinatesForPlots(detId, sensorsize);
				double iux = (CellCentreXY.first < 0 ) ? (CellCentreXY.first + 0.0001) : (CellCentreXY.first - 0.0001) ;
				double iyy = (CellCentreXY.second < 0 ) ? (CellCentreXY.second + 0.0001) : (CellCentreXY.second - 0.0001);

				int nsample = 0;
				h_digi_layer[nsample][detId.layer() - 1]->Fill(iux , iyy, SKI[nsample].adcLow());
				h_digi_layer_profile[nsample][detId.layer() - 1]->Fill(counter1++, SKI[nsample].adcLow(), 1);
				if(eid.iskiroc() > 0)	h_digi_layer_channel[eid.iskiroc() - 1][eid.ichan()][nsample]->Fill(SKI[nsample].adcLow());
				nsample = 1;
				h_digi_layer[nsample][detId.layer() - 1]->Fill(iux , iyy, SKI[nsample - 1].adcHigh());
				h_digi_layer_profile[nsample][detId.layer() - 1]->Fill(counter2++, SKI[nsample - 1].adcHigh(), 1);
				if(eid.iskiroc() > 0)	h_digi_layer_channel[eid.iskiroc() - 1][eid.ichan()][nsample]->Fill(SKI[nsample - 1].adcHigh());
			}
			
		}
	} else {
		edm::LogWarning("DQM") << "No SKIROC2 Digis";
	}

}//analyze method ends here


// ------------ method called once each job just before starting event loop  ------------
void
DigiPlotter::beginJob()
{
	HGCalCondObjectTextIO io(0);
	edm::FileInPath fip(mapfile_);
	if (!io.load(fip.fullPath(), essource_.emap_)) {
		throw cms::Exception("Unable to load electronics map");
	}
}

// ------------ method called once each job just after ending the event loop  ------------
void
DigiPlotter::endJob()
{
	int Code = 0;
	int SENSOR_IX = 0;
	int SENSOR_IV = 0;
	std::ofstream f_highGain, f_lowGain;

	f_highGain.open(m_pedestalsHighGain.c_str());
	f_highGain << "SCHEME_CODE 0" << std::endl;
	f_highGain << "# CODE  LAYER SENSOR_IX SENSOR_IV  IX  IV TYPE  VALUE" << std::endl;

	f_lowGain.open(m_pedestalsLowGain.c_str());
	f_lowGain << "SCHEME_CODE 0" << std::endl;
	f_lowGain << "# CODE  LAYER SENSOR_IX SENSOR_IV  IX  IV TYPE  VALUE" << std::endl;

	for(int ISkiroc = 1; ISkiroc <= MAXSKIROCS; ISkiroc++) {
		for(int Channel = 0; Channel < SKIROC::NCHANNELS; Channel++) {
			HGCalTBElectronicsId ElId(ISkiroc, Channel);
			HGCalTBDetId DetId = essource_.emap_.eid2detId(ElId);
			if(DetId.layer() != 0) {
				f_highGain << " " << Code << " " << DetId.layer() << " " << SENSOR_IX << " " << SENSOR_IV << " " << DetId.iu() << " " << DetId.iv() << " " << " " << DetId.cellType() << " " << h_digi_layer_channel[ISkiroc - 1][Channel][1]->GetMean() << std::endl;
				f_lowGain << " " << Code << " " << DetId.layer() << " " << SENSOR_IX << " " << SENSOR_IV << " " << DetId.iu() << " " << DetId.iv() << " " << " " << DetId.cellType() << " " << h_digi_layer_channel[ISkiroc - 1][Channel][0]->GetMean() << std::endl;
			}
		}
	}
}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
DigiPlotter::fillDescriptions(edm::ConfigurationDescriptions& descriptions)
{
	//The following says we do not know what parameters are allowed so do no validation
	// Please change this to state exactly what you do use, even if it is no parameters
	// not sure what the above means or if it is still valid after my changes -- Tanmay
	edm::ParameterSetDescription desc;
	// desc.setUnknown();
	// descriptions.addDefault(desc);
	desc.addUntracked<std::string>("pedestalsHighGain", "");
	desc.addUntracked<std::string>("pedestalsLowGain", "");
	descriptions.add("hgcaltbdigisplotter", desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(DigiPlotter);
