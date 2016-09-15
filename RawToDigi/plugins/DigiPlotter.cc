// -*- C++ -*-
//
// Package:    HGCal/DigiPlotter
// Class:      DigiPlotter
//
/**\class DigiPlotter DigiPlotter.cc HGCal/DigiPlotter/plugins/DigiPlotter.cc

 Description: [one line class summary]

 Implementation:
     [Notes on implementation]
*/
//
// Original Author:  Rajdeep Mohan Chatterjee
//         Created:  Mon, 15 Feb 2016 09:47:43 GMT
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
#include "HGCal/DataFormats/interface/HGCalTBRecHitCollections.h"
#include "HGCal/DataFormats/interface/HGCalTBDetId.h"
#include "HGCal/DataFormats/interface/HGCalTBRecHit.h"
#include "HGCal/Geometry/interface/HGCalTBCellVertices.h"
#include "HGCal/Geometry/interface/HGCalTBTopology.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "HGCal/CondObjects/interface/HGCalElectronicsMap.h"
#include "HGCal/CondObjects/interface/HGCalCondObjectTextIO.h"
#include "HGCal/DataFormats/interface/HGCalTBElectronicsId.h"
#include "HGCal/DataFormats/interface/HGCalTBDataFrameContainers.h"
#include "HGCal/Geometry/interface/HGCalTBGeometryParameters.h"

using namespace std;


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
	bool DEBUG = 0;
	HGCalTBTopology IsCellValid;
	HGCalTBCellVertices TheCell;
	std::string mapfile_ = "HGCal/CondObjects/data/map_CERN_8Layers_Sept2016.txt";
	struct {
		HGCalElectronicsMap emap_;
	} essource_;
	int sensorsize = 128;// The geometry for a 256 cell sensor hasnt been implemted yet. Need a picture to do this.
	std::vector<std::pair<double, double>> CellXY;
	std::pair<double, double> CellCentreXY;
	std::vector<std::pair<double, double>>::const_iterator it;
	const static int NSAMPLES = 2;
	TH2Poly *h_digi_layer[NSAMPLES][MAXLAYERS];
	TH1F    *h_digi_layer_summed[NSAMPLES][MAXLAYERS];
	TProfile    *h_digi_layer_profile[NSAMPLES][MAXLAYERS];
	const static int cellx = 15;
	const static int celly = 15;
	int Sensor_Iu = 0;
	int Sensor_Iv = 0;
	TH2F* Noise_2D_Profile[NSAMPLES][MAXLAYERS];
	TH1F  *h_digi_layer_channel[MAXSKIROCS][64][NSAMPLES];
//        TH1F  *h_digi_layer_cell_event[NSAMPLES][MAXLAYERS][cellx][celly][512];
	char name[50], title[50];
	double ADC_Sum_SKI_Layer[2][MAXLAYERS][2]; // 2 SKIROCs per layer, High gain and low gain ADC HARD CODED
	int Cell_Count_SKI_Layer[2][4]; // 2 SKIROCs per layer, High gain and low gain ADC HARD CODED
	string m_pedestalsHighGain;
	string m_pedestalsLowGain;
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
	const int HalfHexVertices = 4;
	double HalfHexX[HalfHexVertices] = {0.};
	double HalfHexY[HalfHexVertices] = {0.};
	const int FullHexVertices = 6;
	double FullHexX[FullHexVertices] = {0.};
	double FullHexY[FullHexVertices] = {0.};
	for(int nsample = 0; nsample < NSAMPLES; nsample++) {
		for(int nlayers = 0; nlayers < MAXLAYERS; nlayers++) {
			sprintf(name, "Noise_2D_Profile_ADC%i_Layer%i", nsample, nlayers);
			sprintf(title, "Noise 2D Profile ADC%i Layer%i", nsample, nlayers);
			Noise_2D_Profile[nsample][nlayers] = fs->make<TH2F>(name, title, 128, 0, 127, 2000, -1000, 1000);
		}
	}
	for(int ISkiroc = 1; ISkiroc <= MAXSKIROCS; ISkiroc++) {
		for(int Channel = 0; Channel < 64; Channel++) {
			for(int iii = 0; iii < NSAMPLES; iii++) {
				sprintf(name, "Ski_%i_Channel_%i_ADC%i", ISkiroc, Channel, iii);
				sprintf(title, "Ski %i Channel %i ADC%i", ISkiroc, Channel, iii);
				h_digi_layer_channel[ISkiroc - 1][Channel][iii] = fs->make<TH1F>(name, title, 4096, 0., 4095.);
			}
		}
	}
	int iii = 0;
	for(int nsample = 0; nsample < NSAMPLES; nsample++) {
		for(int nlayers = 0; nlayers < MAXLAYERS; nlayers++) {
//Booking a "hexagonal" histograms to display the sum of Digis for NSAMPLES, in 1 SKIROC in 1 layer. To include all layers soon. Also the 1D Digis per cell in a sensor is booked here for NSAMPLES.
			sprintf(name, "FullLayer_ADC%i_Layer%i", nsample, nlayers + 1);
			sprintf(title, "Sum of adc counts per cell for ADC%i Layer%i", nsample, nlayers + 1);
			h_digi_layer[nsample][nlayers] = fs->make<TH2Poly>();
			h_digi_layer[nsample][nlayers]->SetName(name);
			h_digi_layer[nsample][nlayers]->SetTitle(title);
			sprintf(name, "FullLayer_ADC%i_Layer%i_summed", nsample, nlayers + 1);
			sprintf(title, "Sum of adc counts for all cells in ADC%i Layer%i", nsample, nlayers + 1);
			h_digi_layer_summed[nsample][nlayers] = fs->make<TH1F>(name, title, 4096, 0., 4095.);
			h_digi_layer_summed[nsample][nlayers]->GetXaxis()->SetTitle("Digis[adc counts]");
			sprintf(name, "FullLayer_ADC%i_Layer%i_profile", nsample, nlayers + 1);
			sprintf(title, "profile of adc counts for all cells in ADC%i Layer%i", nsample, nlayers + 1);
			h_digi_layer_profile[nsample][nlayers] = fs->make<TProfile>(name, title, 128, 0, 127, 0., 4095.);
			h_digi_layer_profile[nsample][nlayers]->GetXaxis()->SetTitle("Channel #");
			h_digi_layer_profile[nsample][nlayers]->GetYaxis()->SetTitle("ADC counts");


			for(int iv = -7; iv < 8; iv++) {
				for(int iu = -7; iu < 8; iu++) {
					if(!IsCellValid.iu_iv_valid(nlayers, Sensor_Iu, Sensor_Iv, iu, iv, sensorsize)) continue;
					CellXY = TheCell.GetCellCoordinatesForPlots(nlayers, Sensor_Iu, Sensor_Iv, iu, iv, sensorsize);
					int NumberOfCellVertices = CellXY.size();
					iii = 0;
					if(NumberOfCellVertices == 4) {
						for(it = CellXY.begin(); it != CellXY.end(); it++) {
							HalfHexX[iii] =  it->first;
							HalfHexY[iii++] =  it->second;
						}
//Somehow cloning of the TH2Poly was not working. Need to look at it. Currently physically booked another one.
						h_digi_layer[nsample][nlayers]->AddBin(NumberOfCellVertices, HalfHexX, HalfHexY);
					} else if(NumberOfCellVertices == 6) {
						iii = 0;
						for(it = CellXY.begin(); it != CellXY.end(); it++) {
							FullHexX[iii] =  it->first;
							FullHexY[iii++] =  it->second;
						}
						h_digi_layer[nsample][nlayers]->AddBin(NumberOfCellVertices, FullHexX, FullHexY);
					}

				}//loop over iu
			}//loop over iv
		}//loop over nlayers
	}//loop over nsamples

	m_pedestalsHighGain = iConfig.getUntrackedParameter<string>("pedestalsHighGain", "");
	m_pedestalsLowGain = iConfig.getUntrackedParameter<string>("pedestalsLowGain", "");
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
				const SKIROC2DataFrame& SKI_1 = *k ;
				int n_layer = (SKI_1.detid()).layer();
				int n_sensor_IU = (SKI_1.detid()).sensorIU();
				int n_sensor_IV = (SKI_1.detid()).sensorIV();
				int n_cell_iu = (SKI_1.detid()).iu();
				int n_cell_iv = (SKI_1.detid()).iv();
				uint32_t EID = essource_.emap_.detId2eid(SKI_1.detid());
				HGCalTBElectronicsId eid(EID);
				if(DEBUG) cout << endl << " Layer = " << n_layer << " Sensor IU = " << n_sensor_IU << " Sensor IV = " << n_sensor_IV << " Cell iu = " << n_cell_iu << " Cell iu = " << n_cell_iv << endl;
				if(!IsCellValid.iu_iv_valid(n_layer, n_sensor_IU, n_sensor_IV, n_cell_iu, n_cell_iv, sensorsize))  continue;
//                                ADC_Sum_SKI_Layer[eid.iskiroc() - 2*(n_layer - 1) - 1][n_layer - 1][1] += SKI_1[0].adcHigh();
//                                ADC_Sum_SKI_Layer[eid.iskiroc() - 2*(n_layer - 1) - 1][n_layer - 1][0] += SKI_1[0].adcLow();
//                                Cell_Count_SKI_Layer[eid.iskiroc() - 2*(n_layer - 1) - 1][n_layer - 1] += 1;
			}



//			cout << "SKIROC2 Digis: " << i->provenance()->branchName() << endl;
			for(SKIROC2DigiCollection::const_iterator j = Coll.begin(); j != Coll.end(); j++) {
				const SKIROC2DataFrame& SKI = *j ;
				int n_layer = (SKI.detid()).layer();
				int n_sensor_IU = (SKI.detid()).sensorIU();
				int n_sensor_IV = (SKI.detid()).sensorIV();
				int n_cell_iu = (SKI.detid()).iu();
				int n_cell_iv = (SKI.detid()).iv();
				uint32_t EID = essource_.emap_.detId2eid(SKI.detid());
				HGCalTBElectronicsId eid(EID);
				if(DEBUG) cout << endl << " Layer = " << n_layer << " Sensor IU = " << n_sensor_IU << " Sensor IV = " << n_sensor_IV << " Cell iu = " << n_cell_iu << " Cell iu = " << n_cell_iv << endl;
				if(!IsCellValid.iu_iv_valid(n_layer, n_sensor_IU, n_sensor_IV, n_cell_iu, n_cell_iv, sensorsize))  continue;
				CellCentreXY = TheCell.GetCellCentreCoordinatesForPlots(n_layer, n_sensor_IU, n_sensor_IV, n_cell_iu, n_cell_iv, sensorsize);
				double iux = (CellCentreXY.first < 0 ) ? (CellCentreXY.first + 0.0001) : (CellCentreXY.first - 0.0001) ;
				double iyy = (CellCentreXY.second < 0 ) ? (CellCentreXY.second + 0.0001) : (CellCentreXY.second - 0.0001);
				int nsample = 0;
				h_digi_layer[nsample][n_layer - 1]->Fill(iux , iyy, SKI[nsample].adcLow());
				h_digi_layer_profile[nsample][n_layer - 1]->Fill(counter1++, SKI[nsample].adcLow(), 1);
//				h_digi_layer_summed[nsample][n_layer - 1]->Fill(ADC_Sum_SKI_Layer[eid.iskiroc() - 2*(n_layer - 1) - 1][n_layer - 1][0]);
				if(eid.iskiroc() > 0)	h_digi_layer_channel[eid.iskiroc() - 1][eid.ichan()][nsample]->Fill(SKI[nsample].adcLow());
				nsample = 1;
				h_digi_layer[nsample][n_layer - 1]->Fill(iux , iyy, SKI[nsample - 1].adcHigh());
//				Noise_2D_Profile[nsample][n_layer - 1]->Fill();
				h_digi_layer_profile[nsample][n_layer - 1]->Fill(counter2++, SKI[nsample - 1].adcHigh(), 1);
//				h_digi_layer_summed[nsample][n_layer - 1]->Fill(ADC_Sum_SKI_Layer[eid.iskiroc() - 2*(n_layer - 1) - 1][n_layer - 1][1]);
//                                        if(((SKI.detid()).cellType() != 4) && (eid.ichan() == 0) ) cout<<endl<<"SKIROC=  "<<eid.iskiroc()<<" Chan= "<<eid.ichan()<<" u= "<<n_cell_iu<<" v = "<<n_cell_iv<<endl;
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
	ofstream fs1, fs2;
	fs1.open(m_pedestalsHighGain.c_str());
	fs1 << "SCHEME_CODE 0" << endl;
	fs1 << "# CODE  LAYER SENSOR_IX SENSOR_IV  IX  IV TYPE  VALUE" << endl;
	fs2.open(m_pedestalsLowGain.c_str());
	fs2 << "SCHEME_CODE 0" << endl;
	fs2 << "# CODE  LAYER SENSOR_IX SENSOR_IV  IX  IV TYPE  VALUE" << endl;

	for(int ISkiroc = 1; ISkiroc <= MAXSKIROCS; ISkiroc++) {
		for(int Channel = 0; Channel < 64; Channel++) {
			HGCalTBElectronicsId ElId(ISkiroc, Channel);
			HGCalTBDetId DetId = essource_.emap_.eid2detId(ElId);
			if(DetId.layer() != 0) {
				fs1 << " " << Code << " " << DetId.layer() << " " << SENSOR_IX << " " << SENSOR_IV << " " << DetId.iu() << " " << DetId.iv() << " " << " " << DetId.cellType() << " " << h_digi_layer_channel[ISkiroc - 1][Channel][1]->GetMean() << endl;
				fs2 << " " << Code << " " << DetId.layer() << " " << SENSOR_IX << " " << SENSOR_IV << " " << DetId.iu() << " " << DetId.iv() << " " << " " << DetId.cellType() << " " << h_digi_layer_channel[ISkiroc - 1][Channel][0]->GetMean() << endl;
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
	desc.addUntracked<string>("pedestalsHighGain", "");
	desc.addUntracked<string>("pedestalsLowGain", "");
	descriptions.add("hgcaltbdigisplotter", desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(DigiPlotter);
