// -*- C++ -*-
//
// Package:    HGCal/DigiPlotter_New
// Class:      DigiPlotter_New
//
/**\class DigiPlotter_New DigiPlotter_New.cc HGCal/DigiPlotter_New/plugins/DigiPlotter_New.cc

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
#include "HGCal/DataFormats/interface/HGCalTBDataFrameContainers.h"
#include "HGCal/CondObjects/interface/HGCalElectronicsMap.h"
#include "HGCal/DataFormats/interface/HGCalTBElectronicsId.h"

#define MAXVERTICES 6

class Pedestal
{
	std::vector<std::pair<double, unsigned int> > _pedestals;

	Pedestal();
	void Sum(HGCalTBDetId detId, double pedestal)
	{
		_pedestals[detId.cellType()].first += pedestal;
		_pedestals[detId.cellType()].second++;
	}

	float GetPedestal(HGCalTBDetId detId)
	{
		/// \todo here the logic can be adapted if one want to use MouseBites and HalfCells
		// for example:
		// if(detId.cellType()==HGCalTBDetId::kCellHalf || detId.cellType()==HGCalTBDetId::kCellMouseBite)
		// return _pedestals[HGCalTBDetId::kCellMouseBite].first+_pedestals[HGCalTBDetId::kCellHalf].first/(_pedestals[HGCalTBDetId::kCellMouseBite].second+_pedestals[HGCalTBDetId::kCellHalf].second);
		return _pedestals[detId.cellType()].first / _pedestals[detId.cellType()].second;
	}
};

//
// class declaration
//

// If the analyzer does not use TFileService, please remove
// the template argument to the base class so the class inherits
// from  edm::one::EDAnalyzer<> and also remove the line from
// constructor "usesResource("TFileService");"
// This will improve performance in multithreaded jobs.


class DigiPlotter_New : public edm::one::EDAnalyzer<edm::one::SharedResources>
{

public:
	explicit DigiPlotter_New(const edm::ParameterSet&);
	~DigiPlotter_New();
	static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);
private:
	virtual void beginJob() override;
	void analyze(const edm::Event& , const edm::EventSetup&) override;
	virtual void endJob() override;
	// ----------member data ---------------------------
	Pedestal _pedestals[2]; // one set of pedestals per skiroc
	HGCalTBTopology IsCellValid;
	HGCalTBCellVertices TheCell;
	int sensorsize = 128;// The geometry for a 256 cell sensor hasnt been implemted yet. Need a picture to do this.
	std::vector<HGCalTB::Point2D> CellXY;
	HGCalTB::Point2D CellCentreXY;
	std::vector<HGCalTB::Point2D>::const_iterator it;
	const static int NSAMPLES = 2;
	const static int NLAYERS  = 1;

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
DigiPlotter_New::DigiPlotter_New(const edm::ParameterSet& iConfig)
{
	//now do what ever initialization is needed
	usesResource("TFileService");
	edm::Service<TFileService> fs;
	consumesMany<SKIROC2DigiCollection>();
}//contructor ends here


DigiPlotter_New::~DigiPlotter_New()
{

	// do anything here that needs to be done at desctruction time
	// (e.g. close files, deallocate resources etc.)

}


void DigiPlotter_New::InitTH2Poly(TH2Poly& poly, const HGCalTBDetId detId)
{
	double HexX[MAXVERTICES] = {0.};
	double HexY[MAXVERTICES] = {0.};

	for(int iv = -7; iv < 8; iv++) {
		for(int iu = -7; iu < 8; iu++) {
			CellXY = TheCell.GetCellCoordinatesForPlots(detId, sensorsize);
			assert(CellXY.size() == 4 || CellXY.size() == 6);
			unsigned int iVertex = 0;
			for(it = CellXY.begin(); it != CellXY.end(); it++) {
				HexX[iVertex] =  it->first;
				HexY[iVertex] =  it->second;
				++iVertex;
			}
//Somehow cloning of the TH2Poly was not working. Need to look at it. Currently physically booked another one.
			poly.AddBin(CellXY.size(), HexX, HexY);
		}//loop over iu
	}//loop over iv
}

//
// member functions
//

// ------------ method called for each event  ------------
void
DigiPlotter_New::analyze(const edm::Event& event, const edm::EventSetup& setup)
{
	int evId = event.id().event() - 1;

	using namespace edm;
	using namespace std;
	std::vector<edm::Handle<SKIROC2DigiCollection> > ski;
	event.getManyByType(ski);
//        int Event = event.id().event();
	if(!ski.empty()) {

		std::vector<edm::Handle<SKIROC2DigiCollection> >::iterator i;
		double Average_Pedestal_Per_Event1 = 0, Average_Pedestal_Per_Event2 = 0, Average_Pedestal_Half_Cell_Event1 = 0, Average_Pedestal_Half_Cell_Event2 = 0, Average_Pedestal_Calib_Cell_Event1 = 0;
		int Cell_counter1 = 0, Cell_counter2 = 0, Cell_counter1_Half_Cell1 = 0, Cell_counter1_Half_Cell2 = 0, Cell_counter1_Calib_Cell1 = 0;
//                int counter1=0, counter2=0;
		for(i = ski.begin(); i != ski.end(); i++) {
			const SKIROC2DigiCollection& Coll = *(*i);
#ifdef DEBUG
			if(DEBUG) cout << "SKIROC2 Digis: " << i->provenance()->branchName() << endl;
#endif
//////////////////////////////////Evaluate average pedestal per event to subtract out//////////////////////////////////
			for(SKIROC2DigiCollection::const_iterator k = Coll.begin(); k != Coll.end(); k++) {
				const SKIROC2DataFrame& skiFrame = *k ;
				HGCalTBDetId detId = skiFrame.detid();

#ifdef DEBUG
				if(DEBUG) cout << endl << " Layer = " << n_layer << " Sensor IU = " << n_sensor_IU << " Sensor IV = " << n_sensor_IV << " Cell iu = " << n_cell_iu << " Cell iu = " << n_cell_iv << endl;
#endif
				if(!IsCellValid.iu_iv_valid(detId))  continue;

				CellCentreXY = TheCell.GetCellCentreCoordinatesForPlots(detId, sensorsize);
				double iux = (CellCentreXY.first < 0 ) ? (CellCentreXY.first + 0.0001) : (CellCentreXY.first - 0.0001) ;
				double iyy = (CellCentreXY.second < 0 ) ? (CellCentreXY.second + 0.0001) : (CellCentreXY.second - 0.0001);
				int iSample = 0;

				_pedestals[detId.skiIndexInSensor()].Sum(detId, skiFrame.adcLow());


				// if(detId().cellType() == kCellTypeCalibInner || detId.cellType() == kCellTypeCalibOuter) {
				// 	Average_Pedestal_Calib_Cell_Event1 += skiFrame[iSample].adcLow();
				// 	Cell_counter1_Calib_Cell1++;
				// };


				// if(((iux <= 0.25 && iyy >= -0.25 ) || (iux < -0.5)) && ((detId.cellType() == kCellHalf) || detId.cellType() == kCellMouseBite) ) {
				// 	Average_Pedestal_Half_Cell_Event1 += skiFrame[iSample].adcLow();
				// 	Cell_counter1_Half_Cell1++;
				// };

				// if(((iux > -0.25 && iyy < -0.50 ) || (iux > 0.50)) && ((detId.cellType() == kCellHalf) || detId.cellType() == kCellMouseBite )  ) {
				// 	Average_Pedestal_Half_Cell_Event2 += skiFrame[iSample].adcLow();
				// 	Cell_counter1_Half_Cell2++;
				// };
				if(((iux <= 0.25 && iyy >= -0.25 ) || (iux < -0.5 && iyy < 0)) && detId.cellType() == kCellTypeStandard) {
					// 	Cell_counter1++;
					// 	Average_Pedestal_Per_Event1 += skiFrame[iSample].adcLow();
					Sum_Hist_cells_SKI1->Fill(skiFrame.adcLow());
				}

				// if(((iux > -0.25 && iyy < -0.50 ) || (iux > 0.50 && iyy < 0) ) && detId.cellType() == kCellTypeStandard) {
				// 	Cell_counter2++;
				// 	Average_Pedestal_Per_Event2 += skiFrame[iSample].adcLow();
				// }


			}
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//

			char name[300], title[300];
			TH2Poly *h_digi_layer[MAXLAYERS], *highGain_hpoly[MAXLAYERS];

			for(unsigned int iLayer = 0; iLayer < MAXLAYERS; ++iLayer) {

				h_digi_layer[iLayer] = fs->make<TH2Poly>();
				sprintf(name, "FullLayer_ADC%i_Layer%i_Event%i", 0, iLayer + 1, evId);
				sprintf(title, "Sum of adc counts per cell for ADC%i Layer%i Event%i", 0, iLayer + 1, evId);
				h_digi_layer->SetName(name);
				h_digi_layer->SetTitle(title);
				InitTH2Poly(*h_digi_layer[iLayer], detId);

				highGain_hpoly[iLayer] = fs->make<TH2Poly>();
				sprintf(name, "FullLayer_ADC%i_Layer%i_Event%i", 1, iLayer + 1, evId);
				sprintf(title, "Sum of adc counts per cell for ADC%i Layer%i Event%i", 1, iLayer + 1, evId);
				highGain_hpoly->SetName(name);
				highGain_hpoly->SetTitle(title);
				InitTH2Poly(*highGain_hpoly[iLayer], detId);
			}

			sprintf(name, "SumCellHist_SKI1_Event%i", evId);
			sprintf(title, "Sum of adc counts of cells  Event%i", evId);
			Sum_Hist_cells_SKI1 = fs->make<TH1F>(name, title, 4096, 0 , 4095);

			for(SKIROC2DigiCollection::const_iterator j = Coll.begin(); j != Coll.end(); j++) {
				bool flag_calib = false;
				const SKIROC2DataFrame& SKI = *j ;

				int iSample = 0; // using only 1 sample ///\todo to be generalized

				HGCalTBDetId detId = SKI[iSample].detid();
				if(!IsCellValid.iu_iv_valid(detId, sensorsize)) continue;


#ifdef DEBUG
				if(DEBUG) cout << endl << " Layer = " << n_layer << " Sensor IU = " << n_sensor_IU << " Sensor IV = " << n_sensor_IV << " Cell iu = " << n_cell_iu << " Cell iu = " << n_cell_iv << endl;
#endif


				CellCentreXY = TheCell.GetCellCentreCoordinatesForPlots(detId, sensorsize);

				double iux = (CellCentreXY.first < 0 ) ? (CellCentreXY.first + 0.0001) : (CellCentreXY.first - 0.0001) ;
				double iyy = (CellCentreXY.second < 0 ) ? (CellCentreXY.second + 0.0001) : (CellCentreXY.second - 0.0001);

				if(detId.cellType() == HGCalTBDetId::kCellTypeCalibOuter || detId.cellType() == HGCalTBDetId::kCellTypeCalibInner) {
					h_digi_layer->Fill(iux, iyy, 0.5 * (SKI[iSample].adcLow() - _pedestals[detId.skiIndexInSensor()].GetPedestal(detId)));
				} else {
					h_digi_layer->Fill(iux, iyy, (SKI[iSample].adcLow() - _pedestals[detId.skiIndexInSensor()].GetPedestal(detId)));
				}
#ifdef DEBUG					if(SKI.detid().cellType() == HGCalDetId::kCellMouseBite && event.id().event() == 34) cout << endl << " iux= " << iux << " iyy= " << iyy << " ADC = " << SKI[iSample].adcLow() << " Ped= " << Average_Pedestal_Half_Cell_Event1 / (Cell_counter1_Half_Cell1) << endl;
#endif

				highGain_hpoly->Fill(iux , iyy, (SKI[iSample].adcHigh() - Sum_Hist_cells_SKI1->GetMean()));
			}
		}
	} else {
		edm::LogWarning("DQM") << "No SKIROC2 Digis";
	}

}//analyze method ends here


// ------------ method called once each job just before starting event loop  ------------
void
DigiPlotter_New::beginJob()
{

}

// ------------ method called once each job just after ending the event loop  ------------
void
DigiPlotter_New::endJob()
{
}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
DigiPlotter_New::fillDescriptions(edm::ConfigurationDescriptions& descriptions)
{
	//The following says we do not know what parameters are allowed so do no validation
	// Please change this to state exactly what you do use, even if it is no parameters
	edm::ParameterSetDescription desc;
	desc.setUnknown();
	descriptions.addDefault(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(DigiPlotter_New);
