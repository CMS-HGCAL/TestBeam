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
#include "HGCal/Geometry/interface/HGCalTBGeometryParameters.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "HGCal/DataFormats/interface/HGCalTBDataFrameContainers.h"
#include "HGCal/CondObjects/interface/HGCalElectronicsMap.h"
#include "HGCal/DataFormats/interface/HGCalTBElectronicsId.h"
#include "HGCal/DataFormats/interface/SKIROCParameters.h"
#define MAXVERTICES 6

class Pedestal
{
public:
	std::vector<std::pair<double, unsigned int> > _pedestals;
 	Pedestal(){};

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


//usesResource("TFileService");
edm::Service<TFileService> fs;

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
        void InitTH2Poly(TH2Poly& poly);
	// ----------member data ---------------------------
	Pedestal _pedestals[2]; // one set of pedestals per skiroc
	HGCalTBTopology IsCellValid;
	HGCalTBCellVertices TheCell;
	int sensorsize = 128;// The geometry for a 256 cell sensor hasnt been implemted yet. Need a picture to do this.
	std::vector<std::pair<double, double>> CellXY;
	std::pair<double, double> CellCentreXY;
	std::vector<std::pair<double, double>>::const_iterator it;
	const static int NSAMPLES = 2;
	const static int NLAYERS  = 1;
        int Sensor_Iu = 0;
        int Sensor_Iv = 0;
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
	consumesMany<SKIROC2DigiCollection>();
}//contructor ends here


DigiPlotter_New::~DigiPlotter_New()
{

	// do anything here that needs to be done at desctruction time
	// (e.g. close files, deallocate resources etc.)

}


void DigiPlotter_New::InitTH2Poly(TH2Poly& poly)
{
	double HexX[MAXVERTICES] = {0.};
	double HexY[MAXVERTICES] = {0.};

	for(int iv = -7; iv < 8; iv++) {
		for(int iu = -7; iu < 8; iu++) {
                        if(!IsCellValid.iu_iv_valid(NLAYERS, Sensor_Iu, Sensor_Iv, iu, iv, sensorsize)) continue;
			CellXY = TheCell.GetCellCoordinatesForPlots(NLAYERS, Sensor_Iu, Sensor_Iv, iu, iv, sensorsize);
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
        char name[300], title[300];
	std::vector<edm::Handle<SKIROC2DigiCollection> > ski;
	event.getManyByType(ski);
//        int Event = event.id().event();
        
	if(!ski.empty()) {

		std::vector<edm::Handle<SKIROC2DigiCollection> >::iterator i;
//		double Average_Pedestal_Per_Event1 = 0, Average_Pedestal_Per_Event2 = 0, Average_Pedestal_Half_Cell_Event1 = 0, Average_Pedestal_Half_Cell_Event2 = 0, Average_Pedestal_Calib_Cell_Event1 = 0;
//		int Cell_counter1 = 0, Cell_counter2 = 0, Cell_counter1_Half_Cell1 = 0, Cell_counter1_Half_Cell2 = 0, Cell_counter1_Calib_Cell1 = 0;
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
//                                 int iSample = 0; // using only 1 sample ///\todo to be generalized
                                 int n_layer = (detId).layer();
                                 int n_sensor_IU = (detId).sensorIU();
                                 int n_sensor_IV = (detId).sensorIV();
                                 int n_cell_iu = (detId).iu();
                                 int n_cell_iv = (detId).iv();
         

#ifdef DEBUG
				if(DEBUG) cout << endl << " Layer = " << n_layer << " Sensor IU = " << n_sensor_IU << " Sensor IV = " << n_sensor_IV << " Cell iu = " << n_cell_iu << " Cell iu = " << n_cell_iv << endl;
#endif
				if(!IsCellValid.iu_iv_valid(n_layer, n_sensor_IU, n_sensor_IV, n_cell_iu, n_cell_iv, sensorsize))  continue;

				CellCentreXY = TheCell.GetCellCentreCoordinatesForPlots(n_layer, n_sensor_IU, n_sensor_IV, n_cell_iu, n_cell_iv, sensorsize);
				double iux = (CellCentreXY.first < 0 ) ? (CellCentreXY.first + 0.0001) : (CellCentreXY.first - 0.0001) ;
				double iyy = (CellCentreXY.second < 0 ) ? (CellCentreXY.second + 0.0001) : (CellCentreXY.second - 0.0001);

//				_pedestals[detId.skiIndexInSensor()].Sum(detId, skiFrame[iSample].adcLow());

				if(((iux <= 0.25 && iyy >= -0.25 ) || (iux < -0.5 && iyy < 0)) && detId.cellType() == 0) {
				   cout<<endl<<iux*iyy<<endl;	
				}



			}
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//

			TH2Poly *h_digi_layer[MAXLAYERS], *highGain_hpoly[MAXLAYERS];
                        if(evId%10 == 0){
			for(unsigned int iLayer = 0; iLayer < MAXLAYERS; ++iLayer) {

				h_digi_layer[iLayer] = fs->make<TH2Poly>();
				sprintf(name, "FullLayer_ADC%i_Layer%i_Event%i", 0, iLayer + 1, evId);
				sprintf(title, "Sum of adc counts per cell for ADC%i Layer%i Event%i", 0, iLayer + 1, evId);
				h_digi_layer[iLayer]->SetName(name);
				h_digi_layer[iLayer]->SetTitle(title);
				InitTH2Poly(*h_digi_layer[iLayer]);

				highGain_hpoly[iLayer] = fs->make<TH2Poly>();
				sprintf(name, "FullLayer_ADC%i_Layer%i_Event%i", 1, iLayer + 1, evId);
				sprintf(title, "Sum of adc counts per cell for ADC%i Layer%i Event%i", 1, iLayer + 1, evId);
				highGain_hpoly[iLayer]->SetName(name);
				highGain_hpoly[iLayer]->SetTitle(title);
				InitTH2Poly(*highGain_hpoly[iLayer]);
			}
                       }

			for(SKIROC2DigiCollection::const_iterator j = Coll.begin(); j != Coll.end(); j++) {
//				bool flag_calib = false;
				const SKIROC2DataFrame& SKI = *j ;

				int iSample = 0; // using only 1 sample ///\todo to be generalized
                                int n_layer = (SKI.detid()).layer();
                                int n_sensor_IU = (SKI.detid()).sensorIU();
                                int n_sensor_IV = (SKI.detid()).sensorIV();
                                int n_cell_iu = (SKI.detid()).iu();
                                int n_cell_iv = (SKI.detid()).iv();
				HGCalTBDetId detId = SKI.detid();
				if(!IsCellValid.iu_iv_valid(n_layer, n_sensor_IU, n_sensor_IV, n_cell_iu, n_cell_iv, sensorsize)) continue;


#ifdef DEBUG
				if(DEBUG) cout << endl << " Layer = " << n_layer << " Sensor IU = " << n_sensor_IU << " Sensor IV = " << n_sensor_IV << " Cell iu = " << n_cell_iu << " Cell iu = " << n_cell_iv << endl;
#endif


				CellCentreXY = TheCell.GetCellCentreCoordinatesForPlots(n_layer, n_sensor_IU, n_sensor_IV, n_cell_iu, n_cell_iv, sensorsize);

				double iux = (CellCentreXY.first < 0 ) ? (CellCentreXY.first + 0.0001) : (CellCentreXY.first - 0.0001) ;
				double iyy = (CellCentreXY.second < 0 ) ? (CellCentreXY.second + 0.0001) : (CellCentreXY.second - 0.0001);

				if(detId.cellType() == 1 || detId.cellType() == 4) {
					if(evId%10 == 0) h_digi_layer[n_layer-1]->Fill(iux, iyy, 0.5 * (SKI[iSample].adcLow()));
				} else {
					if(evId%10 == 0) h_digi_layer[n_layer-1]->Fill(iux, iyy, (SKI[iSample].adcLow()));
				}

				if(evId%10 == 0) highGain_hpoly[n_layer-1]->Fill(iux , iyy, (SKI[iSample].adcHigh()));
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
