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
#include "HGCal/Geometry/interface/HGCalTBGeometryParameters.h"

#define MAXVERTICES 6

class Pedestal
{
public:
	std::vector<std::pair<double, unsigned int> > _pedestals;
	Pedestal() {};

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
	void InitTH2Poly(TH2Poly& poly);

	// ----------member data ---------------------------
	Pedestal _pedestals[MAXLAYERS][2]; // one set of pedestals per skiroc
	HGCalTBTopology IsCellValid;
	HGCalTBCellVertices TheCell;
	int sensorsize = 128;// The geometry for a 256 cell sensor hasnt been implemted yet. Need a picture to do this.
	std::vector<HGCalTB::Point2D> CellXY;
	HGCalTB::Point2D CellCentreXY;
	std::vector<HGCalTB::Point2D>::const_iterator it;
	const static int NSAMPLES = 2;
	edm::Service<TFileService> fs;

	const static int NLAYERS  = 1;
	int Sensor_Iu = 0;
	int Sensor_Iv = 0;

	unsigned int _DQMprescale;
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
	consumesMany<SKIROC2DigiCollection>();

	_DQMprescale = iConfig.getParameter<unsigned int>("DQMprescale");
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

///\todo Rajdeep still using loop over iv and iu
	std::vector<HGCalTBDetId> detIds = IsCellValid.validDetIds(0, 128); // hard coding that should be fixed
	for( auto detId : detIds) {
		CellXY = TheCell.GetCellCoordinatesForPlots(detId, sensorsize);
		assert(CellXY.size() == 4 || CellXY.size() == 6);
		unsigned int iVertex = 0;
		for(it = CellXY.begin(); it != CellXY.end(); it++) {
			HexX[iVertex] =  it->first;
			HexY[iVertex] =  it->second;
			++iVertex;
		}
	}
//Somehow cloning of the TH2Poly was not working. Need to look at it. Currently physically booked another one.
	poly.AddBin(CellXY.size(), HexX, HexY);

}

//
// member functions
//

// ------------ method called for each event  ------------
void
DigiPlotter_New::analyze(const edm::Event& event, const edm::EventSetup& setup)
{
	using namespace edm;

	int evId = event.id().event() - 1; // starts from 0 because is used as index in DQM plots....
	if(evId % _DQMprescale != 0) return;

	std::vector<edm::Handle<SKIROC2DigiCollection> > ski;
	event.getManyByType(ski);

	if(ski.empty()) return;
	assert(ski.size() == 1); // the histograms have the same name for every digi collection... it should be fixed

	std::vector<edm::Handle<SKIROC2DigiCollection> >::iterator i;
	for(i = ski.begin(); i != ski.end(); i++) {
		const SKIROC2DigiCollection& Coll = *(*i);
#ifdef DEBUG
		if(DEBUG) std::cout << "SKIROC2 Digis: " << i->provenance()->branchName() << std::endl;
#endif

		char name[300], title[300];
		TH2Poly *lowGain_hpoly[MAXLAYERS], *highGain_hpoly[MAXLAYERS];

		for(unsigned int iLayer = 0; iLayer < MAXLAYERS; ++iLayer) {

			lowGain_hpoly[iLayer] = fs->make<TH2Poly>();
			sprintf(name, "FullLayer_ADC%i_Layer%i_Event%i", 0, iLayer + 1, evId);
			sprintf(title, "Sum of adc counts per cell for ADC%i Layer%i Event%i", 0, iLayer + 1, evId);
			lowGain_hpoly[iLayer]->SetName(name);
			lowGain_hpoly[iLayer]->SetTitle(title);
			InitTH2Poly(*lowGain_hpoly[iLayer]);

			highGain_hpoly[iLayer] = fs->make<TH2Poly>();
			sprintf(name, "FullLayer_ADC%i_Layer%i_Event%i", 1, iLayer + 1, evId);
			sprintf(title, "Sum of adc counts per cell for ADC%i Layer%i Event%i", 1, iLayer + 1, evId);
			highGain_hpoly[iLayer]->SetName(name);
			highGain_hpoly[iLayer]->SetTitle(title);
			InitTH2Poly(*highGain_hpoly[iLayer]);
		}

		sprintf(name, "SumCellHist_SKI1_Event%i", evId);
		sprintf(title, "Sum of adc counts of cells  Event%i", evId);
		TH1F* Sum_Hist_cells_SKI1 = fs->make<TH1F>(name, title, 4096, 0 , 4095);

//////////////////////////////////Evaluate average pedestal per event to subtract out//////////////////////////////////
		for(SKIROC2DigiCollection::const_iterator skiFrame_itr = Coll.begin(); skiFrame_itr != Coll.end(); skiFrame_itr++) {
			const SKIROC2DataFrame& skiFrame = *skiFrame_itr ;
			HGCalTBDetId detId = skiFrame.detid();
			unsigned int iSample = 0; /// \todo now only using one sample, make it more general
			if(!IsCellValid.isValidDetId(detId, sensorsize))  continue;

			CellCentreXY = TheCell.GetCellCentreCoordinatesForPlots(detId, sensorsize);
			double iux = (CellCentreXY.first < 0 ) ? (CellCentreXY.first + 0.0001) : (CellCentreXY.first - 0.0001) ;
			double iyy = (CellCentreXY.second < 0 ) ? (CellCentreXY.second + 0.0001) : (CellCentreXY.second - 0.0001);

			_pedestals[detId.layer()][detId.skiIndexInSensor()].Sum(detId, skiFrame[iSample].adcLow());

			if(((iux <= 0.25 && iyy >= -0.25 ) || (iux < -0.5 && iyy < 0)) && detId.cellType() == HGCalTBDetId::kCellTypeStandard) {
//						std::cout << std::endl << iux*iyy << std::endl;
				Sum_Hist_cells_SKI1->Fill(skiFrame[iSample].adcLow());// the histogram is defined only according to the DQM prescale
			}
		}


///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//

		for(SKIROC2DigiCollection::const_iterator skiFrame_itr = Coll.begin(); skiFrame_itr != Coll.end(); skiFrame_itr++) {
			const SKIROC2DataFrame& SKI = *skiFrame_itr ;

			int iSample = 0; // using only 1 sample ///\todo to be generalized

			HGCalTBDetId detId = SKI.detid();
			if(!IsCellValid.isValidDetId(detId, sensorsize)) continue;

			CellCentreXY = TheCell.GetCellCentreCoordinatesForPlots(detId, sensorsize);

			double iux = (CellCentreXY.first < 0 ) ? (CellCentreXY.first + 0.0001) : (CellCentreXY.first - 0.0001) ;
			double iyy = (CellCentreXY.second < 0 ) ? (CellCentreXY.second + 0.0001) : (CellCentreXY.second - 0.0001);

			if(detId.cellType() == HGCalTBDetId::kCellTypeCalibOuter || detId.cellType() == HGCalTBDetId::kCellTypeCalibInner) { //1 and 4
				lowGain_hpoly[detId.layer()]->Fill(iux, iyy, 0.5 * (SKI[iSample].adcLow()));
//lowGain_hpoly[detId.layer()]->Fill(iux, iyy, 0.5 * (SKI[iSample].adcLow() - _pedestals[detId.skiIndexInSensor()].GetPedestal(detId)));
			} else {
				lowGain_hpoly[detId.layer()]->Fill(iux, iyy, (SKI[iSample].adcLow() ));
//					lowGain_hpoly[detId.layer()]->Fill(iux, iyy, (SKI[iSample].adcLow() - _pedestals[detId.skiIndexInSensor()].GetPedestal(detId)));
			}
#ifdef DEBUG
			if(SKI.detid().cellType() == HGCalDetId::kCellMouseBite && event.id().event() == 34) std::cout << std::endl << " iux= " << iux << " iyy= " << iyy << " ADC = " << SKI[iSample].adcLow() << " Ped= " << Average_Pedestal_Half_Cell_Event1 / (Cell_counter1_Half_Cell1) << std::endl;
#endif

			highGain_hpoly[detId.layer()]->Fill(iux , iyy, (SKI[iSample].adcHigh() - _pedestals[detId.layer()][detId.skiIndexInSensor()].GetPedestal(detId)));
		}
	}
}	//analyze method ends here


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
DigiPlotter_New::fillDescriptions(edm::ConfigurationDescriptions & descriptions)
{
	//The following says we do not know what parameters are allowed so do no validation
	// Please change this to state exactly what you do use, even if it is no parameters
	edm::ParameterSetDescription desc;
	desc.setUnknown();
	descriptions.addDefault(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(DigiPlotter_New);
