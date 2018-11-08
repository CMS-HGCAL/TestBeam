/* 
 * Computation of variables.
 */

/**
	@Author: Thorben Quast <tquast>
		22 November 2017
		thorben.quast@cern.ch / thorben.quast@rwth-aachen.de
*/
//source /cvmfs/cms.cern.ch/cmsset_default.sh
//source before compilation: source /cvmfs/cms.cern.ch/slc7_amd64_gcc630/external/tensorflow-c/1.1.0-cms/etc/profile.d/init.sh;

// system include files
#include <iostream>
#include <string>
#include <vector>
#include <algorithm>
#include <map>
#include <math.h>
// user include files
#include "FWCore/Framework/interface/EDProducer.h"
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "HGCal/CondObjects/interface/HGCalCondObjectTextIO.h"
#include "HGCal/CondObjects/interface/HGCalElectronicsMap.h"
#include "HGCal/DataFormats/interface/HGCalTBElectronicsId.h"
#include "HGCal/DataFormats/interface/HGCalTBRunData.h"	//for the runData type definition
#include "HGCal/DataFormats/interface/HGCalTBWireChamberData.h"
#include "HGCal/DataFormats/interface/HGCalTBRecHitCollections.h"
#include "HGCal/DataFormats/interface/HGCalTBClusterCollection.h"
#include "HGCal/DataFormats/interface/HGCalTBRecHit.h"
#include "HGCal/DataFormats/interface/HGCalTBCommonModeNoise.h"
#include "HGCal/DataFormats/interface/HGCalTBDWCTrack.h"
#include "HGCal/CondObjects/interface/HGCalTBDetectorLayout.h"

#include "HGCal/Geometry/interface/HGCalTBCellVertices.h"
#include "HGCal/Geometry/interface/HGCalTBTopology.h"
#include "HGCal/Geometry/interface/HGCalTBGeometryParameters.h"

#include "HGCal/Reco/interface/PositionResolutionHelpers.h"
#include "HGCal/Reco/interface/Sensors.h"

#include "TFile.h"
#include "TH2F.h"
#include "TH1F.h"  
#include <sstream>
#include <fstream>
#include <iomanip>
#include <set>

typedef std::map<int, std::vector<double> > WindowMap;

//#define DEBUG
//see presentation at CHEF 06 Oct 2017 by Thorben Quast
double X0PosJuly2017[7] = {6.3, 16.8, 25.3, 32.7, 41.1, 48.4, 48.4};	//includes upstream material, last value is a dummy 
double Lambda0PosJuly2017[7] = {0.35, 0.89, 1.6, 2.4, 3.3, 4.0, 4.0};
double weightsJuly2017[7] = {1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0};
double MIP2GeVJuly2017 = 84.9e-6;

double X0PosSeptember2017[18] = {2.76674,4.37279,5.97883,7.58488,9.19093,13.0331,16.1884,23.9158,26.6938,29.4718,32.4454,35.2234,38.0014,41.9177,44.8912,47.8647,50.8382,55.9605};
double Lambda0PosSeptember2017[18] = {0.15558,0.237188,0.318797,0.400405,0.482013,0.665655,0.85385,1.4237,1.70922,1.99473,2.28343,2.56895,2.85447,3.25925,3.54794,3.83664,4.12533,4.662};
double weightsSeptember2017[18] = {24.523, 17.461, 17.461, 17.461, 27.285, 38.737, 75.867, 83.382, 55.394, 55.823, 55.823, 55.394, 66.824, 67.253, 56.252, 56.252, 79.871, 103.49};
double MIP2GeVSeptember2017 = 84.9e-6;


double X0PosMarch2018[3] = {0.08, 4.7+0.28, 4.7+0.28+0.01};
double Lambda0PosMarch2018[3] = {1., 2., 3.};
double weightsMarch2018[3] = {1., 1., 1.};
double MIP2GeVMarch2018 = 84.9e-6;


double X0PosJune2018[28] = {0.932753,1.90872,2.81787,3.79384,4.70298,5.67895,6.5881,7.56407,8.47322,9.44918,10.3583,11.3343,12.2434,13.2194,14.1286,15.2717,16.1808,17.1568,18.0659,19.209,20.1182,21.0941,22.0033,23.0628,23.972,24.9479,25.8571,26.8331};
double Lambda0PosJune2018[28] = {0.0367242,0.0979071,0.129342,0.190524,0.221959,0.283142,0.314576,0.375759,0.407194,0.468377,0.499811,0.560994,0.592429,0.653611,0.685046,0.761895,0.793329,0.854512,0.885947,0.962795,0.99423,1.05541,1.08685,1.15586,1.1873,1.24848,1.27991,1.3411};
double weightsJune2018[28] = {7.84088,12.27451,6.92551,12.27451,6.92551,12.27451,6.92551,12.27451,6.92551,12.27451,6.92551,12.27451,6.92551,12.27451,6.92551,15.29152,6.92551,12.27451,6.92551,15.29152,6.92551,12.27451,6.92551,13.78302,6.92551,12.27451,6.92551,12.27451};
double MIP2GeVJune2018 = 84.9e-6;


double X0PosOctober2018_setup1[40] = {0.932753, 1.90872, 2.81787, 3.79384, 4.70298, 5.67895, 6.5881, 7.56407, 8.47322, 9.44918, 10.3583, 11.3343, 12.2434, 13.2194, 14.1286, 15.1045, 16.0137, 16.9896, 17.8988, 18.8748, 19.7839, 20.927, 21.8362, 22.9793, 23.8884, 25.1341, 26.0433, 27.3262, 30.1386, 33.0571, 35.9756, 38.8941, 41.8126, 44.6475, 45.7447, 48.6632, 51.5816, 54.6956, 57.726, 60.561};
double Lambda0PosOctober2018_setup1[40] = {0.0367242, 0.0979071, 0.129342, 0.190524, 0.221959, 0.283142, 0.314576, 0.375759, 0.407194, 0.468377, 0.499811, 0.560994, 0.592429, 0.653611, 0.685046, 0.746229, 0.777663, 0.838846, 0.870281, 0.931464, 0.962898, 1.03975, 1.07118, 1.14803, 1.17946, 1.25129, 1.28272, 1.35602, 1.6535, 1.95281, 2.25213, 2.55144, 2.85076, 3.14224, 3.25073, 3.55005, 3.84936, 4.15185, 4.44651, 4.73799};
double weightsOctober2018_setup1[40] = {1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 2.5, 2.5, 2.5, 2.5, 2.5, 2.5, 2.5, 2.5, 2.5, 2.5, 2.5, 2.5};
double MIP2GeVOctober2018_setup1 = 84.9e-6;

double X0PosOctober2018_setup2[39] = {0.932753, 1.90872, 2.81787, 3.79384, 4.70298, 5.67895, 6.5881, 7.56407, 8.47322, 9.44918, 10.3583, 11.3343, 12.2434, 13.2194, 14.1286, 15.1045, 16.0137, 16.9896, 17.8988, 18.8748, 19.7839, 20.927, 21.8362, 22.9793, 23.8884, 25.1341, 26.0433, 27.3262, 27.862, 28.5589, 31.3938, 34.3123, 37.1472, 42.9115, 45.83, 48.7485, 51.8624, 54.7809, 57.6994 };
double Lambda0PosOctober2018_setup2[39] = {0.0367242, 0.0979071, 0.129342, 0.190524, 0.221959, 0.283142, 0.314576, 0.375759, 0.407194, 0.468377, 0.499811, 0.560994, 0.592429, 0.653611, 0.685046, 0.746229, 0.777663, 0.838846, 0.870281, 0.931464, 0.962898, 1.03975, 1.07118, 1.14803, 1.17946, 1.25129, 1.28272, 1.35602, 1.41498, 1.46515, 1.75663, 2.05595, 2.34743, 2.94489, 3.24421, 3.54352, 3.84601, 4.14533, 4.44464};
double weightsOctober2018_setup2[39] = {1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 2.5, 2.5, 2.5, 2.5, 2.5, 2.5, 2.5, 2.5, 2.5, 2.5, 2.5};
double MIP2GeVOctober2018_setup2 = 84.9e-6;

//todo (23 October 2018): update values
double X0PosOctober2018_setup3[19] = {1.90872, 3.79384, 5.67895, 13.2194, 18.8748, 20.927, 22.9793, 27.862, 28.5589, 31.3938, 34.3123, 37.1472, 42.9115, 45.83, 48.7485, 51.8624, 54.7809, 57.6994 , 60.561};
double Lambda0PosOctober2018_setup3[19] = {0.129342, 0.221959, 0.314576, 0.407194, 0.685046, 0.838846, 1.17946, 1.41498, 1.46515, 1.75663, 2.05595, 2.34743, 2.94489, 3.24421, 3.54352, 3.84601, 4.14533, 4.44464, 4.73799};
double weightsOctober2018_setup3[19] = {1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 2.5, 2.5, 2.5, 2.5, 2.5, 2.5, 2.5, 2.5, 2.5, 2.5, 2.5, 2.5};
double MIP2GeVOctober2018_setup3 = 84.9e-6;



std::vector<double> getEigenValuesOfSymmetrix3x3(double A11, double A22, double A33, double A12, double A13, double A23) {
	//eigenvalue computation: https://arxiv.org/pdf/1306.6291.pdf
	double b = A11 + A22 + A33;
	double c = A11*A22 + A11*A33 + A22*A33 - pow(A12,2) - pow(A13,2) - pow(A23,2);
	double d = A11*pow(A23,2) + A22*pow(A13,2) + A33*pow(A12,2) - A11*A22*A33 - 2*A12*A13*A23;

	double p = pow(b, 2) - 3*c;
	double q = 2*pow(b, 3) - 9*b*c - 27*d;

	double delta = acos(q/sqrt(4*pow(p, 3)));
	double lambda1 = (b+2*sqrt(p)*cos(delta/3))/3;
	double lambda2 = (b+2*sqrt(p)*cos((delta+2*M_PI)/3))/3;
	double lambda3 = (b+2*sqrt(p)*cos((delta-2*M_PI)/3))/3;

	std::vector<double> EVs;
	EVs.push_back(lambda1);
	EVs.push_back(lambda2);
	EVs.push_back(lambda3);
	std::sort(EVs.begin(), EVs.end());

	return EVs;
}

//#define DEBUG

class VariableComputation : public edm::EDProducer {
	public:
		explicit VariableComputation(const edm::ParameterSet&);
		~VariableComputation();
		static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);

	private:
		virtual void produce(edm::Event& , const edm::EventSetup&);
		void ReadDWCWindows();
		void ReadCurrentDWCWindows(int);

		// ----------member data ---------------------------
		edm::EDGetTokenT<HGCalTBRecHitCollection> HGCalTBRecHitCollection_Token;	 		
		edm::EDGetTokenT<RunData> RunDataToken;	
		edm::EDGetTokenT<std::map<int, WireChamberData> > DWCToken;		
		edm::EDGetTokenT<HGCalTBDWCTrack> DWCTrackToken;		

		std::string m_UserRecordCollectionName;

		std::string m_electronicMap;
		std::string m_detectorLayoutFile;
		struct {
			HGCalElectronicsMap emap_;
			HGCalTBDetectorLayout layout_;
		} essource_;

		std::string m_layerPositionFile;

		int m_NHexaBoards;
		int m_NLayers;
		int N_layers_EE;		
		int N_layers_FH;
		int N_layers_BH;


		std::map<int, double> layerPositions;

		//energy sums:
		double MIP_cut_for_energy;	
		std::map<int, SensorHitMap*> Sensors;

		double energyAll_tot, energyE1_tot, energyE7_tot, energyE19_tot, energyE37_tot, energyE61_tot;
		double energyAllHG_tot, energyAllLG_tot, energyAllTOT_tot;
		double energyAll_weight, energyE1_weight, energyE7_weight, energyE19_weight, energyE37_weight, energyE61_weight;
		std::vector<double> energyAll_layer, energyE1_layer, energyE7_layer, energyE19_layer, energyE37_layer, energyE61_layer;
		std::vector<int> NAll_layer, NE1_layer, NE7_layer, NE19_layer, NE37_layer, NE61_layer;

		double layer_10Percent, layer_90Percent;
		double depthX0, depthLambda0, showerStartDepth;

		//distance information
		std::vector<double> mainCoreWidth;
		std::vector<double> cellDistance_layer;

		double Ixx, Iyy, Izz, Ixy, Ixz, Iyz;


		std::vector<std::string> pathsToMIPWindowFiles;
	  	std::map<std::pair<int, int> ,WindowMap  >loadedDWCWindows;
		WindowMap currentDWCWindows;
		


  		//coordinate system transformation
     	int x_max;
    	int x_min;
    	uint range_x;
    	int y_max;
    	int y_min;
    	uint range_y;

};

VariableComputation::VariableComputation(const edm::ParameterSet& iConfig) {	
	
	HGCalTBRecHitCollection_Token = consumes<HGCalTBRecHitCollection>(iConfig.getParameter<edm::InputTag>("HGCALTBRECHITS"));
	RunDataToken= consumes<RunData>(iConfig.getParameter<edm::InputTag>("RUNDATA"));
	DWCToken= consumes<std::map<int, WireChamberData> >(iConfig.getParameter<edm::InputTag>("MWCHAMBERS"));
	DWCTrackToken= consumes<HGCalTBDWCTrack>(iConfig.getParameter<edm::InputTag>("DWCTRACKS"));
	
	m_UserRecordCollectionName = iConfig.getUntrackedParameter<std::string>("UserRecordCollectionName","DoubleUserRecords");
	
	m_electronicMap = iConfig.getUntrackedParameter<std::string>("ElectronicMap","HGCal/CondObjects/data/map_CERN_Hexaboard_28Layers_AllFlipped.txt");
	m_detectorLayoutFile = iConfig.getUntrackedParameter<std::string>("DetectorLayout","HGCal/CondObjects/data/layerGeom_oct2017_h2_17layers.txt");
	m_layerPositionFile = iConfig.getParameter<std::string>("layerPositionFile");
	m_NHexaBoards= iConfig.getUntrackedParameter<int>("NHexaBoards", 17);
	m_NLayers= iConfig.getUntrackedParameter<int>("NLayers", 17);

	MIP_cut_for_energy = iConfig.getUntrackedParameter<double>("CellEnergyCut", 4.);
	
	produces <UserRecords<double> >(m_UserRecordCollectionName);




	for( int ilayer=0; ilayer<m_NLayers; ilayer++ ){
		energyAll_layer.push_back(0.);
		energyE1_layer.push_back(0.);
		energyE7_layer.push_back(0.);
		energyE19_layer.push_back(0.);
		energyE37_layer.push_back(0.);
		energyE61_layer.push_back(0.);
		NAll_layer.push_back(0);
		NE1_layer.push_back(0);
		NE7_layer.push_back(0);
		NE19_layer.push_back(0);
		NE37_layer.push_back(0);
		NE61_layer.push_back(0);
		mainCoreWidth.push_back(0.);
		cellDistance_layer.push_back(0.);
	}


	std::fstream file; 
	char fragment[100];
	int readCounter = -1;

	file.open(m_layerPositionFile.c_str(), std::fstream::in);

	std::cout<<"Reading file "<<m_layerPositionFile<<" -open: "<<file.is_open()<<std::endl;
	int layer=0;
	while (file.is_open() && !file.eof()) {
		readCounter++;
		file >> fragment;
		if (readCounter==0) layer=atoi(fragment);
		if (readCounter==1) {
			layerPositions[layer]=atof(fragment);
			readCounter=-1;
		}
	}

	HGCalCondObjectTextIO io(0);
	edm::FileInPath fip(m_electronicMap);
	if (!io.load(fip.fullPath(), essource_.emap_)) {
		throw cms::Exception("Unable to load electronics map");
	};

	//for single cell spectra
	//ReadDWCWindows();

}


VariableComputation::~VariableComputation() {
	return;
}

// ------------ method called for each event  ------------
void VariableComputation::produce(edm::Event& event, const edm::EventSetup& setup) {
	edm::Handle<RunData> rd;
	event.getByToken(RunDataToken, rd);
	
	edm::Handle<HGCalTBDWCTrack> dwctrack;
	edm::Handle<std::map<int, WireChamberData> > dwcs;
	if (rd->booleanUserRecords.has("hasValidDWCMeasurement")&&rd->booleanUserRecords.get("hasValidDWCMeasurement")) {	
		try {
			event.getByToken(DWCTrackToken, dwctrack);
			event.getByToken(DWCToken, dwcs);
		} catch(const std::exception& e) {
		}
	}

	edm::Handle<HGCalTBRecHitCollection> Rechits;
	event.getByToken(HGCalTBRecHitCollection_Token, Rechits);

	std::unique_ptr<UserRecords<double> > UR(new UserRecords<double>);
	
	/**********                                 ****************/
	


	//first some basic information
	UR->add("eventID", rd->event);
	UR->add("run", rd->run);
	UR->add("pdgID", rd->pdgID);
	UR->add("beamEnergy", rd->energy);
	UR->add("configuration", rd->configuration);
	UR->add("runType", rd->runType);
	
	//



	/**********                                 ****************/
	//selecting rechits certain layers
	switch(rd->configuration) {
		case 1:
			N_layers_EE = 2;		//June 2017, H2
			N_layers_FH = 4;
			N_layers_BH = 12;
			break;
  	case 2:							//September 2017, H2
			N_layers_EE = 7;
			N_layers_FH = 10;
			N_layers_BH = 12;
			break;
  	case 3:
  	case 4:							//October 2017, H6
			N_layers_EE = 4;		
			N_layers_FH = 6;		
			N_layers_BH = 12;		
			break;
	case 5:							//March 2018, DESY 
			N_layers_EE = 1;		
			N_layers_FH = 0;		
			N_layers_BH = 0;		
			break;			
  	case 6:							
  	case 7:
  	case 8:
  	case 9:
  	case 10:
  	case 13:
  	case 14:
  	case 15:
  	case 16:
			N_layers_EE = 3;		
			N_layers_FH = 0;		
			N_layers_BH = 0;
			break;
  	case 11:
  	case 12:					//March 2018, DESY 
			N_layers_EE = 2;		
			N_layers_FH = 0;		
			N_layers_BH = 0;
			break;	
  	case 17:
  	case 18:					//June 2018, H2
  	case 19:					//June 2018, H2
  	case 20:					//June 2018, H2
  	case 21:					//June 2018, H2
			N_layers_EE = 28;		
			N_layers_FH = 0;		
			N_layers_BH = 0;
			break;			
  	case 22:					//October 2018 - setup 1, H2
			N_layers_EE = 28;		
			N_layers_FH = 12;		
			N_layers_BH = 12;
			break;	
  	case 23:					//October 2018 - setup 2, H2
			N_layers_EE = 28;		
			N_layers_FH = 11;		
			N_layers_BH = 12;
			break;			
  	case 24:					//October 2018 - setup 3, H2
			N_layers_EE = 7;		
			N_layers_FH = 12;		
			N_layers_BH = 12;
			break;				
	}
	std::vector<HGCalTBRecHit> rechits_selected;
	for(auto Rechit : *Rechits) {
		if ((abs(rd->pdgID)==11) && (Rechit.id()).layer()>N_layers_EE) continue;
		rechits_selected.push_back(Rechit);
	}
		
	/**********                                 ****************/
	//filling and sorting of the rechits:

	std::vector<double> rechit_energies;
	std::vector<HGCalTBRecHit> rechits;
	std::vector<HGCalTBRecHit> MIPhits;
	std::vector<HGCalTBRecHit> noisehits;

	double M=0, xmean=0, ymean=0, zmean=0;

	for (int layer=1; layer<=m_NLayers; layer++) {
		Sensors[layer] = new SensorHitMap(layer);
		Sensors[layer]->setSensorSize(133);
		
	}

	for(auto Rechit : rechits_selected) {	
		int layer = (Rechit.id()).layer();

		if (Rechit.energy() > MIP_cut_for_energy) {
			Sensors[layer]->addHit(Rechit, 1.);
			if ((Rechit.id()).cellType() == 0) {
				double x = Rechit.getCellCenterCartesianCoordinate(0)*10.;		//conversion to mm
				double y = Rechit.getCellCenterCartesianCoordinate(1)*10.;		//conversion to mm
				double z = layerPositions[layer];
				double m = Rechit.energy();
				M += m;
				xmean += m*x;
				ymean += m*y;
				zmean += m*z;
				
				rechit_energies.push_back(m);
				rechits.push_back(Rechit);
			}
		} else if (Rechit.energy() > 0.5) {
			if ((Rechit.id()).cellType() == 0) MIPhits.push_back(Rechit);
		} else  {
			if ((Rechit.id()).cellType() == 0) noisehits.push_back(Rechit);
		}
	}

	/**********                                 ****************/
	//mean positions

	xmean /= M;
	ymean /= M;
	zmean /= M;

	UR->add("xmean", xmean);
	UR->add("ymean", ymean);
	UR->add("zmean", zmean);


	/**********                                 ****************/
	//rechit spectra positions

	int NRechits = rechit_energies.size();
	int N25PercentsRechits = NRechits*1./4.;
	int N50PercentsRechits = NRechits*2./4.;
	int N75PercentsRechits = NRechits*3./4.;
	std::sort(rechit_energies.begin(), rechit_energies.end());

	UR->add("NRechits", NRechits);
	UR->add("NMIPHits", MIPhits.size());
	UR->add("NNoisehits", noisehits.size());
	UR->add("25PercentQuantileRechitSpectrum", (NRechits>4) ? rechit_energies[N25PercentsRechits-1] : -1.);
	UR->add("50PercentQuantileRechitSpectrum", (NRechits>4) ? rechit_energies[N50PercentsRechits-1] : -1.);
	UR->add("75PercentQuantileRechitSpectrum", (NRechits>4) ? rechit_energies[N75PercentsRechits-1] : -1.);

	/**********                                 ****************/
	//determine inertia tensor:
	Ixx=Iyy=Izz=Ixy=Ixz=Iyz=0;
	for(auto Rechit : rechits_selected) {	
		int layer = (Rechit.id()).layer();
		if (Rechit.energy() > MIP_cut_for_energy) {
			if ((Rechit.id()).cellType() == 0) {
				double x = Rechit.getCellCenterCartesianCoordinate(0)*10.;
				double y = Rechit.getCellCenterCartesianCoordinate(1)*10.;
				double z = layerPositions[layer];
				double m = Rechit.energy();

				Ixx += m*(pow(y-ymean, 2)+pow(z-zmean, 2));
				Iyy += m*(pow(x-xmean, 2)+pow(z-zmean, 2));
				Izz += m*(pow(x-xmean, 2)+pow(y-ymean, 2));
				Ixy -= m*(x-xmean)*(y-ymean);
				Ixz -= m*(x-xmean)*(z-zmean);
				Iyz -= m*(y-ymean)*(z-zmean);

			}
		}
	}

	Ixx /= M;
	Iyy /= M;
	Izz /= M;
	Ixy /= M;
	Ixz /= M;
	Iyz /= M;
	std::vector<double> I_EV = getEigenValuesOfSymmetrix3x3(Ixx, Iyy, Izz, Ixy, Ixz, Iyz);
	UR->add("Ixx", Ixx);
	UR->add("Iyy", Iyy);
	UR->add("Izz", Izz);
	UR->add("Ixy", Ixy);
	UR->add("Ixz", Ixz);
	UR->add("Iyz", Iyz);
	UR->add("I_EV1", I_EV[0]);
	UR->add("I_EV2", I_EV[1]);
	UR->add("I_EV3", I_EV[2]);


	/**********                                 ****************/
	
	//Energy information
	energyAllHG_tot = energyAllLG_tot = energyAllTOT_tot = energyAll_tot = 0;
	for(auto Rechit : rechits_selected) {
		energyAllHG_tot+=Rechit.energyHigh();
		energyAllLG_tot+=Rechit.energyLow();
		energyAllTOT_tot+=Rechit.energyTot();
		energyAll_tot+=Rechit.energy();
	}
	UR->add("EAllHG_tot", energyAllHG_tot);
	UR->add("EAllLG_tot", energyAllLG_tot);
	UR->add("EAllTOT_tot", energyAllTOT_tot);
	UR->add("EAll_tot", energyAll_tot);

	energyE1_tot = energyE7_tot = energyE19_tot = energyE37_tot = energyE61_tot = 0.;	
	energyE1_weight = energyE7_weight = energyE19_weight = energyE37_weight = energyE61_weight = energyAll_weight = 0.;
	
	depthX0 = 0, depthLambda0 = 0;
	showerStartDepth = -1.;
	layer_10Percent = layer_90Percent = -1;

	double energySum_layers = 0, energySumAll_layers = 0;
	std::vector<std::pair<double, double> > relevantHitPositions;
	for (std::map<int, SensorHitMap*>::iterator it=Sensors.begin(); it!=Sensors.end(); it++) {
		it->second->calculateCenterPosition(CONSIDERALL, LINEARWEIGHTING);
		energySumAll_layers+= it->second->getTotalWeight();
	}

	
	for (std::map<int, SensorHitMap*>::iterator it=Sensors.begin(); it!=Sensors.end(); it++) {
		//most intensive cell
		it->second->calculateCenterPosition(CONSIDERALL, MOSTINTENSIVE);
		energyE1_tot += it->second->getTotalWeight();
		
		energyE1_layer[it->first-1] = it->second->getTotalWeight();
		relevantHitPositions = it->second->getHitPositionsForPositioning();
		NE1_layer[it->first-1] = (int)relevantHitPositions.size();
		relevantHitPositions.clear();

		//one ring around
		it->second->calculateCenterPosition(CONSIDERSEVEN, LINEARWEIGHTING);
		energyE7_tot += it->second->getTotalWeight();
		energyE7_layer[it->first-1] = it->second->getTotalWeight();
		relevantHitPositions = it->second->getHitPositionsForPositioning();
		NE7_layer[it->first-1] = (int)relevantHitPositions.size();
		relevantHitPositions.clear();

		//two rings around
		it->second->calculateCenterPosition(CONSIDERNINETEEN, LINEARWEIGHTING);
		energyE19_tot += it->second->getTotalWeight();
		energyE19_layer[it->first-1] = it->second->getTotalWeight();
		relevantHitPositions = it->second->getHitPositionsForPositioning();
		NE19_layer[it->first-1] = (int)relevantHitPositions.size();
		
		relevantHitPositions.clear();
	
		//three rings around
 		it->second->calculateCenterPosition(CONSIDERTHIRTYSEVEN, LINEARWEIGHTING);
		energyE37_tot += it->second->getTotalWeight();
		energyE37_layer[it->first-1] = it->second->getTotalWeight();
		relevantHitPositions = it->second->getHitPositionsForPositioning();
		NE37_layer[it->first-1] = (int)relevantHitPositions.size();
		relevantHitPositions.clear();

     	//four rings around
 		it->second->calculateCenterPosition(CONSIDERSIXTYONE, LINEARWEIGHTING);
		energyE61_tot += it->second->getTotalWeight();
		energyE61_layer[it->first-1] = it->second->getTotalWeight();
		relevantHitPositions = it->second->getHitPositionsForPositioning();
		NE61_layer[it->first-1] = (int)relevantHitPositions.size();
		relevantHitPositions.clear();
		

		//sum of all
		it->second->calculateCenterPosition(CONSIDERALL, LINEARWEIGHTING);
		energyAll_layer[it->first-1] = it->second->getTotalWeight();
		energySum_layers += it->second->getTotalWeight();
		if ((layer_10Percent==-1) && (energySum_layers>0.1 * energySumAll_layers)) layer_10Percent = ((it->first) * energyAll_layer[it->first-1] + (it->first-1) * energyAll_layer[it->first-2]) / (energyAll_layer[it->first-1] + energyAll_layer[it->first-2]);
		if ((layer_90Percent==-1) && (energySum_layers>0.9 * energySumAll_layers)) layer_90Percent = ((it->first) * energyAll_layer[it->first-1] + (it->first-1) * energyAll_layer[it->first-2]) / (energyAll_layer[it->first-1] + energyAll_layer[it->first-2]);

		relevantHitPositions = it->second->getHitPositionsForPositioning();
		NAll_layer[it->first-1] = (int)relevantHitPositions.size();
		UR->add("NAll_layer"+std::to_string(it->first), NAll_layer[it->first-1]);
		relevantHitPositions.clear();
		UR->add("EAll_layer"+std::to_string(it->first), energyAll_layer[it->first-1]);
	
		UR->add("E1PerE7_layer"+std::to_string(it->first), energyE1_layer[it->first-1]/energyE7_layer[it->first-1]);
		UR->add("E1PerE19_layer"+std::to_string(it->first), energyE1_layer[it->first-1]/energyE19_layer[it->first-1]);
		UR->add("E7PerE19_layer"+std::to_string(it->first), energyE7_layer[it->first-1]/energyE19_layer[it->first-1]);
		UR->add("E19PerE37_layer"+std::to_string(it->first), energyE19_layer[it->first-1]/energyE37_layer[it->first-1]);
		UR->add("E37PerE61_layer"+std::to_string(it->first), energyE37_layer[it->first-1]/energyE61_layer[it->first-1]);

		double weight = 0., MIP2GeV=1., X0 = 0., lambda0 = 0.;
		if (rd->configuration==1) {
			weight = weightsJuly2017[it->first-1]*1e-3;
			MIP2GeV = MIP2GeVJuly2017;
			X0 = X0PosJuly2017[it->first-1];
			lambda0 = Lambda0PosJuly2017[it->first-1];
		}
		else if (rd->configuration==2) {
			weight = weightsSeptember2017[it->first-1]*1e-3;
			MIP2GeV = MIP2GeVSeptember2017;
			X0 = X0PosSeptember2017[it->first-1];
			lambda0 = Lambda0PosSeptember2017[it->first-1];
		}
		else if (rd->configuration>=5 && rd->configuration<=17) {
			weight = weightsMarch2018[it->first-1]*1e-3;
			MIP2GeV = MIP2GeVMarch2018;
			X0 = X0PosMarch2018[it->first-1];
			lambda0 = Lambda0PosMarch2018[it->first-1];
		}
		else if ((rd->configuration>=17) && (rd->configuration<=21)) {
			weight = weightsJune2018[it->first-1]*1e-3;
			MIP2GeV = MIP2GeVJune2018;
			X0 = X0PosJune2018[it->first-1];
			lambda0 = Lambda0PosJune2018[it->first-1];
		}
		else if ((rd->configuration==22)) {
			weight = weightsOctober2018_setup1[it->first-1]*1e-3;
			MIP2GeV = MIP2GeVOctober2018_setup1;
			X0 = X0PosOctober2018_setup1[it->first-1];
			lambda0 = Lambda0PosOctober2018_setup1[it->first-1];
		}	
		else if ((rd->configuration==23)) {
			weight = weightsOctober2018_setup2[it->first-1]*1e-3;
			MIP2GeV = MIP2GeVOctober2018_setup2;
			X0 = X0PosOctober2018_setup2[it->first-1];
			lambda0 = Lambda0PosOctober2018_setup2[it->first-1];
		}		
		else if ((rd->configuration>=24)) {
			weight = weightsOctober2018_setup3[it->first-1]*1e-3;
			MIP2GeV = MIP2GeVOctober2018_setup3;
			X0 = X0PosOctober2018_setup3[it->first-1];
			lambda0 = Lambda0PosOctober2018_setup3[it->first-1];
		}				
		depthX0 += X0*energyAll_layer[it->first-1]*(MIP2GeV+weight);
		depthLambda0 += lambda0*energyAll_layer[it->first-1]*(MIP2GeV+weight);

		if ((showerStartDepth==-1.)&&(energyE19_layer[it->first-1] > 20.)&&(NE19_layer[it->first-1]>1)) {
			showerStartDepth = lambda0;
		}

		energyE1_weight += energyE1_layer[it->first-1]*(MIP2GeV+weight); 
		energyE7_weight += energyE7_layer[it->first-1]*(MIP2GeV+weight); 
		energyE19_weight += energyE19_layer[it->first-1]*(MIP2GeV+weight); 
		energyE37_weight += energyE37_layer[it->first-1]*(MIP2GeV+weight); 
		energyE61_weight += energyE61_layer[it->first-1]*(MIP2GeV+weight); 
		energyAll_weight += energyAll_layer[it->first-1]*(MIP2GeV+weight); 

		//position resolution
		
		it->second->calculateCenterPosition(CONSIDERNINETEEN, LOGWEIGHTING_35_10);
		//investigate here
		//x = -x in DWC coordinate system
		UR->add("RecoPosX_layer"+std::to_string(it->first),(it->second->getLabHitPosition().first));
		UR->add("RecoPosY_layer"+std::to_string(it->first),(it->second->getLabHitPosition().second));
	
	}
	UR->add("layer_10Percent", layer_10Percent);
	UR->add("layer_90Percent", layer_90Percent);

	if (rd->booleanUserRecords.get("hasValidDWCMeasurement")) {	
		UR->add("dwc1_multiplicity", dwcs->at(0).averageHitMultiplicty);
		UR->add("dwc2_multiplicity", dwcs->at(1).averageHitMultiplicty);
		UR->add("dwc3_multiplicity", dwcs->at(2).averageHitMultiplicty);
		UR->add("dwc4_multiplicity", dwcs->at(3).averageHitMultiplicty);
		UR->add("dwctrack_type", dwctrack->referenceType);
		UR->add("dwctrack_chi2x", dwctrack->chi2_x);
		UR->add("dwctrack_chi2y", dwctrack->chi2_y);
	}

	depthX0 /= energyAll_weight;
	depthLambda0 /= energyAll_weight;

	UR->add("E1_tot", energyE1_tot);
	UR->add("E7_tot", energyE7_tot);
	UR->add("E19_tot", energyE19_tot);
	UR->add("E37_tot", energyE37_tot);
	UR->add("E61_tot", energyE61_tot);
	
	UR->add("E1_weight", energyE1_weight);
	UR->add("E7_weight", energyE7_weight);
	UR->add("E19_weight", energyE19_weight);
	UR->add("E37_weight", energyE37_weight);
	UR->add("E61_weight", energyE61_weight);
	UR->add("EAll_weight", energyAll_weight);

	UR->add("depthX0", depthX0);
	UR->add("depthLambda0", depthLambda0);
	UR->add("showerStartDepth", showerStartDepth);


	double E_EE = 0, E_FH = 0;
	for (int l=0; l<N_layers_EE; l++) E_EE += energyAll_layer[l];
	for (int l=N_layers_EE; l<N_layers_EE+N_layers_FH; l++) E_FH += energyAll_layer[l];

	UR->add("E_EE", E_EE);
	UR->add("E_FH", E_FH);
	UR->add("E_EEperE_tot", E_EE/(E_EE+E_FH));
	UR->add("E_FHperE_tot", E_FH/(E_EE+E_FH));

	

	/**********                                 ****************/
	
	//Spacing information

	for (std::map<int, SensorHitMap*>::iterator it=Sensors.begin(); it!=Sensors.end(); it++) {
		//two rings around
		it->second->calculateCenterPosition(CONSIDERNINETEEN, LINEARWEIGHTING);		
		std::vector<double> gaussianFitParameters = it->second->fit2DGaussian();
		if (gaussianFitParameters[0] != 0) mainCoreWidth[it->first-1] = -1.;
		else mainCoreWidth[it->first-1] = sqrt(pow(gaussianFitParameters[3], 2) + pow(gaussianFitParameters[5],2));
		UR->add("width_E19_layer"+std::to_string(it->first), mainCoreWidth[it->first-1]);

		//consideration of all
		it->second->calculateCenterPosition(CONSIDERALL, LINEARWEIGHTING);	
		cellDistance_layer[it->first-1] = it->second->getDistanceBetweenMostIntenseCells();		

		UR->add("d2_maxE_layer"+std::to_string(it->first), cellDistance_layer[it->first-1]);
	}

	/**********                                 ****************/
		

	//spectra of selected cells
	ReadCurrentDWCWindows(1680);
	std::vector<double> cell_chip1_ch36_energySpectra; for (size_t l=0; (int)l<m_NLayers; l++) cell_chip1_ch36_energySpectra.push_back(-1.);
	
	if (rd->booleanUserRecords.get("hasValidDWCMeasurement")&&dwctrack->valid&&(dwctrack->referenceType>=7) && (dwctrack->chi2_x<=5.) && (dwctrack->chi2_y<=5.)) {

		for(auto Rechit : rechits_selected) {	

			HGCalTBElectronicsId eid( essource_.emap_.detId2eid( Rechit.id().rawId() ) );
			int board = eid.iskiroc_rawhit() / 4;
			int skiroc = eid.iskiroc_rawhit() % 4;
			int channel = eid.ichan();
			

			if (skiroc!=1)	continue;	//specific for September 2017 setup (7EE, 10FH layers)
			if ((channel!=36)&&(channel!=38)&&(channel!=44)&&(channel!=46)&&(channel!=52)&&(channel!=54)&&(channel!=56)) continue;

			int key = board*1000+skiroc*100+channel;

			int layer = (Rechit.id()).layer();
			double dwc_x_layer = dwctrack->DWCExtrapolation_XY(layer).first;
			double dwc_y_layer = dwctrack->DWCExtrapolation_XY(layer).second;

			if (rd->runType==HGCAL_TB_BEAM) {
				//if (!Rechit.checkFlag(HGCalTBRecHit::kLowGainSaturated)) continue;
		 		if (currentDWCWindows.find(key) == currentDWCWindows.end())	continue;	
				if (-dwc_x_layer < currentDWCWindows[key][0] || -dwc_x_layer > currentDWCWindows[key][1]) continue;
				if (dwc_y_layer < currentDWCWindows[key][2] || dwc_y_layer > currentDWCWindows[key][3]) continue;
			}
			
			cell_chip1_ch36_energySpectra[layer-1] = Rechit.energy();
			
		}
	}

	for (size_t l=0; (int)l<m_NLayers; l++) UR->add("energy_chip1_layer"+std::to_string(l), cell_chip1_ch36_energySpectra[l]);

	/**********                                 ****************/
	event.put(std::move(UR), m_UserRecordCollectionName);

	for (std::map<int, SensorHitMap*>::iterator it=Sensors.begin(); it!=Sensors.end(); it++) {
		delete (*it).second;
	};	Sensors.clear();	

}// analyze ends here



void VariableComputation::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
	edm::ParameterSetDescription desc;
	desc.setUnknown();
	descriptions.addDefault(desc);
}


void VariableComputation::ReadDWCWindows() {
  	
	std::fstream file; 
	char fragment[100];
	int readCounter = -2;

	WindowMap _parameters;
	  	
	std::cout<<"Opening: "<<"/afs/cern.ch/user/t/tquast/CMS_HGCal_Upgrade/workflows/tbAnalysis2017/configurations/MIP_DWC_Windows_1658_1683.txt"<<std::endl;
	file.open("/afs/cern.ch/user/t/tquast/CMS_HGCal_Upgrade/workflows/tbAnalysis2017/configurations/MIP_DWC_Windows_1658_1683.txt", std::fstream::in);

	int minRun, maxRun;
	if (file.is_open()) {
		file >> fragment;
		minRun = atoi(fragment);
		file >> fragment;
		maxRun = atoi(fragment);
	}

	int iboard = -1, iskiroc = -1, ichannel = -1;
	double DWC_x_min, DWC_x_max, DWC_y_min, DWC_y_max;

	while (file.is_open() && !file.eof()) {		
		if (readCounter!=-2) readCounter++;
			file >> fragment;

		if (std::string(fragment)=="y_max" ) readCounter = -1;  //first parameter is read out

		if (readCounter==0) iboard = atoi(fragment);
		if (readCounter==1) iskiroc = atoi(fragment); 
		if (readCounter==2) ichannel = atoi(fragment);
		if (readCounter==3) DWC_x_min = atof(fragment);
		if (readCounter==4) DWC_x_max = atof(fragment);
		if (readCounter==5) DWC_y_min = atof(fragment);
		if (readCounter==6) { DWC_y_max = atof(fragment);
			int key = iboard*1000+iskiroc*100+ichannel;
			_parameters[key].push_back(DWC_x_min);
			_parameters[key].push_back(DWC_x_max);
			_parameters[key].push_back(DWC_y_min);
			_parameters[key].push_back(DWC_y_max);
			readCounter=-1;
		}
	}
	
	WindowMap ::iterator it;
	#ifdef DEBUG
		for (it=_parameters.begin(); it!=_parameters.end(); it++) {
			std::cout<<"key: "<<it->first;
			for (int i=0; i<4; i++) {
				std::cout<<"  "<<it->second[i];
			}
			std::cout<<std::endl;
		}
	#endif

	loadedDWCWindows[std::make_pair(minRun, maxRun)] = _parameters;
}


void VariableComputation::ReadCurrentDWCWindows(int this_run) {
	std::map<std::pair<int, int> ,WindowMap  >::iterator it;
	for (it=loadedDWCWindows.begin(); it!=loadedDWCWindows.end(); it++) {
		int run_min = it->first.first;
		int run_max = it->first.second;

		if (this_run>=run_min && (this_run<=run_max || run_max==-1) ) {
			currentDWCWindows = it->second;
			break;
		}
	}
}

//define this as a plug-in
DEFINE_FWK_MODULE(VariableComputation);