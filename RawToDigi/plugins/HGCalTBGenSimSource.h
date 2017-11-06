#include "FWCore/Framework/interface/Event.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Sources/interface/ProducerSourceFromFiles.h"
#include "FWCore/ParameterSet/interface/ConfigurationDescriptions.h"
#include "FWCore/ParameterSet/interface/ParameterSetDescription.h"
#include "HGCal/CondObjects/interface/HGCalCondObjectTextIO.h"
#include "HGCal/CondObjects/interface/HGCalElectronicsMap.h"
#include "HGCal/DataFormats/interface/HGCalTBElectronicsId.h"
#include "HGCal/DataFormats/interface/HGCalTBDetId.h"
#include "HGCal/DataFormats/interface/HGCalTBRunData.h"
#include "HGCal/DataFormats/interface/HGCalTBMultiWireChamberData.h"
#include "HGCal/DataFormats/interface/HGCalTBRecHitCollections.h"
#include "HGCal/Geometry/interface/HGCalTBCellVertices.h"
#include "HGCal/Geometry/interface/HGCalWaferGeometry.h"
#include <fstream>
#include <iostream>
#include <stdio.h>
#include <string>
#include <sstream>
#include "stdlib.h"
#include <vector>
#include <cmath>

#include "TFile.h"
#include "TBranch.h"
#include "TTree.h"
#include "TDirectory.h"
#include "TRandom.h"
/**
 *
 *
 *
 *
 *
**/
double ADCtoMIP_CERN[16] =  {17.31, 17.12, 16.37, 17.45, 17.31, 16.98, 16.45, 16.19, 17.55, 17.19, 16.99, 17.92, 15.95, 16.64, 16.79, 15.66};
double MIP2GeV_sim = 51.91e-06; //mpv muon EMM pysics list

//must be a globally available variable, otherwise a segmentation fault is thrown  
struct {
  HGCalElectronicsMap emap_;
} essource_;


struct FileInfo {
	int index;
	double energy;
	std::string runType;
	int config;
	std::string name;
};

//to the EDM::Event via auxiliary information

class HGCalTBGenSimSource : public edm::ProducerSourceFromFiles
{
private:
	void fillConfiguredRuns(std::fstream& map_file);
	bool setRunAndEventInfo(edm::EventID& id, edm::TimeValue_t& time, edm::EventAuxiliary::ExperimentType&);
	virtual void produce(edm::Event & e);
  	virtual void endJob() override;
	
	std::string outputCollectionName;

	std::string runEnergyMapFile;
	std::string inputPathFormat;
	
	std::vector<FileInfo> _fileNames;
	std::vector<FileInfo>::iterator fileIterator;

	int currentRun;
	int currentEvent;
	int eventCounter;

	TFile *rootFile;
  	TTree *tree;   //!pointer to the analyzed TTree or TChain
  	TDirectory *dir;

  	std::vector<unsigned int> *simHitCellIdE;
  	std::vector<float>        *simHitCellEnE;
  	double					  beamX;
  	double					  beamY;
  	double					  beamP;
  	TBranch                   *b_simHitCellIdE;   
  	TBranch                   *b_simHitCellEnE;   
  	TBranch                   *b_beamX;   
  	TBranch                   *b_beamY;   
  	TBranch                   *b_beamP;   

	//getting the required electronic mapping
	std::string _e_mapFile;

	HGCalTBCellVertices TheCell;
	HexGeometry* geomc;

	TRandom* randgen;
	double energyNoise;
	double energyNoiseResolution;
	double smearingResolution;

public:
	explicit HGCalTBGenSimSource(const edm::ParameterSet & pset, edm::InputSourceDescription const& desc);
	virtual ~HGCalTBGenSimSource()
	{
		delete rootFile;
			
	}
	
};
