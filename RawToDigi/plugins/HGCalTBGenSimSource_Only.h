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

//must be a globally available variable, otherwise a segmentation fault is thrown  
struct {
  HGCalElectronicsMap emap_;
} essource_;


//to the EDM::Event via auxiliary information

class HGCalTBGenSimSource_Only : public edm::ProducerSourceFromFiles
{
private:
	bool setRunAndEventInfo(edm::EventID& id, edm::TimeValue_t& time, edm::EventAuxiliary::ExperimentType&);
	virtual void produce(edm::Event & e);
  	virtual void endJob() override;
	
	std::string outputCollectionName;
	std::vector<std::string> fileNames_;
	

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

public:
	explicit HGCalTBGenSimSource_Only(const edm::ParameterSet & pset, edm::InputSourceDescription const& desc);
	virtual ~HGCalTBGenSimSource_Only()
	{
		delete rootFile;
			
	}
	
};
