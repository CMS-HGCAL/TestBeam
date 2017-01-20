#include "FWCore/Framework/interface/Event.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Sources/interface/ProducerSourceFromFiles.h"
#include "FWCore/ParameterSet/interface/ConfigurationDescriptions.h"
#include "FWCore/ParameterSet/interface/ParameterSetDescription.h"
#include "HGCal/DataFormats/interface/HGCalTBDetId.h"
#include "HGCal/DataFormats/interface/HGCalTBRunData.h"
#include "HGCal/DataFormats/interface/HGCalTBRecHitCollections.h"
#include "HGCal/Geometry/interface/HGCalTBCellVertices.h"
#include "HGCal/Reco/interface/WaferGeometry.h"
#include <fstream>
#include <iostream>
#include <stdio.h>
#include <string>
#include <sstream>
#include "stdlib.h"
#include <vector>

#include "TFile.h"
#include "TBranch.h"
#include "TTree.h"
#include "TDirectory.h"
/**
 *
 *
 *
 *
 *
**/
double ADCtoMIP_CERN[16] =  {17.24, 16.92, 17.51, 16.4, 17.35, 17.49, 16.29, 16.32, 1., 1., 1., 1., 1., 1., 1., 1.};
double MIP2GeV_sim = 51.91e-06; //mpv muon EMM pysics list


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
	
	std::string outputCollectionName;

	std::string runEnergyMapFile;
	std::string inputPathFormat;
	
	std::vector<FileInfo> _fileNames;
	std::vector<FileInfo>::iterator fileIterator;

	int currentRun;
	int currentEvent;

	TFile *rootFile;
  TTree *tree;   //!pointer to the analyzed TTree or TChain
  TDirectory *dir;

  std::vector<unsigned int> *simHitCellIdE;
  std::vector<float>        *simHitCellEnE;
  TBranch                   *b_simHitCellIdE;   
  TBranch                   *b_simHitCellEnE;   

  HGCalTBCellVertices TheCell;

public:
	explicit HGCalTBGenSimSource(const edm::ParameterSet & pset, edm::InputSourceDescription const& desc);
	virtual ~HGCalTBGenSimSource()
	{
		delete rootFile;
			
	}
	
};
