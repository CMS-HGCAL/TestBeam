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
#include "HGCal/DataFormats/interface/HGCalTBDATURATelescopeData.h"
#include "HGCal/DataFormats/interface/HGCalTBWireChamberData.h"
#include "HGCal/DataFormats/interface/HGCalTBRecHitCollections.h"
#include "HGCal/Geometry/interface/HGCalTBCellVertices.h"
#include "HGCal/Geometry/interface/HGCalWaferGeometry.h"
#include "HGCal/Geometry/interface/HGCalTBGeometryParameters.h"
#include "HGCal/Reco/interface/PositionResolutionHelpers.h"
#include "SimDataFormats/CaloTest/interface/HGCalTestNumbering.h"
#include <fstream>
#include <iostream>
#include <stdio.h>
#include <string>
#include <sstream>
#include "stdlib.h"
#include <vector>
#include <cmath>
#include <unistd.h>

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


struct FileInfo {
	int index;
	double energy;
	int pdgID;
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
	void makeRecHit(int, int, bool, int, double, std::unique_ptr<HGCalTBRecHitCollection>&);
	
	std::string RechitOutputCollectionName;
	std::string RunDataOutputCollectionName;
	
	bool produceDATURATracksInsteadOfDWCs;
	std::string DWCOutputCollectionName;
	std::string DATURAOutputCollectionName;

	bool m_maskNoisyChannels;
	std::string m_channelsToMask_filename;
	std::vector<int> m_noisyChannels;

	double energyNoise;
	double energyNoiseResolution;

	//options related to the dwc
	std::vector<double> wc_resolutions;
	std::vector<double> referenceTracking_zPositions;		//filled by area specification

	//options related to the datura
	std::vector<double> datura_resolutions;
    std::string m_layerPositionFile;
    std::map<int, double> layerPositions;	

  	double beamEnergy;
  	int beamParticlePDGID;
  	unsigned int setupConfiguration;
  	double GeVToMip;
	std::string areaSpecification;

	std::string physicsListUsed;
	RUNTYPES _enumPhysicsListUsed;


	std::vector<FileInfo> _fileNames;
	std::vector<FileInfo>::iterator fileIterator;

	int currentRun;
	int currentEvent;
	int eventCounter;

	TFile *rootFile;
  	TTree *tree;   //!pointer to the analyzed TTree or TChain
  	TDirectory *dir;

  	std::vector<unsigned int> *simHitCellIdEE;
  	std::vector<unsigned int> *simHitCellIdFH;
  	std::vector<unsigned int> *simHitCellIdBH;
  	std::vector<float>        *simHitCellEnEE;
  	std::vector<float>        *simHitCellEnFH;
  	std::vector<float>        *simHitCellEnBH;
  	double					  beamX;
  	double					  beamY;
  	double					  beamP;

  	TBranch                   *b_simHitCellIdEE;   
  	TBranch                   *b_simHitCellEnEE;   
   	TBranch                   *b_simHitCellIdFH;   
  	TBranch                   *b_simHitCellEnFH;   
   	TBranch                   *b_simHitCellIdBH;   
  	TBranch                   *b_simHitCellEnBH;   
  	TBranch                   *b_beamX;   
  	TBranch                   *b_beamY;   
  	TBranch                   *b_beamP;   

	//getting the required electronic mapping
	std::string _e_mapFile;

	HGCalTBCellVertices TheCell;
	HexGeometry* geomc;

	TRandom* randgen;

	unsigned short N_layers_EE, N_layers_FH, N_layers_BH;
	unsigned short firstNLayersFH_asDaisies;


public:
	explicit HGCalTBGenSimSource(const edm::ParameterSet & pset, edm::InputSourceDescription const& desc);
	virtual ~HGCalTBGenSimSource()
	{
		delete rootFile;
			
	}
	
};
