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
#include "HGCal/DataFormats/interface/HGCalTBWireChamberData.h"
#include "HGCal/DataFormats/interface/HGCalTBRecHitCollections.h"
#include "HGCal/Geometry/interface/HGCalTBCellVertices.h"
#include "HGCal/Geometry/interface/HGCalWaferGeometry.h"
#include "HGCal/Geometry/interface/HGCalTBGeometryParameters.h"
#include "HGCal/Geometry/interface/HGCalTBTopology.h"
#include <iostream>
#include <stdio.h>
#include <string>
#include <sstream>
#include "stdlib.h"
#include <vector>
#include <cmath>
#include <unistd.h>

#include "DNN/TensorFlow/interface/TensorFlow.h"
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

class HGCalTBGANSimSource : public edm::ProducerSourceFromFiles
{
private:
	void fillConfiguredRuns(std::fstream& map_file);
	bool setRunAndEventInfo(edm::EventID& id, edm::TimeValue_t& time, edm::EventAuxiliary::ExperimentType&);
	virtual void produce(edm::Event & e);
  virtual void endJob() override;
	void makeRecHit(int, int, int, float, std::unique_ptr<HGCalTBRecHitCollection>&);
	
	std::string RechitOutputCollectionName;
	std::string DWCOutputCollectionName;
	std::string RunDataOutputCollectionName;

	bool m_maskNoisyChannels;
	std::string m_channelsToMask_filename;
	std::vector<int> m_noisyChannels;

	std::vector<double> wc_resolutions;
  int sensorSize;
  int u_max;
  int u_min;
  int v_max;
  int v_min;

  std::string areaSpecification;
  unsigned int beamEnergy;
  int beamParticlePDGID;
  unsigned int setupConfiguration;

	std::string GANModelIndex;
	RUNTYPES _enumPhysicsListUsed;

	std::vector<double> dwc_zPositions;		//filled by area specification

  int NEvents;
  int zDim;	
	int currentEvent;

	//getting the required electronic mapping
	std::string _e_mapFile;

	HGCalTBCellVertices TheCell;
  std::pair<double, double> CellCentreXY;
  HGCalTBTopology HGCalDetectorTopology;

	TRandom* randgen;
  tf::Tensor* energy_tensor;
  tf::Tensor* z_tensor;
  tf::Tensor* simImage;
  tf::Graph* GAN_graph;
  tf::Session* GAN_session;


	int N_layers_HGCal, N_layers_BH;

public:
	explicit HGCalTBGANSimSource(const edm::ParameterSet & pset, edm::InputSourceDescription const& desc);
};
