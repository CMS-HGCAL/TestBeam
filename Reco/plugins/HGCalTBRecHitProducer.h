#ifndef HGCALTBRECHITPRODUCER_H
#define HGCALTBRECHITPRODUCER_H

#include "FWCore/Framework/interface/EDProducer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "FWCore/Framework/interface/ESHandle.h"

#include "DataFormats/CaloRecHit/interface/CaloRecHit.h"
#include "HGCal/DataFormats/interface/HGCalTBCommonModeNoise.h"
#include "HGCal/DataFormats/interface/HGCalTBRecHit.h"
#include "HGCal/DataFormats/interface/HGCalTBRecHitCollections.h"
#include "HGCal/DataFormats/interface/HGCalTBDetId.h"
#include "HGCal/DataFormats/interface/HGCalTBRawHitCollection.h"
#include "HGCal/CondObjects/interface/HGCalElectronicsMap.h"
#include "HGCal/CondObjects/interface/HGCalTBDetectorLayout.h"
#include "HGCal/CondObjects/interface/HGCalTBADCConversionsMap.h"
#include "HGCal/CondObjects/interface/HGCalCondObjectTextIO.h"
#include "HGCal/DataFormats/interface/HGCalTBElectronicsId.h"
#include "HGCal/DataFormats/interface/HGCalTBRunData.h" //for the runData type definition
#include "HGCal/DataFormats/interface/HGCalTBGlobalTimestamps.h" //for the runData type definition

#include "HGCal/Geometry/interface/HGCalTBCellVertices.h"

//#include "DNN/Tensorflow/interface/Graph.h"
//#include "DNN/Tensorflow/interface/Tensor.h"

#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "TH2F.h"
#include "TTree.h"
#include <sstream>
#include <fstream>
#include <vector>
#include <map>


//#define DEBUG

class HGCalTBRecHitProducer : public edm::EDProducer
{
 public:
  HGCalTBRecHitProducer(const edm::ParameterSet&);
  virtual void produce(edm::Event&, const edm::EventSetup&);
 private:
  virtual void beginJob() override;
  std::string m_CommonModeNoiseCollectionName;
  std::string m_outputCollectionName;
  std::string m_electronicMap;
  std::string m_detectorLayoutFile;
  std::string m_adcCalibrationsFile;
  double m_timeSample3ADCCut;

  edm::EDGetTokenT<HGCalTBRawHitCollection> m_HGCalTBRawHitCollection;
  edm::EDGetTokenT<RunData> RunDataToken; 
  edm::EDGetTokenT<HGCalTBGlobalTimestamps> HGCalTBGlobalTimestampsToken; 

  bool m_maskNoisyChannels;
  std::string m_channelsToMask_filename;
  int m_NHexaBoards;

  bool investigatePulseShape;
  std::map<int, TH2F*> shapesLG;
  std::map<int, TH2F*> shapesHG;
  std::map<int, TH2F*> ToARisevsTMaxLG;
  std::map<int, TH2F*> ToARisevsTMaxHG;
  std::map<int, TH2F*> ToAFallvsTMaxLG;
  std::map<int, TH2F*> ToAFallvsTMaxHG;
  std::map<int, TH2F*> TMaxHGvsTMaxLG;

  std::vector<int> m_noisyChannels;

  std::map<int, TH2F*> m_h_HighVsLowGainAmpl;
  std::map<int, TH2F*> m_h_LowGainVsTOTAmpl;

  TTree* outtree;
  int event_for_tree;
  int run_for_tree;
  int pdgID_for_tree;
  double beamEnergy_for_tree;
  int board_for_tree;
  int skiroc_for_tree;
  int channel_for_tree;
  uint32_t globalTimestamp_for_tree;
  double MIPEnergy_for_tree;
  double HG_max_for_tree;
  double LG_max_for_tree;
  double TMax_HG_for_tree;
  double TMax_LG_for_tree;
  double TOT_for_tree;
  double TOA_rise_for_tree;
  double TOA_fall_for_tree;
  double chi2_HG_for_tree;
  double chi2_LG_for_tree;
  double trise_HG_for_tree;
  double trise_LG_for_tree;
  double errortmax_HG_for_tree;
  double errortmax_LG_for_tree;
  double erroramplitude_HG_for_tree;
  double erroramplitude_LG_for_tree;


  std::pair<double, double> CellCentreXY;
  HGCalTBCellVertices TheCell;

  //https://gitlab.cern.ch/mrieger/CMSSW-DNN/tree/80X
  std::string timingNNFilePath;
  std::map<int, std::string> timingNNFilePaths;
  //dnn::tf::Graph* timingCalibrationNN;
  //dnn::tf::Tensor* xNNInput;
  //dnn::tf::Tensor* xIN;
  //dnn::tf::Tensor* yNNOutput;
  //dnn::tf::Tensor* yIN;


  #ifdef DEBUG
    int eventCounter;
  #endif

  struct {
    HGCalElectronicsMap emap_;
    HGCalTBDetectorLayout layout_;
    HGCalTBADCConversionsMap adccalibmap_;
  } essource_;
  
  void setupTimingNNs(std::string);
  
};

#endif
