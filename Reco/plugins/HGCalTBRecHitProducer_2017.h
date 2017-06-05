#ifndef HGCALTBRECHITPRODUCER_2017_H
#define HGCALTBRECHITPRODUCER_2017_H
/** \class Reco/plugins/HGCalTBRecHitProducer_2017.h HGCalTBRecHitProducer_2017 HGCalTBRecHitProducer_2017
  \brief
  \author Thorben Quast
 */

#include "FWCore/Framework/interface/EDProducer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "FWCore/Framework/interface/ESHandle.h"

#include "HGCal/DataFormats/interface/HGCalTBRecHitCollections.h"
//#include "HGCal/DataFormats/interface/HGCalTBDetId.h"
#include "HGCal/DataFormats/interface/HGCalTBRawHitCollection.h"


// here are defined objects containing pedestals and ADCtoGeV factors
#include "HGCal/CondObjects/interface/HGCalCondObjects.h"
#include "HGCal/CondObjects/interface/HGCalCondObjectTextIO.h"
#include "HGCal/CondObjects/interface/HGCalTBNumberingScheme.h"
#include "HGCal/CondObjects/interface/HGCalElectronicsMap.h"
#include "HGCal/DataFormats/interface/HGCalTBElectronicsId.h"
#include "HGCal/Geometry/interface/HGCalTBCellVertices.h"
//#define DEBUG


#ifdef DEBUG
#include <iostream>
#endif

class HGCalTBRecHitProducer_2017 : public edm::EDProducer
{

public:
  HGCalTBRecHitProducer_2017(const edm::ParameterSet&);
  virtual void produce(edm::Event&, const edm::EventSetup&);
private:
  std::string outputCollectionName;     ///<label name of collection made by this producer
  edm::EDGetTokenT<HGCalTBRawHitCollection> _rawHitsToken;
  std::string _pedestalLow_filename, _pedestalHigh_filename;
  std::string _gainsLow_filename, _gainsHigh_filename;
  int _adcSaturation;
  std::vector<double> _LG2HG_value;
  std::string _mapFile;
  int _layers_config;
  int _N_timestamps;
  struct {
    HGCalElectronicsMap emap_;
        } essource_;
  HGCalTBCellVertices TheCell;
  std::pair<double, double> CellCenterXY;
};


#endif