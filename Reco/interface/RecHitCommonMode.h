#ifndef HGCAL_RECO_RECHITCOMMONMODE_H
#define HGCAL_RECO_RECHITCOMMONMODE_H
// user include files
#include <TH1F.h>
#include <TF1.h>
#include <TCanvas.h>

#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "HGCal/DataFormats/interface/HGCalTBRecHitCollections.h"
#include "HGCal/DataFormats/interface/HGCalTBDetId.h"
#include "HGCal/DataFormats/interface/HGCalTBRecHit.h"
#include "HGCal/Geometry/interface/HGCalTBCellVertices.h"
#include "HGCal/Geometry/interface/HGCalTBTopology.h"
#include "HGCal/CondObjects/interface/HGCalElectronicsMap.h"
#include "HGCal/DataFormats/interface/HGCalTBElectronicsId.h"
#include "HGCal/Geometry/interface/HGCalTBGeometryParameters.h"

using namespace std;
//
// class declaration
//



class RecHitCommonMode
{

public:
	RecHitCommonMode(edm::Handle<HGCalTBRecHitCollection> Rechits);
	~RecHitCommonMode();
   
        void evaluate();
        float getCommonModeNoise(HGCalTBDetId id);
        float getCommonModeNoise(int layer, HGCalTBDetId::CellType type, std::string const& = "");

private:

        edm::Handle<HGCalTBRecHitCollection> Rechits_;
	HGCalTBTopology IsCellValid;
	HGCalTBCellVertices TheCell;
        std::string mapfile_ = "HGCal/CondObjects/data/map_FNAL_16Layers.txt";
        struct {
                HGCalElectronicsMap emap_;
        } essource_;
	int sensorsize = 128;// The geometry for a 256 cell sensor hasnt been implemted yet. Need a picture to do this.
	std::pair<double, double> CellCentreXY;

        TH1F* Full_Cell[MAXLAYERS];
	char name[50], title[50];
        float FullCell_CMNoise[MAXLAYERS], HalfCell_CMNoise[MAXLAYERS], CalibPad_CMNoise[MAXLAYERS], MB_CMNoise[MAXLAYERS], MergedCell_CMNoise[MAXLAYERS];
};

#endif
