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
	RecHitCommonMode(edm::Handle<HGCalTBRecHitCollection> Rechits, HGCalElectronicsMap& emap);
	~RecHitCommonMode();
   
        void evaluate(float maxEcut = 1000.0);
        float getGaussCommonModeNoise(HGCalTBDetId id);
        float getGaussCommonModeNoise(int idSKIROC, int type);
        float getMeanCommonModeNoise(HGCalTBDetId id);
        float getMeanCommonModeNoise(int idSKIROC, int type);
        

private:

        edm::Handle<HGCalTBRecHitCollection> Rechits_;
        struct {
                HGCalElectronicsMap emap_;
        } essource_;

        TH1F* Full_Cell[MAXSKIROCS];
	char name[50], title[50];
        float FullCell_CMNoise_Fit[MAXSKIROCS];
        float FullCell_CMNoise_Mean[MAXSKIROCS];
        float CalibPad_CMNoise[MAXSKIROCS];
        float HalfCell_CMNoise[MAXSKIROCS];
        float MB_CMNoise[MAXSKIROCS];
        float OuterCalibPad_CMNoise[MAXSKIROCS];
        float MergedCell_CMNoise[MAXSKIROCS];
};

#endif
