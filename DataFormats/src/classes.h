#include "HGCal/DataFormats/interface/HGCalTBRecHit.h"
#include "HGCal/DataFormats/interface/HGCalTBRecHitCollections.h"
#include "HGCal/DataFormats/interface/HGCalTBCluster.h"
#include "HGCal/DataFormats/interface/HGCalTBClusterCollection.h"
#include "HGCal/DataFormats/interface/HGCalTBCaloTrack.h"
#include "HGCal/DataFormats/interface/HGCalTBCaloTrackCollection.h"
#include "HGCal/DataFormats/interface/HGCalTBDataFrameContainers.h"
#include "HGCal/DataFormats/interface/HGCalTBTrackCollection.h"

#include "DataFormats/Common/interface/RefProd.h"
#include "DataFormats/Common/interface/Wrapper.h"
#include "DataFormats/Common/interface/RefToBase.h"
#include "DataFormats/Common/interface/Holder.h"
#include <vector>

namespace DataFormats_HGCTB
{
struct dictionary {
	HGCalTBRecHit _aRecHit;
	std::vector<HGCalTBRecHit> _HGCTBRHitVect;
//	edm::SortedCollection<HGCalTBRecHit> _theHGCTBRsc;
//	edm::Wrapper< HGCalTBRecHitCollection > _HGCTBeeRHitProd;
	
	reco::HGCalTBCluster _aCluster;
	std::vector<reco::HGCalTBCluster> _HGCTBClusterVect;

	reco::HGCalTBCaloTrack _aCaloTrack;
	std::vector<reco::HGCalTBCaloTrack> _HGCTBCaloTrackVect;
	
	SKIROC2DigiCollection _SR2DC;
	edm::Wrapper<SKIROC2DigiCollection> _theSR2DC;

	HGCalTBTrack _aTrack;
	std::vector<HGCalTBTrack> _HGCTBTRackVect;
//	edm::Wrapper< HGCalTBTrackCollection > _HGCTBTrackProd;

};
}

