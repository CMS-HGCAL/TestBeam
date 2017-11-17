#include "HGCal/DataFormats/interface/HGCalTBRecHit.h"
#include "HGCal/DataFormats/interface/HGCalTBRecHitCollections.h"
#include "HGCal/DataFormats/interface/HGCalTBCluster.h"
#include "HGCal/DataFormats/interface/HGCalTBClusterCollection.h"
#include "HGCal/DataFormats/interface/HGCalTBSkiroc2CMS.h"
#include "HGCal/DataFormats/interface/HGCalTBSkiroc2CMSCollection.h"
#include "HGCal/DataFormats/interface/HGCalTBRawHit.h"
#include "HGCal/DataFormats/interface/HGCalTBRawHitCollection.h"
#include "HGCal/DataFormats/interface/HGCalTBDetId.h"
#include "HGCal/DataFormats/interface/HGCalTBDataFrameContainers.h"
#include "HGCal/DataFormats/interface/HGCalTBTrackCollection.h"
#include "HGCal/DataFormats/interface/HGCalTBRunData.h"
#include "HGCal/DataFormats/interface/HGCalTBDWCTrack.h"
#include "HGCal/DataFormats/interface/HGCalTBWireChamberData.h"

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

	HGCalTBSkiroc2CMS _aSki;
	std::vector<HGCalTBSkiroc2CMS> _skiVect;
	
	HGCalTBRawHit _aRawHit;
	std::vector<HGCalTBRawHit> _rawHitVec;

	SKIROC2DigiCollection _SR2DC;
	edm::Wrapper<SKIROC2DigiCollection> _theSR2DC;

	UserRecords<bool> _aBooleanUserRecord;
	UserRecords<double> _aDoubleUserRecord;

	RunData _aRunData;
	edm::Wrapper<RunData> _aRunDataWrapper;

	HGCalTBDWCTrack _aHGCalTBDWCTrack;
	edm::Wrapper<HGCalTBDWCTrack> _aHGCalTBDWCTrackWrapper;

	WireChamberData _aWireChamberData;
	std::vector<WireChamberData> _aWireChamberDataVector;
	edm::Wrapper<WireChamberData> __aWireChamberDataWrapper;

	HGCalTBTrack _aTrack;
	std::vector<HGCalTBTrack> _HGCTBTRackVect;
//	edm::Wrapper< HGCalTBTrackCollection > _HGCTBTrackProd;

};
}

