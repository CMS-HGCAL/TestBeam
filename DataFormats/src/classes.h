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
#include "HGCal/DataFormats/interface/HGCalTBCommonModeNoise.h"
#include "HGCal/DataFormats/interface/HGCalTBRunData.h"
#include "HGCal/DataFormats/interface/HGCalTBGlobalTimestamps.h"
#include "HGCal/DataFormats/interface/HGCalTBDWCTrack.h"
#include "HGCal/DataFormats/interface/HGCalTBWireChamberData.h"
#include "HGCal/DataFormats/interface/HGCalTBDATURATelescopeData.h"

#include "DataFormats/Common/interface/RefProd.h"
#include "DataFormats/Common/interface/Wrapper.h"
#include "DataFormats/Common/interface/RefToBase.h"
#include "DataFormats/Common/interface/Holder.h"
#include "FWCore/Utilities/interface/TypeWithDict.h"
#include <vector>

namespace DataFormats_HGCTB
{
struct dictionary {
	HGCalTBRecHit _aRecHit;
	std::vector<HGCalTBRecHit> _HGCTBRHitVect;
	std::pair<Float16_t, Float16_t> _aFloat16_tPair;
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

  	commonModeNoise _acommonModeNoise;
  	std::map<int, commonModeNoise> _acommonModeNoiseMap;
  	edm::Wrapper<std::map<int, commonModeNoise> > _acommonModeNoiseMapWrapper;

	UserRecords<bool> _aBooleanUserRecord;
	UserRecords<float> _aFloatUserRecord;
	UserRecords<double> _aDoubleUserRecord;
	UserRecords<int> _aIntUserRecord;
  	edm::Wrapper<UserRecords<double> > _aDoubleUserRecordWrapper;
	UserRecords<std::vector<short> > _aShortVectorUserRecord;
	UserRecords<std::vector<float> > _aFloatVectorUserRecord;

	RunData _aRunData;
	edm::Wrapper<RunData> _aRunDataWrapper;

	HGCalTBGlobalTimestamps _aHGCalTBGlobalTimestamps;
	std::map<int, uint32_t> askiroc_to_timestamps;
	edm::Wrapper<HGCalTBGlobalTimestamps> _aHGCalTBGlobalTimestampsWrapper;

	HGCalTBDWCTrack _aHGCalTBDWCTrack;
	edm::Wrapper<HGCalTBDWCTrack> _aHGCalTBDWCTrackWrapper;

  	WireChamberData _aWireChamberData;
  	std::map<int, WireChamberData> _aWireChamberDataMap;
  	edm::Wrapper<std::map<int, WireChamberData> > _aWireChamberDataWrapper;

  	HGCalTBDATURATelescopeData _aHGCalTBDATURATelescopeData;
  	std::vector<HGCalTBDATURATelescopeData> _HGCTBHGCalTBDATURATelescopeDataVect;
  	edm::Wrapper<std::vector<HGCalTBDATURATelescopeData> > _aHGCalTBDATURATelescopeDataVecWrapper;

	HGCalTBTrack _aTrack;
	std::vector<HGCalTBTrack> _HGCTBTRackVect;
//	edm::Wrapper< HGCalTBTrackCollection > _HGCTBTrackProd;

};
}

