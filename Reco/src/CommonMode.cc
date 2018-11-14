#include <HGCal/Reco/interface/CommonMode.h>
#include <HGCal/Geometry/interface/HGCalTBTopology.h>
#include <HGCal/Geometry/interface/HGCalTBGeometryParameters.h>

CommonMode::CommonMode( HGCalElectronicsMap &emap, bool useMedian, bool cmPerChip, float threshold, int expectedMaxTS ) : _useMedian(useMedian),
															  _cmPerChip(cmPerChip),
															  _threshold(threshold),
															  _expectedMaxTS(expectedMaxTS)
{
  _emap = emap;
}

void CommonMode::Evaluate( edm::Handle<HGCalTBRawHitCollection> hits )
{
  if ( _useMedian && _cmPerChip ) EvaluateMedianPerChip( hits );
  else if ( _useMedian  ) EvaluateMedianPerBoard( hits );
  else if ( _cmPerChip  ) EvaluateAveragePerChip( hits );
  else EvaluateMedianPerBoardWithThr( hits );
}

void CommonMode::EvaluateMedianPerChip( edm::Handle<HGCalTBRawHitCollection> hits )
{
  std::vector<float> cmMapHG[HGCAL_TB_GEOMETRY::NUMBER_OF_HEXABOARD * HGCAL_TB_GEOMETRY::N_SKIROC_PER_HEXA * NUMBER_OF_TIME_SAMPLES];
  std::vector<float> cmMapLG[HGCAL_TB_GEOMETRY::NUMBER_OF_HEXABOARD * HGCAL_TB_GEOMETRY::N_SKIROC_PER_HEXA * NUMBER_OF_TIME_SAMPLES];
  std::vector<float> cmMapHGHalf[HGCAL_TB_GEOMETRY::NUMBER_OF_HEXABOARD * HGCAL_TB_GEOMETRY::N_SKIROC_PER_HEXA * NUMBER_OF_TIME_SAMPLES];
  std::vector<float> cmMapLGHalf[HGCAL_TB_GEOMETRY::NUMBER_OF_HEXABOARD * HGCAL_TB_GEOMETRY::N_SKIROC_PER_HEXA * NUMBER_OF_TIME_SAMPLES];
  for ( int iski = 0; iski < HGCAL_TB_GEOMETRY::NUMBER_OF_HEXABOARD * HGCAL_TB_GEOMETRY::N_SKIROC_PER_HEXA; iski++ ) {
    for ( int it = 0; it < NUMBER_OF_TIME_SAMPLES; it++ ) {
      cmMapHG[NUMBER_OF_TIME_SAMPLES * iski + it].reserve(32);
      cmMapLG[NUMBER_OF_TIME_SAMPLES * iski + it].reserve(32);
      cmMapHGHalf[NUMBER_OF_TIME_SAMPLES * iski + it].reserve(32);
      cmMapLGHalf[NUMBER_OF_TIME_SAMPLES * iski + it].reserve(32);
    }
  }
  for ( auto hit : *hits ) {
    HGCalTBElectronicsId eid( _emap.detId2eid(hit.detid().rawId()) );
    if ( !_emap.existsEId(eid) ) continue;
    int iski = hit.skiroc();
    for ( size_t it = 0; it < NUMBER_OF_TIME_SAMPLES; it++ ) {
      if ( hit.detid().cellType() == 0 || hit.detid().cellType() == 4 ) {
        cmMapHG[NUMBER_OF_TIME_SAMPLES * iski + it].push_back(hit.highGainADC(it));
        cmMapLG[NUMBER_OF_TIME_SAMPLES * iski + it].push_back(hit.lowGainADC(it));
      }
      else if ( hit.detid().cellType() != 5 ) {
        cmMapHGHalf[NUMBER_OF_TIME_SAMPLES * iski + it].push_back(hit.highGainADC(it));
        cmMapLGHalf[NUMBER_OF_TIME_SAMPLES * iski + it].push_back(hit.lowGainADC(it));
      }
    }
  }

  HGCalTBTopology topo;
  for ( int iski = 0; iski < HGCAL_TB_GEOMETRY::NUMBER_OF_HEXABOARD * HGCAL_TB_GEOMETRY::N_SKIROC_PER_HEXA; iski++ ) {
    for ( int it = 0; it < NUMBER_OF_TIME_SAMPLES; it++ ) {
      int index = NUMBER_OF_TIME_SAMPLES * iski + it;
      std::sort( cmMapHG[index].begin(), cmMapHG[index].end() );
      std::sort( cmMapLG[index].begin(), cmMapLG[index].end() );
      std::sort( cmMapHGHalf[index].begin(), cmMapHGHalf[index].end() );
      std::sort( cmMapLGHalf[index].begin(), cmMapLGHalf[index].end() );
      if ( _cmMap.find( iski ) == _cmMap.end() ) {
        commonModeNoise cm;
        uint16_t size = cmMapHG[index].size();
        cm.fullHG[it] = size % 2 == 0 ? cmMapHG[index][size / 2 - 1] : cmMapHG[index][size / 2] ;
        cm.outerHG[it] = cm.fullHG[it];
        cm.mergedHG[it] = cm.fullHG[it] * 1.5;
        cm.fullLG[it] = size % 2 == 0 ? cmMapLG[index][size / 2 - 1] : cmMapLG[index][size / 2] ;
        cm.outerLG[it] = cm.fullLG[it];
        cm.mergedLG[it] = cm.fullLG[it] * 1.5;

        size = cmMapHGHalf[index].size();
        cm.halfHG[it] = size % 2 == 0 ? cmMapHGHalf[index][size / 2 - 1] : cmMapHGHalf[index][size / 2] ;
        cm.mouseBiteHG[it] = cm.halfHG[it] * topo.Cell_Area(3) / topo.Cell_Area(2);
        cm.halfLG[it] = size % 2 == 0 ? cmMapLGHalf[index][size / 2 - 1] : cmMapLGHalf[index][size / 2] ;
        cm.mouseBiteLG[it] = cm.halfLG[it] * topo.Cell_Area(3) / topo.Cell_Area(2);
        std::pair<int, commonModeNoise> p(iski, cm);
        _cmMap.insert(p);
      }
      else {
        uint16_t size = cmMapHG[index].size();
        _cmMap[iski].fullHG[it] = size % 2 == 0 ? cmMapHG[index][size / 2 - 1] : cmMapHG[index][size / 2] ;
        _cmMap[iski].outerHG[it] = _cmMap[iski].fullHG[it];
        _cmMap[iski].mergedHG[it] = _cmMap[iski].fullHG[it] * 1.5;
        _cmMap[iski].fullLG[it] = size % 2 == 0 ? cmMapLG[index][size / 2 - 1] : cmMapLG[index][size / 2] ;
        _cmMap[iski].outerLG[it] = _cmMap[iski].fullLG[it];
        _cmMap[iski].mergedLG[it] = _cmMap[iski].fullLG[it] * 1.5;
        size = cmMapHGHalf[index].size();
        _cmMap[iski].halfHG[it] = size % 2 == 0 ? cmMapHGHalf[index][size / 2 - 1] : cmMapHGHalf[index][size / 2] ;
        _cmMap[iski].mouseBiteHG[it] = _cmMap[iski].halfHG[it] * topo.Cell_Area(3) / topo.Cell_Area(2);
        _cmMap[iski].halfLG[it] = size % 2 == 0 ? cmMapLGHalf[index][size / 2 - 1] : cmMapLGHalf[index][size / 2] ;
        _cmMap[iski].mouseBiteLG[it] = _cmMap[iski].halfLG[it] * topo.Cell_Area(3) / topo.Cell_Area(2);
      }
    }
  }
}

void CommonMode::EvaluateMedianPerBoard( edm::Handle<HGCalTBRawHitCollection> hits )
{
  std::vector<float> cmMapHG[HGCAL_TB_GEOMETRY::NUMBER_OF_HEXABOARD * NUMBER_OF_TIME_SAMPLES];
  std::vector<float> cmMapLG[HGCAL_TB_GEOMETRY::NUMBER_OF_HEXABOARD * NUMBER_OF_TIME_SAMPLES];
  std::vector<float> cmMapHGHalf[HGCAL_TB_GEOMETRY::NUMBER_OF_HEXABOARD * NUMBER_OF_TIME_SAMPLES];
  std::vector<float> cmMapLGHalf[HGCAL_TB_GEOMETRY::NUMBER_OF_HEXABOARD * NUMBER_OF_TIME_SAMPLES];
  for ( int iboard = 0; iboard < HGCAL_TB_GEOMETRY::NUMBER_OF_HEXABOARD; iboard++ ) {
    for ( int it = 0; it < NUMBER_OF_TIME_SAMPLES; it++ ) {
      cmMapHG[NUMBER_OF_TIME_SAMPLES * iboard + it].reserve(128);
      cmMapLG[NUMBER_OF_TIME_SAMPLES * iboard + it].reserve(128);
      cmMapHGHalf[NUMBER_OF_TIME_SAMPLES * iboard + it].reserve(128);
      cmMapLGHalf[NUMBER_OF_TIME_SAMPLES * iboard + it].reserve(128);
    }
  }
  for ( auto hit : *hits ) {
    HGCalTBElectronicsId eid( _emap.detId2eid(hit.detid().rawId()) );
    if ( !_emap.existsEId(eid) ) continue;
    int iboard = hit.skiroc()/HGCAL_TB_GEOMETRY::N_SKIROC_PER_HEXA;
    for ( size_t it = 0; it < NUMBER_OF_TIME_SAMPLES; it++ ) {
      if ( hit.detid().cellType() == 0 || hit.detid().cellType() == 4 ) {
        cmMapHG[NUMBER_OF_TIME_SAMPLES * iboard + it].push_back(hit.highGainADC(it));
        cmMapLG[NUMBER_OF_TIME_SAMPLES * iboard + it].push_back(hit.lowGainADC(it));
      }
      else if ( hit.detid().cellType() != 5 ) {
        cmMapHGHalf[NUMBER_OF_TIME_SAMPLES * iboard + it].push_back(hit.highGainADC(it));
        cmMapLGHalf[NUMBER_OF_TIME_SAMPLES * iboard + it].push_back(hit.lowGainADC(it));
      }
    }
  }

  HGCalTBTopology topo;
  for ( int iboard = 0; iboard < HGCAL_TB_GEOMETRY::NUMBER_OF_HEXABOARD; iboard++ ) {
    for ( int it = 0; it < NUMBER_OF_TIME_SAMPLES; it++ ) {
      int index = NUMBER_OF_TIME_SAMPLES * iboard + it;
      std::sort( cmMapHG[index].begin(), cmMapHG[index].end() );
      std::sort( cmMapLG[index].begin(), cmMapLG[index].end() );
      std::sort( cmMapHGHalf[index].begin(), cmMapHGHalf[index].end() );
      std::sort( cmMapLGHalf[index].begin(), cmMapLGHalf[index].end() );
      if ( _cmMap.find( iboard ) == _cmMap.end() ) {
        commonModeNoise cm;
        uint16_t size = cmMapHG[index].size();
        cm.fullHG[it] = size % 2 == 0 ? cmMapHG[index][size / 2 - 1] : cmMapHG[index][size / 2] ;
        cm.outerHG[it] = cm.fullHG[it];
        cm.mergedHG[it] = cm.fullHG[it] * 1.5;
        cm.fullLG[it] = size % 2 == 0 ? cmMapLG[index][size / 2 - 1] : cmMapLG[index][size / 2] ;
        cm.outerLG[it] = cm.fullLG[it];
        cm.mergedLG[it] = cm.fullLG[it] * 1.5;

        size = cmMapHGHalf[index].size();
        cm.halfHG[it] = size % 2 == 0 ? cmMapHGHalf[index][size / 2 - 1] : cmMapHGHalf[index][size / 2] ;
        cm.mouseBiteHG[it] = cm.halfHG[it] * topo.Cell_Area(3) / topo.Cell_Area(2);
        cm.halfLG[it] = size % 2 == 0 ? cmMapLGHalf[index][size / 2 - 1] : cmMapLGHalf[index][size / 2] ;
        cm.mouseBiteLG[it] = cm.halfLG[it] * topo.Cell_Area(3) / topo.Cell_Area(2);
        std::pair<int, commonModeNoise> p(iboard, cm);
        _cmMap.insert(p);
      }
      else {
        uint16_t size = cmMapHG[index].size();
        _cmMap[iboard].fullHG[it] = size % 2 == 0 ? cmMapHG[index][size / 2 - 1] : cmMapHG[index][size / 2] ;
        _cmMap[iboard].outerHG[it] = _cmMap[iboard].fullHG[it];
        _cmMap[iboard].mergedHG[it] = _cmMap[iboard].fullHG[it] * 1.5;
        _cmMap[iboard].fullLG[it] = size % 2 == 0 ? cmMapLG[index][size / 2 - 1] : cmMapLG[index][size / 2] ;
        _cmMap[iboard].outerLG[it] = _cmMap[iboard].fullLG[it];
        _cmMap[iboard].mergedLG[it] = _cmMap[iboard].fullLG[it] * 1.5;
        size = cmMapHGHalf[index].size();
        _cmMap[iboard].halfHG[it] = size % 2 == 0 ? cmMapHGHalf[index][size / 2 - 1] : cmMapHGHalf[index][size / 2] ;
        _cmMap[iboard].mouseBiteHG[it] = _cmMap[iboard].halfHG[it] * topo.Cell_Area(3) / topo.Cell_Area(2);
        _cmMap[iboard].halfLG[it] = size % 2 == 0 ? cmMapLGHalf[index][size / 2 - 1] : cmMapLGHalf[index][size / 2] ;
        _cmMap[iboard].mouseBiteLG[it] = _cmMap[iboard].halfLG[it] * topo.Cell_Area(3) / topo.Cell_Area(2);
      }
    }
  }
}

void CommonMode::EvaluateAveragePerChip( edm::Handle<HGCalTBRawHitCollection> hits )
{
}

void CommonMode::EvaluateMedianPerBoardWithThr( edm::Handle<HGCalTBRawHitCollection> hits )
{
  std::vector<float> cmMapHG[HGCAL_TB_GEOMETRY::NUMBER_OF_HEXABOARD * NUMBER_OF_TIME_SAMPLES];
  std::vector<float> cmMapLG[HGCAL_TB_GEOMETRY::NUMBER_OF_HEXABOARD * NUMBER_OF_TIME_SAMPLES];
  for ( int iboard = 0; iboard < HGCAL_TB_GEOMETRY::NUMBER_OF_HEXABOARD; iboard++ ) {
    for ( int it = 0; it < NUMBER_OF_TIME_SAMPLES; it++ ) {
      cmMapHG[NUMBER_OF_TIME_SAMPLES * iboard + it].reserve(128);
      cmMapLG[NUMBER_OF_TIME_SAMPLES * iboard + it].reserve(128);
    }
  }
  for ( auto hit : *hits ) {
    HGCalTBElectronicsId eid( _emap.detId2eid(hit.detid().rawId()) );
    if ( !_emap.existsEId(eid) ) continue;
    int iboard = hit.skiroc()/HGCAL_TB_GEOMETRY::N_SKIROC_PER_HEXA;
    if ( hit.highGainADC(_expectedMaxTS)>_threshold ) continue;
    for ( size_t it = 0; it < NUMBER_OF_TIME_SAMPLES; it++ ) {
      cmMapHG[NUMBER_OF_TIME_SAMPLES * iboard + it].push_back(hit.highGainADC(it));
      cmMapLG[NUMBER_OF_TIME_SAMPLES * iboard + it].push_back(hit.lowGainADC(it));
    }
  }
  
  HGCalTBTopology topo;
  for ( int iboard = 0; iboard < HGCAL_TB_GEOMETRY::NUMBER_OF_HEXABOARD; iboard++ ) {
    for ( int it = 0; it < NUMBER_OF_TIME_SAMPLES; it++ ) {
      int index = NUMBER_OF_TIME_SAMPLES * iboard + it;
      std::sort( cmMapHG[index].begin(), cmMapHG[index].end() );
      std::sort( cmMapLG[index].begin(), cmMapLG[index].end() );
      if ( _cmMap.find( iboard ) == _cmMap.end() ) {
        commonModeNoise cm;
        uint16_t size = cmMapHG[index].size();
	if( size>1 ){
	  cm.fullHG[it] = size % 2 == 0 ? cmMapHG[index][size / 2 - 1] : cmMapHG[index][size / 2] ;
	  cm.outerHG[it] = cm.fullHG[it];
	  cm.mergedHG[it] = cm.fullHG[it] * 1.5;
	  cm.fullLG[it] = size % 2 == 0 ? cmMapLG[index][size / 2 - 1] : cmMapLG[index][size / 2] ;
	  cm.outerLG[it] = cm.fullLG[it];
	  cm.mergedLG[it] = cm.fullLG[it] * 1.5;
	  cm.halfHG[it] = cm.fullHG[it] * topo.Cell_Area(2) / topo.Cell_Area(0);
	  cm.mouseBiteHG[it] = cm.fullHG[it] * topo.Cell_Area(3) / topo.Cell_Area(0);
	  cm.halfLG[it] = cm.fullLG[it] * topo.Cell_Area(2) / topo.Cell_Area(0);
	  cm.mouseBiteLG[it] = cm.fullLG[it] * topo.Cell_Area(3) / topo.Cell_Area(0);
	}
	else{
	  cm.fullHG[it] = 0;
	  cm.outerHG[it] = 0;
	  cm.mergedHG[it] = 0;
	  cm.fullLG[it] = 0;
	  cm.outerLG[it] = 0;
	  cm.mergedLG[it] = 0;
	  cm.halfHG[it] = 0;
	  cm.mouseBiteHG[it] = 0;
	  cm.halfLG[it] = 0;
	  cm.mouseBiteLG[it] = 0;
	}
        std::pair<int, commonModeNoise> p(iboard, cm);
        _cmMap.insert(p);
      }
      else {
        uint16_t size = cmMapHG[index].size();
	if( size!= 0 ){
	  _cmMap[iboard].fullHG[it] = size % 2 == 0 ? cmMapHG[index][size / 2 - 1] : cmMapHG[index][size / 2] ;
	  _cmMap[iboard].outerHG[it] = _cmMap[iboard].fullHG[it];
	  _cmMap[iboard].mergedHG[it] = _cmMap[iboard].fullHG[it] * 1.5;
	  _cmMap[iboard].fullLG[it] = size % 2 == 0 ? cmMapLG[index][size / 2 - 1] : cmMapLG[index][size / 2] ;
	  _cmMap[iboard].outerLG[it] = _cmMap[iboard].fullLG[it];
	  _cmMap[iboard].mergedLG[it] = _cmMap[iboard].fullLG[it] * 1.5;
	  _cmMap[iboard].halfHG[it] = _cmMap[iboard].fullHG[it] * topo.Cell_Area(2) / topo.Cell_Area(0);
	  _cmMap[iboard].mouseBiteHG[it] = _cmMap[iboard].fullHG[it] * topo.Cell_Area(3) / topo.Cell_Area(0);
	  _cmMap[iboard].halfLG[it] = _cmMap[iboard].fullLG[it] * topo.Cell_Area(2) / topo.Cell_Area(0);
	  _cmMap[iboard].mouseBiteLG[it] = _cmMap[iboard].fullLG[it] * topo.Cell_Area(3) / topo.Cell_Area(0);
	}
	else{
	  _cmMap[iboard].fullHG[it] = 0;
	  _cmMap[iboard].outerHG[it] = 0;
	  _cmMap[iboard].mergedHG[it] = 0;
	  _cmMap[iboard].fullLG[it] = 0;
	  _cmMap[iboard].outerLG[it] = 0;
	  _cmMap[iboard].mergedLG[it] = 0;
	  _cmMap[iboard].halfHG[it] = 0;
	  _cmMap[iboard].mouseBiteHG[it] = 0;
	  _cmMap[iboard].halfLG[it] = 0;
	  _cmMap[iboard].mouseBiteLG[it] = 0;
	}
      }
    }
  }
}
