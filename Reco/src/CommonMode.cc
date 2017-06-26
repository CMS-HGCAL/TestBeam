#include <HGCal/Reco/interface/CommonMode.h>

CommonMode::CommonMode( HGCalElectronicsMap &emap, bool useMedian, bool cmPerChip, float threshold ) : _useMedian(useMedian),
												      _cmPerChip(cmPerChip),
												      _threshold(threshold)
{
  _emap=emap;
}

void CommonMode::Evaluate( edm::Handle<HGCalTBRawHitCollection> hits )
{
  if( _useMedian && _cmPerChip ) EvaluateMedianPerChip( hits );
  else if( _useMedian  ) EvaluateMedianPerLayer( hits );
  else if( _cmPerChip  ) EvaluateAveragePerChip( hits );
  else EvaluateAveragePerLayer( hits );
}

void CommonMode::EvaluateMedianPerChip( edm::Handle<HGCalTBRawHitCollection> hits )
{
  std::map< int,std::vector<float> > cmMapHG;
  std::map< int,std::vector<float> > cmMapLG;
  for( auto hit : *hits ){
    HGCalTBElectronicsId eid( _emap.detId2eid(hit.detid().rawId()) );
    if( !_emap.existsEId(eid) ) continue;
    int iski=eid.iskiroc();
    for( size_t it=0; it<NUMBER_OF_TIME_SAMPLES; it++ ){
      if( cmMapHG.find(100*iski+it)==cmMapHG.end() ){
	std::vector<float> vech; vech.push_back( hit.highGainADC(it) );
	std::vector<float> vecl; vecl.push_back( hit.lowGainADC(it) );
	std::pair< int,std::vector<float> > ph(100*iski+it,vech);
	std::pair< int,std::vector<float> > pl(100*iski+it,vecl);
	cmMapHG.insert(ph);
	cmMapLG.insert(pl);
      }
      else{
	cmMapHG[100*iski+it].push_back( hit.highGainADC(it) );
	cmMapLG[100*iski+it].push_back( hit.lowGainADC(it) );
      }
    }
  }

  for( std::map< int,std::vector<float> >::iterator it=cmMapHG.begin(); it!=cmMapHG.end(); ++it ){
    std::sort( it->second.begin(),it->second.end() );
    uint16_t size = it->second.size();
    uint16_t ts=it->first%100;
    uint16_t iski=it->first/100;
    if( _cmMap.find( iski )==_cmMap.end() ){
      commonModeNoise cm;
      cm.fullHG[ts] = size%2==0 ? it->second.at(size/2-1) : it->second.at(size/2) ;
      cm.outerHG[ts] = cm.fullHG[ts];
      cm.mouseBiteHG[ts] = cm.fullHG[ts]/2;
      cm.halfHG[ts] = cm.fullHG[ts]/2;
      std::pair<int,commonModeNoise> p(iski,cm);
      _cmMap.insert(p);
    }
    else{
      _cmMap[iski].fullHG[ts] = size%2==0 ? it->second.at(size/2-1) : it->second.at(size/2) ;
      _cmMap[iski].outerHG[ts] = _cmMap[iski].fullHG[ts];
      _cmMap[iski].mouseBiteHG[ts] = _cmMap[iski].fullHG[ts]/2;
      _cmMap[iski].halfHG[ts] = _cmMap[iski].fullHG[ts]/2;
    }
  }
  for( std::map< int,std::vector<float> >::iterator it=cmMapLG.begin(); it!=cmMapLG.end(); ++it ){
    std::sort( it->second.begin(),it->second.end() );
    uint16_t size = it->second.size();
    uint16_t ts=it->first%100;
    uint16_t iski=it->first/100;
    if( _cmMap.find( iski )==_cmMap.end() ){
      std::cout << "the _cmMap should already contain element with key (skirocId) = " << iski << " -> exit(1) " << std::endl;
      exit(1);
    }
    _cmMap[iski].fullLG[ts] = size%2==0 ? it->second.at(size/2-1) : it->second.at(size/2) ;
    _cmMap[iski].outerLG[ts] = _cmMap[iski].fullLG[ts];
    _cmMap[iski].mouseBiteLG[ts] = _cmMap[iski].fullLG[ts]/2;
    _cmMap[iski].halfLG[ts] = _cmMap[iski].fullLG[ts]/2;
  }
}

void CommonMode::EvaluateMedianPerLayer( edm::Handle<HGCalTBRawHitCollection> hits )
{
  std::map< int,std::vector<float> > cmMapHG;
  std::map< int,std::vector<float> > cmMapLG;
  for( auto hit : *hits ){
    HGCalTBElectronicsId eid( _emap.detId2eid(hit.detid().rawId()) );
    if( !_emap.existsEId(eid) ) continue;
    int ilayer=hit.detid().layer();
    for( size_t it=0; it<NUMBER_OF_TIME_SAMPLES; it++ ){
      if( cmMapHG.find(100*ilayer+it)==cmMapHG.end() ){
	std::vector<float> vech; vech.push_back( hit.highGainADC(it) );
	std::vector<float> vecl; vecl.push_back( hit.lowGainADC(it) );
	std::pair< int,std::vector<float> > ph(100*ilayer+it,vech);
	std::pair< int,std::vector<float> > pl(100*ilayer+it,vecl);
	cmMapHG.insert(ph);
	cmMapLG.insert(pl);
      }
      else{
	cmMapHG[100*ilayer+it].push_back( hit.highGainADC(it) );
	cmMapLG[100*ilayer+it].push_back( hit.lowGainADC(it) );
      }
    }
  }
  for( std::map< int,std::vector<float> >::iterator it=cmMapHG.begin(); it!=cmMapHG.end(); ++it ){
    std::sort( it->second.begin(),it->second.end() );
    uint16_t size = it->second.size();
    uint16_t ts=it->first%100;
    uint16_t ilayer=it->first/100;
    if( _cmMap.find( ilayer )==_cmMap.end() ){
      commonModeNoise cm;
      cm.fullHG[ts] = size%2==0 ? it->second.at(size/2-1) : it->second.at(size/2) ;
      cm.outerHG[ts] = cm.fullHG[ts];
      cm.mouseBiteHG[ts] = cm.fullHG[ts]/2;
      cm.halfHG[ts] = cm.fullHG[ts]/2;
      std::pair<int,commonModeNoise> p(ilayer,cm);
      _cmMap.insert(p);
    }
    else{
      _cmMap[ilayer].fullHG[ts] = size%2==0 ? it->second.at(size/2-1) : it->second.at(size/2) ;
      _cmMap[ilayer].outerHG[ts] = _cmMap[ilayer].fullHG[ts];
      _cmMap[ilayer].mouseBiteHG[ts] = _cmMap[ilayer].fullHG[ts]/2;
      _cmMap[ilayer].halfHG[ts] = _cmMap[ilayer].fullHG[ts]/2;
    }
  }
  for( std::map< int,std::vector<float> >::iterator it=cmMapLG.begin(); it!=cmMapLG.end(); ++it ){
    std::sort( it->second.begin(),it->second.end() );
    uint16_t size = it->second.size();
    uint16_t ts=it->first%100;
    uint16_t ilayer=it->first/100;
    if( _cmMap.find( ilayer )==_cmMap.end() ){
      std::cout << "the _cmMap should already contain element with key (layerId) = " << ilayer << " -> exit(1) " << std::endl;
      exit(1);
    }
    _cmMap[ilayer].fullLG[ts] = size%2==0 ? it->second.at(size/2-1) : it->second.at(size/2) ;
    _cmMap[ilayer].outerLG[ts] = _cmMap[ilayer].fullLG[ts];
    _cmMap[ilayer].mouseBiteLG[ts] = _cmMap[ilayer].fullLG[ts]/2;
    _cmMap[ilayer].halfLG[ts] = _cmMap[ilayer].fullLG[ts]/2;
  }
}

void CommonMode::EvaluateAveragePerChip( edm::Handle<HGCalTBRawHitCollection> hits )
{
}

void CommonMode::EvaluateAveragePerLayer( edm::Handle<HGCalTBRawHitCollection> hits )
{
}

