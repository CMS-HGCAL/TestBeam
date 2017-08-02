#include <HGCal/Reco/interface/CommonMode.h>
#include <HGCal/Geometry/interface/HGCalTBTopology.h>

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
  std::map< int,std::vector<float> > cmMapHGHalf;
  std::map< int,std::vector<float> > cmMapLGHalf;
  for( auto hit : *hits ){
    HGCalTBElectronicsId eid( _emap.detId2eid(hit.detid().rawId()) );
    if( !_emap.existsEId(eid) ) continue;
    int iski=hit.skiroc();
    if( hit.detid().cellType()==0 || hit.detid().cellType()==4 ) {
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
    else{
      for( size_t it=0; it<NUMBER_OF_TIME_SAMPLES; it++ ){
    	if( cmMapHGHalf.find(100*iski+it)==cmMapHGHalf.end() ){
    	  std::vector<float> vech; vech.push_back( hit.highGainADC(it) );
    	  std::vector<float> vecl; vecl.push_back( hit.lowGainADC(it) );
    	  std::pair< int,std::vector<float> > ph(100*iski+it,vech);
    	  std::pair< int,std::vector<float> > pl(100*iski+it,vecl);
    	  cmMapHGHalf.insert(ph);
    	  cmMapLGHalf.insert(pl);
    	}
    	else{
    	  cmMapHGHalf[100*iski+it].push_back( hit.highGainADC(it) );
    	  cmMapLGHalf[100*iski+it].push_back( hit.lowGainADC(it) );
    	}
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
      // cm.mouseBiteHG[ts] = cm.fullHG[ts]/2;
      // cm.halfHG[ts] = cm.fullHG[ts]/2;
      std::pair<int,commonModeNoise> p(iski,cm);
      _cmMap.insert(p);
    }
    else{
      _cmMap[iski].fullHG[ts] = size%2==0 ? it->second.at(size/2-1) : it->second.at(size/2) ;
      _cmMap[iski].outerHG[ts] = _cmMap[iski].fullHG[ts];
      //_cmMap[iski].mouseBiteHG[ts] = _cmMap[iski].fullHG[ts]/2;
      //_cmMap[iski].halfHG[ts] = _cmMap[iski].fullHG[ts]/2;
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
    // _cmMap[iski].mouseBiteLG[ts] = _cmMap[iski].fullLG[ts]/2;
    // _cmMap[iski].halfLG[ts] = _cmMap[iski].fullLG[ts]/2;
  }
  HGCalTBTopology topo;
  for( std::map< int,std::vector<float> >::iterator it=cmMapHGHalf.begin(); it!=cmMapHGHalf.end(); ++it ){
    std::sort( it->second.begin(),it->second.end() );
    uint16_t size = it->second.size();
    uint16_t ts=it->first%100;
    uint16_t iski=it->first/100;
    if( _cmMap.find( iski )==_cmMap.end() ){
      std::cout << "the _cmMap should already contain element with key (skirocId) = " << iski << " -> exit(1) " << std::endl;
      exit(1);
    }
    _cmMap[iski].halfHG[ts] = size%2==0 ? it->second.at(size/2-1) : it->second.at(size/2) ;
    _cmMap[iski].mouseBiteHG[ts] = _cmMap[iski].halfHG[ts]*topo.Cell_Area(3)/topo.Cell_Area(2);
  }
  for( std::map< int,std::vector<float> >::iterator it=cmMapLGHalf.begin(); it!=cmMapLGHalf.end(); ++it ){
    std::sort( it->second.begin(),it->second.end() );
    uint16_t size = it->second.size();
    uint16_t ts=it->first%100;
    uint16_t iski=it->first/100;
    if( _cmMap.find( iski )==_cmMap.end() ){
      std::cout << "the _cmMap should already contain element with key (skirocId) = " << iski << " -> exit(1) " << std::endl;
      exit(1);
    }
    _cmMap[iski].halfLG[ts] = size%2==0 ? it->second.at(size/2-1) : it->second.at(size/2) ;
    _cmMap[iski].mouseBiteLG[ts] = _cmMap[iski].halfLG[ts]*topo.Cell_Area(3)/topo.Cell_Area(2);
  }
}

void CommonMode::EvaluateMedianPerLayer( edm::Handle<HGCalTBRawHitCollection> hits )
{
  std::map< int,std::vector<float> > cmMapHG;
  std::map< int,std::vector<float> > cmMapLG;
  std::map< int,std::vector<float> > cmMapHGHalf;
  std::map< int,std::vector<float> > cmMapLGHalf;
  for( auto hit : *hits ){
    HGCalTBElectronicsId eid( _emap.detId2eid(hit.detid().rawId()) );
    if( !_emap.existsEId(eid) ) continue;
    int ilayer=hit.detid().layer();
    if( hit.detid().cellType()==0 || hit.detid().cellType()==4 ) {
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
    else{
      for( size_t it=0; it<NUMBER_OF_TIME_SAMPLES; it++ ){
	if( cmMapHGHalf.find(100*ilayer+it)==cmMapHGHalf.end() ){
	  std::vector<float> vech; vech.push_back( hit.highGainADC(it) );
	  std::vector<float> vecl; vecl.push_back( hit.lowGainADC(it) );
	  std::pair< int,std::vector<float> > ph(100*ilayer+it,vech);
	  std::pair< int,std::vector<float> > pl(100*ilayer+it,vecl);
	  cmMapHGHalf.insert(ph);
	  cmMapLGHalf.insert(pl);
	}
	else{
	  cmMapHGHalf[100*ilayer+it].push_back( hit.highGainADC(it) );
	  cmMapLGHalf[100*ilayer+it].push_back( hit.lowGainADC(it) );
	}
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
      std::pair<int,commonModeNoise> p(ilayer,cm);
      _cmMap.insert(p);
    }
    else{
      _cmMap[ilayer].fullHG[ts] = size%2==0 ? it->second.at(size/2-1) : it->second.at(size/2) ;
      _cmMap[ilayer].outerHG[ts] = _cmMap[ilayer].fullHG[ts];
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
  }
  for( std::map< int,std::vector<float> >::iterator it=cmMapHGHalf.begin(); it!=cmMapHGHalf.end(); ++it ){
    std::sort( it->second.begin(),it->second.end() );
    uint16_t size = it->second.size();
    uint16_t ts=it->first%100;
    uint16_t ilayer=it->first/100;
    if( _cmMap.find( ilayer )==_cmMap.end() ){
      std::cout << "the _cmMap should already contain element with key (layerId) = " << ilayer << " -> exit(1) " << std::endl;
      exit(1);
    }
    else{
      _cmMap[ilayer].mouseBiteHG[ts] = size%2==0 ? it->second.at(size/2-1) : it->second.at(size/2) ;
      _cmMap[ilayer].halfHG[ts] = _cmMap[ilayer].mouseBiteHG[ts];
    }
  }
  for( std::map< int,std::vector<float> >::iterator it=cmMapLGHalf.begin(); it!=cmMapLGHalf.end(); ++it ){
    std::sort( it->second.begin(),it->second.end() );
    uint16_t size = it->second.size();
    uint16_t ts=it->first%100;
    uint16_t ilayer=it->first/100;
    if( _cmMap.find( ilayer )==_cmMap.end() ){
      std::cout << "the _cmMap should already contain element with key (layerId) = " << ilayer << " -> exit(1) " << std::endl;
      exit(1);
    }
    else{
      _cmMap[ilayer].mouseBiteLG[ts] = size%2==0 ? it->second.at(size/2-1) : it->second.at(size/2) ;
      _cmMap[ilayer].halfLG[ts] = _cmMap[ilayer].mouseBiteLG[ts];
    }
  }
}

void CommonMode::EvaluateAveragePerChip( edm::Handle<HGCalTBRawHitCollection> hits )
{
}

void CommonMode::EvaluateAveragePerLayer( edm::Handle<HGCalTBRawHitCollection> hits )
{
}

