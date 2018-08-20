#ifndef HGCAL_TB_DETECTOR_LAYOUT_H
#define HGCAL_TB_DETECTOR_LAYOUT_H

#include <HGCal/DataFormats/interface/HGCalTBModule.h>
#include <HGCal/DataFormats/interface/HGCalTBLayer.h>
#include <iostream>
#include <set>
#include <algorithm>

class HGCalTBDetectorLayout{
 public:
  HGCalTBDetectorLayout(){;}
  ~HGCalTBDetectorLayout(){;}
    
  void add( HGCalTBLayer layer ){ _layers.push_back(layer); }
  HGCalTBLayer at(int i){ return _layers.at(i); }
  // need to be sure  how the vector _layers has been filled (should be ok if geometry file filled with logical order) : layout.at( hit.detid().layer()-1 )
  // -1 needed if elec map starts from 1
  HGCalTBLayer getLayerWithModuleIndex( int index )
  {
    for( auto layer : _layers ){
      for( auto module: layer.modules() )
	if( module.moduleID()==index )
	  return layer;
    }
    std::cout << "Major problems in HGCalTBLayer getLayerWithModuleIndex( int index ) -> should have already returned a HGCalTBLayer" << std::endl;
    return HGCalTBLayer(0,0,0);
  }
  std::vector<HGCalTBLayer>& layers(){return _layers;}
  int nlayers(){
    std::set<int> l;
    for( auto il : _layers )
      l.insert(il.layerID());
    return l.size();
  }
  bool layerExists( HGCalTBLayer layer )
  {
    if( std::find(_layers.begin(), _layers.end(), layer)!=_layers.end() )
      return true;
    else return false;
  }
 private:
  std::vector<HGCalTBLayer> _layers;
};

#endif
