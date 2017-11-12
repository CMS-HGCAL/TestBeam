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
  void add( HGCalTBModule module ){ _modules.push_back(module); }
  HGCalTBLayer at(int i){ return _layers.at(i); }
  HGCalTBModule atBoard(int i) {return _modules.at(i); }
  // need to be sure  how the vector _layers has been filled (should be ok if geometry file filled with logical order) : layout.at( hit.detid().layer()-1 )
  // -1 needed if elec map starts from 1
  std::vector<HGCalTBLayer>& layers(){return _layers;}
  std::vector<HGCalTBModule>& modules(){return _modules;}
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
  std::vector<HGCalTBModule> _modules;
};

#endif
