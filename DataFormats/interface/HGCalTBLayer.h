#ifndef HGCAL_TB_LAYER_H
#define HGCAL_TB_LAYER_H

#include <HGCal/DataFormats/interface/HGCalTBModule.h>
#include <iostream>
#include <set>

class HGCalTBLayer
{
 public:
  HGCalTBLayer(int layerId, float z, std::string subdet)
    {
      _layerId=layerId;
      _z=z;
      _subdet=subdet==std::string("EE") ? 0 : 1;
    }
  HGCalTBLayer(int layerId, float z, std::string subdet, std::vector<HGCalTBModule> modules)
    {
      _layerId=layerId;
      _z=z;
      _modules=modules;
      _subdet=subdet==std::string("EE") ? 0 : 1;
    }
  ~HGCalTBLayer(){;}

  std::vector<HGCalTBModule> &modules(){return _modules;}
  void add( HGCalTBModule module ){ _modules.push_back(module); }
  HGCalTBModule at(int i){ return _modules.at(i); }
  HGCalTBModule at(int IU, int IV){
    for(std::vector<HGCalTBModule>::iterator it=_modules.begin(); it!=_modules.end(); ++it)
      if( (*it).sensorIU()==IU && (*it).sensorIV()==IV )
	return (*it);
    for(std::vector<HGCalTBModule>::iterator it=_modules.begin(); it!=_modules.end(); ++it)
      (*it).print();
    std::cout << "Major problems in method at(int IU, int IV) from HGCalTBLayer.h (" << IU << "," << IV << ") is not valid -> exit(1)" << std::endl;
    exit(1);
    return _modules.at(0);
  }
  const float z(){return _z;}
  int layerID() const {return _layerId;}
  const int subdet(){return _subdet;}
  void print()
  {
    for( auto module : _modules ){
      std::cout << "z = " << _z << "\t";
      module.print();
    }
  }
 private:
  std::vector<HGCalTBModule> _modules;
  int _layerId;
  int _subdet; // 0 for EE, 1 for FH
  float _z;
};

inline bool operator==(const HGCalTBLayer& rhs0, const HGCalTBLayer& rhs1) { 
  return (rhs0.layerID() == rhs1.layerID()); 
} 

#endif
