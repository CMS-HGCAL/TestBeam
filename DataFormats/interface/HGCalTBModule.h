#ifndef HGCAL_TB_MODULE_H
#define HGCAL_TB_MODULE_H

#include <iostream>
#include <set>

class HGCalTBModule
{
 public:
  HGCalTBModule(int layer, std::string det, int sensor_iu, int sensor_iv, int moduleId)
    {
      int detid = det==std::string("EE") ? 0 : 1;
      int iuid,ivid;
      if( sensor_iu>=0 ) iuid=sensor_iu*2;
      else iuid=-2*sensor_iu-1;
      if( sensor_iv>=0 ) ivid=sensor_iv*2;
      else ivid=-2*sensor_iv-1;
      _id = layer; //layer<52 (28 CE-E + 24 CE-H) -> 6bits
      _id |= detid << 0x6; // EE or FH -> 1  bit
      _id |= iuid << 0x7; // 4 bits is far enough
      _id |= ivid << 0xb; // 4 bits is far enough
      _id |= moduleId << 0x10; //
    }
  ~HGCalTBModule(){;}
  void print()
  {
    std::cout << "id = " << id() << "\t"
	      << "layerID = " << layerID() << "\t"
	      << "subDetId = " << subDetID() << "\t"
	      << "sensorIU = " << sensorIU() << "\t"
	      << "sensorIV = " << sensorIV() << "\t"
	      << "moduleID = " << moduleID() << std::endl;
  }
  const int id(){return _id;}
  const int layerID(){return _id&0x3f;}
  const int subDetID(){return (_id&0x40)>>0x6;} // 0 for CE-E, 1 for CE-H
  const int sensorIU(){int tmp=(_id&0x780)>>0x7; return tmp%2==0 ? tmp/2 : -tmp/2-1;}
  const int sensorIV(){int tmp=(_id&0xfff0)>>0xb; return tmp%2==0 ? tmp/2 : -tmp/2-1;}
  const int moduleID(){return _id>>0x10;}
 private:
  unsigned int _id;
};

#endif

