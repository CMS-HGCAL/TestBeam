#ifndef HGCalTBGlobalTimestamps_H
#define HGCalTBGlobalTimestamps_H
#include <string>
#include <map>


struct HGCalTBGlobalTimestamps {
  std::map<int, uint32_t> skiroc_to_timestamps;
};

#endif