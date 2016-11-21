#ifndef HGCalTBRunData_H
#define HGCalTBRunData_H
#include <string>

struct RunData {
  explicit RunData(int r, double e, double lt, std::string rt): run(r), energy(e), layerThickness(lt), runType(rt) {};
  RunData(): run(0), energy(0), layerThickness(0), runType(""){};
  int run;
  double energy;
  double layerThickness;
  std::string runType;
};

typedef std::map<int, RunData> runMap;


#endif