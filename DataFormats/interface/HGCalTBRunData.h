#ifndef HGCalTBRunData_H
#define HGCalTBRunData_H
#include <string>

struct RunData {
  explicit RunData(int config, int r, double e, std::string rt): configuration(config), run(r), energy(e), runType(rt) {};
  RunData(): configuration(0), run(0), energy(0), runType(""), trueEnergy(-1){};
  int configuration;
  int run;
  int event;
  double energy;
  std::string runType;
  bool hasValidMWCMeasurement;
  bool hasDanger;
  double trueEnergy;
};

typedef std::map<int, RunData> runMap;


#endif