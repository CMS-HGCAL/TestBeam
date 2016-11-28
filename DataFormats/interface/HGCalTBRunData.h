#ifndef HGCalTBRunData_H
#define HGCalTBRunData_H
#include <string>

struct RunData {
  explicit RunData(int config, int r, double e, std::string rt): configuration(config), run(r), energy(e), runType(rt) {};
  RunData(): configuration(0), run(0), energy(0), runType(""){};
  int configuration;
  int run;
  double energy;
  std::string runType;
};

typedef std::map<int, RunData> runMap;


#endif