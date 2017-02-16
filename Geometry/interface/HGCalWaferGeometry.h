#include <algorithm>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <fstream>
#include <string>
#include <vector>

#include "HGCal/Geometry/interface/HGCalTBCellParameters.h"
  

//this code patches the fact that the digitization step of the simulation is not implemented yet
//it serves as a map of simulated hits to the actual x-y positions

class HexGeometry {

public :
  HexGeometry(bool fine);
  ~HexGeometry();

  std::pair<double,double> position(const int cell);
  int cellType(const int cell);


private :
  std::vector<std::pair<double,double> > xypos;
  std::vector<int> types;

};

void testGeometry(); 
