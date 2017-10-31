#include "HGCal/Geometry/interface/HGCalTBCellParameters.h"
  

#include <algorithm>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <fstream>
#include <string>
#include <vector>
  
class HexGeometry {

public :
  HexGeometry(bool fine);
  virtual ~HexGeometry() {}

  void initCellGeom(bool fine);
  void initWaferGeom();
  
  std::pair<double,double> position_cell(const int cell);
  int cellType(const int cell);
  std::pair<double,double> position_wafer(const int wafer);

private :
  std::vector<std::pair<double,double> > xypos_cell;
  std::vector<int> types_cell;
  std::vector<std::pair<double,double> > xypos_wafer;

};
