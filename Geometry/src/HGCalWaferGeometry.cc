#include "HGCal/Geometry/interface/HGCalWaferGeometry.h"
  

HexGeometry::HexGeometry(bool fine) {
  const int nC(15), nF(20);
  int nCoarse(11), nyCoarse(-42), nFine(15), nyFine(-56);
  int cellCoarse[nC] = {2,5,8,11,12,11,12,11,12,11,12,11,8,5,2};
  int cellFine[nF] = {3,6,9,12,15,16,15,16,15,16,15,16,15,16,15,14,11,8,5,2};
  double wafer(123.7);

  int    rows = (fine) ? nF : nC;
  double cell = (fine) ? wafer/nFine : wafer/nCoarse;
  double dx   = 0.5*cell;
  double dy   = 0.5*dx*tan(30.0*M_PI/180.0);
  int    ny   = (fine) ? nyFine : nyCoarse;
  for (int ir = 0; ir < rows; ++ir) {
    int    column = (fine) ? cellFine[ir] : cellCoarse[ir];
    int    nx     = 1 - column;
    double ypos   = dy*ny;
    for (int ic = 0; ic<column; ++ic) {
      double xpos = dx*nx;
      nx += 2;
      double x_rot = cos(90.0*M_PI/180.0) * xpos + sin(90.0*M_PI/180.0) * ypos; 
      double y_rot = -sin(90.0*M_PI/180.0) * xpos + cos(90.0*M_PI/180.0) * ypos;
      xypos.push_back(std::pair<double,double>(x_rot,y_rot));

      int type = 0;
      if (ir==0) type = 2;
      else if (ir==1 && (ic==0 || ic==4)) type = 2;
      else if (ir==2 && (ic==0 || ic==7)) type = 2;
      else if (ir==3 && (ic==0 || ic==10)) type = 2;
      else if (ir==4 && (ic==0 || ic==11)) type = 2;
      else if (ir==6 && (ic==0 || ic==11)) type = 2;
      else if (ir==8 && (ic==0 || ic==11)) type = 2;
      else if (ir==10 && (ic==0 || ic==11)) type = 2;
      else if (ir==11 && (ic==0 || ic==10)) type = 2;
      else if (ir==12 && (ic==0 || ic==7)) type = 2;
      else if (ir==13 && (ic==0 || ic==4)) type = 2;
      else if (ir==14) type = 2;


      types.push_back(type);
    }
    ny += 6;
  }
}


HexGeometry::~HexGeometry() {}

std::pair<double,double> HexGeometry::position(const int cell) {
  std::pair<double,double> xy;
  if (cell >= 0 && cell < (int)(xypos.size())) {
    xy = xypos[cell];
  } else {
    xy = std::pair<double,double>(0,0);
  }
  return xy;
}

int HexGeometry::cellType(const int cell) {
  int t = -1;
  if (cell >= 0 && cell < (int)(types.size())) {
    t = types[cell];
  } 
  return t; 
  
}

void testGeometry() {

  HexGeometry geomc(false);
  for (int k = 0; k < 133; ++k) {
    std::pair<double,double> xy = geomc.position(k);
    std::cout << "Coarse Cell[" << k << "] " << xy.first << ":" << xy.second
        << std::endl;
  }

  HexGeometry geomf(true);
  for (int k = 0; k < 240; ++k) {
    std::pair<double,double> xy = geomf.position(k);
    std::cout << "Fine Cell[" << k << "] " << xy.first << ":" << xy.second
        << std::endl;
  }
}
