#include "HGCal/Geometry/interface/HGCalWaferGeometry.h"
  

HexGeometry::HexGeometry(bool fine) {

  initCellGeom(fine);
  initWaferGeom();

}


/////cell geometry
void HexGeometry::initCellGeom(bool fine){
  const int nC(15), nF(20);
  //int nCoarse(11), nyCoarse(-42), nFine(15), nyFine(-56);
  int nCoarse(11), nyCoarse(21), nFine(15), nyFine(28);
  int cellCoarse[nC] = {2,5,8,11,12,11,12,11,12,11,12,11,8,5,2};
  int cellFine[nF] = {3,6,9,12,15,16,15,16,15,16,15,16,15,16,15,14,11,8,5,2};
  double wafer(19*HGCAL_TB_CELL::FULL_CELL_SIDE);

  int    rows = (fine) ? nF : nC;
  double cell = (fine) ? wafer/nFine : wafer/nCoarse;
  double dx   = 0.5*cell;
  //double dy   = 0.5*dx*tan(30.0*M_PI/180.0);
  double dy   = 0.5*cell*tan(30.0*M_PI/180.0);
  int    ny   = (fine) ? nyFine : nyCoarse;
  for (int ir = 0; ir < rows; ++ir) {
    int    column = (fine) ? cellFine[ir] : cellCoarse[ir];
    int    nx     = 1 - column;
    double ypos   = dy*ny;
    for (int ic = 0; ic<column; ++ic) {
      double xpos = dx*nx;
      nx += 2;
      //xypos.push_back(std::pair<double,double>(xpos,ypos));
      xypos_cell.push_back(std::pair<double,double>(ypos,xpos));  ///currently we have rotated the geometry by 90 degrees so x becomes y and y becomes x.

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

      types_cell.push_back(type);

    }
    //ny += 6;
    ny -= 3;
  }
  std::cout << "Initialize HexGeometry for a cell " << xypos_cell.size() << " cells"
      << std::endl;

}


/////wafer geometry
void HexGeometry::initWaferGeom(){
  const int nC(3);
  int nCoarse(3), nyCoarse(3);
  int cellCoarse[nC] = {2,3,2};
  double wafer(19*HGCAL_TB_CELL::FULL_CELL_SIDE*3);

  int    rows = nC;
  double cell = wafer/nCoarse;
  double dx   = 0.5*cell;
  double dy   = 0.5*cell*tan(30.0*M_PI/180.0);
  int    ny   = nyCoarse;
  for (int ir = 0; ir < rows; ++ir) {
    int    column = cellCoarse[ir];
    int    nx     = 1 - column;
    double ypos   = dy*ny;
    for (int ic = 0; ic<column; ++ic) {
      double xpos = dx*nx;
      nx += 2;
      //xypos.push_back(std::pair<double,double>(xpos,ypos));
      xypos_wafer.push_back(std::pair<double,double>(ypos,xpos));  ///currently we have rotated the geometry by 90 degrees so x becomes y and y becomes x.
    }
    //ny += 6;
    ny -= 3;
  }
  std::cout << "Initialize HexGeometry for a wafer " << xypos_wafer.size() << " wafer"
      << std::endl;

}

int HexGeometry::cellType(const int cell) {
  int t = -1;
  if (cell >= 0 && cell < (int)(types_cell.size())) {
    t = types_cell[cell];
  } 
  return t; 
  
}

std::pair<double,double> HexGeometry::position_cell(const int cell) {
  std::pair<double,double> xy;
  if (cell >= 0 && cell < (int)(xypos_cell.size())) {
    xy = xypos_cell[cell];
  } else {
    xy = std::pair<double,double>(0,0);
  }
  return xy;
}


std::pair<double,double> HexGeometry::position_wafer(const int wafer) {
  std::pair<double,double> xy;
  if (wafer >= 0 && wafer < (int)(xypos_wafer.size())) {
    xy = xypos_wafer[wafer];
  } else {
    xy = std::pair<double,double>(0,0);
  }
  return xy;
}


void testGeometry() {

  HexGeometry geomc(false);
  for (int k = 0; k < 133; ++k) {
    std::pair<double,double> xy = geomc.position_cell(k);
    std::cout << "Coarse Cell[" << k << "] " << xy.first << ":" << xy.second
        << std::endl;
  }


  for (int k = 0; k < 7; ++k) {
    std::pair<double,double> xy_w = geomc.position_wafer(k);
    std::cout << "Wafer[" << k << "] " << xy_w.first << ":" << xy_w.second
        << std::endl;
  }

  /*HexGeometry geomf(true);
  for (int k = 0; k < 240; ++k) {
    std::pair<double,double> xy = geomf.position(k);
    std::cout << "Fine Cell[" << k << "] " << xy.first << ":" << xy.second
        << std::endl;
  }
  */


}