#include "HGCal/Reco/interface/WaferGeometry.h"
  

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
      std::cout<<xypos.size()<<"  "<<xpos<<"  "<<ypos<<std::endl;
      xypos.push_back(std::pair<double,double>(xpos,ypos));
    }
    ny += 6;
  }
}

/*
HexGeometry::HexGeometry(bool fine) {
  double alpha = 30.0*M_PI/180.0;

  const int nRows(23);      //row is along the y-axis, we start iterating from negative values
                            //http://cms-hgcal.github.io/TestBeam/coordinates_.html
  int nColumns[nRows] = {4, 5, 4, 5, 6, 5, 6, 7, 6, 7, 8, 7 ,8, 7, 6, 7, 6, 5, 6, 5, 4, 5, 4};


  double dx = 2*HGCAL_TB_CELL::FULL_CELL_SIDE * 10;        //use the correct dimension of the 133 cell sensor in cm
  double dy = dx*cos(alpha);
  

  std::cout<<"DX = "<<dx<<std::endl<<std::endl<<std::endl;

  for (int ir = 0; ir < nRows; ++ir) {
    double ypos = (ir-11)*dy/2.;

    int nColumn = nColumns[ir];
    for (int ic = 0; ic<nColumn; ++ic) {
      double xpos;
      
      if (nColumn % 2) {
       int center_index = (nColumn/2);    //automatic rounding to the lower integer value
       xpos = (ic - center_index) * dx * (1+sin(alpha));
      } else {
        xpos = (ic - (nColumn/2. - 0.5)) * dx * (1+sin(alpha));
      }
      
      xypos.push_back(std::pair<double,double>(xpos,ypos));
      std::cout<<xypos.size()<<"  "<<xpos<<"  "<<ypos<<std::endl;
    }
  }
}
*/


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
