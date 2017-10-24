#ifndef HGCalTBDWCTrack_H
#define HGCalTBDWCTrack_H
#include <utility>


class HGCalTBDWCTrack {    //simple straight lines
  public:
    HGCalTBDWCTrack(): valid(false), b_x(0), m_x(0), b_y(0), m_y(0){};
    bool valid;

    double b_x;
    double m_x;
    double b_y;
    double m_y;

    int referenceType;
    int N_points;
    double chi2_x;
    double chi2_y;

    std::pair<double, double> DWCExtrapolation_XY(double z) const{
      return std::make_pair(b_x+z*m_x, b_y+z*m_y);
    }

    int NDWCTrackPoints() const {
      return N_points;
    }
};


#endif