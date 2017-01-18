#ifndef HGCAL_RECO_RESHELPERS_H
#define HGCAL_RECO_RESHELPERS_H

#include <iostream>
#include <utility>
#include <algorithm>
#include <cmath>
#include <fstream>

#include "TF1.h"
#include "TGraphErrors.h"
#include "TFitResult.h"


void parseAlignmentFile(std::map<int, double> &alignmentParameters, std::string path);

//
// class declarations
//

//class that performs an analytical straight line fit
class LineFitter {
  private:
    std::vector<double> _x;
    std::vector<double> _y;
    std::vector<double> _sigma_y;
    double _S, _S_x, _S_xx, _S_xy, _S_y, _Delta;
  public:
    LineFitter(std::vector<double> x, std::vector<double> y, std::vector<double> sigma_y);
    void addPoint(double x, double y, double sigma_y);
    void fit();
    bool converged();
    double getM();
    double getMError();
    double getB();
    double getBError();
    double getMBCovariance();

    double eval(double x);
    double evalError(double x);
};

#endif