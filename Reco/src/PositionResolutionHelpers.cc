#include "HGCal/Reco/interface/PositionResolutionHelpers.h"


//****   Parsing of alignment values    ****//


    
    
AlignmentParameters::AlignmentParameters(std::vector<std::string> files, double _defaultRun) : defaultRun(_defaultRun) {
  boost::regex run_regex("^.*RUN_([0-9]+)_.*$");
  boost::smatch run_matches;

  for(size_t i = 0; i<files.size(); i++) {
    std::string filePath = files[i];

    std::cout<<filePath<<std::endl;
    if (boost::regex_match(filePath, run_matches, run_regex)) {
      std::string run_str(run_matches[1].first, run_matches[1].second);

      _params[atof(run_str.c_str())] = parseFile(filePath);
    }

  }
};

AlignmentParameters::AlignmentParameters(std::vector<std::string> files) : AlignmentParameters(files, 887){};

std::map<int, double> AlignmentParameters::parseFile(std::string path) {
  std::map<int, double> fileParams;
  std::fstream file;
  
  char fragment[100];
  int readCounter = -2, currentParameter = 0;
  file.open(path.c_str(), std::fstream::in);

  while (file.is_open() && !file.eof()) {
    if (readCounter!=-2) readCounter++;
    file >> fragment;
    if (std::string(fragment)=="111") readCounter = 0;  //first parameter is read out

    if (readCounter==0) currentParameter = atoi(fragment);
    if (readCounter==1) currentParameter = fileParams[currentParameter] = atof(fragment);
    if (readCounter==4) readCounter = -1;
  }

  if (readCounter==-2) {
    for (int i=1; i<= 8; i++) {
      fileParams[i*100+11] = 0;
      fileParams[i*100+12] = 0;
      fileParams[i*100+13] = 0;
      fileParams[i*100+21] = 0;
      fileParams[i*100+22] = 0;
      fileParams[i*100+23] = 0;
    }
  }

  return fileParams;
}

double AlignmentParameters::getValue(double run, int paramId, bool tryDefault) {
  if (_params.find(run) == _params.end() && tryDefault) return getValue(defaultRun, paramId, false);
  else if (_params[run].find(paramId) == _params[run].end()) return 0.;
  else return _params[run][paramId];
};



//****   Line Fiting Class    ****//

LineFitter::LineFitter(std::vector<double> x, std::vector<double> y, std::vector<double> sigma_y) {
  if (x.size() != y.size()) {
    std::cout<<"LineFitter class: x and y vectors have different dimension!"<<std::endl;
    return;
  }
  if (x.size() != sigma_y.size()) {
    std::cout<<"LineFitter class: x and sigma_y vectors have different dimension!"<<std::endl;
    return;
  }
  if (y.size() != sigma_y.size()) {
    std::cout<<"LineFitter class: y and sigma_y vectors have different dimension!"<<std::endl;
    return;
  }
  _x = x; _y = y; _sigma_y = sigma_y;
  _S_x = _S_xx = _S_y = _S_xy = 0.;
  _S = 1.;
}

void LineFitter::addPoint(double x, double y, double sigma_y) {
  _x.push_back(x); _y.push_back(y); _sigma_y.push_back(sigma_y);
}

void LineFitter::fit() {
  if (_x.size()==0 || _y.size()==0 || _sigma_y.size()==0) {
    std::cout<<"One of the input vectors has zero dimension!"<<std::endl;
  }
  _S_x = _S_xx = _S_y = _S_xy = _S = _Delta = 0.;

  for (size_t i=0; i<_x.size(); i++) {
    if (_sigma_y[i] == 0) continue;
    _S += 1.0/pow(_sigma_y[i], 2);
    _S_x += _x[i]/pow(_sigma_y[i], 2);
    _S_y += _y[i]/pow(_sigma_y[i], 2);
    _S_xx += _x[i]*_x[i]/pow(_sigma_y[i], 2);
    _S_xy += _x[i]*_y[i]/pow(_sigma_y[i], 2);
  }

  _Delta = _S*_S_xx-_S_x*_S_x;
}
bool LineFitter::converged() {
  return (_Delta > 0.);
}
double LineFitter::getM() {
  return (_S*_S_xy-_S_x*_S_y)/_Delta;
}
double LineFitter::getMError() {
  return sqrt(_S/_Delta);
}
double LineFitter::getB() {
  return (_S_xx*_S_y-_S_x*_S_xy)/_Delta;
}
double LineFitter::getBError() {
  return sqrt(_S_xx/_Delta);
}
double LineFitter::getMBCovariance() {
  return - _S_x/_Delta;
}

double LineFitter::eval(double x) {
  if (!converged()) return 0;
  else 
    return this->getB() + x * this->getM();
};
double LineFitter::evalError(double x) {
  if (!converged()) return 0;
  else
    return sqrt(pow(this->getBError(), 2) + pow(x*this->getMError(),2) + 2*fabs(x)*this->getMBCovariance());
};