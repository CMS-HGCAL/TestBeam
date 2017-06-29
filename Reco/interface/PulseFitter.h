#ifndef HGCAL_PULSEFITTER
#define HGCAL_PULSEFITTER

#include <vector>
#include <cmath>
#include <iostream>

struct PulseFitterResult{
PulseFitterResult() : tmax(0.), amplitude(0.),
    errortmax(0.), erroramplitude(0.) {;}
  double tmax;
  double amplitude;
  double chi2;
  double trise;
  double errortmax;
  double erroramplitude;
  double errorchi2;
  int status;
  int ncalls;
};

class PulseFitter{
 public:
  PulseFitter(  int printLevel=1 , double maxTime=225. , double alpha=10 , double trise=50 );
  ~PulseFitter(){;}
  void run(std::vector<double> &time, std::vector<double> &energy, PulseFitterResult &fit, double noise=-1);
 private:
  int _printLevel;
};

#endif
