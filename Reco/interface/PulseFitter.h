#ifndef HGCAL_PULSEFITTER
#define HGCAL_PULSEFITTER

#include <vector>
#include <cmath>
#include <iostream>

#include <Math/Minimizer.h>
#include <Math/Factory.h>
#include <Math/Functor.h>

struct PulseFitterResult{
PulseFitterResult() : tmax(0.), trise(0.), alpha(0.), amplitude(0.),
    errortmax(0.), errortrise(0.), erroralpha(0.), erroramplitude(0.) {;}
  double tmax;
  double trise;
  double alpha;
  double amplitude;
  double chi2;
  double errortmax;
  double errortrise;
  double erroralpha;
  double erroramplitude;
  double errorchi2;
  int status;
  int ncalls;
};

class PulseFitter{
 public:
  PulseFitter(  int printLevel=1 , double maxTime=225.);
  ~PulseFitter();
  void run(std::vector<double> &time, std::vector<double> &energy, PulseFitterResult &fit);
 private:
  int _printLevel;
  ROOT::Math::Minimizer* m;
  const double *xm;
  const double *errors;
};

#endif
