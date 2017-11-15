#ifndef HGCAL_PULSEFITTER
#define HGCAL_PULSEFITTER

#include <vector>
#include <cmath>
#include <iostream>

#include <Math/Minimizer.h>
#include <Math/Factory.h>
#include <Math/Functor.h>

struct PulseFitterResult{
PulseFitterResult() : tmax(0.), amplitude(0.),
    errortmax(0.), erroramplitude(0.), status(-1) {;}
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

struct fitterParameter{
fitterParameter():tmax0(90),
    tmaxRangeUp(120),
    tmaxRangeDown(65),
    nMaxIterations(100)
  {;}
  double tmax0;
  double tmaxRangeUp;
  double tmaxRangeDown;
  int nMaxIterations;
};

class PulseFitter{
 public:
  PulseFitter(  int printLevel=1 , double maxTime=225. , double alpha=10 , double trise=50 );
  ~PulseFitter(){;}
  void run(std::vector<double> &time, std::vector<double> &energy, PulseFitterResult &fit, double noise=-1);
  void setFitterParameter( fitterParameter params ){ m_fitterParameter=params; }
 private:
  int m_printLevel;
  fitterParameter m_fitterParameter;

};

double parabolicFit(std::vector<double>, std::vector<double>);


#endif
