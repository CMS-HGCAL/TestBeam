#ifndef HGCAL_PULSEFITTER
#define HGCAL_PULSEFITTER

#include <vector>
#include <cmath>
#include <iostream>

struct PulseFitterResult{
PulseFitterResult() : tmax(0.), amplitude(0.), chi2(1e6), 
    errortmax(0.), erroramplitude(0.), status(-1), ncalls(10000) {;}
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
fitterParameter():tmax0(140),
    tmaxRangeUp(200),
    tmaxRangeDown(30),
    nMaxIterations(100)
  {;}
  double tmax0;
  double tmaxRangeUp;
  double tmaxRangeDown;
  int nMaxIterations;
};

class PulseFitter{
 public:
  PulseFitter( int printLevel, double maxTime=225 , double trise=45 , double ampl_norm=1.85 , double tau=25 , int n_ord=3 );
  ~PulseFitter(){;}
  void run(std::vector<double> &time, std::vector<double> &energy, PulseFitterResult &fit, double noise=-1);
  void setFitterParameter( fitterParameter params ){ m_fitterParameter=params; }
 private:
  int m_printLevel;
  fitterParameter m_fitterParameter;
};

#endif
