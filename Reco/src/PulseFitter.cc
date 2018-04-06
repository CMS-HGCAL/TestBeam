#include <HGCal/Reco/interface/PulseFitter.h>

#include <limits>

#include <Math/Minimizer.h>
#include <Math/Factory.h>
#include <Math/Functor.h>

const int numTS = 8;
double _time[numTS],_energy[numTS];double _maxTime=225.; //seems to be mandatory since we need static function
double _alpha=10.;
double _trise=50.;
double _noise=8.;

double pulseShape_fcn(double t, double tmax, double amp, double amp0 = 0., double tau = 19., int n_ord = 3){

    const double ampl_norm = 1.85; // amplitude normalization factor for tau = 19, n = 3

    if( t>tmax-_trise )
	return (amp*ampl_norm * (1 - ((t-(tmax-_trise))/tau)/(n_ord+1)) * std::pow((t-(tmax-_trise))/tau, n_ord) * std::exp(-(t-(tmax-_trise))/tau)) + amp0;
    else return 0;

}

double pulseShape_chi2(const double *x)
{
  double sum = 0.0;
  for(size_t i=0; i<numTS; i++){
    //if( _energy[i]<0 || _time[i]>_maxTime ) continue;
    if( _energy[i] < -150 || _time[i]>_maxTime ) continue;
    double zero = _energy[i]-pulseShape_fcn( _time[i],
	       x[0],x[1] );
    sum += zero * zero / _noise / _noise;
  }
  return sum;
}

/*
/// Old pulse shape
double pulseShape_fcn(double t, double tmax, double amp)
{
  if( t>tmax-_trise ) return amp*std::pow( (t-(tmax-_trise))/_trise,_alpha )*std::exp(-_alpha*(t-tmax)/_trise);
  else return 0;
}
double pulseShape_chi2(const double *x)
{
  double sum = 0.0;
  for(size_t i=0; i<numTS; i++){
    if( _energy[i]<0 || _time[i]>_maxTime ) continue;
    double zero = _energy[i]-pulseShape_fcn( _time[i],
	       x[0],x[1] );
    sum += zero * zero / _noise / _noise;
  }
  return sum;
}
*/

PulseFitter::PulseFitter( int printLevel, double maxTime , double alpha , double trise ) : m_printLevel(printLevel)
{
  _maxTime=maxTime;
  _alpha=alpha;
  _trise=trise;
}

void PulseFitter::run(std::vector<double> &time, std::vector<double> &energy, PulseFitterResult &fit, double noise)
{
  if( time.size()!=energy.size() ){
    std::cout << "ERROR : we should have the same vector size in PulseFitter::run(std::vector<double> time, std::vector<double> energy, PulseFitterResult fit) -> return without fitting" << std::endl;
    return;
  }
  if( time.size()!=11 ){
    std::cout << "ERROR : we should have less than 13 time sample in PulseFitter::run(std::vector<double> time, std::vector<double> energy, PulseFitterResult fit) -> return without fitting" << std::endl;
    return;
  }
  for( uint16_t i=0; i<numTS; i++ ){
    _time[i] = time[i];
    _energy[i] = energy[i];
  }

  if( noise>0 )
    _noise=noise;

  ROOT::Math::Minimizer* m = ROOT::Math::Factory::CreateMinimizer("Minuit", "Migrad");
  m->SetMaxFunctionCalls(m_fitterParameter.nMaxIterations);
  m->SetMaxIterations(m_fitterParameter.nMaxIterations);
  m->SetTolerance(0.001);
  m->SetPrintLevel(m_printLevel);

  ROOT::Math::Functor f(&pulseShape_chi2, 2);

  m->SetFunction(f);

  m->Clear(); // just a precaution

  m->SetVariable(0, "tmax", m_fitterParameter.tmax0, 0.001);
  m->SetVariableLimits(0,
	   m_fitterParameter.tmaxRangeDown,
	   m_fitterParameter.tmaxRangeUp);
  m->SetVariable(1, "amp", _energy[3], 0.001);
  m->SetVariableLimits(1,0,10000);

  m->Minimize();

  const double *xm = m->X();
  const double *errors = m->Errors();
  fit.tmax=xm[0];
  fit.trise=_trise;
  fit.amplitude=xm[1];
  fit.errortmax=errors[0];
  fit.erroramplitude=errors[1];
  fit.chi2=m->MinValue();
  fit.status=(fabs(xm[0]-m_fitterParameter.tmaxRangeDown)>std::numeric_limits<double>::epsilon()&&
	      fabs(xm[0]-m_fitterParameter.tmaxRangeUp)>std::numeric_limits<double>::epsilon()) ? m->Status() : 6;

  fit.ncalls=m->NCalls();

  delete m;

}


double parabolicFit(std::vector<double> x, std::vector<double> y) {
      if (x.size()!=3) return -1;

      //energy reconstruction from parabolic function a*x^2 + b*x + c
      double _a = ( (y[2]-y[1])/(x[2]-x[1]) - (y[1]-y[0])/(x[1]-x[0]) )/( (pow(x[2],2)-pow(x[1],2))/(x[2]-x[1]) - (pow(x[1],2)-pow(x[0],2))/(x[1]-x[0]) );
      double _b = (y[2]-y[1]+_a*(pow(x[1],2)-pow(x[2],2)))/(x[2]-x[1]);
      //double _c = y[0]-_a*pow(x[0],2)-_b*x[0];

      double max_x = (_a < 0) ? -_b/(2*_a) : 0;   //require maximum <--> a<0, unit is ns
      double f_x = (max_x<=150. && max_x>=25.) ? _a*pow(max_x,2)+_b*max_x : -1000.;

      return f_x;
}
