#include "HGCal/Reco/interface/MCPwaveform.h"


MCPwaveform::MCPwaveform(){
  verbosity_ = 0;
  nSamples_ = 0;
  polarity_ = "neg";
  sigWinMin_ = 0;
  sigWinMax_ = 1000;
  peak_ = -1; 
  amp_ = -1;
  base_ = -1;
  noise_ = -1;
  charge5nsS_ = -1;
  charge5nsB_ = -1;
  fQuality_ = 0;
  Aparam_ = 0.;
  Bparam_ = 0.;
  fitAmp_ = -1;
  fitTimeMax_ = -1;
  fitTimeCF_ = -1;
  nFitSamples_ = 0;
  minForAmpFit_ = 0;
  maxForAmpFit_ = 0;
}



MCPwaveform::MCPwaveform(int verbosity, int nsamples, short* pulse, std::string polarity):
  verbosity_(verbosity), 
  nSamples_(nsamples),   polarity_(polarity), 
  sigWinMin_(0), sigWinMax_(1000), 
  peak_(-1), amp_(-1), base_(-1), noise_(-1), 
  charge5nsS_(-1), charge5nsB_(-1),
  fQuality_(0), Aparam_(-1), Bparam_(-1),
  fitAmp_(-1), fitTimeMax_(-1), fitTimeCF_(-1),
  nFitSamples_(5), minForAmpFit_(5), maxForAmpFit_(5)
{
  for(int ij=0; ij<nSamples_; ++ij) wf_.push_back(pulse[ij]);

  if(nSamples_ < sigWinMax_) {
    std::cout << " MCP WF non valid " << std::endl;
  }

  findAbsolutePeak();
}


MCPwaveform::MCPwaveform(int nsamples, short* pulse, std::string polarity):
  nSamples_(nsamples),   polarity_(polarity), 
  sigWinMin_(0), sigWinMax_(1000), 
  peak_(-1), amp_(-1), base_(-1), noise_(-1), 
  charge5nsS_(-1), charge5nsB_(-1),
  fQuality_(0), Aparam_(-1), Bparam_(-1),
  fitAmp_(-1), fitTimeMax_(-1), fitTimeCF_(-1),
  nFitSamples_(5), minForAmpFit_(5), maxForAmpFit_(5)
{
  verbosity_ = 1;

  for(int ij=0; ij<nSamples_; ++ij) wf_.push_back(pulse[ij]);

  if(nSamples_ < sigWinMax_) {
    std::cout << " MCP WF non valid " << std::endl;
  }

  findAbsolutePeak();
}


void MCPwaveform::analysePeak(){
  if(verbosity_)
    std::cout << " in analysePeak " << std::endl;

  //compute integral of baseline and signal in a 5ns window (25samples * 0.2)
  //noise_ = getIntegral(30, 55);
  getModIntegral(25, 50);
  getSignalIntegral(13, 12);

  //make criteria as loose as possible to get an estimate of the timeCF
  //if( charge5nsS_ > 3. * charge5nsB_ && amp_ > 200.) fQuality_ = 1;
  if(amp_ > 3. * noise_ && amp_ > 100.) fQuality_ = 1;

  if(verbosity_ == 1)
    std::cout << " charge5nsS_ =  " << charge5nsS_ << " charge5nsB_ = " << charge5nsB_  << std::endl;
}


//Get the max/min amplitude wrt polarity
//maybe add criteria to identify double genuine signal (multiple particles)
int MCPwaveform::findAbsolutePeak(std::string polarity){
  if(verbosity_ == 1)
    std::cout << " in findAbsolutePeak polarity_ = " << polarity_ << std::endl;

  //loop over all samples looking for min value
  //taken by Florian and Thorben version: exclude some first samples
  int startSample = 5;
  amp_ = wf_[startSample-1];

  for(int i = startSample; i < nSamples_; ++i){
    if((polarity_ == "neg" && wf_[i] < amp_) ||
       (polarity_ == "pos" && wf_[i] > amp_)){
      amp_ = wf_[i];
      peak_ = i;
    }
    else if(polarity_ != "pos" && polarity_ != "neg"){
      if(wf_[i] > amp_){
	amp_ = wf_[i];
	peak_ = i;
      }
    }
  }

  if(verbosity_ == 1)
    std::cout << " amp_ = " << amp_ << " peak_ = " << peak_ << std::endl;

  return peak_;
}


// calculate mean and variance of all samples outside of peak region
// keep original call with (90, -1) => consider range [10, peak_-90]
float MCPwaveform::getBaseline(int nbinsExcludedLeftOfPeak, int nbinsExcludedRightOfPeak) {
  if(verbosity_ == 1) std::cout << " in getBaseline fQuality_ = " << fQuality_ << std::endl;

  float sum = 0;
  float cnt = 0;
  float var = 0;
  
  if(peak_ == -1) findAbsolutePeak(polarity_);

  int startBin = sigWinMin_;
  int lastBin = sigWinMax_;

  //treat "excluded" as good range interval
  if(nbinsExcludedRightOfPeak != -1 && nbinsExcludedLeftOfPeak != -1){
    startBin = nbinsExcludedLeftOfPeak;
    lastBin = nbinsExcludedRightOfPeak;
  }
  else if(nbinsExcludedRightOfPeak != -1 && peak_ + nbinsExcludedRightOfPeak < 1000){
    startBin = peak_ + nbinsExcludedRightOfPeak;
    lastBin = 1000;
  }
  else if(nbinsExcludedLeftOfPeak != -1 && peak_ - nbinsExcludedLeftOfPeak > 0){
    startBin = 10;
    lastBin = peak_ - nbinsExcludedLeftOfPeak;
  }
  else {
    startBin = nbinsExcludedLeftOfPeak;
    lastBin = nbinsExcludedRightOfPeak;
  }
  
  for(int i = startBin; i < lastBin; ++i){
    sum += wf_[i];
    var += pow(wf_[i], 2);
    cnt += 1;
  }
  
  // break if no cnts
  if (cnt == 0) {
    base_ = -1;
    noise_ = -1;
    return -1;
  }
  
  // calculate variance
  base_ = sum / cnt;
  var = (cnt * var - pow(sum, 2)) / pow(cnt, 2);
  noise_ = sqrt(var);

  return base_;
}


void MCPwaveform::substractBaseline(int nSamples, short* pulse){
  float sum = 0.;
  float cnt = 0.;

  for(int i = 5; i < 50; i++) {
    sum += pulse[i];
    cnt += 1;
  }
  
  float baseline = sum/cnt;
  if(verbosity_ == 1)
    std::cout << " >>> baseline = " << baseline << std::endl;

  for(int i = 0; i < nSamples; i++) {
    pulse[i] = pulse[i] - baseline;
  }
  return;
}


// Get the waveform integral in the given range
float MCPwaveform::getIntegral(int min, int max){
  if(verbosity_ == 1) std::cout << " in getIntegral fQuality_ = " << fQuality_ << std::endl;

  float integral = 0.;
  for(int iSample=min; iSample<max; ++iSample){
    if(iSample < 0) continue;
    if(iSample >= nSamples_) break;
    integral += wf_.at(iSample);
  }

  charge5nsB_ = integral;
  return integral;
}


//Get the signal integral around the the max
float MCPwaveform::getSignalIntegral(int riseWin, int fallWin){
  if(verbosity_ == 1) std::cout << " in getSignalIntegral fQuality_ = " << fQuality_ << std::endl;

  //compute position of the max
  if(peak_ == -1) findAbsolutePeak();
  if(peak_-riseWin < 0) {
    riseWin = 0;
    fallWin = riseWin+25;
  }

  //compute integral
  float integral = 0.;
  for(int iSample=peak_-riseWin; iSample<peak_+fallWin; ++iSample){
    //if signal window goes out of bound return a bad value
    if(iSample < 0) continue;
    if(iSample >= nSamples_) break;
    integral += wf_.at(iSample);
  }

  charge5nsS_ = integral;
  return integral;
}


//Get the integral of Abs(WF) over the given range
float MCPwaveform::getModIntegral(int min, int max){   
  float integral = 0.;
  for(int iSample=min; iSample<max; ++iSample){
    if(iSample < 0) continue;
    if(iSample >= nSamples_) break;

    if(wf_.at(iSample) < 0) integral -= wf_.at(iSample);
    else integral += wf_.at(iSample);
  }

  charge5nsB_ = integral;
  return integral;
}


// Get CF time for a given fraction and in a given range
// note rise time is about 1ns = 5samples => nFitSamples to adapt
float MCPwaveform::getTimeCF(float frac, int nFitSamples, int min, int max){
  if(verbosity_ == 1)   
    std::cout << " in getTimeCF fQuality_ = " << fQuality_ << " fitAmp_ = " << fitAmp_ << " frac = " << frac <<  std::endl;

  if(!fQuality_) return -1;

  if(frac == 1) return fitTimeMax_;
  if(Aparam_ != -1 && Bparam_ != -1) return (fitAmp_ * frac - Aparam_) / Bparam_;

  if(nFitSamples != 0) nFitSamples_ = nFitSamples;
  if(min != 0 && max != 0){
    minForAmpFit_ = min;
    maxForAmpFit_ = max;
  }

  //good with min, max == 3
  if(fitAmp_ == -1) getInterpolatedAmpMax(minForAmpFit_, maxForAmpFit_);

  int cfSample = (sigWinMin_ == -1) ? 0 : sigWinMin_;        
  //find first sample above Amax*frac
  for(int iSample=peak_; iSample>min; --iSample){
    if(wf_.at(iSample) < fitAmp_*frac){
      cfSample = iSample;
      break;
    }
  }
  
  if(verbosity_ == 1)
    std::cout << " in getTimeCF cfSample = " << cfSample << " peak_ = " << peak_ << std::endl;
  
  //rescale nFitSamples in case bigger than useful interval
  //maybe restrict to 5 samples anyway
  if(std::abs(cfSample - peak_) < (nFitSamples_-1)/2){
    nFitSamples_ = 2*std::abs(cfSample - peak_)+1;
  }

  //interpolate with A+Bx = frac * amp
  float chi2cf = linearInterpolation(cfSample - (nFitSamples_-1)/2, cfSample + (nFitSamples_-1)/2);
  fitTimeCF_ = (fitAmp_ * frac - Aparam_) / Bparam_;

  if(verbosity_ == 1)
    std::cout << " fit time = " << fitTimeCF_ << " peak_ = " << peak_  << std::endl;

  return fitTimeCF_;
}


//Linear interpolation with  A+Bx
float MCPwaveform::linearInterpolation(const int& min, const int& max, const int& sampleToSkip){
  if(verbosity_ == 1) 
  std::cout << " in linearInterpolation " << std::endl;

  // assume y = A + B x
  // B = (NSxy Sy - Sx Sy)/(NSxx - SxSx)
  // A = Sy/N - B (Sx/N) = (Sxx Sy - Sx Sxy)/(NSxx - SxSx)
  float xx = 0.;
  float xy = 0.;
  float Sx = 0.;
  float Sy = 0.;
  float Sxx = 0.;
  float Sxy = 0.;

  int usedSamples=0;
  for(int iSample=min; iSample<=max; ++iSample){

    if(iSample<0 || iSample >= int(wf_.size())) continue;

    xx = iSample * iSample;
    xy = iSample * wf_[iSample];
    Sx += iSample;
    Sy += wf_[iSample];
    Sxx += xx;
    Sxy += xy;
    ++usedSamples;
  }
    
  float Delta = usedSamples * Sxx - Sx * Sx;
  Aparam_ = (Sxx * Sy - Sx * Sxy) / Delta;
  Bparam_ = (usedSamples * Sxy - Sx * Sy) / Delta;

  //---compute chi2---
  float chi2 = 0.;
  float sigma2 = pow(noise_, 2);
  for(int iSample=min; iSample<=max; ++iSample){
    if(iSample<0 || iSample >= int(wf_.size())) continue;

    chi2 += pow(wf_[iSample] - Aparam_ - Bparam_ * iSample, 2)/sigma2;
  } 

  return chi2/(usedSamples-2);
}


//Linear interpolation with A+Bx+Cxx
float MCPwaveform::pol2Interpolation(const int& min, const int& max, const int& nFitSamples){
  if(verbosity_ == 1) 
    std::cout << " in pol2Interpolation fQuality_ = " << fQuality_  
	      << " min = " << min << "max = " << max << std::endl;

  TH1F h_max("h_max", "", nFitSamples, min, max+1);
  TF1 f_max("f_max", "pol2", min, max+1);

  int bin = 1;
  for(int iSample=min; iSample<=max; ++iSample){
    h_max.SetBinContent(bin, wf_[iSample]);
    h_max.SetBinError(bin, noise_);
    ++bin;
    if(verbosity_ == 1)
      std::cout << "tg->SetPoint(" << bin-1 << ", " << iSample << ", " << wf_[iSample] << ");" << std::endl;
  }

  auto fit_result = h_max.Fit(&f_max, "QRS");
  fitTimeMax_ = -1. * f_max.GetParameter(1) / (2.*f_max.GetParameter(2));
  fitAmp_ = f_max.Eval(fitTimeMax_);

  float chi2 = f_max.GetChisquare()/(nFitSamples-3);

  return chi2;
}


//Get the interpolated max/min amplitude 
float MCPwaveform::getInterpolatedAmpMax(int min, int max){
  if(verbosity_ == 1) std::cout << " in getInterpolatedAmpMax fQuality_ = " << fQuality_ << std::endl;

  if(!fQuality_) return -1;

  //get maximum sample
  if(peak_ == -1) 
    findAbsolutePeak();

  if(peak_ - min < sigWinMin_) min = minForAmpFit_;
  if(peak_ + max > sigWinMax_) max = maxForAmpFit_;

  pol2Interpolation(peak_- min, peak_+ max+1, max+min+1);

  return fitAmp_;
}
