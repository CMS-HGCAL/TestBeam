#ifndef MCPwaveform_h
#define MCPwaveform_h

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <iostream>
#include <sstream>
#include <iomanip>

#include "TTree.h"
#include "TBranch.h"
#include "TROOT.h"
#include "TStyle.h"
#include "TFile.h"
#include "TTree.h"
#include "TBranch.h"
#include "TCanvas.h"
#include "TMath.h"
#include "TH1.h"
#include "TH2.h"
#include "TH3.h"
#include "TF1.h"
#include "TLinearFitter.h"
#include "TGraph.h"
#include "TGraphErrors.h"



class MCPwaveform{
 public:

  MCPwaveform();
  MCPwaveform(int nsamples, short* pulse, std::string);
  MCPwaveform(int verbosity, int nsamples, short* pulse, std::string);

  void analysePeak();

  int findAbsolutePeak(std::string polarity="pos");

  void substractBaseline(int nSamples, short* pulse);
  float getBaseline(int nbinsExcludedLeftOfPeak, int nbinsExcludedRightOfPeak);

  float getIntegral(int min, int max);
  float getSignalIntegral(int riseWin, int fallWin);
  float getModIntegral(int min, int max);

  float getTimeCF(float frac, int nFitSamples=0, int min=0, int max=0);
  float linearInterpolation(const int& min, const int& max, const int& sampleToSkip=0);
  float pol2Interpolation(const int& min, const int& max, const int& sampleToSkip=0);
  float getInterpolatedAmpMax(int min, int max);


  inline float getNoise(){ return noise_; };
  inline int getPeak(){ return peak_; };
  inline float getAmp(){ return amp_; };
  inline float getFitAmp(){ return fitAmp_; };
  inline float getFitPeak(){ return fitTimeMax_; };
  inline int getQuality(){ return fQuality_; };
  inline float getCharge5nsS(){ return charge5nsS_; };
  inline float getCharge5nsB(){ return charge5nsB_; };



 private:

  int verbosity_;
  std::vector<short> wf_;
  int nSamples_;
  std::string polarity_;
  int sigWinMin_;
  int sigWinMax_;

  Int_t peak_;
  Float_t amp_;
  Float_t base_;
  Float_t noise_;

  Float_t charge5nsS_;
  Float_t charge5nsB_;

  UInt_t fQuality_;
  Float_t Aparam_;
  Float_t Bparam_;
  Float_t fitAmp_;
  Float_t fitTimeMax_;
  Float_t fitTimeCF_;

  int nFitSamples_;
  int minForAmpFit_;
  int maxForAmpFit_;
};

#endif
