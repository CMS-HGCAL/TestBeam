/********************************************************************
* File: waveform.h
* ------------------------
*
* Description:
* Tools to analyse raw waveforms.
*
*
* Version:
* Original Author: Florian Pitters,
* Extended by Thorben Quast at the indicated positions, 25 October 2018
*
*******************************************************************/


#ifndef __WAVEFORM_H_INCLUDED__
#define __WAVEFORM_H_INCLUDED__

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
#include "TGraph.h"
#include "TGraphErrors.h"


using namespace std;




// --------------------------------------------------------
// Defintions
// --------------------------------------------------------

struct peakValues{
	Int_t id;
	Int_t peak;
	Float_t amp;
	Float_t charge;
	Float_t base;
	Float_t noise;
	Float_t amppeak;
	Float_t tpeak;
	Float_t tlinear15;
	Float_t tlinear30;
	Float_t tlinear45;
	Float_t tlinear60;

	UInt_t fQuality;
	UInt_t fSaturated;
};


float getBaseline(int peak, short *pulse, int nbinsExcludedLeftOfPeak, int nbinsExcludedRightOfPeak);
int substractBaseline(int nsamples, short* pulse, int baseline);		//added by T.Q.
int findAbsolutePeak(int nsamples, short *pulse, std::string polarity);
int findRealPeak(int nsamples, short *pulse, int rank, std::string option);
int findFirstMinAboveNoise(int nsamples, short *pulse, short noise);
float getPulseIntegral(int peak, short *pulse, int roi_left, int roi_right, std::string option);
float getNoise(int peak, short *pulse, int nbinsExcludedLeftOfPeak, int nbinsExcludedRightOfPeak);
int qualityPulse(int nsamples, short *pulse);

peakValues *analysePeak(int nsamples, short *pulse);


#endif
