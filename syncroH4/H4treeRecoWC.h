#ifndef H4treeRecoWC_h
#define H4treeRecoWC_h

#include "H4treeWC.h"

#include "TFile.h"
#include "TRandom.h"
#include "TString.h"
#include "TChain.h"

#include <set>

class H4treeRecoWC : public H4treeWC
{

 public:

  H4treeRecoWC(TChain *, TString outUrl="H4treeRecoWCOut.root");
  void FillEvent();   // (UInt_t evN, ULong64_t evT); 
  void FillTDC();
  void InitDigi();

  inline float timeSampleUnit(int drs4Freq)
  {
    if (drs4Freq == 0)
      return 0.2E-9;
    else if (drs4Freq == 1)
      return 0.4E-9;
    else if (drs4Freq == 2)
      return 1.E-9;    
    return -999.;
  }


  ~H4treeRecoWC();  

  //TDC readings
  UInt_t MaxTdcChannels_, MaxTdcReadings_;
  std::vector< std::vector<Float_t> > tdc_readings_;
  Float_t wc_recox_[16], wc_recoy_[16];
  UInt_t wc_xl_hits_[16], wc_xr_hits_[16], wc_yu_hits_[16], wc_yd_hits_[16]; 
  UInt_t nwc_,wcXl_[4],wcXr_[4],wcYu_[4],wcYd_[4];


 protected:
  TTree *recoT_;
  TFile *fOut_;

  // private:
};

#endif
