#include "H4treeRecoWC.h"

#include <iostream>

#include "TChain.h"
#include "TString.h"

//#define DEBUG_VERBOSE

//
H4treeRecoWC::H4treeRecoWC(TChain *tree, TString outUrl):
  H4treeWC(tree), MaxTdcChannels_(16), MaxTdcReadings_(20)
{
  fOut_  = TFile::Open(outUrl,"RECREATE");
  recoT_ = new TTree("H4treeRecoWC","H4treeRecoWC");
  recoT_->SetDirectory(fOut_);

  //event header
  recoT_->Branch("runNumber",    &runNumber,    "runNumber/i");
  recoT_->Branch("spillNumber",  &spillNumber,  "spillNumber/i");
  recoT_->Branch("evtNumber",    &evtNumber,    "evtNumber/i");
  recoT_->Branch("evtTimeDist",  &evtTimeDist,  "evtTimeDist/i");
  recoT_->Branch("evtTimeStart", &evtTimeStart, "evtTimeStart/i");
  recoT_->Branch("nEvtTimes",    &nEvtTimes,    "nEvtTimes/i");

  //TDC 
  tdc_readings_.resize(MaxTdcChannels_);
  recoT_->Branch("nwc",       &nwc_,      "nwc/i");
  recoT_->Branch("wc_recox",   wc_recox_, "wc_recox[nwc]/F");
  recoT_->Branch("wc_recoy",   wc_recoy_, "wc_recoy[nwc]/F");
  recoT_->Branch("wc_xl_hits",   wc_xl_hits_, "wc_xl_hits[nwc]/i");
  recoT_->Branch("wc_xr_hits",   wc_xr_hits_, "wc_xr_hits[nwc]/i");
  recoT_->Branch("wc_yu_hits",   wc_yu_hits_, "wc_yu_hits[nwc]/i");
  recoT_->Branch("wc_yd_hits",   wc_yd_hits_, "wc_yd_hits[nwc]/i");

  InitDigi();
}

void H4treeRecoWC::FillEvent()
{
  recoT_->Fill();
}

//
void H4treeRecoWC::InitDigi()
{
  //WCD map                                                                                                                                                                 
  std::map<std::string, int> WCD_channel;
  WCD_channel["l"] = 0;
  WCD_channel["r"] = 1;
  WCD_channel["d"] = 2;
  WCD_channel["u"] = 3;

  //WCE map                                                                                                                                                                 
  std::map<std::string, int> WCE_channel;
  WCE_channel["l"] = 4;
  WCE_channel["r"] = 5;
  WCE_channel["d"] = 6;
  WCE_channel["u"] = 7;

  std::vector< std::map<std::string, int>  > WC;
  // 0 == WCD                                                                                                                                                               
  WC.push_back(WCD_channel);
  WC.push_back(WCE_channel);


  nwc_ = WC.size();
  for(size_t i=0; i<WC.size(); ++i){
    wcXl_[i] = WC.at(i)["l"];
    wcXr_[i] = WC.at(i)["r"];
    wcYd_[i] = WC.at(i)["d"];
    wcYu_[i] = WC.at(i)["u"];
  }
}



void H4treeRecoWC::FillTDC()
{
  //reset data  
  for (uint j=0; j<MaxTdcChannels_; j++){ tdc_readings_[j].clear();}

  //fill with new data 
  for (uint i=0; i<nTdcChannels; i++)
    {
      if (tdcChannel[i]<MaxTdcChannels_)
        {
          tdc_readings_[tdcChannel[i]].push_back((float)tdcData[i]);
        }
    }
  
  //compute average positions
  for(UInt_t iwc=0; iwc<nwc_; iwc++)
    {
      //default values in case no tdc data available  
      wc_recox_[iwc]=-999;
      wc_recoy_[iwc]=-999;
      
      wc_xl_hits_[iwc]=tdc_readings_[wcXl_[iwc]].size();
      wc_xr_hits_[iwc]=tdc_readings_[wcXr_[iwc]].size();
      wc_yu_hits_[iwc]=tdc_readings_[wcYu_[iwc]].size();
      wc_yd_hits_[iwc]=tdc_readings_[wcYd_[iwc]].size();

      if (tdc_readings_[wcXl_[iwc]].size()!=0 && tdc_readings_[wcXr_[iwc]].size()!=0){
	float TXl = *std::min_element(tdc_readings_[wcXl_[iwc]].begin(),tdc_readings_[wcXl_[iwc]].begin()+tdc_readings_[wcXl_[iwc]].size());
	float TXr = *std::min_element(tdc_readings_[wcXr_[iwc]].begin(),tdc_readings_[wcXr_[iwc]].begin()+tdc_readings_[wcXr_[iwc]].size());
	wc_recox_[iwc] = (TXr-TXl)*0.005; // = /40./5./10. //position in cm 0.2mm/ns with 25ps LSB TDC                                                 
      }
      if (tdc_readings_[wcYd_[iwc]].size()!=0 && tdc_readings_[wcYu_[iwc]].size()!=0){
	float TYd = *std::min_element(tdc_readings_[wcYd_[iwc]].begin(),tdc_readings_[wcYd_[iwc]].begin()+tdc_readings_[wcYd_[iwc]].size());
	float TYu = *std::min_element(tdc_readings_[wcYu_[iwc]].begin(),tdc_readings_[wcYu_[iwc]].begin()+tdc_readings_[wcYu_[iwc]].size());
	wc_recoy_[iwc] = (TYu-TYd)*0.005; // = /40./5./10. //position in cm 0.2mm/ns with 25ps LSB TDC 
      }
    }
}



H4treeRecoWC::~H4treeRecoWC()
{
  fOut_->cd();
  std::cout << " recoT_->GetEntries() = " << recoT_->GetEntries() << std::endl;
  recoT_->Write();
  fOut_->Close();
}
