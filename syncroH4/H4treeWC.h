#ifndef H4treeWC_h
#define H4treeWC_h

#include <TROOT.h>
#include <TChain.h>
#include <TString.h>
#include <TFile.h>

// Header file for the classes stored in the TTree if any.
// Fixed size dimensions of array or collections stored in the TTree if any.

class H4treeWC {
 public :
  TChain         *fChain;   //!pointer to the analyzed TTree or TChain
  Int_t           fCurrent; //!current Tree number in a TChain

  // Declaration of leaf types
  UInt_t          runNumber;
  UInt_t          spillNumber;
  UInt_t          evtNumber;
  UInt_t          evtTimeDist;
  UInt_t          evtTimeStart;
  UInt_t          nEvtTimes;
  ULong64_t       evtTime[1];   //[nEvtTimes]
  UInt_t          evtTimeBoard[1];   //[nEvtTimes]

  UInt_t          nTdcChannels;
  UInt_t          tdcBoard[22];   //[nTdcChannels]
  UInt_t          tdcChannel[22];   //[nTdcChannels]
  UInt_t          tdcData[22];   //[nTdcChannels]

  // List of branches
  TBranch        *b_runNumber;   //!
  TBranch        *b_spillNumber;   //!
  TBranch        *b_evtNumber;   //!
  TBranch        *b_evtTimeDist;   //!
  TBranch        *b_evtTimeStart;   //!
  TBranch        *b_nEvtTimes;   //!
  TBranch        *b_evtTime;   //!
  TBranch        *b_evtTimeBoard;   //!

  TBranch        *b_nTdcChannels;   //!
  TBranch        *b_tdcBoard;   //!
  TBranch        *b_tdcChannel;   //!
  TBranch        *b_tdcData;   //!

  H4treeWC(TChain *chain=0);
  virtual ~H4treeWC();

  virtual Int_t    GetEntry(Long64_t entry);
  virtual Int_t    TotEntries();
  virtual Long64_t LoadTree(Long64_t entry);
  virtual void     Init(TChain *tree);
  /* virtual void     Loop(); */
  /* virtual Bool_t   Notify(); */
  /* virtual void     Show(Long64_t entry = -1); */
};

#endif
