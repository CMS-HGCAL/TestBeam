#include "H4treeWC.h"

#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>

// void H4treeWC::Loop()
// {
//   if (fChain == 0) return;

//   Long64_t nentries = fChain->GetEntries();

//   Long64_t nbytes = 0, nb = 0;
//   for (Long64_t jentry=0; jentry<nentries;jentry++) {
//     Long64_t ientry = LoadTree(jentry);
//     if (ientry < 0) break;
//     nb = fChain->GetEntry(jentry);   nbytes += nb;
//     // if (Cut(ientry) < 0) continue;
//   }
// }

H4treeWC::H4treeWC(TChain *tree) : fChain(0) 
{
  // if parameter tree is not specified (or zero), connect the file
  // used to generate this class and read the Tree.
  if (tree == 0) {
    TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject("root://eoscms//eos/cms/store/group/dpg_ecal/alca_ecalcalib/TimingTB_H2_Jul2015/raw/DataTree/3374/14.root");
    if (!f || !f->IsOpen()) {
      f = new TFile("root://eoscms//eos/cms/store/group/dpg_ecal/alca_ecalcalib/TimingTB_H2_Jul2015/raw/DataTree/3374/14.root");
    }
    f->GetObject("H4treeWC",tree);

  }
  Init(tree);
}

H4treeWC::~H4treeWC()
{
  if (!fChain) return;
  delete fChain->GetCurrentFile();
}

Int_t H4treeWC::GetEntry(Long64_t entry)
{
  // Read contents of entry.
  if (!fChain) return 0;
  return fChain->GetEntry(entry);
}

Int_t H4treeWC::TotEntries()
{
  // Read contents of entry.
  if (!fChain) return 0;
  return fChain->GetEntries();
}

Long64_t H4treeWC::LoadTree(Long64_t entry)
{
  // Set the environment to read one entry
  if (!fChain) return -5;
  Long64_t centry = fChain->LoadTree(entry);
  if (centry < 0) return centry;
  if (fChain->GetTreeNumber() != fCurrent) {
    fCurrent = fChain->GetTreeNumber();
    //    Notify();
  }
  return centry;
}

void H4treeWC::Init(TChain *tree)
{
  // The Init() function is called when the selector needs to initialize
  // a new tree or chain. Typically here the branch addresses and branch
  // pointers of the tree will be set.
  // It is normally not necessary to make changes to the generated
  // code, but the routine can be extended by the user if needed.
  // Init() will be called many times when running on PROOF
  // (once per file to be processed).

  // Set branch addresses and branch pointers
  if (!tree) return;
  fChain = tree;
  fCurrent = -1;
  fChain->SetMakeClass(1);

  fChain->SetBranchAddress("runNumber", &runNumber, &b_runNumber);
  fChain->SetBranchAddress("spillNumber", &spillNumber, &b_spillNumber);
  fChain->SetBranchAddress("evtNumber", &evtNumber, &b_evtNumber);
  fChain->SetBranchAddress("evtTimeDist", &evtTimeDist, &b_evtTimeDist);
  fChain->SetBranchAddress("evtTimeStart", &evtTimeStart, &b_evtTimeStart);
  fChain->SetBranchAddress("nEvtTimes", &nEvtTimes, &b_nEvtTimes);
  fChain->SetBranchAddress("evtTime", evtTime, &b_evtTime);
  fChain->SetBranchAddress("evtTimeBoard", evtTimeBoard, &b_evtTimeBoard);
  fChain->SetBranchAddress("nTdcChannels", &nTdcChannels, &b_nTdcChannels);
  fChain->SetBranchAddress("tdcBoard", tdcBoard, &b_tdcBoard);
  fChain->SetBranchAddress("tdcChannel", tdcChannel, &b_tdcChannel);
  fChain->SetBranchAddress("tdcData", tdcData, &b_tdcData);
  //  Notify();
}

// Bool_t H4treeWC::Notify()
// {
//   // The Notify() function is called when a new file is opened. This
//   // can be either for a new TTree in a TChain or when when a new TTree
//   // is started when using PROOF. It is normally not necessary to make changes
//   // to the generated code, but the routine can be extended by the
//   // user if needed. The return value is currently not used.

//   return kTRUE;
// }

// void H4treeWC::Show(Long64_t entry)
// {
//   // Print contents of entry.
//   // If entry is not specified, print current entry
//   if (!fChain) return;
//   fChain->Show(entry);
// }
// Int_t ee::Cut(Long64_t entry)
// {
//   // This function may be called from Loop.
//   // returns  1 if entry is accepted.
//   // returns -1 otherwise.
//   return 1;
// }
