#include <stdio.h>
#include <string.h>
#include <sstream>
#include <TFile.h>
#include <TTree.h>
#include "HGCal/CondObjects/interface/HGCalCondObjectROOTIO.h"
#include "HGCal/DataFormats/interface/HGCalTBElectronicsId.h"
#include "HGCal/Geometry/interface/HGCalTBGeometryParameters.h"

bool HGCalCondObjectROOTIO::loadADCConversion(const std::string& filename, HGCalTBADCConversionsMap &adcConvMap)
{
  TFile* f=new TFile(filename.c_str());
  if( !f->IsOpen() ){
    fprintf(stderr, "Unable to open '%s'\n", filename.c_str());
    return false;
  }
  int module,chip,channel;
  double hgcoeff,lgcoeff;
  int hgsat,lgsat;
  double totcoeff,totped,totthr,totnorm,totpower;
  TTree* tree=(TTree*)f->Get("calib");
  tree->SetBranchAddress("module",&module);
  tree->SetBranchAddress("chip",&chip);
  tree->SetBranchAddress("channel",&channel);
  tree->SetBranchAddress("hgcoeff",&hgcoeff);
  tree->SetBranchAddress("lgcoeff",&lgcoeff);
  tree->SetBranchAddress("hgsat",&hgsat);
  tree->SetBranchAddress("lgsat",&lgsat);
  tree->SetBranchAddress("totcoeff",&totcoeff);
  tree->SetBranchAddress("totthr",&totthr);
  tree->SetBranchAddress("totped",&totped);
  tree->SetBranchAddress("totpower",&totpower);
  tree->SetBranchAddress("totnorm",&totnorm);
  
  for( int ientry=0; ientry<tree->GetEntries(); ientry++ ){
    tree->GetEntry(ientry);
    HGCalTBDetId detid;
    int layerId=m_layout.getLayerWithModuleIndex(module).layerID();
    int skiId=HGCAL_TB_GEOMETRY::N_SKIROC_PER_HEXA*layerId+(HGCAL_TB_GEOMETRY::N_SKIROC_PER_HEXA-chip)%HGCAL_TB_GEOMETRY::N_SKIROC_PER_HEXA+1; //skiId from 1 to nlayers*4
    HGCalTBElectronicsId eid(skiId,channel);      
    if(m_emap.existsEId(eid.rawId()))
      detid = m_emap.eid2detId(eid);
    else{
      std::cout << "module layerId chip skiId channel " << module << " " << layerId << " " << chip << " " << skiId << " " << channel << " not found in emap" << std::endl;
      continue;
    }
    ADCConversions adcConv(detid,1/hgcoeff,1/lgcoeff,totcoeff,totped,totthr,totnorm,totpower);
    adcConvMap.addEntry( adcConv );
  }
  return true;
}
