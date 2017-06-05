/*
  Simplified version of the original ShowerShape class for the 2016 test-beams by Shilpi Jain
  @Author: Thorben Quast
  2nd June 2016
  thorben.quast@cern.ch

  >> specified only for the variables of interest, no common mode subtraction
*/



// system include files
#include <memory>
#include <iostream>
#include "TH2Poly.h"
#include "TH1F.h"
#include "TF1.h"
#include <sstream>
#include <fstream>
#include <math.h>
// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/one/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "HGCal/DataFormats/interface/HGCalTBRecHitCollections.h"
#include "HGCal/DataFormats/interface/HGCalTBDetId.h"
#include "HGCal/DataFormats/interface/HGCalTBRecHit.h"
#include "HGCal/Geometry/interface/HGCalTBCellVertices.h"
#include "HGCal/Geometry/interface/HGCalTBTopology.h"
#include "HGCal/Geometry/interface/HGCalTBGeometryParameters.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "HGCal/CondObjects/interface/HGCalElectronicsMap.h"
#include "HGCal/CondObjects/interface/HGCalCondObjectTextIO.h"
#include "HGCal/DataFormats/interface/HGCalTBElectronicsId.h"
#include "HGCal/DataFormats/interface/HGCalTBDataFrameContainers.h"
#include "HGCal/Geometry/interface/HGCalTBCellVertices.h"
#include "HGCal/Geometry/interface/HGCalTBCellParameters.h"
#include "HGCal/Geometry/interface/HGCalTBSpillParameters.h"


#include "TMatrixDSym.h"
#include "TMatrixDSymEigen.h"

using namespace std;

enum ENERGY_TYPE {
  ADC_BASELINE_SUBTRACTED = 0,
  ADC_HIGH,
  ADC_LOW,
  ADC_HIGH_CMSubtracted,
  MIP
};

class ShowerShape2017 {
  public:
    void init();
    double getE1(int layer) { return e1[layer]; }
    double getE7(int layer) { return e7[layer]; }
    double getE19(int layer) { return e19[layer]; }
    double getEAll(int layer) { return eAll[layer]; }

    void getAllEnergy(int layer, double& en1, double& en7, double& en19, double& enAll)  {
      en1 = e1[layer]; 
      en7 = e7[layer]; 
      en19 = e19[layer]; 
      enAll = eAll[layer]; 
      return;
    }
    ShowerShape2017(std::string _mapFile, int _sensorsize, edm::Handle<HGCalTBRecHitCollection> Rechits_, ENERGY_TYPE _energytype, double ADCtoMIP_[MAXSKIROCS], double maxX[MAXLAYERS], double maxY[MAXLAYERS]);
    ~ShowerShape2017();

  private:

    // ----------member data ---------------------------
    struct {
      HGCalElectronicsMap emap_;
    } essource_;


    edm::Handle<HGCalTBRecHitCollection> Rechits;

    int sensorsize;
    std::vector<std::pair<double, double>> CellXY;
    std::pair<double, double> CellCentreXY;
    HGCalTBCellVertices TheCell;

    double maxdist = (sqrt(3) / 2.) * HGCAL_TB_CELL::FULL_CELL_SIDE; 
    double maxdist_secNB = (0.5 + sqrt(3)) * HGCAL_TB_CELL::FULL_CELL_SIDE;
    double maxdist_thirdNB = (0.5 + sqrt(3) * 2.) * HGCAL_TB_CELL::FULL_CELL_SIDE;
    double maxdist_fourthNB = (0.5 + sqrt(3) * 3.) * HGCAL_TB_CELL::FULL_CELL_SIDE;

    ENERGY_TYPE energytype;

    double ADCtoMIP[MAXSKIROCS];

    double maxX[MAXLAYERS], maxY[MAXLAYERS];

    double e1[MAXLAYERS];
    double e7[MAXLAYERS];
    double e19[MAXLAYERS];
    double eAll[MAXLAYERS];


  };



ShowerShape2017::ShowerShape2017(std::string _mapFile, int _sensorsize, edm::Handle<HGCalTBRecHitCollection> Rechits_, ENERGY_TYPE _energytype, double ADCtoMIP_[MAXSKIROCS], double maxX_[MAXLAYERS], double maxY_[MAXLAYERS]) {
  HGCalCondObjectTextIO io(0);
  edm::FileInPath fip(_mapFile);
  if (!io.load(fip.fullPath(), essource_.emap_)) {
    throw cms::Exception("Unable to load electronics map");
  };
  
  sensorsize = _sensorsize;

  Rechits =  Rechits_;

  energytype = _energytype;

  for(int i=0; i<MAXSKIROCS; ++i)
    ADCtoMIP[i] = ADCtoMIP_[i];
  for(int i=0; i<MAXLAYERS; ++i){
    maxX[i] = maxX_[i];
    maxY[i] = maxY_[i];
  }

  init();

}//constructor ends here


ShowerShape2017::~ShowerShape2017() {

}


void ShowerShape2017::init(){
  double ncell1[MAXLAYERS], ncell7[MAXLAYERS], ncell19[MAXLAYERS];
  

  for(int iL=0; iL<MAXLAYERS; ++iL){
    e1[iL] = e7[iL] = e19[iL] = eAll[iL] = 0.;
    ncell1[iL] = ncell7[iL] = ncell19[iL] = 0;
  } 

  double radius = 0.;

  for(auto Rechit : *Rechits){
    // looping over each rechit to fill histogram
    CellCentreXY.first = Rechit.getCellCenterCartesianCoordinate(0);
    CellCentreXY.second = Rechit.getCellCenterCartesianCoordinate(1);


    
    int eLayer = (Rechit.id()).layer()-1;
    int eCellType = (Rechit.id()).cellType();
    uint32_t EID = essource_.emap_.detId2eid(Rechit.id());
    HGCalTBElectronicsId eid(EID);
    int nSkiroc = (eid.iskiroc() - 1);

    if(eCellType != 0 && eCellType != 1 && eCellType != 4) continue;

    double energyCMsub;
    if (energytype == ADC_BASELINE_SUBTRACTED) {
      energyCMsub = (Rechit.energy()) / ADCtoMIP[nSkiroc];    
    } else if (energytype == ADC_HIGH) {
      energyCMsub = (Rechit.energyHigh()) / ADCtoMIP[nSkiroc];    
    } else if (energytype == ADC_LOW) {
      energyCMsub = (Rechit.energyLow()) / ADCtoMIP[nSkiroc];    
    } else {
      throw cms::Exception("Not-implemented energy type.");
    }
    
    radius = sqrt( pow(CellCentreXY.first - maxX[eLayer], 2) + pow(CellCentreXY.second - maxY[eLayer], 2) );

    if((radius < maxdist) && (energyCMsub > 0)){
      e1[eLayer] += energyCMsub;
      if(eCellType != 1 && eCellType != 2) ncell1[eLayer] += 1;
    }    

 
    if((radius < maxdist_secNB) && (energyCMsub > 0)){
      e7[eLayer] += energyCMsub;
      if(eCellType != 1 && eCellType != 2) ncell7[eLayer] += 1.;
    }    

    if((radius < maxdist_thirdNB) && (energyCMsub > 0)){
      e19[eLayer] += energyCMsub;
      if(eCellType != 1 && eCellType != 2) ncell19[eLayer] += 1.;
    }    
  
    if(energyCMsub > 0) eAll[eLayer] += energyCMsub;
  }


}
