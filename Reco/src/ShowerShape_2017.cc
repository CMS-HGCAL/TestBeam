/*
  Simplified version of the original ShowerShape class for the 2016 test-beams by Shilpi Jain
  @Author: Thorben Quast
  2nd June 2016
  thorben.quast@cern.ch

  >> specified only for the variables of interest, no common mode subtraction
*/



// system include files
#include <iostream>
#include <fstream>
#include <math.h>
// user include files
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "HGCal/DataFormats/interface/HGCalTBRecHitCollections.h"
#include "HGCal/DataFormats/interface/HGCalTBDetId.h"
#include "HGCal/DataFormats/interface/HGCalTBRecHit.h"
#include "HGCal/Geometry/interface/HGCalTBCellVertices.h"
#include "HGCal/Geometry/interface/HGCalTBTopology.h"
#include "HGCal/Geometry/interface/HGCalTBGeometryParameters.h"
#include "HGCal/CondObjects/interface/HGCalElectronicsMap.h"
#include "HGCal/CondObjects/interface/HGCalCondObjectTextIO.h"
#include "HGCal/DataFormats/interface/HGCalTBElectronicsId.h"
#include "HGCal/Geometry/interface/HGCalTBCellVertices.h"
#include "HGCal/Geometry/interface/HGCalTBCellParameters.h"


using namespace std;

bool reverse_sort (double a, double b) { return (a>b); }


enum ENERGY_TYPE {
  TOT = 0,
  ADC_HIGH,
  ADC_LOW,
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

    double getCellIntensity(int layer, size_t index);

    ShowerShape2017(std::string _mapFile, edm::Handle<HGCalTBRecHitCollection> Rechits_, ENERGY_TYPE _energytype, int _n_layers, int _n_skirocs, double ADCtoMIP_[], double maxX_[], double maxY_[]);
    ~ShowerShape2017();

  private:

    // ----------member data ---------------------------
    struct {
      HGCalElectronicsMap emap_;
    } essource_;


    edm::Handle<HGCalTBRecHitCollection> Rechits;

    std::vector<std::pair<double, double>> CellXY;
    std::pair<double, double> CellCentreXY;
    HGCalTBCellVertices TheCell;

    double maxdist = (sqrt(3) / 2.) * HGCAL_TB_CELL::FULL_CELL_SIDE; 
    double maxdist_secNB = (0.5 + sqrt(3)) * HGCAL_TB_CELL::FULL_CELL_SIDE;
    double maxdist_thirdNB = (0.5 + sqrt(3) * 2.) * HGCAL_TB_CELL::FULL_CELL_SIDE;
    double maxdist_fourthNB = (0.5 + sqrt(3) * 3.) * HGCAL_TB_CELL::FULL_CELL_SIDE;

    ENERGY_TYPE energytype;

    int n_layers;
    int n_skirocs;

    double* ADCtoMIP;

    double* maxX;
    double* maxY;

    double* e1;
    double* e7;
    double* e19;
    double* eAll;

    std::map<int, std::vector<double> > cell_intensities;

  };



ShowerShape2017::ShowerShape2017(std::string _mapFile, edm::Handle<HGCalTBRecHitCollection> Rechits_, ENERGY_TYPE _energytype, int _n_layers, int _n_skirocs, double ADCtoMIP_[], double maxX_[], double maxY_[]) {
  HGCalCondObjectTextIO io(0);
  edm::FileInPath fip(_mapFile);
  if (!io.load(fip.fullPath(), essource_.emap_)) {
    throw cms::Exception("Unable to load electronics map");
  };
  
  Rechits =  Rechits_;
  energytype = _energytype;

  n_layers = _n_layers;
  n_skirocs = _n_skirocs;

  ADCtoMIP = new double[n_skirocs];

  maxX = new double[n_layers];
  maxY = new double[n_layers];

  e1 = new double[n_layers];
  e7 = new double[n_layers];
  e19 = new double[n_layers];
  eAll = new double[n_layers];

  for(int i=0; i<n_skirocs; ++i)
    ADCtoMIP[i] = ADCtoMIP_[i];
  for(int i=0; i<n_layers; ++i){
    maxX[i] = maxX_[i];
    maxY[i] = maxY_[i];
  }

  init();

}//constructor ends here


ShowerShape2017::~ShowerShape2017() {
  delete ADCtoMIP;

  delete maxX;
  delete maxY;

  delete e1;
  delete e7;
  delete e19;
  delete eAll;
}


void ShowerShape2017::init(){
  double ncell1[n_layers], ncell7[n_layers], ncell19[n_layers];
  
  for(int iL=0; iL<n_layers; ++iL){
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

    if(eCellType != 0 && eCellType != 4) continue;

    double geometricScaleFactor = 1.;
    if (eCellType == 4) //outer calibration pad cell
      geometricScaleFactor = 9./8.;

    double energyScaled;
    if (energytype == TOT) {
      energyScaled = geometricScaleFactor * (Rechit.energy()) / ADCtoMIP[nSkiroc];    
    } else if (energytype == ADC_HIGH) {
      energyScaled = geometricScaleFactor * (Rechit.energyHigh()) / ADCtoMIP[nSkiroc];    
    } else if (energytype == ADC_LOW) {
      energyScaled = geometricScaleFactor * (Rechit.energyLow()) / ADCtoMIP[nSkiroc];    
    } else {
      throw cms::Exception("Not-implemented energy type.");
    }
    
    radius = sqrt( pow(CellCentreXY.first - maxX[eLayer], 2) + pow(CellCentreXY.second - maxY[eLayer], 2) );

    if((radius < maxdist) && (energyScaled > 0)){
      e1[eLayer] += energyScaled;
      if(eCellType != 1 && eCellType != 2) ncell1[eLayer] += 1;
    }    

 
    if((radius < maxdist_secNB) && (energyScaled > 0)){
      e7[eLayer] += energyScaled;
      if(eCellType != 1 && eCellType != 2) ncell7[eLayer] += 1.;
    }    

    if((radius < maxdist_thirdNB) && (energyScaled > 0)){
      e19[eLayer] += energyScaled;
      if(eCellType != 1 && eCellType != 2) ncell19[eLayer] += 1.;
    }    
  
    if(energyScaled > 0) eAll[eLayer] += energyScaled;
  
    cell_intensities[eLayer].push_back(energyScaled);
  }

  //sort the cell intensities
  for(int iL=0; iL<n_layers; ++iL) std::sort(cell_intensities[iL].begin(), cell_intensities[iL].end()+cell_intensities[iL].size(), reverse_sort);

}

double ShowerShape2017::getCellIntensity(int layer, size_t index) {
  return cell_intensities[layer][index];
}
