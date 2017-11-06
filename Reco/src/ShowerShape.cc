/*
  originally copied from Layer_Sum_Analyzer and then modified accordingly
  @Author: Shilpi Jain
  28th August 2016
  dhilpi.jain@cern.ch

  >> updated for nLayers and simplifying the structure (amartell)
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


class ShowerShape 
{

public:
  void init();

  float e1By7(int layer) {return e1over7[layer]; };
  float e1By19(int layer) {return e1over19[layer]; };
  //  float eMax(int layer) {return  &maxE, float &maxX, float &maxY);
  void getAllEnergy(int layer, float& en1, float& en7, float& en19, float& enAll) 
  {en1 = e1[layer]; en7 = e7[layer]; en19 = e19[layer]; enAll = eAll[layer]; return;}
  void getPos(int layer, float& xtmp, float& ytmp)
  {xtmp = xEnWeig19[layer]; ytmp = yEnWeig19[layer]; return;};

  void logWeightedSpread(int layer, float &pox, float &poy) {pox = spreadx43[layer]; poy = spready43[layer]; return;};
  void logWeightedSqSpread19(int layer, float &spreadx, float &spready, float &spreadrad) 
  {spreadx = spreadx19[layer]; spready = spready19[layer]; spreadrad = spreadrad19[layer]; return;};
  void logWeightedPosition(int layer, float &posx, float &posy)
  {posx = posx43[layer]; posy = posy43[layer]; return;};
  void logWeightedPosition19(int layer, float &posx, float &posy)
  {posx = posx19[layer]; posy = posy19[layer]; return;};

  void F7F19(int layer, float &F7, float &F19) {F7 = log(e6[layer]/e7[layer]); F19 = log(e7[layer]/e19[layer]); return;};

  void eccentricity(int layer, float &ratio) {ratio = eccentr[layer]; return;};
  
  ShowerShape(std::string _mapFile, edm::Handle<HGCalTBRecHitCollection> Rechits_, double ADCtoMIP_[MAXSKIROCS], double commonmode_[MAXLAYERS], int CMTHRESHOLD, double maxE[MAXLAYERS], double maxX[MAXLAYERS], double maxY[MAXLAYERS]);
  ~ShowerShape();
  
  
private:

  // ----------member data ---------------------------
  struct {
    HGCalElectronicsMap emap_;
  } essource_;


  edm::Handle<HGCalTBRecHitCollection> Rechits;
  
  int sensorsize = 128;
  std::vector<std::pair<double, double>> CellXY;
  std::pair<double, double> CellCentreXY;
  HGCalTBCellVertices TheCell;
  //float maxdist = (1 + sqrt (3) / 2) * HGCAL_TB_CELL::FULL_CELL_SIDE;
  //  float maxdist = (0.5 + sqrt (3)) * HGCAL_TB_CELL::FULL_CELL_SIDE;
  double maxdist = (sqrt(3) / 2.) * HGCAL_TB_CELL::FULL_CELL_SIDE;  // <<< FIXME maxdist > HGCAL_TB_CELL::FULL_CELL_SIDE !!                                           
  double maxdist_secNB = (0.5 + sqrt(3)) * HGCAL_TB_CELL::FULL_CELL_SIDE;
  double maxdist_thirdNB = (0.5 + sqrt(3) * 2.) * HGCAL_TB_CELL::FULL_CELL_SIDE;
  double maxdist_fourthNB = (0.5 + sqrt(3) * 3.) * HGCAL_TB_CELL::FULL_CELL_SIDE;

  
  double ADCtoMIP[MAXSKIROCS];
  
  double commonmode[MAXLAYERS];
  double maxE[MAXLAYERS], maxX[MAXLAYERS], maxY[MAXLAYERS];
  int CMTHRESHOLD;

  float e1over7[MAXLAYERS];
  float e1over19[MAXLAYERS];

  float e1[MAXLAYERS];
  float e6[MAXLAYERS];
  float e7[MAXLAYERS];
  float e12[MAXLAYERS];
  float e19[MAXLAYERS];
  float e43[MAXLAYERS];
  float eAll[MAXLAYERS];

  //energy weighted position
  float xEnWeig19[MAXLAYERS];
  float yEnWeig19[MAXLAYERS];

  //transvers profile log weighted spread and position
  float spreadx43[MAXLAYERS];
  float spready43[MAXLAYERS];
  float posx43[MAXLAYERS];
  float posy43[MAXLAYERS];

  ///in 19 crystals log weighted spread and position
  float posx19[MAXLAYERS]; 
  float posy19[MAXLAYERS];
  float spreadx19[MAXLAYERS]; 
  float spready19[MAXLAYERS];
  float spreadrad19[MAXLAYERS];

  //eccentricity
  float eccentr[MAXLAYERS];


};



ShowerShape::ShowerShape(std::string _mapFile, edm::Handle<HGCalTBRecHitCollection> Rechits_, double ADCtoMIP_[MAXSKIROCS], double commonmode_[MAXLAYERS], int CMTHRESHOLD_, 
			 double maxE_[MAXLAYERS], double maxX_[MAXLAYERS], double maxY_[MAXLAYERS])
{
  
  Rechits =  Rechits_;
  CMTHRESHOLD = CMTHRESHOLD_;
  for(int i=0; i<MAXSKIROCS; ++i)
    ADCtoMIP[i] = ADCtoMIP_[i];
  for(int i=0; i<MAXLAYERS; ++i){
    commonmode[i] = commonmode_[i];
    maxE[i] = maxE_[i];
    maxX[i] = maxX_[i];
    maxY[i] = maxY_[i];
  }

  HGCalCondObjectTextIO io(0);
  edm::FileInPath fip(_mapFile);
  if (!io.load(fip.fullPath(), essource_.emap_)) {
    throw cms::Exception("Unable to load electronics map");
  };

  init();

}//constructor ends here


ShowerShape::~ShowerShape()
{

	// do anything here that needs to be done at desctruction time
	// (e.g. close files, deallocate resources etc.)

}


void ShowerShape::init(){
  float ncell1[MAXLAYERS], ncell6[MAXLAYERS], ncell12[MAXLAYERS];
  float ncell7[MAXLAYERS], ncell19[MAXLAYERS], ncell43[MAXLAYERS];

  //  float ncell1S[MAXLAYERS], ncell6S[MAXLAYERS], ncell12S[MAXLAYERS];
  //float ncell7S[MAXLAYERS], 
  float ncell19S[MAXLAYERS], ncell43S[MAXLAYERS];

  float deno43[MAXLAYERS], lposx43[MAXLAYERS], lposy43[MAXLAYERS], Dposx43[MAXLAYERS], Dposy43[MAXLAYERS], eccx43[MAXLAYERS], eccy43[MAXLAYERS], eccxy43[MAXLAYERS];
  float deno19[MAXLAYERS], lposx19[MAXLAYERS], lposy19[MAXLAYERS], Dposx19[MAXLAYERS], Dposy19[MAXLAYERS];

  float radius = 0.;



  for(int iL=0; iL<MAXLAYERS; ++iL){
    e1[iL] = e6[iL] = e7[iL] = e12[iL] = e19[iL] = e43[iL] = eAll[iL] = eccentr[iL] = 0.;
    ncell1[iL] = ncell6[iL] = ncell7[iL] = ncell12[iL] = ncell19[iL] = 0;
    //ncell1S[iL] = ncell6S[iL] = ncell7S[iL] = 
    ncell43S[iL] = ncell19S[iL] = 0;
    deno43[iL] = lposx43[iL] = lposy43[iL] = Dposx43[iL] = Dposy43[iL] = eccx43[iL] = eccy43[iL] = eccxy43[iL] = 0.;
    deno19[iL] = lposx19[iL] = lposy19[iL] = Dposx19[iL] = Dposy19[iL] =  xEnWeig19[iL] =  yEnWeig19[iL] = 0.;
  } 

  //  std::cout << "loop recHits size = " << Rechits->size() << std::endl;
  for(auto Rechit : *Rechits){
    // looping over each rechit to fill histogram
    CellCentreXY = TheCell.GetCellCentreCoordinatesForPlots((Rechit.id()).layer(), (Rechit.id()).sensorIU(), (Rechit.id()).sensorIV(), (Rechit.id()).iu(), (Rechit.id()).iv(), sensorsize);
    
    int eLayer = (Rechit.id()).layer()-1;
    int eCellType = (Rechit.id()).cellType();
    uint32_t EID = essource_.emap_.detId2eid(Rechit.id());
    HGCalTBElectronicsId eid(EID);
    int eSkiroc = (eid.iskiroc() - 1);

    if(eCellType != 0 && eCellType != 1 && eCellType != 4) continue;

    float energyCMsub = (Rechit.energy() - commonmode[eLayer]) / ADCtoMIP[eSkiroc];
    if(eCellType == 1) energyCMsub = (Rechit.energy()) / ADCtoMIP[eSkiroc];

    radius = sqrt( pow(CellCentreXY.first - maxX[eLayer], 2) + pow(CellCentreXY.second - maxY[eLayer], 2) );

    if((radius < maxdist) && (energyCMsub > CMTHRESHOLD)){
      //      if(ncell1[eLayer] == 1) continue;
      //      std::cout << " x cell = " << CellCentreXY.first << " y cell = " << CellCentreXY.second << " layer = " << eLayer << std::endl;
      //      std::cout << " x MAX = " << maxX[eLayer] << " y MAX = " << maxY[eLayer] << std::endl;

      //      std::cout << " id = " << eCellType << " ncell1[eLayer] = " << ncell1[eLayer] << std::endl;
      e1[eLayer] += energyCMsub;
      if(eCellType == 1 || eCellType == 2) ncell1[eLayer] += 0.;
      else ncell1[eLayer] += 1;
    }    
    if((radius > maxdist && radius < maxdist_secNB) && (energyCMsub > CMTHRESHOLD)){
      //      if(ncell6[eLayer] == 6) continue;
      e6[eLayer] += energyCMsub;
      if(eCellType == 1 || eCellType == 2) ncell6[eLayer] += 0.;
      else ncell6[eLayer] += 1;
    }    
    if((radius < maxdist_secNB) && (energyCMsub > CMTHRESHOLD)){
      //      if(ncell7[eLayer] == 7) continue;
      e7[eLayer] += energyCMsub;
      if(eCellType == 1 || eCellType == 2) ncell7[eLayer] += 0.;
      else ncell7[eLayer] += 1;
    }    
    if((radius > maxdist_secNB && radius < maxdist_thirdNB) && (energyCMsub > CMTHRESHOLD)){
      //      if(ncell12[eLayer] == 12) continue;
      e12[eLayer] += energyCMsub;
      if(eCellType == 1 || eCellType == 2) ncell12[eLayer] += 0.;
      else ncell12[eLayer] += 1;
    }    
    if((radius < maxdist_thirdNB) && (energyCMsub > CMTHRESHOLD)){
      //      if(ncell19[eLayer] == 19) continue;
      e19[eLayer] += energyCMsub;
      if(eCellType == 1 || eCellType == 2) ncell19[eLayer] += 0.;
      else ncell19[eLayer] += 1;
    }    
    if((radius < maxdist_fourthNB) && (energyCMsub > CMTHRESHOLD)){
      //      if(ncell43[eLayer] == 43) continue;
      e43[eLayer] += energyCMsub;
      if(eCellType == 1 || eCellType == 2) ncell43[eLayer] += 0.;
      else ncell43[eLayer] += 1;
    }    
    if(energyCMsub > CMTHRESHOLD) eAll[eLayer] += energyCMsub;
  }//recHits


  double w0 = 4.7;
  for(auto Rechit1 : *Rechits){
    //getting X and Y coordinates
    CellCentreXY = TheCell.GetCellCentreCoordinatesForPlots((Rechit1.id()).layer(), (Rechit1.id()).sensorIU(), (Rechit1.id()).sensorIV(), (Rechit1.id()).iu(), (Rechit1.id()).iv(), sensorsize);
    
    int eLayer = (Rechit1.id()).layer()-1;
    int eCellType = (Rechit1.id()).cellType();
    uint32_t EID = essource_.emap_.detId2eid(Rechit1.id());
    HGCalTBElectronicsId eid(EID);
    int eSkiroc = (eid.iskiroc() - 1);

    if(eCellType != 0 && eCellType != 1 && eCellType != 4) continue;

    float energyCMsub = (Rechit1.energy() - commonmode[eLayer]) / ADCtoMIP[eSkiroc];
    if(eCellType == 1) energyCMsub = (Rechit1.energy()) / ADCtoMIP[eSkiroc];

    radius = sqrt( pow(CellCentreXY.first - maxX[eLayer], 2) + pow(CellCentreXY.second - maxY[eLayer], 2) );


    double x = CellCentreXY.first;
    double y = CellCentreXY.second;
    if((radius < maxdist_fourthNB) && energyCMsub > CMTHRESHOLD){
      //      if(ncell43S[eLayer] == 43) continue;
      if(eCellType == 1 || eCellType == 2) ncell43S[eLayer] += 0.;
      else ncell43S[eLayer] += 1;
      lposx43[eLayer] += (x) * TMath::Max(0., (log(energyCMsub/e43[eLayer]) + w0) ); 
      lposy43[eLayer] += (y) * TMath::Max(0., (log(energyCMsub/e43[eLayer]) + w0) );       
      Dposx43[eLayer] += (x-maxX[eLayer]) * TMath::Max(0., (log(energyCMsub/e43[eLayer]) + w0) ); 
      Dposy43[eLayer] += (y-maxY[eLayer]) * TMath::Max(0., (log(energyCMsub/e43[eLayer]) + w0) );       
      
      eccx43[eLayer] += pow( (x-posx43[eLayer]),2 )*TMath::Max(0., (log(energyCMsub/e43[eLayer]) + w0) ); 
      eccy43[eLayer] += pow( (y-posy43[eLayer]),2 )*TMath::Max(0., (log(energyCMsub/e43[eLayer]) + w0) ); 
      eccxy43[eLayer] += (x-posx43[eLayer])*(y-posy43[eLayer])*TMath::Max(0., (log(energyCMsub/e43[eLayer]) + w0) ); 
      
      deno43[eLayer] += TMath::Max(0., (log(energyCMsub/e43[eLayer]) + w0) );
    } // radius 4th ring
    
    if((radius < maxdist_thirdNB) && energyCMsub > CMTHRESHOLD){
      //      if(ncell19S[eLayer] == 19) continue;
      if(eCellType == 1 || eCellType == 2) ncell19S[eLayer] += 0.;
      else ncell19S[eLayer] += 1;

      xEnWeig19[eLayer] += CellCentreXY.first * energyCMsub;
      yEnWeig19[eLayer] += CellCentreXY.second * energyCMsub;

      lposx19[eLayer] += (x)*TMath::Max(0., (log(energyCMsub/e19[eLayer]) + w0) ); 
      lposy19[eLayer] += (y)*TMath::Max(0., (log(energyCMsub/e19[eLayer]) + w0) );       
      Dposx19[eLayer] += (x-maxX[eLayer])*TMath::Max(0., (log(energyCMsub/e19[eLayer]) + w0) ); 
      Dposy19[eLayer] += (y-maxY[eLayer])*TMath::Max(0., (log(energyCMsub/e19[eLayer]) + w0) );       
      deno19[eLayer] += TMath::Max(0., (log(energyCMsub/e19[eLayer]) + w0) );
    } // radius 3rd ring

  }//for(auto Rechit1 : *Rechits)


  for(int iL=0; iL<MAXLAYERS; ++iL){
    if(ncell1[iL] > 1 ) std::cout << " Problem expected 1 cells found => " << ncell1[iL] << std::endl;
    if(ncell6[iL] > 6) std::cout << " Problem expected 6 cells found => " << ncell6[iL] << std::endl; 
    if(ncell7[iL] > 7) std::cout << " Problem expected 7 cells found => " << ncell7[iL] << std::endl;
    if(ncell12[iL] > 12) std::cout << " Problem expected 12 cells found => " << ncell12[iL] << std::endl;
    if(ncell19[iL] > 19) std::cout << " Problem expected 19 cells found => " << ncell19[iL] << " in layer = " << iL << std::endl;
    
  
    if(ncell7[iL] > 0) e1over7[iL] = maxE[iL]/e7[iL];
    else e1over7[iL] = -99.;
    if(ncell19[iL] > 0) e1over19[iL] = maxE[iL]/e19[iL];
    else e1over19[iL] = -99.;

    posx19[iL] = lposx19[iL]/deno19[iL];
    posy19[iL] = lposy19[iL]/deno19[iL];
    posx43[iL] = lposx43[iL]/deno43[iL];
    posy43[iL] = lposy43[iL]/deno43[iL];
    
    spreadx43[iL] = Dposx43[iL]/deno43[iL];
    spready43[iL] = Dposy43[iL]/deno43[iL];

    spreadx19[iL] = Dposx19[iL]/deno19[iL];
    spready19[iL] = Dposy19[iL]/deno19[iL];
    spreadrad19[iL] = sqrt(spreadx19[iL] * spreadx19[iL] + spready19[iL] * spready19[iL]);


    //eccentricy
    const int N = 2;
    double spread[N*N] = {
      eccx43[iL], eccxy43[iL],
      eccxy43[iL], eccy43[iL] };

    TMatrixDSym m(N, spread);
    TMatrixDSymEigen me(m);
    
    TVectorD eigenval = me.GetEigenValues();
    TMatrixD eigenvec = me.GetEigenVectors();
    
    float e1 = eigenval[0];
    float e2 = eigenval[1];
  
    eccentr[iL] = e1/e2;
    //eigenval.Print();  
  }
}
