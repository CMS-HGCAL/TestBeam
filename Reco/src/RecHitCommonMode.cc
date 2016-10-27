#include "HGCal/Reco/interface/RecHitCommonMode.h"
using namespace std;
//
// class declaration
//

RecHitCommonMode::RecHitCommonMode(edm::Handle<HGCalTBRecHitCollection> Rechits, HGCalElectronicsMap& emap)
{
	using namespace edm;
        Rechits_ = Rechits;
        essource_.emap_ = emap;
        if(essource_.emap_.size() == 0)
          throw cms::Exception("InvalidElectronicsMap") << "HGCalElectronicsMap has zero size"; 

        for(int iSKI(0);iSKI<MAXSKIROCS;iSKI++){
          FullCell_CMNoise_Fit[iSKI] = 0;
          FullCell_CMNoise_Mean[iSKI] = 0;
          CalibPad_CMNoise[iSKI] = 0;
          HalfCell_CMNoise[iSKI] = 0;
          MB_CMNoise[iSKI] = 0;
          OuterCalibPad_CMNoise[iSKI] = 0;
          MergedCell_CMNoise[iSKI] = 0;
        }

}


RecHitCommonMode::~RecHitCommonMode()
{
  for(int iSKI(0);iSKI<MAXSKIROCS;iSKI++){
    delete Full_Cell[iSKI];
  }
}

void
RecHitCommonMode::evaluate(float maxEcut/*=0*/)
{

        for(int iSKI(0);iSKI<MAXSKIROCS;iSKI++){ 
          sprintf(name, "FullCell_Layer_%i",iSKI);
          sprintf(title, "FullCell Layer %i",iSKI);
          Full_Cell[iSKI] = new TH1F(name, title, 400,-200., 200.);
        }

        double Average_Pedestal_Per_Event_Full[MAXSKIROCS] = {0};
        int Cell_counter_Full[MAXSKIROCS] = {0};
        double Average_Pedestal_Per_Event_Calib_Pad[MAXSKIROCS] = {0};
        int Cell_counter_Calib_Pad[MAXSKIROCS] = {0};
        double Average_Pedestal_Per_Event_Half[MAXSKIROCS] = {0};
        int Cell_counter_Half[MAXSKIROCS] = {0};
        double Average_Pedestal_Per_Event_MB[MAXSKIROCS] = {0};
        int Cell_counter_MB[MAXSKIROCS] = {0};
        double Average_Pedestal_Per_Event_Merged_Cell[MAXSKIROCS] = {0};
        int Cell_counter_Merged_Cell[MAXSKIROCS] = {0};

        for(HGCalTBRecHitCollection::const_iterator rechit = Rechits_->begin(); rechit != Rechits_->end(); rechit++) {
          HGCalTBRecHit hit = (*rechit);
          if(hit.energyHigh() > maxEcut)continue;
          uint32_t EID = essource_.emap_.detId2eid(rechit->id());
          HGCalTBElectronicsId eid(EID);
          int idSKIROC = eid.iskiroc()-1;
	  int type = (rechit->id()).cellType();
          switch(type){
            case 0:     
              Full_Cell[idSKIROC]->Fill(hit.energyHigh());
              Cell_counter_Full[idSKIROC]+=1; 
              Average_Pedestal_Per_Event_Full[idSKIROC] += hit.energyHigh();
              break; 
            case 1:
	      Cell_counter_Calib_Pad[idSKIROC]+=1;
	      Average_Pedestal_Per_Event_Calib_Pad[idSKIROC] += 0;
              break; 
            case 2:
	      Cell_counter_Half[idSKIROC]+= 1;
	      Average_Pedestal_Per_Event_Half[idSKIROC] += hit.energyHigh();
              break; 
	    case 3:
	      Cell_counter_MB[idSKIROC]+=1;
	      Average_Pedestal_Per_Event_MB[idSKIROC] += hit.energyHigh();
              break; 
            case 4:
              Full_Cell[idSKIROC]->Fill(hit.energyHigh()); // Full cell + outer calib pad  for Gaussian fit
              Cell_counter_Full[idSKIROC]+=1; 
              Average_Pedestal_Per_Event_Full[idSKIROC] += hit.energyHigh();
              break;
            case 5:
              Cell_counter_Merged_Cell[idSKIROC]+=1;
              Average_Pedestal_Per_Event_Merged_Cell[idSKIROC] += hit.energyHigh();    
              break;
            default:
              throw cms::Exception("InvalidCellType") << "rechit celltype if out of range";
              break;
          }               
	}
	
        for(int iSKI(0);iSKI<MAXSKIROCS;iSKI++){
          if(Full_Cell[iSKI]->GetEntries() > 0){
            int fitstatus = Full_Cell[iSKI]->Fit("gaus", "Q");
            if(fitstatus ==0)FullCell_CMNoise_Fit[iSKI] = Full_Cell[iSKI]->GetFunction("gaus")->GetParameter(1);
            else FullCell_CMNoise_Fit[iSKI] = 0;
          } else {
            FullCell_CMNoise_Fit[iSKI] = 0;
          }
          
          if(Cell_counter_Full[iSKI] > 0)FullCell_CMNoise_Mean[iSKI] = Average_Pedestal_Per_Event_Full[iSKI]/Cell_counter_Full[iSKI];
          CalibPad_CMNoise[iSKI] = 0.0; 
          if(Cell_counter_Half[iSKI] > 0)HalfCell_CMNoise[iSKI] = Average_Pedestal_Per_Event_Half[iSKI]/Cell_counter_Half[iSKI];
          if(Cell_counter_MB[iSKI] > 0)MB_CMNoise[iSKI] = Average_Pedestal_Per_Event_MB[iSKI]/Cell_counter_MB[iSKI];
          OuterCalibPad_CMNoise[iSKI] = FullCell_CMNoise_Mean[iSKI];
          if(Cell_counter_Merged_Cell[iSKI] > 0)MergedCell_CMNoise[iSKI] = Average_Pedestal_Per_Event_Merged_Cell[iSKI]/Cell_counter_Merged_Cell[iSKI];
        }
}

float 
RecHitCommonMode::getGaussCommonModeNoise(HGCalTBDetId id)
{
      uint32_t EID = essource_.emap_.detId2eid(id);
      HGCalTBElectronicsId eid(EID);
      int idSKIROC = eid.iskiroc()-1;
      int type = id.cellType();
      float CMNoise(0);
      switch(type){
        case 0: CMNoise = FullCell_CMNoise_Fit[idSKIROC]; break;
        case 1: CMNoise = CalibPad_CMNoise[idSKIROC]; break;
        case 2: CMNoise = HalfCell_CMNoise[idSKIROC]; break;
        case 3: CMNoise = MB_CMNoise[idSKIROC]; break;
        case 4: CMNoise = FullCell_CMNoise_Fit[idSKIROC]; break;
	case 5: CMNoise = MergedCell_CMNoise[idSKIROC]; break;
        default: CMNoise = 0; break;
      };
      
      return CMNoise;
}

float 
RecHitCommonMode::getMeanCommonModeNoise(HGCalTBDetId id)
{
      uint32_t EID = essource_.emap_.detId2eid(id);
      HGCalTBElectronicsId eid(EID);
      int idSKIROC = eid.iskiroc()-1;
      int type = id.cellType();
      float CMNoise(0);
      switch(type){
        case 0: CMNoise = FullCell_CMNoise_Mean[idSKIROC]; break;
        case 1: CMNoise = CalibPad_CMNoise[idSKIROC]; break;
        case 2: CMNoise = HalfCell_CMNoise[idSKIROC]; break;
        case 3: CMNoise = MB_CMNoise[idSKIROC]; break;
        case 4: CMNoise = OuterCalibPad_CMNoise[idSKIROC]; break;
	case 5: CMNoise = MergedCell_CMNoise[idSKIROC]; break;
        default: CMNoise = 0; break;
      };
      
      return CMNoise;
}

float 
RecHitCommonMode::getGaussCommonModeNoise(int idSKIROC, int type)
{
  if( type < 0 || type > 5){ 
    throw cms::Exception("InvalidCellType") << "Invalid cell type: " <<  type;
  } 

  float CMNoise(0);
  switch(type){
    case 0: CMNoise = FullCell_CMNoise_Fit[idSKIROC]; break;
    case 1: CMNoise = CalibPad_CMNoise[idSKIROC]; break;
    case 2: CMNoise = HalfCell_CMNoise[idSKIROC]; break;
    case 3: CMNoise = MB_CMNoise[idSKIROC]; break;
    case 4: CMNoise = FullCell_CMNoise_Fit[idSKIROC]; break;
    case 5: CMNoise = MergedCell_CMNoise[idSKIROC]; break;
    default: CMNoise = 0; break;
  };

  return CMNoise;
}

float
RecHitCommonMode::getMeanCommonModeNoise(int idSKIROC, int type)
{
  if( type < 0 || type > 5){
    throw cms::Exception("InvalidCellType") << "Invalid cell type: " <<  type;
  }
      float CMNoise(0);
      switch(type){
        case 0: CMNoise = FullCell_CMNoise_Mean[idSKIROC]; break;
        case 1: CMNoise = CalibPad_CMNoise[idSKIROC]; break;
        case 2: CMNoise = HalfCell_CMNoise[idSKIROC]; break;
        case 3: CMNoise = MB_CMNoise[idSKIROC]; break;
        case 4: CMNoise = OuterCalibPad_CMNoise[idSKIROC]; break;
        case 5: CMNoise = MergedCell_CMNoise[idSKIROC]; break;
        default: CMNoise = 0; break;
      };

      return CMNoise;
}






 
