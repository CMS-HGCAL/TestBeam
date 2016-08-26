#include "HGCal/Reco/interface/RecHitCommonMode.h"
using namespace std;
//
// class declaration
//

RecHitCommonMode::RecHitCommonMode(edm::Handle<HGCalTBRecHitCollection> Rechits)
{
	using namespace edm;
        Rechits_ = Rechits;
       
        for(int iLayer(0);iLayer<MAXLAYERS;iLayer++){
          FullCell_CMNoise[iLayer] = 0;
          HalfCell_CMNoise[iLayer] = 0;
          CalibPad_CMNoise[iLayer] = 0;
          MB_CMNoise[iLayer] = 0;
          MergedCell_CMNoise[iLayer] = 0;
        }

}


RecHitCommonMode::~RecHitCommonMode()
{
  for(int iLayer(0);iLayer<MAXLAYERS;iLayer++){
    delete Full_Cell[iLayer];
  }
}

void
RecHitCommonMode::evaluate()
{

        for(int iLayer(0);iLayer<MAXLAYERS;iLayer++){ 
          sprintf(name, "FullCell_Layer_%i",iLayer);
          sprintf(title, "FullCell Layer %i",iLayer);
          Full_Cell[iLayer] = new TH1F(name, title, 400,-200., 200.);
        }

        double Average_Pedestal_Per_Event_Half[MAXLAYERS] = {0};
        int Cell_counter_Half[MAXLAYERS] = {0};
        double Average_Pedestal_Per_Event_MB[MAXLAYERS] = {0};
        int Cell_counter_MB[MAXLAYERS] = {0};
        double Average_Pedestal_Per_Event_Calib_Pad[MAXLAYERS] = {0};
        int Cell_counter_Calib_Pad[MAXLAYERS] = {0};
        double Average_Pedestal_Per_Event_Merged_Cell[MAXLAYERS] = {0};
        int Cell_counter_Merged_Cell[MAXLAYERS] = {0};

        for(HGCalTBRecHitCollection::const_iterator rechit = Rechits_->begin(); rechit != Rechits_->end(); rechit++) {
          CellCentreXY = TheCell.GetCellCentreCoordinatesForPlots((rechit->id()).layer(), (rechit->id()).sensorIU(), (rechit->id()).sensorIV(), (rechit->id()).iu(), (rechit->id()).iv(), sensorsize);
          HGCalTBRecHit hit = (*rechit);
          if(hit.energyHigh() > 100)continue;
          CellType type = (CellType)((rechit->id()).cellType());
          if(type == FullCell){     
            Full_Cell[(rechit->id()).layer() - 1]->Fill(hit.energyHigh()); 
          }
	  else if(type == HalfCell){
	    Cell_counter_Half[(rechit->id()).layer() - 1]+= 1;
	    Average_Pedestal_Per_Event_Half[(rechit->id()).layer() - 1] += hit.energyHigh(); 
	  }
	  else if(type == CalibPad){
	    Cell_counter_Calib_Pad[(rechit->id()).layer() - 1]++;
	    Average_Pedestal_Per_Event_Calib_Pad[(rechit->id()).layer() - 1] += hit.energyHigh(); 
	  } 
	  else if(type == MBandMerged && ( ((rechit->id()).iu() == 4 && (rechit->id()).iv() == 3) || ((rechit->id()).iu() == -7 && (rechit->id()).iv() == 4) || ((rechit->id()).iu() == 7 && (rechit->id()).iv() == -3) || ((rechit->id()).iu() == -4 && (rechit->id()).iv() == -3) ) ){
	    Cell_counter_MB[(rechit->id()).layer() - 1]+=1;
	    Average_Pedestal_Per_Event_MB[(rechit->id()).layer() - 1] += hit.energyHigh(); 
	  }
	  else if(type == MBandMerged && ( ((rechit->id()).iu() == -4 && (rechit->id()).iv() == 6) || ((rechit->id()).iu() == -2 && (rechit->id()).iv() == 6) || ((rechit->id()).iu() == 4 && (rechit->id()).iv() == -7) || ((rechit->id()).iu() == 2 && (rechit->id()).iv() == -6) ) ){
	    Cell_counter_Merged_Cell[(rechit->id()).layer() - 1]+=1;
	    Average_Pedestal_Per_Event_Merged_Cell[(rechit->id()).layer() - 1] += hit.energyHigh(); 
	  }
      }

        for(int iLayer(0);iLayer<MAXLAYERS;iLayer++){
          if(Full_Cell[iLayer]->GetEntries() > 0){
            Full_Cell[iLayer]->Fit("gaus", "Q");
            FullCell_CMNoise[iLayer] = Full_Cell[iLayer]->GetFunction("gaus")->GetParameter(1);
          } else {
            FullCell_CMNoise[iLayer] = 0;
          }

          if(Cell_counter_Calib_Pad[iLayer] > 0)CalibPad_CMNoise[iLayer] = Average_Pedestal_Per_Event_Calib_Pad[iLayer]/Cell_counter_Calib_Pad[iLayer];
          if(Cell_counter_Half[iLayer] > 0)HalfCell_CMNoise[iLayer] = Average_Pedestal_Per_Event_Half[iLayer]/Cell_counter_Half[iLayer];
          if(Cell_counter_MB[iLayer] > 0)MB_CMNoise[iLayer] = Average_Pedestal_Per_Event_MB[iLayer]/Cell_counter_MB[iLayer];
          if(Cell_counter_Merged_Cell[iLayer] > 0)MergedCell_CMNoise[iLayer] = Average_Pedestal_Per_Event_Merged_Cell[iLayer]/Cell_counter_Merged_Cell[iLayer];
        }
}

float 
RecHitCommonMode::getCommonModeNoise(int layer, CellType type, int iu, int iv)
{

      float CMNoise(0);
      switch(type){
        case FullCell: CMNoise = FullCell_CMNoise[layer]; break;
        case CalibPad: CMNoise = CalibPad_CMNoise[layer]; break;
        case HalfCell: CMNoise = HalfCell_CMNoise[layer]; break;
        case MBandMerged: if( (iu == 4 && iv == 3) || ( iu == -7 && iv == 4) || ( iu == 7 && iv == -3) || ( iu == -4 && iv == -3) )CMNoise = MB_CMNoise[layer];
                          else if ( (iu == -4 && iv == 6) || ( iu == -2 && iv == 6) || ( iu == 4 && iv == -7) || ( iu == 2 && iv == -6) )CMNoise = MergedCell_CMNoise[layer];

                          break;
        case OuterCalib: CMNoise = FullCell_CMNoise[layer]; break;
        default: CMNoise = 0; break;
      };

   return CMNoise;
}

      
float 
RecHitCommonMode::getCommonModeNoise(int layer, CellType type, std::string const& subtype/*=""*/)
{
      if( (subtype != "") && (subtype != "MB") && (subtype != "MergedCell")){
        throw cms::Exception("InvalidCellType") << "Invalid cell subtype: " <<  subtype;
      } 

      float CMNoise(0);
      switch(type){
        case FullCell: CMNoise = FullCell_CMNoise[layer]; break;
        case CalibPad: CMNoise = CalibPad_CMNoise[layer]; break;
        case HalfCell: CMNoise = HalfCell_CMNoise[layer]; break;
        case MBandMerged: if(subtype == "MB")CMNoise = MB_CMNoise[layer];
                          else if(subtype == "MergedCell")CMNoise = MergedCell_CMNoise[layer];
                          break;
        case OuterCalib: CMNoise = FullCell_CMNoise[layer]; break;
        default: CMNoise = 0; break;
      };

   return CMNoise;
}








 
