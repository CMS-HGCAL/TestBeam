#include "HGCal/RawToDigi/plugins/HGCalTBGenSimSource_Only.h"

using namespace std;

HGCalTBGenSimSource_Only::HGCalTBGenSimSource_Only(const edm::ParameterSet & pset, edm::InputSourceDescription const& desc) :  edm::ProducerSourceFromFiles(pset, desc, true),
	currentRun(-1),
	currentEvent(-1),
	rootFile(NULL)
{
	eventCounter = 0;

	//find and fill the configured runs
	outputCollectionName = pset.getParameter<std::string>("OutputCollectionName");
	fileNames_ = pset.getUntrackedParameter<std::vector<std::string>>("fileNames");
	produces <HGCalTBRecHitCollection>(outputCollectionName);

	_e_mapFile = pset.getParameter<std::string>("e_mapFile_CERN");	
	HGCalCondObjectTextIO io(0);
	edm::FileInPath fip(_e_mapFile);
 	
  if (!io.load(fip.fullPath(), essource_.emap_)) {
	  throw cms::Exception("Unable to load electronics map");
	};
	
	geomc = new HexGeometry(false);

	tree = 0;
  	simHitCellIdE = 0;
  	simHitCellEnE = 0;
  	beamX	 = 0;
  	beamY 	 = 0;

}

bool HGCalTBGenSimSource_Only::setRunAndEventInfo(edm::EventID& id, edm::TimeValue_t& time, edm::EventAuxiliary::ExperimentType& evType)
{	

	if (currentRun == -1) {		//initial loading of a file
		currentRun = 1;
		currentEvent = 0;

		if (!(fileNames_.size())) return false; // need a file...
		std::string fileName = fileNames_[0].c_str(); // currently works with only one file
		if (fileName.find("file:") == 0) fileName = fileName.substr(5);

  		rootFile = new TFile((fileName).c_str());	
  		dir  = (TDirectory*)rootFile->FindObjectAny("HGCalTBAnalyzer");
  		if (dir != NULL) {
  			tree = (TTree*)dir->Get("HGCTB");
  		} else {
  			tree = (TTree*)rootFile->Get("HGCTB");
  		}

   		tree->SetBranchAddress("simHitCellIdE", &simHitCellIdE, &b_simHitCellIdE);
  		tree->SetBranchAddress("simHitCellEnE", &simHitCellEnE, &b_simHitCellEnE);
  		tree->SetBranchAddress("xBeam", &beamX, &b_beamX);
  		tree->SetBranchAddress("yBeam", &beamY, &b_beamY);
  		tree->SetBranchAddress("pBeam", &beamP, &b_beamP);
	}

/*
	if (currentEvent == tree->GetEntries()) {
		currentRun = -1;
		setRunAndEventInfo(id, time, evType);
	}
*/

	if(currentEvent <= tree->GetEntries()){	
		tree->GetEntry(currentEvent);
		currentEvent++;
		return true;
	}
	else return false;

}


void HGCalTBGenSimSource_Only::produce(edm::Event & event)
{	

	eventCounter++;	//indexes each event chronologically passing this plugin 
	//first: fill the rechits
	std::auto_ptr<HGCalTBRecHitCollection> rechits(new HGCalTBRecHitCollection);

	for(unsigned int icell=0; icell<simHitCellIdE->size(); icell++){
		int layer = ((simHitCellIdE->at(icell)>>19)&0x7F);
		
		int cellno = (simHitCellIdE->at(icell)>>0)&0xFF;
		std::pair<double,double> xy = geomc->position(cellno);
		double x = xy.first / 10.;		//values are converted from mm to cm
		double y =  xy.second / 10.;	//values are converted from mm to cm
		int cellType = geomc->cellType(cellno);
		std::pair<int, int> iuiv = TheCell.GetCellIUIVCoordinates(x, y);

		HGCalTBRecHit recHit(HGCalTBDetId(layer, 0, 0, iuiv.first, iuiv.second, cellType), 0., 0., 0., 0., 0); 
			
		recHit.setCellCenterCoordinate(x, y);
	
		uint32_t EID = essource_.emap_.detId2eid(recHit.id());
		HGCalTBElectronicsId eid(EID);	 

		double energy = simHitCellEnE->at(icell);

	 	recHit.setEnergy(energy);
	 	recHit._energyLow = energy;
	 	recHit._energyHigh = energy;
	 	recHit._energyTot = energy;

		rechits->push_back(recHit);
	}	
	event.put(rechits, outputCollectionName);
	
}

void HGCalTBGenSimSource_Only::endJob() {
}


#include "FWCore/Framework/interface/InputSourceMacros.h"

DEFINE_FWK_INPUT_SOURCE(HGCalTBGenSimSource_Only);
