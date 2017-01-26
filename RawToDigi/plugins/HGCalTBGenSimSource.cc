#include "HGCal/RawToDigi/plugins/HGCalTBGenSimSource.h"

using namespace std;

HGCalTBGenSimSource::HGCalTBGenSimSource(const edm::ParameterSet & pset, edm::InputSourceDescription const& desc) :  edm::ProducerSourceFromFiles(pset, desc, true),
	inputPathFormat(""),
	currentRun(-1),
	currentEvent(-1),
	rootFile(NULL)
{

	//find and fill the configured runs
	outputCollectionName = pset.getParameter<std::string>("OutputCollectionName");
	runEnergyMapFile = pset.getUntrackedParameter<std::string>("runEnergyMapFile"); 
	inputPathFormat = pset.getUntrackedParameter<std::string>("inputPathFormat");
	
	produces <HGCalTBRecHitCollection>(outputCollectionName);
	produces<RunData>("RunData");


	if (fileNames()[0] != "file:DUMMY") {
		for (int i = 0; i<(int)(fileNames().size()); i++) {
			FileInfo fInfo;
			fInfo.index = -1;
			fInfo.energy = -1;
			fInfo.runType = "";
			fInfo.config = -1;
			fInfo.name = fileNames()[i];

			_fileNames.push_back(fInfo);
		}
	}
	

	std::fstream map_file;
	map_file.open(runEnergyMapFile.c_str(), std::fstream::in);
	fillConfiguredRuns(map_file);
	map_file.close();

	geomc = new HexGeometry(false);

	tree = 0;
}

void HGCalTBGenSimSource::fillConfiguredRuns(std::fstream& map_file) {
	std::string filePath;

	//perform the loop and fill configuredRuns
	char fragment[100];
	int readCounter = 0;
	int _run = 0, _configuration = 0; double _energy = 0; std::string _runType = ""; 

	while (map_file.is_open() && !map_file.eof()) {
		readCounter++;
		map_file >> fragment;
		if (readCounter <= 4) continue; 	//skip the header
		else if (readCounter % 4 == 1) {
			if (((std::string)fragment).find("//") == std::string::npos)	//skip comments of form //
				_run = atoi(fragment); 
			else
				readCounter = 1;
		}
		else if (readCounter % 4 == 2) _energy = atof(fragment); 
		else if (readCounter % 4 == 3) _runType = (std::string)fragment; 
		else if (readCounter % 4 == 0) {
			_configuration = atoi(fragment); 
			//store
		
			FileInfo fInfo;
			fInfo.index = _run;
			fInfo.energy = _energy;
			fInfo.runType = _runType;
			fInfo.config = _configuration;
			
			
			filePath = inputPathFormat;		
			
			filePath.replace(filePath.find("<RUN>"), 5, std::to_string(_run));
			filePath.replace(filePath.find("<ENERGY>"), 8, std::to_string((int)_energy));

			fInfo.name = filePath;
			
			_fileNames.push_back(fInfo);			
		}
	}

	fileIterator = _fileNames.begin();
}

bool HGCalTBGenSimSource::setRunAndEventInfo(edm::EventID& id, edm::TimeValue_t& time, edm::EventAuxiliary::ExperimentType& evType)
{	

	if (fileIterator ==_fileNames.end()) {
		return false; 		//end of files is reached
	}

	if (currentRun == -1) {		//initial loading of a file
		currentRun = (*fileIterator).index;
		currentEvent = 0;
		/*
		if (rootFile != NULL)
			delete rootFile;
  	*/
  
  	rootFile = new TFile(((*fileIterator).name).c_str());	
  	dir  = (TDirectory*)rootFile->FindObjectAny("HGCalTBAnalyzer");
  	tree = (TTree*)dir->Get("HGCTB");

   	tree->SetBranchAddress("simHitCellIdE", &simHitCellIdE, &b_simHitCellIdE);
  	tree->SetBranchAddress("simHitCellEnE", &simHitCellEnE, &b_simHitCellEnE);
	} 

	if (currentEvent == tree->GetEntries()) {
		fileIterator++;
		currentRun = -1;
		setRunAndEventInfo(id, time, evType);
	}
	
	tree->GetEntry(currentEvent);

	currentEvent++;

	return true;

}


void HGCalTBGenSimSource::produce(edm::Event & event)
{	

	if (fileIterator ==_fileNames.end()) {
		std::cout<<"End of the files in the producer is reached..."<<std::endl;
		return;
	}
	std::auto_ptr<RunData> rd(new RunData);
	rd->energy = (*fileIterator).energy;
	rd->configuration = (*fileIterator).config;
	rd->runType = (*fileIterator).runType;
	rd->run = (*fileIterator).index;

	event.put(std::move(rd), "RunData");	
	std::auto_ptr<HGCalTBRecHitCollection> rechits(new HGCalTBRecHitCollection);
	std::map<int, TH2D*> histograms;
	for (int layer=1; layer<=8; layer++) {
		histograms[layer] = new TH2D(("event_"+std::to_string(currentEvent)+"__layer_"+std::to_string(layer)).c_str(),("event_"+std::to_string(currentEvent)+"__layer_"+std::to_string(layer)).c_str(), 13, -6, 6, 13, -6, 6);
	}

	//given code fragment
	for(unsigned int icell=0; icell<simHitCellIdE->size(); icell++){
	 int layer = ((simHitCellIdE->at(icell)>>19)&0x7F);
	 //int idcell = (simHitCellIdE->at(icell))&0xFF;
	 
	 double energy = simHitCellEnE->at(icell) / MIP2GeV_sim * ADCtoMIP_CERN[layer-1];

	 int cellno = (simHitCellIdE->at(icell)>>0)&0xFF;
	 std::pair<double,double> xy = geomc->position(cellno);
	 double x = xy.first / 10.;		//values are converted from mm to cm
	 double y =  xy.second / 10.;	//values are converted from mm to cm
	 int cellType = geomc->cellType(cellno);

	 std::pair<int, int> iuiv = TheCell.GetCellIUIVCoordinates(x, y);
	 
	 //std::cout<<"cell number: "<<cellno<<"  x: "<<x<<"  y: "<<y<<"   iu: "<<iuiv.first<<"   iv: "<<iuiv.second<<"   type: "<<cellType<<std::endl;

	 HGCalTBRecHit recHit(HGCalTBDetId(layer, 0, 0, iuiv.first, iuiv.second, cellType), energy, energy, energy, 0); 
	 	
	 /* Back and forth computation, if correct: Numbers should be identical!
	 std::pair<double, double> CellCenterXY = TheCell.GetCellCentreCoordinatesForPlots((recHit.id()).layer(), (recHit.id()).sensorIU(), (recHit.id()).sensorIV(), (recHit.id()).iu(), (recHit.id()).iv(), 128); 	//TODO: Hard Coded Number!
	 std::pair<int, int> iuiv2 = TheCell.GetCellIUIVCoordinates(CellCenterXY.first, CellCenterXY.second);
	 std::cout<<"x: "<<CellCenterXY.first<<"  y: "<<CellCenterXY.second<<"   iu: "<<iuiv2.first<<"   iv: "<<iuiv2.second<<std::endl;
	 std::cout<<std::endl;
	 */ 

	 recHit.setCellCenterCoordinate(x, y);
	 rechits->push_back(recHit);
	}	

	event.put(rechits, outputCollectionName);

}

void HGCalTBGenSimSource::endJob() {
}


#include "FWCore/Framework/interface/InputSourceMacros.h"

DEFINE_FWK_INPUT_SOURCE(HGCalTBGenSimSource);
