#include "HGCal/RawToDigi/plugins/HGCalTBGenSimSource.h"

using namespace std;

HGCalTBGenSimSource::HGCalTBGenSimSource(const edm::ParameterSet & pset, edm::InputSourceDescription const& desc) :  edm::ProducerSourceFromFiles(pset, desc, true),
	inputPathFormat(""),
	currentRun(-1),
	currentEvent(-1),
	rootFile(0)
{

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

	//find and fill the configured runs
	outputCollectionName = pset.getUntrackedParameter<std::string>("HGCalTBRecHitCollectionName");
	runEnergyMapFile = pset.getUntrackedParameter<std::string>("runEnergyMapFile"); 
	inputPathFormat = pset.getUntrackedParameter<std::string>("inputPathFormat");
	
	std::fstream map_file;
	map_file.open(runEnergyMapFile.c_str(), std::fstream::in);
	fillConfiguredRuns(map_file);
	map_file.close();

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
		std::cout<<"Fragment: "<<fragment<<std::endl;
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
			std::cout<<"Adding "<<filePath<<std::endl;
			fInfo.name = filePath;
			
			_fileNames.push_back(fInfo);			
		}
	}

	fileIterator = _fileNames.begin();
}

bool HGCalTBGenSimSource::setRunAndEventInfo(edm::EventID& id, edm::TimeValue_t& time, edm::EventAuxiliary::ExperimentType& evType)
{	

	if (fileIterator ==_fileNames.end()) return false; 		//end of files is reached

	if (currentRun == -1) {		//initial loading of a file
		currentRun = (*fileIterator).index;
		currentEvent = 0;
  	std::cout<<"Opening "<<(*fileIterator).name<<std::endl;
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
		//delete rootFile;
	}
	
	tree->GetEntry(currentEvent);

	currentEvent++;

	return true;

}


void HGCalTBGenSimSource::produce(edm::Event & event)
{
	std::auto_ptr<RunData> rd(new RunData);
	rd->energy = (*fileIterator).energy;
	rd->configuration = (*fileIterator).config;
	rd->runType = (*fileIterator).runType;
	rd->run = (*fileIterator).index;

	event.put(std::move(rd), "RunData");	

	std::auto_ptr<HGCalTBRecHitCollection> rechits(new HGCalTBRecHitCollection);

	//given code fragment
	HexGeometry geomc(true);
	for(unsigned int icell=0; icell<simHitCellIdE->size(); icell++){
	 int layerno = ((simHitCellIdE->at(icell)>>19)&0x7F);

	 int cellno = (simHitCellIdE->at(icell)>>0)&0xFF;
	 std::pair<double,double> xy = geomc.position(cellno);
	 double x = xy.first;
	 double y=  xy.second;
	 cout << "testing icell at layer "<< layerno <<":" << icell << " x:"<<x<<" y:"<<y<<endl;
	}

	//loop over all hits in the event:
	//HGCalTBRecHit recHit(digi.detid(), energy, energyLow, energyHigh, digi[iSample].tdc()); 
	HGCalTBRecHit recHit(0, 0, 0, 0, 0); 
	
	recHit.setCellCenterCoordinate(0., 0.);
	recHit.setEnergy(0.0);

	rechits->push_back(recHit);
	event.put(rechits, outputCollectionName);

	std::cout<<"Event "<<currentEvent<<" has "<<simHitCellIdE->size()<<" entries."<<std::endl;
}


#include "FWCore/Framework/interface/InputSourceMacros.h"

DEFINE_FWK_INPUT_SOURCE(HGCalTBGenSimSource);
