#include "HGCal/RawToDigi/plugins/HGCalTBGenSimSource.h"

using namespace std;

HGCalTBGenSimSource::HGCalTBGenSimSource(const edm::ParameterSet & pset, edm::InputSourceDescription const& desc) :  edm::ProducerSourceFromFiles(pset, desc, true),
	inputPathFormat(""),
	currentRun(-1),
	currentEvent(-1),
	rootFile(NULL)
{
	eventCounter = 0;

	//find and fill the configured runs
	outputCollectionName = pset.getParameter<std::string>("OutputCollectionName");
	runEnergyMapFile = pset.getUntrackedParameter<std::string>("runEnergyMapFile"); 
	inputPathFormat = pset.getUntrackedParameter<std::string>("inputPathFormat");
	energyNoise = pset.getParameter<double>("energyNoise");
	energyNoiseResolution = pset.getParameter<double>("energyNoiseResolution");
	smearingResolution = pset.getParameter<double>("MWCSmearingResolution")/1000.;		//configuration is in microns, convert to cm

	produces <HGCalTBRecHitCollection>(outputCollectionName);
	produces<RunData>("RunData");
	produces<MultiWireChambers>("MultiWireChambers");


	if (fileNames()[0] != "file:DUMMY") {
		for (int i = 0; i<(int)(fileNames().size()); i++) {
			FileInfo fInfo;
			fInfo.index = 0;
			fInfo.energy = -1;
			fInfo.runType = "";
			fInfo.config = -1;
			fInfo.name = fileNames()[i];
			_fileNames.push_back(fInfo);
		}
		fileIterator = _fileNames.begin();
	} else {
		std::fstream map_file;
		map_file.open(runEnergyMapFile.c_str(), std::fstream::in);
		fillConfiguredRuns(map_file);
		map_file.close();	
	}
	


	 
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


  	randgen = new TRandom();
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
	if (fileIterator == _fileNames.end()) {
		std::cout<<"End of the files in the producer is reached..."<<std::endl;
		return;
	}

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
			
		/* Back and forth computation, if correct: Numbers should be identical!
		std::pair<double, double> CellCenterXY = TheCell.GetCellCentreCoordinatesForPlots((recHit.id()).layer(), (recHit.id()).sensorIU(), (recHit.id()).sensorIV(), (recHit.id()).iu(), (recHit.id()).iv(), 128); 	//TODO: Hard Coded Number!
		std::pair<int, int> iuiv2 = TheCell.GetCellIUIVCoordinates(CellCenterXY.first, CellCenterXY.second);
		std::cout<<"x: "<<CellCenterXY.first<<"  y: "<<CellCenterXY.second<<"   iu: "<<iuiv2.first<<"   iv: "<<iuiv2.second<<std::endl;
		std::cout<<std::endl;
		*/ 

		recHit.setCellCenterCoordinate(x, y);
	
		uint32_t EID = essource_.emap_.detId2eid(recHit.id());
		HGCalTBElectronicsId eid(EID);	 

		//attention! the electric ID to skiroc mapping takes the cellType information as well which is only 0, 2 in the simulation!
		//any analysis on the simulated data must respect that by computing MIP-ADC factors in the same way and by converting energies into MIP units
		int skiRocIndex = (eid.iskiroc() - 1) > 0 ? eid.iskiroc() - 1 : 0;		

		double energy = simHitCellEnE->at(icell) / MIP2GeV_sim;
	 	
		//additional noise to the energy in MIPs
		energy += randgen->Gaus(energyNoise, energyNoiseResolution);

		energy *= ADCtoMIP_CERN[skiRocIndex];

	 	recHit.setEnergy(energy);
	 	recHit._energyLow = energy;
	 	recHit._energyHigh = energy;
	 	recHit._energyTot = energy;

		rechits->push_back(recHit);
	}	
	event.put(rechits, outputCollectionName);

	
	//second: add fake multi-wire chambers from xBeam, yBeam
	double x1_mc = beamX + randgen->Gaus(0, smearingResolution);
	double y1_mc = beamY + randgen->Gaus(0, smearingResolution);
	double x2_mc = beamX + randgen->Gaus(0, smearingResolution);
	double y2_mc = beamY + randgen->Gaus(0, smearingResolution);
	
	std::auto_ptr<MultiWireChambers> mwcs(new MultiWireChambers);	
	mwcs->push_back(MultiWireChamberData(1, x1_mc*cos(90.0*M_PI/180.0) + sin(90.0*M_PI/180.0)*y1_mc, -x1_mc*sin(90.0*M_PI/180.0) + cos(90.0*M_PI/180.0)*y1_mc, -126.-147.));
	mwcs->push_back(MultiWireChamberData(2, x2_mc*cos(90.0*M_PI/180.0) + sin(90.0*M_PI/180.0)*y2_mc, -x2_mc*sin(90.0*M_PI/180.0) + cos(90.0*M_PI/180.0)*y2_mc, -147.));
	
	event.put(std::move(mwcs), "MultiWireChambers");		
	

	//third: fill the run data
	std::auto_ptr<RunData> rd(new RunData);
	rd->energy = (*fileIterator).energy;		//mean energy of the beam configuration
	rd->configuration = (*fileIterator).config;
	rd->runType = (*fileIterator).runType;
	rd->run = (*fileIterator).index;
	rd->event = eventCounter;
	rd->hasDanger = false;
	rd->trueEnergy = beamP;						//energy as truely simulated

	bool _hasValidMWCMeasurement = true;
	_hasValidMWCMeasurement = (x1_mc != -99.9) && _hasValidMWCMeasurement;
	_hasValidMWCMeasurement = (y1_mc != -99.9) && _hasValidMWCMeasurement;
	_hasValidMWCMeasurement = (x2_mc != -99.9) && _hasValidMWCMeasurement;
	_hasValidMWCMeasurement = (y2_mc != -99.9) && _hasValidMWCMeasurement;
	rd->hasValidMWCMeasurement = _hasValidMWCMeasurement;

	event.put(std::move(rd), "RunData");	


}

void HGCalTBGenSimSource::endJob() {
	delete randgen;
}


#include "FWCore/Framework/interface/InputSourceMacros.h"

DEFINE_FWK_INPUT_SOURCE(HGCalTBGenSimSource);
