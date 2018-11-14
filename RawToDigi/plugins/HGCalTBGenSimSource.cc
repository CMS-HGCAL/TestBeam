#include "HGCal/RawToDigi/plugins/HGCalTBGenSimSource.h"


using namespace std;

//#define DEBUG

HGCalTBGenSimSource::HGCalTBGenSimSource(const edm::ParameterSet & pset, edm::InputSourceDescription const& desc) :  edm::ProducerSourceFromFiles(pset, desc, true),
	currentRun(-1),
	currentEvent(-1),
	rootFile(NULL)
{
	eventCounter = 0;

	//find and fill the configured runs
	RechitOutputCollectionName = pset.getParameter<std::string>("RechitOutputCollectionName");
	RunDataOutputCollectionName = pset.getParameter<std::string>("RunDataOutputCollectionName");

	produceDATURATracksInsteadOfDWCs = pset.getUntrackedParameter<bool>("produceDATURATracksInsteadOfDWCs", false);
	DWCOutputCollectionName = pset.getParameter<std::string>("DWCOutputCollectionName");
	DATURAOutputCollectionName = pset.getParameter<std::string>("DATURAOutputCollectionName");

	m_maskNoisyChannels = pset.getUntrackedParameter<bool>("MaskNoisyChannels", true);
	m_channelsToMask_filename = pset.getUntrackedParameter<std::string>("ChannelsToMaskFileName", "HGCal/CondObjects/data/noisyChannels.txt");

	energyNoise = pset.getUntrackedParameter<double>("energyNoise", 0.);	//in MIPS
	energyNoiseResolution = pset.getUntrackedParameter<double>("energyNoiseResolution", 0.);


	std::vector<double> v1(4, 1.0);
	wc_resolutions = pset.getUntrackedParameter<std::vector<double> >("wc_resolutions", v1);

	std::vector<double> v2(6, 1.0);
	datura_resolutions = pset.getUntrackedParameter<std::vector<double> >("datura_resolutions", v2);
	m_layerPositionFile = pset.getParameter<std::string>("layerPositionFile");

	beamEnergy = pset.getUntrackedParameter<double> ("beamEnergy", 250);
	beamParticlePDGID = pset.getUntrackedParameter<int> ("beamParticlePDGID", 211);
	setupConfiguration = pset.getUntrackedParameter<unsigned int> ("setupConfiguration", 1);

	switch (setupConfiguration) {
	case 1:	//July 2017
		N_layers_EE = 2;		//must shift remove layer 0 artificially
		N_layers_FH = 4;
		firstNLayersFH_asDaisies = 0;
		N_layers_BH = 12;
		break;
	case 2: //September 2017
		N_layers_EE = 7;
		N_layers_FH = 10;
		firstNLayersFH_asDaisies = 0;
		N_layers_BH = 12;
		break;
	case 3:	//October 2017 (29th May 2018: no simulated samples exist yet)
		N_layers_EE = 4;
		N_layers_FH = 6;
		firstNLayersFH_asDaisies = 6;
		N_layers_BH = 12;
		break;
	case 6:	//DESY March 2018, "CaloRuns, i.e. config #4"
	case 7:
	case 8:
	case 9:
	case 10:
		N_layers_EE = 3;
		N_layers_FH = 0;
		firstNLayersFH_asDaisies = 0;
		N_layers_BH = 0;
		break;
	case 11:	//DESY March 2018, configs #5.1/5.2
	case 12:
		N_layers_EE = 2;
		N_layers_FH = 0;
		firstNLayersFH_asDaisies = 0;
		N_layers_BH = 0;
		break;
	case 13:
	default:
		N_layers_EE = 28;
		N_layers_FH = 0;
		firstNLayersFH_asDaisies = 0;
		N_layers_BH = 0;
		break;
	case 22:					//October 2018 - setup 1, H2
		N_layers_EE = 28;
		N_layers_FH = 12;
		firstNLayersFH_asDaisies = 9;
		N_layers_BH = 12;
		break;
	case 23:					//October 2018 - setup 2, H2
		N_layers_EE = 28;
		N_layers_FH = 11;
		firstNLayersFH_asDaisies = 9;
		N_layers_BH = 12;
		break;
	case 24:					//October 2018 - setup 3, H2
		N_layers_EE = 7;
		N_layers_FH = 12;
		firstNLayersFH_asDaisies = 12;
		N_layers_BH = 12;
		break;
	}

	GeVToMip = pset.getUntrackedParameter<double>("GeVToMip", 1. / (86.5e-6));			//value from Shilpi

	areaSpecification = pset.getUntrackedParameter<std::string>("areaSpecification", "H2");
	//given in cm
	if (areaSpecification == "H6A_October2017") {
		referenceTracking_zPositions.push_back(-500.);
	} else if (areaSpecification == "H2_Summer2017") {
		referenceTracking_zPositions.push_back(-109.);
		referenceTracking_zPositions.push_back(-235.);
		referenceTracking_zPositions.push_back(-1509.);
		referenceTracking_zPositions.push_back(-1769.);
	} else if (areaSpecification == "H2") {
		referenceTracking_zPositions.push_back(-120.);
		referenceTracking_zPositions.push_back(-246.);
		referenceTracking_zPositions.push_back(-1520.);
		referenceTracking_zPositions.push_back(-1780.);
	} else if (areaSpecification == "H2_October2018") {
		referenceTracking_zPositions.push_back(-160.);
		referenceTracking_zPositions.push_back(-880.);
		referenceTracking_zPositions.push_back(-2700.);
		referenceTracking_zPositions.push_back(-3200.);
	} else if (areaSpecification == "DESY_T21_March2018") {
		referenceTracking_zPositions.push_back(0.);
		referenceTracking_zPositions.push_back(15.3);
		referenceTracking_zPositions.push_back(30.5);
		referenceTracking_zPositions.push_back(64.8);
		referenceTracking_zPositions.push_back(80.0);
		referenceTracking_zPositions.push_back(95.3);
	}


	//reads the layer positions
	std::fstream file;
	char fragment[100];
	int readCounter = -1;

	file.open(m_layerPositionFile.c_str(), std::fstream::in);
	std::cout << "Reading file " << m_layerPositionFile << " -open: " << file.is_open() << std::endl;
	int layer = 0;
	while (file.is_open() && !file.eof()) {
		readCounter++;
		file >> fragment;
		if (readCounter == 0) layer = atoi(fragment);
		if (readCounter == 1) {
			layerPositions[layer] = atof(fragment) / 10.;   //values are given in mm and should be converted into cm
			readCounter = -1;
		}
	}
	file.close();

	//indicate the physics list
	physicsListUsed = pset.getUntrackedParameter<std::string>("physicsListUsed", "");

	if (physicsListUsed == "FTF_BIC")
		_enumPhysicsListUsed = HGCAL_TB_SIM_FTF_BIC;
	else if (physicsListUsed == "FTFP_BERT")
		_enumPhysicsListUsed = HGCAL_TB_SIM_FTFP_BERT ;
	else if (physicsListUsed == "FTFP_BERT_EML")
		_enumPhysicsListUsed = HGCAL_TB_SIM_FTFP_BERT_EML ;
	else if (physicsListUsed == "FTFP_BERT_EMM")
		_enumPhysicsListUsed = HGCAL_TB_SIM_FTFP_BERT_EMM;
	else if (physicsListUsed == "FTFP_BERT_EMZ")
		_enumPhysicsListUsed = HGCAL_TB_SIM_FTFP_BERT_EMZ;
	else if (physicsListUsed == "FTFP_BERT_HP_EML")
		_enumPhysicsListUsed = HGCAL_TB_SIM_FTFP_BERT_HP_EML;
	else if (physicsListUsed == "FTFP_BERT_EMY")
		_enumPhysicsListUsed = HGCAL_TB_SIM_FTFP_BERT_EMY;
	else if (physicsListUsed == "QGSP_BERT")
		_enumPhysicsListUsed = HGCAL_TB_SIM_QGSP_BERT;
	else if (physicsListUsed == "QGSP_BERT_EML")
		_enumPhysicsListUsed = HGCAL_TB_SIM_QGSP_BERT_EML;
	else if (physicsListUsed == "QGSP_BERT_HP_EML")
		_enumPhysicsListUsed = HGCAL_TB_SIM_QGSP_BERT_HP_EML;
	else if (physicsListUsed == "QGSP_FTFP_BERT")
		_enumPhysicsListUsed = HGCAL_TB_SIM_QGSP_FTFP_BERT;
	else if (physicsListUsed == "QGSP_FTFP_BERT_EML")
		_enumPhysicsListUsed = HGCAL_TB_SIM_QGSP_FTFP_BERT_EML;
	else if (physicsListUsed == "QGSP_FTFP_BERT_EML_New")
		_enumPhysicsListUsed = HGCAL_TB_SIM_QGSP_FTFP_BERT_EML_New;
	else if (physicsListUsed == "QGSP_FTFP_BERT_EMM")
		_enumPhysicsListUsed = HGCAL_TB_SIM_QGSP_FTFP_BERT_EMM;
	else
		_enumPhysicsListUsed = HGCAL_TB_SIM;

	produces <HGCalTBRecHitCollection>(RechitOutputCollectionName);
	produces<RunData>(RunDataOutputCollectionName);

	if (!produceDATURATracksInsteadOfDWCs) produces<std::map<int, WireChamberData> >(DWCOutputCollectionName);
	else produces<std::vector<HGCalTBDATURATelescopeData> >(DATURAOutputCollectionName);

	if (fileNames()[0] != "file:DUMMY") {
		for (int i = 0; i < (int)(fileNames().size()); i++) {
			FileInfo fInfo;
			fInfo.index = i;
			fInfo.energy = beamEnergy;
			fInfo.pdgID = beamParticlePDGID;
			fInfo.config = setupConfiguration;
			fInfo.name = fileNames()[i];
			_fileNames.push_back(fInfo);
		}
		fileIterator = _fileNames.begin();
	}

	_e_mapFile = pset.getUntrackedParameter<std::string>("e_mapFile_CERN");
	HGCalCondObjectTextIO io(0);
	edm::FileInPath fip(_e_mapFile);

	if (!io.load(fip.fullPath(), essource_.emap_)) {
		throw cms::Exception("Unable to load electronics map");
	};

	geomc = new HexGeometry(false);

	tree = 0;
	simHitCellIdEE = 0;
	simHitCellEnEE = 0;
	simHitCellIdFH = 0;
	simHitCellEnFH = 0;
	simHitCellIdBH = 0;
	simHitCellEnBH = 0;
	beamX	 = 0;
	beamY 	 = 0;
	beamP 	 = 0;

	randgen = new TRandom();

	if ( m_maskNoisyChannels ) {
		FILE* file;
		char buffer[300];
		//edm::FileInPath fip();
		file = fopen (m_channelsToMask_filename.c_str() , "r");
		if (file == NULL) {
			perror ("Error opening noisy channels file"); exit(1);
		} else {

			while ( ! feof (file) ) {
				if ( fgets (buffer , 300 , file) == NULL ) break;
				const char* index = buffer;
				int layer, skiroc, channel, ptr, nval;
				nval = sscanf( index, "%d %d %d %n", &layer, &skiroc, &channel, &ptr );
				int skiId = HGCAL_TB_GEOMETRY::N_SKIROC_PER_HEXA * layer + (HGCAL_TB_GEOMETRY::N_SKIROC_PER_HEXA - skiroc) % HGCAL_TB_GEOMETRY::N_SKIROC_PER_HEXA + 1;
				if ( nval == 3 ) {
					HGCalTBElectronicsId eid(skiId, channel);
					if (essource_.emap_.existsEId(eid.rawId()))
						m_noisyChannels.push_back(eid.rawId());
				} else continue;
			}
		}
		fclose (file);
	}

}


bool HGCalTBGenSimSource::setRunAndEventInfo(edm::EventID& id, edm::TimeValue_t& time, edm::EventAuxiliary::ExperimentType& evType)
{
	if (fileIterator == _fileNames.end()) {
		std::cout << "End of the files in the setRundAndEvent is reached..." << std::endl;
		delete randgen;
		return false; 		//end of files is reached
	}

	if (currentRun == -1) {		//initial loading of a file
		currentRun = (*fileIterator).index;
		currentEvent = 0;

		std::string fileName = (*fileIterator).name;
		std::string filePrefix = "file:";
		fileName.replace(fileName.find(filePrefix), filePrefix.size(), "");
		if (access( (fileName).c_str(), F_OK ) == -1) {
			std::cout << fileName << " does not exist and is skipped" << std::endl;
			fileIterator++;
			currentRun = -1;
			if (! setRunAndEventInfo(id, time, evType)) return false;
		}
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

		tree->SetBranchAddress("simHitCellIdEE", &simHitCellIdEE, &b_simHitCellIdEE);
		tree->SetBranchAddress("simHitCellEnEE", &simHitCellEnEE, &b_simHitCellEnEE);
		tree->SetBranchAddress("simHitCellIdFH", &simHitCellIdFH, &b_simHitCellIdFH);
		tree->SetBranchAddress("simHitCellEnFH", &simHitCellEnFH, &b_simHitCellEnFH);
		tree->SetBranchAddress("simHitCellIdBH", &simHitCellIdBH, &b_simHitCellIdBH);
		tree->SetBranchAddress("simHitCellEnBH", &simHitCellEnBH, &b_simHitCellEnBH);
		tree->SetBranchAddress("xBeam", &beamX, &b_beamX);
		tree->SetBranchAddress("yBeam", &beamY, &b_beamY);
		tree->SetBranchAddress("pBeam", &beamP, &b_beamP);
	}

	if (currentEvent == tree->GetEntries()) {
		fileIterator++;
		currentRun = -1;
		if (! setRunAndEventInfo(id, time, evType)) return false;
	}

	tree->GetEntry(currentEvent);

	currentEvent++;

	return true;

}



void HGCalTBGenSimSource::produce(edm::Event & event)
{
	if (fileIterator == _fileNames.end()) {
		std::cout << "End of the files in the producer is reached..." << std::endl;
		return;
	}

	eventCounter++;	//indexes each event chronologically passing this plugin

	//first: fill the rechits
	std::unique_ptr<HGCalTBRecHitCollection> rechits(new HGCalTBRecHitCollection);

	//EE part
#ifdef DEBUG
	std::cout << "Hits in EE: " << simHitCellIdEE->size() << std::endl;
#endif
	for (unsigned int icell = 0; icell < simHitCellIdEE->size(); icell++) {
		uint32_t id = simHitCellIdEE->at(icell);
		int subdet, z, layer, wafer, celltype, cellno;
		HGCalTestNumbering::unpackHexagonIndex(id, subdet, z, layer, wafer, celltype, cellno);
		bool onDaisyStructure = false;

		if (setupConfiguration == 1)
			layer -= 1;		//remove the first layer as it is not present in the data
		if (layer == 0) continue;		//no layers with index 0 allowed
		double energy = simHitCellEnEE->at(icell);
#ifdef DEBUG
		std::cout << icell << "  layer: " << layer << "   wafer:  " << wafer << "   cellno:  " << cellno << "  energy: " << energy << std::endl;
#endif

		makeRecHit(layer, wafer, onDaisyStructure, cellno, energy, rechits);

	}

	//FH part
#ifdef DEBUG
	std::cout << "Hits in FH: " << simHitCellIdFH->size() << std::endl;
#endif
	for (unsigned int icell = 0; icell < simHitCellIdFH->size(); icell++) {

		uint32_t id = simHitCellIdFH->at(icell);
		int subdet, z, layer, wafer, celltype, cellno;
		HGCalTestNumbering::unpackHexagonIndex(id, subdet, z, layer, wafer, celltype, cellno);
		bool onDaisyStructure = (layer<=firstNLayersFH_asDaisies);
		
		layer = layer + N_layers_EE;
		double energy = simHitCellEnFH->at(icell);
#ifdef DEBUG
		std::cout << icell << "  layer: " << layer << "   wafer:  " << wafer << "   cellno:  " << cellno << "  energy: " << energy << std::endl;
#endif
		makeRecHit(layer, wafer, onDaisyStructure, cellno, energy, rechits);

	}

	event.put(std::move(rechits), RechitOutputCollectionName);


	//second: fill the run data
	std::unique_ptr<RunData> rd(new RunData);
	rd->energy = (*fileIterator).energy;		//mean energy of the beam configuration
	rd->configuration = (*fileIterator).config;
	rd->pdgID = (*fileIterator).pdgID;
	rd->runType = _enumPhysicsListUsed;
	rd->run = (*fileIterator).index;
	rd->event = eventCounter;
	rd->booleanUserRecords.add("hasDanger", false);
	rd->doubleUserRecords.add("trueEnergy", beamP);


	if (!produceDATURATracksInsteadOfDWCs) {
		//third option A: add fake dwcs from xBeam, yBeam
		std::unique_ptr<std::map<int, WireChamberData> > dwcs(new std::map<int, WireChamberData>);
		rd->booleanUserRecords.add("hasValidDWCMeasurement", true);
		for (size_t d = 0; d < referenceTracking_zPositions.size(); d++) {
			double dwc_x = 10 * beamX + randgen->Gaus(0, wc_resolutions[d]);		//conversion to cm performed at a different the track computation stage
			double dwc_y = 10 * beamY + randgen->Gaus(0, wc_resolutions[d]);		//conversion to cm performed at a different the track computation stage

			WireChamberData* dwc = new WireChamberData(d + 1, dwc_x , dwc_y, referenceTracking_zPositions[d]);
			dwc->goodMeasurement_X = dwc->goodMeasurement_Y = dwc->goodMeasurement = true;
			dwc->res_x = dwc->res_y = wc_resolutions[d];
			dwc->averageHitMultiplicty = 1.;

			(*dwcs)[d] = *dwc;
		}
		event.put(std::move(dwcs), DWCOutputCollectionName);
	} else {
		std::unique_ptr<std::vector<HGCalTBDATURATelescopeData> > DATURATracks(new std::vector<HGCalTBDATURATelescopeData>);
		rd->booleanUserRecords.add("hasValidDATURAMeasurement", true);
		//in this block: logical copy from HGCalTBDATURATelescopeProducer.cc
		HGCalTBDATURATelescopeData DATURATelescopeTrack(1);		//just one track by default
		LineFitter TripletTrack1X;
		LineFitter TripletTrack1Y;
		LineFitter TripletTrack2X;
		LineFitter TripletTrack2Y;
		for (int MIMOSA_index = 1; MIMOSA_index <= 6; MIMOSA_index++) {
			double x_measured = beamX + randgen->Gaus(0, datura_resolutions[MIMOSA_index - 1] / 10);		//value in cm
			double y_measured = beamY + randgen->Gaus(0, datura_resolutions[MIMOSA_index - 1] / 10);		//value in cm
			double z_measured = referenceTracking_zPositions[MIMOSA_index - 1];

			DATURATelescopeTrack.addPointForTracking(x_measured, y_measured, z_measured, datura_resolutions[MIMOSA_index - 1], datura_resolutions[MIMOSA_index - 1]);
			if (MIMOSA_index <= 3) {
				TripletTrack1X.addPoint(z_measured, x_measured, datura_resolutions[MIMOSA_index - 1]);
				TripletTrack1Y.addPoint(z_measured, y_measured, datura_resolutions[MIMOSA_index - 1]);
			} else {
				TripletTrack2X.addPoint(z_measured, x_measured, datura_resolutions[MIMOSA_index - 1]);
				TripletTrack2Y.addPoint(z_measured, y_measured, datura_resolutions[MIMOSA_index - 1]);
			}
		}

		TripletTrack1X.fit();       TripletTrack1Y.fit();
		TripletTrack2X.fit();       TripletTrack2Y.fit();

		for (std::map<int, double>::iterator layerIt = layerPositions.begin(); layerIt != layerPositions.end(); layerIt++) {
			double layer_ref_x = TripletTrack2X.eval(layerIt->second);
			double layer_ref_y = TripletTrack2Y.eval(layerIt->second);
			double layer_ref_x_chi2 = TripletTrack2X.GetChisquare();
			double layer_ref_y_chi2 = TripletTrack2Y.GetChisquare();

			//mean between two triplets for the intermediate layer
			if (layerIt->first <= 1 ) {
				layer_ref_x += TripletTrack1X.eval(layerIt->second);
				layer_ref_y += TripletTrack1Y.eval(layerIt->second);
				layer_ref_x_chi2 += TripletTrack1X.GetChisquare();
				layer_ref_y_chi2 += TripletTrack1Y.GetChisquare();

				layer_ref_x /= 2.;
				layer_ref_y /= 2.;
				layer_ref_x_chi2 /= 2.;
				layer_ref_y_chi2 /= 2.;
			}
			DATURATelescopeTrack.addLayerReference(layerIt->first, layer_ref_x, layer_ref_y, layer_ref_x_chi2, layer_ref_y_chi2);
		}

		float kinkAngleX_DUT1 = atan(TripletTrack1X.getM()) - atan(TripletTrack2X.getM());
		float kinkAngleY_DUT1 = atan(TripletTrack1Y.getM()) - atan(TripletTrack2Y.getM());
		DATURATelescopeTrack.floatUserRecords.add("kinkAngleX_DUT1", kinkAngleX_DUT1);
		DATURATelescopeTrack.floatUserRecords.add("kinkAngleY_DUT1", kinkAngleY_DUT1);
		DATURATracks->push_back(DATURATelescopeTrack);

		event.put(std::move(DATURATracks), DATURAOutputCollectionName);
	}


	event.put(std::move(rd), RunDataOutputCollectionName);
}

void HGCalTBGenSimSource::makeRecHit(int layer, int wafer, bool onDaisyStructure, int cellno, double energy, std::unique_ptr<HGCalTBRecHitCollection> &rechits) {
	int cellType = geomc->cellType(cellno);

	std::pair<double, double> xy = geomc->position_cell(cellno);
	double x = xy.first;
	double y =  xy.second;

	double X = 0, Y = 0;
	if (onDaisyStructure) {
		std::pair<double, double> XY = geomc->position_wafer(wafer);
		X = XY.first;
		Y =  XY.second;
	}

	//get hexagonal coordinates
	std::pair<int, int> iuiv = TheCell.GetCellIUIVCoordinates(x, y);
	std::pair<int, int> iUiV = TheCell.GetSensorIUIVCoordinates(X, Y);
#ifdef DEBUG
	std::cout << "X: " << X << "  Y: " << Y << "    x: " << x << " y: " << y << std::endl;
	std::cout << "u: " << iuiv.first << "  v: " << iuiv.second << std::endl;
	std::cout << "U: " << iUiV.first << "  V: " << iUiV.second << std::endl;
#endif
	HGCalTBRecHit recHit(HGCalTBDetId(layer, iUiV.first, iUiV.second, iuiv.first, iuiv.second, cellType), 0., 0., 0., 0., 0);
	//attention! the electric ID to skiroc mapping takes the cellType information as well which is only 0, 2 in the simulation!
	//any analysis on the simulated data must respect that by computing MIP-ADC factors in the same way and by converting energies into MIP units
	uint32_t EID = essource_.emap_.detId2eid(recHit.id());
	HGCalTBElectronicsId eid(EID);
	if ( !essource_.emap_.existsEId(eid.rawId()) || std::find(m_noisyChannels.begin(), m_noisyChannels.end(), eid.rawId()) != m_noisyChannels.end() ) {
#ifdef DEBUG
		std::cout << "Not a valid hit" << std::endl;
#endif
		return;
	}

	int skiRocIndex = eid.iskiroc_rawhit();
	skiRocIndex = skiRocIndex % HGCAL_TB_GEOMETRY::N_SKIROC_PER_HEXA;

	recHit.setTime(-1);
	recHit.setCellCenterCoordinate(x + X, y + Y);


	energy *= GeVToMip;
	//additional noise to the energy in MIPs
	energy += randgen->Gaus(energyNoise, energyNoiseResolution);

	//gain assumptions from average estimates, 28th June 2018
	//1 MIP - 40 HG ADC
	//1 LG ADC - 9 HG ADC
	//1 TOT ADC - 1./0.194 LG ADC + 192 offset
	recHit._energyHigh = energy * 40;
	recHit._energyLow = recHit._energyHigh / 9.;
	recHit._energyTot = recHit._energyLow * 0.194 + 192;

	recHit.setEnergy(energy);
	recHit.setFlag(HGCalTBRecHit::kGood);


	rechits->push_back(recHit);
}

#include "FWCore/Framework/interface/InputSourceMacros.h"

DEFINE_FWK_INPUT_SOURCE(HGCalTBGenSimSource);
