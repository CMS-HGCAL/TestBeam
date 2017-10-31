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
	DWCOutputCollectionName = pset.getParameter<std::string>("DWCOutputCollectionName"); 
	RunDataOutputCollectionName = pset.getParameter<std::string>("RunDataOutputCollectionName"); 
	m_maskNoisyChannels = pset.getUntrackedParameter<bool>("MaskNoisyChannels", true);
  	m_channelsToMask_filename = pset.getUntrackedParameter<std::string>("ChannelsToMaskFileName","HGCal/CondObjects/data/noisyChannels.txt");

	energyNoise = pset.getUntrackedParameter<double>("energyNoise", 0.);	//in MIPS
	energyNoiseResolution = pset.getUntrackedParameter<double>("energyNoiseResolution", 0.);
  	std::vector<double> v1(4, 1.0);
	wc_resolutions = pset.getUntrackedParameter<std::vector<double> >("wc_resolutions", v1);


  	beamEnergy = pset.getUntrackedParameter<unsigned int> ("beamEnergy", 250);
  	beamParticlePDGID = pset.getUntrackedParameter<std::string> ("beamParticlePDGID", "211");
  	setupConfiguration = pset.getUntrackedParameter<unsigned int> ("setupConfiguration", 1);

  	switch(setupConfiguration) {
  		case 1:
  			N_layers_EE = 2;		//must shift remove layer 0 artificially
  			N_layers_FH = 4;
  			N_layers_BH = 12;
  			break;
    	case 2:
  			N_layers_EE = 7;
  			N_layers_FH = 10;
  			N_layers_BH = 12;
  			break;
    	case 3:
  			N_layers_EE = 4;		//not sure about this number, 31.10.18 (T. Quast)
  			N_layers_FH = 6;		//not sure about this number, 31.10.18 (T. Quast)
  			N_layers_BH = 12;		//fix
  			break;
  		default:
    	case 4:
  			N_layers_EE = 4;		//not sure about this number, 31.10.18 (T. Quast)
  			N_layers_FH = 6;		//not sure about this number, 31.10.18 (T. Quast)
  			N_layers_BH = 12;
  			break;
  	}

  	GeVToMip = pset.getUntrackedParameter<double>("GeVToMip", 1./(86.5e-6));			//value from Shilpi

	areaSpecification = pset.getUntrackedParameter<std::string>("areaSpecification", "H2");

	if (areaSpecification=="H6A") {
		dwc_zPositions.push_back(-500.);
	} else {
		dwc_zPositions.push_back(-109.);
		dwc_zPositions.push_back(-235.);
		dwc_zPositions.push_back(-1509.);
		dwc_zPositions.push_back(-1769.);
	}

	produces <HGCalTBRecHitCollection>(RechitOutputCollectionName);
	produces<WireChambers>(DWCOutputCollectionName);
	produces<RunData>(RunDataOutputCollectionName);

	if (fileNames()[0] != "file:DUMMY") {
		for (int i = 0; i<(int)(fileNames().size()); i++) {
			FileInfo fInfo;
			fInfo.index = i;
			fInfo.energy = beamEnergy;
			fInfo.runType = beamParticlePDGID;
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
		setRunAndEventInfo(id, time, evType);
	}
	
	tree->GetEntry(currentEvent);

	currentEvent++;

	return true;

}


//TODO: Why are the cell calibration pads ignored?
//Check the orientation of the sensors

void HGCalTBGenSimSource::produce(edm::Event & event)
{	
	if (fileIterator == _fileNames.end()) {
		std::cout<<"End of the files in the producer is reached..."<<std::endl;
		return;
	}

	eventCounter++;	//indexes each event chronologically passing this plugin 

	//first: fill the rechits
	std::auto_ptr<HGCalTBRecHitCollection> rechits(new HGCalTBRecHitCollection);

	//EE part
	#ifdef DEBUG
		std::cout<<"Hits in EE: "<<simHitCellIdEE->size()<<std::endl;
	#endif
	for(unsigned int icell=0; icell<simHitCellIdEE->size(); icell++){
		int layer = ((simHitCellIdEE->at(icell)>>19)&0x7F);
		if (setupConfiguration==1)
			layer-=1;		//remove the first layer as it is not present in the data		
		if (layer==0) continue;		//no layers with index 0 allowed
		int cellno = (simHitCellIdEE->at(icell)>>0)&0xFF;
	
		double energy = simHitCellEnEE->at(icell);
		#ifdef DEBUG
			std::cout<<icell<<"  layer: "<<layer<<"   cellno:  "<<cellno<<"  energy: "<<energy<<std::endl;
		#endif
		makeRecHit(layer, cellno, energy, rechits);
		
	}	

	//FH part
	#ifdef DEBUG
		std::cout<<"Hits in FH: "<<simHitCellIdFH->size()<<std::endl;
	#endif
	for(unsigned int icell=0; icell<simHitCellIdFH->size(); icell++){
		int layer = ((simHitCellIdFH->at(icell)>>19)&0x7F) + N_layers_EE;		
		int cellno = (simHitCellIdFH->at(icell)>>0)&0xFF;
	
		double energy = simHitCellEnFH->at(icell);
		#ifdef DEBUG
			std::cout<<icell<<"  layer: "<<layer<<"   cellno:  "<<cellno<<"  energy: "<<energy<<std::endl;
		#endif
		makeRecHit(layer, cellno, energy, rechits);
		
	}	
	
	event.put(rechits, RechitOutputCollectionName);

	
	//second: add fake dwcs from xBeam, yBeam
	std::auto_ptr<WireChambers> dwcs(new WireChambers);	
	for (size_t d=0; d<dwc_zPositions.size(); d++) {
		double dwc_x = beamX * 10. + randgen->Gaus(0, wc_resolutions[d]);
		double dwc_y = beamY * 10. + randgen->Gaus(0, wc_resolutions[d]);

		WireChamberData* dwc = new WireChamberData(d+1, dwc_x*cos(90.0*M_PI/180.0) + sin(90.0*M_PI/180.0)*dwc_y, -dwc_x*sin(90.0*M_PI/180.0) + cos(90.0*M_PI/180.0)*dwc_y, dwc_zPositions[d]);
		dwc->goodMeasurement_X = dwc->goodMeasurement_Y = dwc->goodMeasurement = true;
		dwc->res_x = dwc->res_y = wc_resolutions[d];
		
		dwcs->push_back(*dwc);
	}
	event.put(std::move(dwcs), DWCOutputCollectionName);		


	//third: fill the run data
	std::auto_ptr<RunData> rd(new RunData);
	rd->energy = (*fileIterator).energy;		//mean energy of the beam configuration
	rd->configuration = (*fileIterator).config;
	rd->runType = (*fileIterator).runType;
	rd->run = (*fileIterator).index;
	rd->event = eventCounter;
	rd->booleanUserRecords.add("hasDanger", false);
	rd->doubleUserRecords.add("trueEnergy", beamP);
	rd->booleanUserRecords.add("hasValidMWCMeasurement", true);
	event.put(std::move(rd), RunDataOutputCollectionName);	

}

void HGCalTBGenSimSource::beginJob() {
  if( m_maskNoisyChannels ){
    FILE* file;
    char buffer[300];
    //edm::FileInPath fip();
    file = fopen (m_channelsToMask_filename.c_str() , "r");
    if (file == NULL){
      perror ("Error opening noisy channels file"); exit(1); 
    } else{
    
      while ( ! feof (file) ){
        if ( fgets (buffer , 300 , file) == NULL ) break;
        const char* index = buffer;
        int layer,skiroc,channel,ptr,nval;
        nval=sscanf( index, "%d %d %d %n",&layer,&skiroc,&channel,&ptr );
        int skiId=HGCAL_TB_GEOMETRY::N_SKIROC_PER_HEXA*layer+(HGCAL_TB_GEOMETRY::N_SKIROC_PER_HEXA-skiroc)%HGCAL_TB_GEOMETRY::N_SKIROC_PER_HEXA+1;
        if( nval==3 ){
          HGCalTBElectronicsId eid(skiId,channel);      
          if (essource_.emap_.existsEId(eid.rawId()))
            m_noisyChannels.push_back(eid.rawId());
        } else continue;
      }
    }
    fclose (file);
  }
}

void HGCalTBGenSimSource::endJob() {
	delete randgen;
}

void HGCalTBGenSimSource::makeRecHit(int layer, int cellno, double energy, std::auto_ptr<HGCalTBRecHitCollection> &rechits) {
		std::pair<double,double> xy = geomc->position_cell(cellno);
		double x = xy.first / 10.;		//values are converted from mm to cm
		double y =  xy.second / 10.;	//values are converted from mm to cm
		int cellType = geomc->cellType(cellno);

		std::pair<int, int> iuiv = TheCell.GetCellIUIVCoordinates(x, y);


		HGCalTBRecHit recHit(HGCalTBDetId(layer, 0, 0, iuiv.first, iuiv.second, cellType), 0., 0., 0., 0., 0); 
		//attention! the electric ID to skiroc mapping takes the cellType information as well which is only 0, 2 in the simulation!
		//any analysis on the simulated data must respect that by computing MIP-ADC factors in the same way and by converting energies into MIP units
		uint32_t EID = essource_.emap_.detId2eid(recHit.id());
		HGCalTBElectronicsId eid(EID);	 
	    if( !essource_.emap_.existsEId(eid.rawId()) || std::find(m_noisyChannels.begin(),m_noisyChannels.end(),eid.rawId())!=m_noisyChannels.end() )
	      return;

		int skiRocIndex = (eid.iskiroc() - 1) > 0 ? eid.iskiroc() - 1 : 0;
		skiRocIndex = skiRocIndex % HGCAL_TB_GEOMETRY::N_SKIROC_PER_HEXA;		
		int channelIndex = eid.ichan();

	    recHit.setBoard(layer);
	    recHit.setSkiroc(skiRocIndex);
	    recHit.setChannel(channelIndex);
	    recHit.setTime(-1);
		recHit.setCellCenterCoordinate(x, y);
	 	
	 	#ifdef DEBUG
	 		std::cout<<"skiRocIndex: "<<skiRocIndex<<"   channelIndex: "<<channelIndex<<std::endl;
	 	#endif

		energy *= GeVToMip;
		//additional noise to the energy in MIPs
		energy += randgen->Gaus(energyNoise, energyNoiseResolution);

		MIP_to_HG = 49.3;			//must be fed into from the calibration
		HG_to_LG = 1./9.;			//must be fed into from the calibration
		LG_to_TOT = 1./3.;			//must be fed into from the calibration
		highGainADCSaturation = 2500.;		//must be fed into from the calibration
		lowGainADCSaturation = 2500.;		//must be fed into from the calibration

	 	recHit._energyHigh = energy * MIP_to_HG;
	 	recHit._energyLow = recHit._energyHigh * HG_to_LG;
	 	recHit._energyTot = recHit._energyLow * LG_to_TOT;
	 	
	 	recHit.setEnergy(recHit._energyHigh);


	    if(recHit._energyHigh < highGainADCSaturation){
	      recHit.setFlag(HGCalTBRecHit::kGood);
	    }     
	    else if(recHit._energyLow < lowGainADCSaturation){
	      recHit.setFlag(HGCalTBRecHit::kHighGainSaturated);
	      recHit.setFlag(HGCalTBRecHit::kGood);
	    }
	    else {
	      recHit.setFlag(HGCalTBRecHit::kLowGainSaturated);
	      recHit.setFlag(HGCalTBRecHit::kGood);
	    }

		rechits->push_back(recHit);
}

#include "FWCore/Framework/interface/InputSourceMacros.h"

DEFINE_FWK_INPUT_SOURCE(HGCalTBGenSimSource);
