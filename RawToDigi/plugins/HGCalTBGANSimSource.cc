#include "HGCal/RawToDigi/plugins/HGCalTBGANSimSource.h"
//source before compilation: source /cvmfs/cms.cern.ch/slc7_amd64_gcc630/external/tensorflow-c/1.1.0/etc/profile.d/init.sh;

using namespace std;


HGCalTBGANSimSource::HGCalTBGANSimSource(const edm::ParameterSet & pset, edm::InputSourceDescription const& desc) :  edm::ProducerSourceFromFiles(pset, desc, true),
	currentEvent(-1) {
	
	//find and fill the configured runs
	RechitOutputCollectionName = pset.getParameter<std::string>("RechitOutputCollectionName"); 
	DWCOutputCollectionName = pset.getParameter<std::string>("DWCOutputCollectionName"); 
	RunDataOutputCollectionName = pset.getParameter<std::string>("RunDataOutputCollectionName"); 
	m_maskNoisyChannels = pset.getUntrackedParameter<bool>("MaskNoisyChannels", true);
  	m_channelsToMask_filename = pset.getUntrackedParameter<std::string>("ChannelsToMaskFileName","HGCal/CondObjects/data/noisyChannels.txt");

  	std::vector<double> v1(4, 1.0);
	wc_resolutions = pset.getUntrackedParameter<std::vector<double> >("wc_resolutions", v1);
	sensorSize = pset.getUntrackedParameter<int> ("sensorSize", 128);

    if (sensorSize==128) {
      	x_max = 7;
    	x_min = -7;
    	y_max = 11;
    	y_min = -11;  
    } else {		//other geometries to be implemented
      	x_max = 7;
    	x_min = -7;
    	y_max = 11;
    	y_min = -11;  	
    }
    range_x = x_max-x_min + 1;
    range_y = ceil((y_max-y_min + 1)/2.);

    NColorsInputImage = pset.getUntrackedParameter<uint>("NColorsInputImage", 17);

	NEvents = pset.getUntrackedParameter<unsigned int> ("NEvents", 20);
	zDim = pset.getUntrackedParameter<unsigned int> ("zDim", 100);
	beamEnergy = pset.getUntrackedParameter<unsigned int> ("beamEnergy", 250);
	beamParticlePDGID = pset.getUntrackedParameter<int> ("beamParticlePDGID", 211);
	//gaussian beam profile, indicated in mm
	beamX_mu = pset.getUntrackedParameter<double>("beamX_mu", 0.0); 
	beamY_mu = pset.getUntrackedParameter<double>("beamY_mu", 0.0); 
	beamX_sigma = pset.getUntrackedParameter<double>("beamX_sigma", 20.0); 
	beamY_sigma = pset.getUntrackedParameter<double>("beamY_sigma", 20.0); 
	
	setupConfiguration = pset.getUntrackedParameter<unsigned int> ("setupConfiguration", 1);

	switch(setupConfiguration) {
		case 1:
			N_layers_HGCal = 6;		//must shift remove layer 0 artificially
			N_layers_BH = 12;
			break;
  		case 2:
			N_layers_HGCal = 17;
			N_layers_BH = 12;
			break;
  		case 3:
			N_layers_HGCal = 10;	
			N_layers_BH = 12;		//fix
			break;
		default:
  		case 4:
			N_layers_HGCal = 10;	
			N_layers_BH = 12;
			break;
	}

	areaSpecification = pset.getUntrackedParameter<std::string>("areaSpecification", "H2");

	if (areaSpecification=="H6A") {
		dwc_zPositions.push_back(-500.);
	} else {
		dwc_zPositions.push_back(-109.);
		dwc_zPositions.push_back(-235.);
		dwc_zPositions.push_back(-1509.);
		dwc_zPositions.push_back(-1769.);
	}


	GANModelIndex = pset.getUntrackedParameter<std::string>("GANModelIndex", "");

	if (GANModelIndex=="1") 
		_enumPhysicsListUsed = HGCAL_TB_SIM_GAN_1;
	else if (GANModelIndex=="2") 
		_enumPhysicsListUsed = HGCAL_TB_SIM_GAN_2;
	else if (GANModelIndex=="3") 
		_enumPhysicsListUsed = HGCAL_TB_SIM_GAN_3;
	else 
		_enumPhysicsListUsed = HGCAL_TB_SIM_GAN;


	produces <HGCalTBRecHitCollection>(RechitOutputCollectionName);
	produces<std::map<int, WireChamberData> >(DWCOutputCollectionName);
	produces<RunData>(RunDataOutputCollectionName);


	_e_mapFile = pset.getUntrackedParameter<std::string>("e_mapFile_CERN");	
	HGCalCondObjectTextIO io(0);
	edm::FileInPath fip(_e_mapFile);
 	
	if (!io.load(fip.fullPath(), essource_.emap_)) {
	  throw cms::Exception("Unable to load electronics map");
	};
	

  randgen = new TRandom();

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


bool HGCalTBGANSimSource::setRunAndEventInfo(edm::EventID& id, edm::TimeValue_t& time, edm::EventAuxiliary::ExperimentType& evType) {	
	if (currentEvent == -1) {		//initial loading of a file
		tf::Shape zShape[] = { 1, zDim }; // 1 = single batch, must ensure that batch size is set to one
		z_tensor = new tf::Tensor(2, zShape);
		tf::Shape eShape[] = { 1, 1 }; // 1 = single batch
		energy_tensor = new tf::Tensor(2, eShape);

  		tf::Shape positionShape[] = { 1, 2 };
  		position_tensor = new tf::Tensor(2, positionShape);

	  	simImage = new tf::Tensor();
	  	std::string fileName = fileNames()[0];
	  	fileName.replace(fileName.find("file:"), std::string("file:").length(), "");
	  	GAN_graph = new tf::Graph(fileName.c_str());
	  	GAN_session = new tf::Session(&(*GAN_graph)); 
			
	  	GAN_session->addInput(z_tensor, "z");
	  	GAN_session->addInput(energy_tensor, "energy_real");
	  	GAN_session->addInput(position_tensor, "impactPoint");
	  	GAN_session->addOutput(simImage, "generator/generator");
		
		currentEvent=0;
	}

	if (currentEvent==NEvents) return false;

	/*
		generate batch here
	*/
	//noise
	std::vector<float> zvalues;
	for (int i=0; i<zDim; i++) zvalues.push_back(randgen->Uniform(-1., 1.));	
	z_tensor->setVector<float>(1, 0, zvalues); // axis 1, batch 0, values

	//energy
	std::vector<float> evalue = {(float)beamEnergy};
	energy_tensor->setVector<float>(1, 0, evalue); // axis 1, batch 0, values
	
	//impact position
	impactX = randgen->Gaus(beamX_mu, beamX_sigma); 
	impactY = randgen->Gaus(beamY_mu, beamY_sigma); 
	std::vector<float> pvalues = {(float)impactX, (float)impactY};
	position_tensor->setVector<float>(1, 0, pvalues); // axis 1, batch 0, values


	GAN_session->run();

	currentEvent++;

	return true;

}


void HGCalTBGANSimSource::produce(edm::Event & event) {	
	//first: fill the rechits
	std::unique_ptr<HGCalTBRecHitCollection> rechits(new HGCalTBRecHitCollection);

	int layer; int u; int v; 
	for(int l=0; l<NColorsInputImage; l++) for(int y_index=0; y_index<(int)range_y; y_index++) for(int x_index=0; x_index<(int)range_x; x_index++) {
		//conversion logic
		layer = l+1;
		//inversion of the coordinate transformation
		int x = x_index;
		v = x+x_min;
		int y = y_index;
		y*=2;
		if (v%2 == 1) y += 1;
		u = y + y_min - v; 
		
		if (HGCalDetectorTopology.iu_iv_valid(layer, 0, 0, u, v, sensorSize)){
			float energy = *(simImage->getPtr<float>(0, u, v, l));
			#ifdef DEBUG
				std::cout<<"Generated: "<<u<<"   "<<v<<"  energy: "<<energy<<std::endl;
			#endif
			makeRecHit(layer, u, v, energy, rechits);
		}
	} 
	event.put(std::move(rechits), RechitOutputCollectionName);

	
	//second: add fake dwcs from xBeam, yBeam
	//dwcs still not part of the GAN
	std::unique_ptr<std::map<int, WireChamberData> > dwcs(new std::map<int, WireChamberData>);	
	for (size_t d=0; d<dwc_zPositions.size(); d++) {
		double dwc_x = impactX + randgen->Gaus(0, wc_resolutions[d]);
		double dwc_y = impactY + randgen->Gaus(0, wc_resolutions[d]);

		WireChamberData* dwc = new WireChamberData(d+1, dwc_x , dwc_y, dwc_zPositions[d]);
		dwc->goodMeasurement_X = dwc->goodMeasurement_Y = dwc->goodMeasurement = true;
		dwc->res_x = dwc->res_y = wc_resolutions[d];
		dwc->averageHitMultiplicty = 1.;

		(*dwcs)[d] = *dwc;
	}
	event.put(std::move(dwcs), DWCOutputCollectionName);		
	


	//third: fill the run data
	std::unique_ptr<RunData> rd(new RunData);
	rd->energy = beamEnergy;		//mean energy of the beam configuration
	rd->configuration = setupConfiguration;
	rd->pdgID = beamParticlePDGID;
	rd->runType = _enumPhysicsListUsed;
	rd->run = 1;
	rd->event = currentEvent;
	rd->booleanUserRecords.add("hasDanger", false);
	rd->doubleUserRecords.add("trueEnergy", beamParticlePDGID);
	rd->booleanUserRecords.add("hasValidDWCMeasurement", true);
	event.put(std::move(rd), RunDataOutputCollectionName);	

}

void HGCalTBGANSimSource::endJob() {
	delete randgen;
	delete z_tensor;
	delete energy_tensor;
	delete simImage;
	delete GAN_graph;
	delete GAN_session;
}

void HGCalTBGANSimSource::makeRecHit(int layer, int u, int v, float energy, std::unique_ptr<HGCalTBRecHitCollection> &rechits) {
	HGCalTBRecHit recHit(HGCalTBDetId(layer, 0, 0, u, v, 0), 0., 0., 0., 0., 0); 
	CellCentreXY = TheCell.GetCellCentreCoordinatesForPlots(layer, 0, 0, u, v, sensorSize );
	double iux = (CellCentreXY.first < 0 ) ? (CellCentreXY.first + HGCAL_TB_GEOMETRY::DELTA) : (CellCentreXY.first - HGCAL_TB_GEOMETRY::DELTA);
	double iuy = (CellCentreXY.second < 0 ) ? (CellCentreXY.second + HGCAL_TB_GEOMETRY::DELTA) : (CellCentreXY.second - HGCAL_TB_GEOMETRY::DELTA);
	recHit.setCellCenterCoordinate(iux, iuy);

	//attention! the electric ID to skiroc mapping takes the cellType information as well which is only 0, 2 in the simulation!
	//any analysis on the simulated data must respect that by computing MIP-ADC factors in the same way and by converting energies into MIP units
	uint32_t EID = essource_.emap_.detId2eid(recHit.id());
	HGCalTBElectronicsId eid(EID);	 
	if( !essource_.emap_.existsEId(eid.rawId()) || std::find(m_noisyChannels.begin(),m_noisyChannels.end(),eid.rawId())!=m_noisyChannels.end() )
		return;

	int skiRocIndex = eid.iskiroc_rawhit();
	skiRocIndex = skiRocIndex % HGCAL_TB_GEOMETRY::N_SKIROC_PER_HEXA;		

	recHit.setTime(-1);

	#ifdef DEBUG
		std::cout<<"skiRocIndex: "<<skiRocIndex<<"   channelIndex: "<<channelIndex<<std::endl;
	#endif


	recHit._energyHigh = -1;
	recHit._energyLow = -1;
	recHit._energyTot = -1;

	recHit.setEnergy(energy);
	recHit.setFlag(HGCalTBRecHit::kGood);

	rechits->push_back(recHit);
		
}

#include "FWCore/Framework/interface/InputSourceMacros.h"

DEFINE_FWK_INPUT_SOURCE(HGCalTBGANSimSource);
