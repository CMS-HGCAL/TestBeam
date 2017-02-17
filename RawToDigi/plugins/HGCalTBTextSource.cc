#include "HGCal/RawToDigi/plugins/HGCalTBTextSource.h"

using namespace std;
//#define DEBUG
/**
RUN=000880	SPILL=01	EVENT=000000	GLOBALTIME=0x01B6EAF2	BOARD=09	RUN=000880	SPILL=01	EVENT=000000	GLOBALTIME=0x01B6EAF2	BOARD=12	RUN=000880	SPILL=01	EVENT=000000	GLOBALTIME=0x01B6EAF2	BOARD=14	RUN=000880	SPILL=01	EVENT=000000	GLOBALTIME=0x01B6EAF2	BOARD=17	RUN=000880	SPILL=01	EVENT=000000	GLOBALTIME=0x01B6EAF2	BOARD=18	RUN=000880	SPILL=01	EVENT=000000	GLOBALTIME=0x01B6EAF2	BOARD=19	RUN=000880	SPILL=01	EVENT=000000	GLOBALTIME=0x01B6EAF2	BOARD=22	RUN=000880	SPILL=01	EVENT=000000	GLOBALTIME=0x01B6EAF2	BOARD=23
T 0x0072D911 0xBEC20B15	T 0x0072D93F 0xBED19F54	T 0x0072D96E 0xBED19F56	T 0x0072D99D 0x3A11F91F	T 0x0072D9CC 0x40DACEB1	T 0x0072D9FA 0xC5EAE85C	T 0x00000000 0xBED19F52	T 0x00000000 0x57B86176
0x118F11BD 0x1081108D	0x11A511B7 0x118A108C	0x119311F2 0x119F1195	0x109E109E 0x1083108D	0x118D11F2 0x11FD1198	0x11881195 0x118D1196	0x11911190 0x1193118A	0x1088109F 0x118C118A
*/


//it is crucial that the format of the configuration file is respected!
void HGCalTBTextSource::fillConfiguredRuns(std::fstream& map_file) {
	std::string run_prefix;
	std::string filePath;

	//perform the loop and fill configuredRuns
	char fragment[100];
	int readCounter = 0;
	int _run = 0, _mwcrun = 0, _configuration = 0; double _energy = 0; std::string _runType = ""; 

	while (map_file.is_open() && !map_file.eof()) {
		readCounter++;
		map_file >> fragment;
		if (readCounter <= 5) continue; 	//skip the header
		else if (readCounter % 5 == 1) {
			if (((std::string)fragment).find("//") == std::string::npos)	//skip comments of form //
				_run = atoi(fragment); 
			else
				readCounter = 1;	//artificially mimics the line as a header
		}
		else if (readCounter % 5 == 2) _mwcrun = atoi(fragment);
		else if (readCounter % 5 == 3) _energy = atof(fragment); 
		else if (readCounter % 5 == 4) _runType = (std::string)fragment; 
		else if (readCounter % 5 == 0) {
			_configuration = atoi(fragment); 
			
			//if readOnlyRuns parameter is set, make sure to only run the analysis for files of exactly that energy
			if (readOnlyRuns.size()>0 && std::find(readOnlyRuns.begin(), readOnlyRuns.end(), _run) == readOnlyRuns.end())
				continue;

			//store
			configuredRuns[_run].energy = _energy;
			configuredRuns[_run].runType = _runType;
			configuredRuns[_run].configuration = _configuration;

			//add the zeros
			if (_run < 10) run_prefix= "00000";			
			else if (_run < 100) run_prefix= "0000";			
			else if (_run < 1000) run_prefix= "000";			
			else if (_run < 10000) run_prefix= "00";			
			else if (_run < 100000) run_prefix= "0";			
			else run_prefix = "0";
			
			
			filePath = inputPathFormat;		
			filePath.replace(filePath.find("<RUN>"), 5, run_prefix+std::to_string(_run));
			std::cout<<"Adding DAQ file "<<filePath<<std::endl;
			_fileNames.push_back(filePath);
			
			filePath = MWCInputPathFormat;		
			filePath.replace(filePath.find("<RUN>"), 5, std::to_string(_mwcrun));
			filePath.replace(filePath.find("file:"), 5, "");
			std::cout<<"Adding MWC file "<<filePath<<std::endl;
			_MWCFileNames.push_back(filePath);
		}
	}
}

void HGCalTBTextSource::readMWCDataFromFile(std::string filepath) {
	EventMultiWireChambers.clear();
	mwcCounter = 0;
	std::fstream mwc_file;
	 mwc_file.open(filepath.c_str(), std::fstream::in);	
	char fragment[100];
	while (mwc_file.is_open()) {
		MultiWireChambers _mwcs;
		mwc_file >> fragment;
		double x1 = atof(fragment);
		mwc_file >> fragment;
		double x2 = atof(fragment);
		mwc_file >> fragment;
		double y1 = atof(fragment);
		mwc_file >> fragment;
		double y2 = atof(fragment);

		if (mwc_file.eof())
			break;

		_mwcs.push_back(MultiWireChamberData(1, x1, y1, -126.-147.));
		_mwcs.push_back(MultiWireChamberData(2, x2, y2, -147.));
		EventMultiWireChambers.push_back(_mwcs);
		mwcCounter++;
	}

	//std::cout<<mwcCounter<<" MWC entries..."<<std::endl;
	mwcCounter = 0;
}

bool HGCalTBTextSource::setRunAndEventInfo(edm::EventID& id, edm::TimeValue_t& time, edm::EventAuxiliary::ExperimentType& evType)
{	
	if (!(_fileNames.size())) return false; // need a file...
	if (newFileIndex==(int)_fileNames.size()) return false;	//do not overshoot in the file name array

	//magic must come here!
	if (currentFileIndex != newFileIndex) {
		currentFileIndex = newFileIndex;
		std::string name = _fileNames[currentFileIndex].c_str(); 

		if (currentFileIndex < (int)_MWCFileNames.size()) {
			std::string MWCname = _MWCFileNames[currentFileIndex]; 
			readMWCDataFromFile(MWCname);
		}

		if (name.find("file:") == 0) name = name.substr(5);
		m_file = fopen(name.c_str(), "r");
		if (m_file == 0) {
			throw cms::Exception("FileNotFound") << "Unable to open file " << name;
		}
		m_event = 0;
	}
	if (feof(m_file)) {
		newFileIndex++;
		return setRunAndEventInfo(id, time, evType);
	}
	// read the header
	while(readHeader() == false) {
		if (feof(m_file)) {
			newFileIndex++;
			return setRunAndEventInfo(id, time, evType);
		}
		// have to skip the data lines and return false
		readLines();
	}

	// now try to read from the file
	if (!readLines()) return false;

	id = edm::EventID(m_run, m_spill, m_event);

	time = (edm::TimeValue_t) m_time;
	if(m_spill <= NSpills) return true;
	else return false;
}

// sets m_run, m_spill, m_even, m_time
///\todo read the number of boards
bool HGCalTBTextSource::readHeader()
{
	char buffer[1024]; // typically 605 chars
	buffer[0] = 0;     //init buffer
	fgets(buffer, 1000, m_file); // read the header line
	std::string b = buffer;
#ifdef DEBUG
	std::cout << "------------------------------\n";
	std::cout << b << std::endl;
#endif
	if(b.find("DANGER=true") != std::string::npos) _hasDanger=true; //return false;
	else _hasDanger=false;
//RUN=000880	SPILL=01	EVENT=000000	GLOBALTIME=0x01B6EAF2	BOARD=09	RUN=000880	SPILL=01	EVENT=000000	GLOBALTIME=0x01B6EAF2	BOARD=12	RUN=000880	SPILL=01	EVENT=000000	GLOBALTIME=0x01B6EAF2	BOARD=14	RUN=000880	SPILL=01	EVENT=000000	GLOBALTIME=0x01B6EAF2	BOARD=17	RUN=000880	SPILL=01	EVENT=000000	GLOBALTIME=0x01B6EAF2	BOARD=18	RUN=000880	SPILL=01	EVENT=000000	GLOBALTIME=0x01B6EAF2	BOARD=19	RUN=000880	SPILL=01	EVENT=000000	GLOBALTIME=0x01B6EAF2	BOARD=22	RUN=000880	SPILL=01	EVENT=000000	GLOBALTIME=0x01B6EAF2	BOARD=23

	if(sscanf(buffer, "RUN=%u\tSPILL=%u\tEVENT=%u\tGLOBALTIME=%x", &m_run, &m_spill, &m_event, &m_time) != 4) {
		return false; // what to do?
	}

	// here count the number of boards for this event

	return true;
}


// each time it returns true, then the event is complete and can be saved in EDM format
bool HGCalTBTextSource::readLines()
{
	max_boards = 0;
	for( auto& board : m_lines) {
		board.clear();
	}

	char buff[1024], buff_SK0[1024], buff_SK1[1024];
	buff[0] = 0;
	unsigned int data_sk0, data_sk1;

	// loop over all the lines
	for(unsigned int i = 0; i < 69 && !feof(m_file); ++i) {
		buff[0] = 0;
		fgets(buff, 1000, m_file);
		// loop over one line of the text file
		std::string b = buff;
		std::istringstream buffer(b);
		unsigned int board_counter = 0;
		if((i < 1) || (i > 64)) continue;// Only these are data words, the rest dont interest us for now
		while( buffer.peek() != EOF && board_counter < MAXLAYERS) { //buffer.good() gives compilation errors
			// read the data of the two skirocs of one board
			buffer >> buff_SK1;
			buffer >> buff_SK0;
			data_sk0 = strtoul(buff_SK0,NULL,0);
			data_sk1 = strtoul(buff_SK1,NULL,0);
			//continue for all the boards
			// extra security
			m_lines[board_counter].push_back(data_sk0);
			m_lines[board_counter].push_back(data_sk1);
			
			++board_counter;
			if(board_counter > max_boards) max_boards = board_counter;
		}
	}
//		if(sscanf(buffer, "Event header for event %x with (200ns) timestamp %x", &tmp_event, &m_time) == 2) {
	return !m_lines.empty();
}

void HGCalTBTextSource::produce(edm::Event & event){
	eventCounter++;	//indexes each event chronologically passing this plugin

	//add the multi-wire chambers only if available
	bool _hasValidMWCMeasurement = true;
	std::auto_ptr<MultiWireChambers> mwcs(new MultiWireChambers);	
	if (mwcCounter < (int)EventMultiWireChambers.size()) {
		for (size_t _imwc=0; _imwc < (size_t)EventMultiWireChambers[mwcCounter].size(); _imwc++) {
			_hasValidMWCMeasurement = (EventMultiWireChambers[mwcCounter][_imwc].x != -999) && _hasValidMWCMeasurement;
			_hasValidMWCMeasurement = (EventMultiWireChambers[mwcCounter][_imwc].y != -999) && _hasValidMWCMeasurement;
			double rotAngle = mwcRotation*M_PI/180.;
			double x_preRot = EventMultiWireChambers[mwcCounter][_imwc].x/10.;		//conersion from mm to cm 
			double y_preRot = EventMultiWireChambers[mwcCounter][_imwc].y/10.; 
			EventMultiWireChambers[mwcCounter][_imwc].x = cos(rotAngle) * x_preRot + sin(rotAngle) * y_preRot; 	//apply the rotation just here because the filtering for -999 must occur first
			EventMultiWireChambers[mwcCounter][_imwc].y = -sin(rotAngle) * x_preRot + cos(rotAngle) * y_preRot; 
			
			if (_imwc==1) {//i.e. the second MWC
				EventMultiWireChambers[mwcCounter][_imwc].x += mwc2DeltaX;
				EventMultiWireChambers[mwcCounter][_imwc].y += mwc2DeltaY;
			}

			mwcs->push_back(EventMultiWireChambers[mwcCounter][_imwc]);
		}
		mwcCounter++;
	} else {
		_hasValidMWCMeasurement = false;
		//push some dummy value for the MWCs, subsequent plugins using this information should check the _hadValidMWCMeasurement flag
		mwcs->push_back(MultiWireChamberData(1, -999, -999, 0));
	}
	event.put(std::move(mwcs), "MultiWireChambers");		


	//add the DAQ data
	std::auto_ptr<FEDRawDataCollection> bare_product(new  FEDRawDataCollection());
	// here we parse the data
	std::vector<uint16_t> skiwords;
	// make sure there are an even number of 32-bit-words (a round number of 64 bit words...

	for (unsigned int i_board = 0 ; i_board < max_boards; ++i_board) {
		auto board = m_lines[i_board];
		for (auto skiword : board) {
                        skiwords.push_back(skiword >> 16);
			skiwords.push_back(skiword);
		}
	}

	FEDRawData& fed = bare_product->FEDData(m_sourceId);
	size_t len = sizeof(uint16_t) * skiwords.size();
	fed.resize(len);
	memcpy(fed.data(), &(skiwords[0]), len);
	event.put(bare_product);



	//add the RunData
	std::auto_ptr<RunData> rd(new RunData);
	if (configuredRuns.find(m_run) != configuredRuns.end()) {
		rd->energy = configuredRuns[m_run].energy;
		rd->configuration = configuredRuns[m_run].configuration;
		rd->runType = configuredRuns[m_run].runType;
		rd->run = m_run;
		rd->event = eventCounter;
	} else {
		rd->energy = -1;
		rd->configuration = -1;
		rd->runType = "-1";
		rd->run = -1;
		rd->event = eventCounter;
	}

	if (eventsPerRun.find(m_run) == eventsPerRun.end()) {
		eventsPerRun[m_run].first = 0;
		eventsPerRun[m_run].second = EventMultiWireChambers.size();
	}
	eventsPerRun[m_run].first++;
	rd->hasDanger = _hasDanger;
	rd->hasValidMWCMeasurement = _hasValidMWCMeasurement;

	event.put(std::move(rd), "RunData");	
}



void HGCalTBTextSource::fillDescriptions(edm::ConfigurationDescriptions& descriptions)
{
	edm::ParameterSetDescription desc;
	desc.setComment("TEST");
	desc.addUntracked<int>("run", 101);
	desc.addUntracked<std::vector<std::string> >("fileNames");
	desc.addUntracked<std::string>("inputPathFormat");
	desc.addUntracked<std::string>("MWCInputPathFormat");
	desc.addUntracked<double>("mwcRotation");
	desc.addUntracked<double>("mwc2DeltaX");
	desc.addUntracked<double>("mwc2DeltaY");
	desc.addUntracked<std::vector<int> >("readOnlyRuns");
	desc.addUntracked<std::string>("runEnergyMapFile");
	desc.addUntracked<unsigned int>("nSpills", 6);
	descriptions.add("source", desc);
}

void HGCalTBTextSource::endJob() {
	std::cout<<"Run 	DAQ events 		MWC events"<<std::endl;
	for (std::map<int, std::pair<int, int> >::iterator it = eventsPerRun.begin(); it != eventsPerRun.end(); it++) {
		std::cout<<it->first<<": "<<it->second.first<<" , "<<it->second.second<<std::endl;
	}
}

#include "FWCore/Framework/interface/InputSourceMacros.h"

DEFINE_FWK_INPUT_SOURCE(HGCalTBTextSource);
