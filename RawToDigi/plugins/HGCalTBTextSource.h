#include "FWCore/Framework/interface/Event.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Sources/interface/ProducerSourceFromFiles.h"
#include "DataFormats/FEDRawData/interface/FEDRawDataCollection.h"

#include "FWCore/ParameterSet/interface/ConfigurationDescriptions.h"
#include "FWCore/ParameterSet/interface/ParameterSetDescription.h"

#include <fstream>
#include <iostream>

//#define DEBUG

/* now we have only 2 layers */
#define MAXSKIROCS 2 

/**
 * \class HGCalTBTextSource HGCal/RawToDigi/plugins/HGCalTBTextSource.h
 *
 * \brief convert data from txt file to FEDRawData
 *
 * \todo efficiency not tested, many improvements can be done
 * FED ID is set to the first chipID, this has to be fixed: maybe the FED ID can be read from the file or the filename
 * The number of skiroc chips is hard coded, but should be set by cfg file based on the setup
 */
class HGCalTBTextSource : public edm::ProducerSourceFromFiles
{

public:
	explicit HGCalTBTextSource(const edm::ParameterSet & pset, edm::InputSourceDescription const& desc) :  edm::ProducerSourceFromFiles(pset, desc, true),
		//m_file(0),
		m_run(0)
	{

		produces<FEDRawDataCollection>();
		if (fileNames().size()<1){
			throw cms::Exception("No input files") << "";
		}
		_fileName_itr = fileNames().begin();
		openFile();
	}

	virtual ~HGCalTBTextSource()
	{
	}

	static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);

private:
	void openFile(){
		std::string name = *_fileName_itr;
		if (name.find("file:") == 0) name = name.substr(5);
#ifdef DEBUG
		std::cout << "[DEBUG] Opening file: " << name << std::endl;
#endif
		
		m_file.open(name);
		if (m_file.fail()) {
			throw cms::Exception("FileNotFound") << "Unable to open file " << *_fileName_itr;
		}
	}

	virtual bool setRunAndEventInfo(edm::EventID& id, edm::TimeValue_t& time, edm::EventAuxiliary::ExperimentType&)
	{

		if(m_file.eof()){
			// reached end of file
			m_file.close();
			_fileName_itr++;
			if(_fileName_itr!=fileNames().end()) openFile();
			else return false;
		} 

		if(!readLines()) return false;
		id = edm::EventID(m_run, 1, m_event);

		// time is a hack
		time = (edm::TimeValue_t) m_time;
		#ifdef DEBUGTIME
		std::cout << m_time << "\t" << time << std::endl;
#endif
		return true;
	}
	virtual void produce(edm::Event & e);
	bool readLines();

	void parseAddSkiword(std::vector<uint16_t>& skiwords, std::string);

//	std::vector<std::pair<int, std::vector<std::string> > > m_lines;
	std::vector<uint16_t> m_skiwords;
	std::vector<std::string>::const_iterator _fileName_itr;
	std::ifstream m_file;
	int m_event, m_run, m_event_tmp, m_run_tmp;
	int m_sourceId, m_sourceId_tmp;
	unsigned int  m_time, m_time_tmp;
};
