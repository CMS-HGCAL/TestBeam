#include "FWCore/Framework/interface/Event.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Sources/interface/ProducerSourceFromFiles.h"
#include "DataFormats/FEDRawData/interface/FEDRawDataCollection.h"

#include "FWCore/ParameterSet/interface/ConfigurationDescriptions.h"
#include "FWCore/ParameterSet/interface/ParameterSetDescription.h"

#include <fstream>
#include <iostream>

//#define DEBUG

/**
 * \class HGCalTBTextSource HGCal/RawToDigi/plugins/HGCalTBTextSource.h
 *
 * \brief convert data from txt file to FEDRawData
 *
 * \todo efficiency not tested, many improvements can be done
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

		return true;
	}
	virtual void produce(edm::Event & e);
	bool readLines();

	std::vector<std::string> m_lines;
	std::vector<std::string>::const_iterator _fileName_itr;
	std::ifstream m_file;
	int m_event, m_run;
	int m_sourceId;
	unsigned int  m_time;
};
