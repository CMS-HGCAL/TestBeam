#include "FWCore/Framework/interface/Event.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Sources/interface/ProducerSourceFromFiles.h"
#include "DataFormats/FEDRawData/interface/FEDRawDataCollection.h"

#include "FWCore/ParameterSet/interface/ConfigurationDescriptions.h"
#include "FWCore/ParameterSet/interface/ParameterSetDescription.h"

#include <fstream>
#include <iostream>


/* now we have only 2 layers */
#define MAXSKIROCS 2

#define _TELESCOPE_FED_ID_ 99

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
		m_run(0),
		_hgcalFiles(pset.getUntrackedParameterSet("hgcalData")),
		_telescopeFiles(pset.getUntrackedParameterSet("telescopeData"))
	{

		produces<FEDRawDataCollection>();

		if(_telescopeFiles.fileNames().size() < 1) {
			//cms::LogWarning("INPUT SOURCE") << "No telescope data";
			std::cerr << "[WARNING] No telescope data" << std::endl;
		}
	}

	virtual ~HGCalTBTextSource()
	{
	}

	static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);

private:

	bool openFile(edm::FromFiles& files, std::ifstream& file); ///< open a new file and update the pointer, it checks if the end of file is reached and increment the file, returns false if the list of files has been completely processed
	bool readLines(); ///< read the hgcal file (two SKIROCS) and return the true if at least one word has been read
	bool readTelescopeLines(); ///< read the telescope file and return the true if at least one word has been read

	virtual bool setRunAndEventInfo(edm::EventID& id, edm::TimeValue_t& time, edm::EventAuxiliary::ExperimentType&);

	virtual void produce(edm::Event & e);

	void parseAddSkiword(std::vector<uint16_t>& skiwords, std::string);
	void parseAddTelescopeWords(std::vector<float>& telescope_words, std::string i);

	std::vector<uint16_t> m_skiwords;
	std::vector<float>    _telescope_words;
	std::ifstream _hgcalFile, _telescopeFile;

	int m_event, m_run, m_event_tmp, m_run_tmp;
	int m_sourceId, m_sourceId_tmp;
	unsigned int  m_time, m_time_tmp;
	unsigned int _triggerID, t_triggerID;

	edm::FromFiles _hgcalFiles;
	edm::FromFiles _telescopeFiles;

};
