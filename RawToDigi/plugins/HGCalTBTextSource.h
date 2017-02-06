#include "FWCore/Framework/interface/Event.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Sources/interface/ProducerSourceFromFiles.h"
#include "DataFormats/FEDRawData/interface/FEDRawDataCollection.h"
#include "HGCal/Geometry/interface/HGCalTBGeometryParameters.h"
#include "HGCal/DataFormats/interface/HGCalTBRunData.h"
#include "HGCal/DataFormats/interface/HGCalTBMultiWireChamberData.h"
#include "FWCore/ParameterSet/interface/ConfigurationDescriptions.h"
#include "FWCore/ParameterSet/interface/ParameterSetDescription.h"
#include <fstream>
#include <iostream>
#include <stdio.h>
#include <string>
#include <sstream>
#include "stdlib.h"
/**
 * \class HGCalTBTextSource HGCal/RawToDigi/plugins/HGCalTBTextSource.h
 *
 * \brief convert data from txt file to FEDRawData
 *
 * \todo replace c-like scanf with c++ versions
 * \todo change run and fed IDs (now are hardcoded)
 */


//to the EDM::Event via auxiliary information


class HGCalTBTextSource : public edm::ProducerSourceFromFiles
{

public:
	explicit HGCalTBTextSource(const edm::ParameterSet & pset, edm::InputSourceDescription const& desc) :  edm::ProducerSourceFromFiles(pset, desc, true),
		currentFileIndex(-1),
		newFileIndex(0),
		inputPathFormat(""),
		MWCInputPathFormat(""),
		m_file(0),
		NSpills(pset.getUntrackedParameter<unsigned int>("nSpills", 6))
	{

		m_sourceId = pset.getUntrackedParameter<int>("fed", 1000); /// \todo check and read from file?
		produces<FEDRawDataCollection>();
		produces<RunData>("RunData");
		produces<MultiWireChambers>("MultiWireChambers");

		if (fileNames()[0] != "file:DUMMY") {
			_fileNames = fileNames();
		}
		//find and fill the configured runs
		runEnergyMapFile = pset.getUntrackedParameter<std::string>("runEnergyMapFile"); 
		inputPathFormat = pset.getUntrackedParameter<std::string>("inputPathFormat");
		MWCInputPathFormat = pset.getUntrackedParameter<std::string>("MWCInputPathFormat");
		
		std::fstream map_file;
		map_file.open(runEnergyMapFile.c_str(), std::fstream::in);
		fillConfiguredRuns(map_file);
		map_file.close();

	}

	virtual ~HGCalTBTextSource()
	{
		if (m_file != 0) fclose(m_file);
		m_file = 0;
	}

	static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);

private:
	int currentFileIndex;	//index of the current file
	int newFileIndex;
	void fillConfiguredRuns(std::fstream& map_file);
	void readMWCDataFromFile(std::string filepath);
	bool setRunAndEventInfo(edm::EventID& id, edm::TimeValue_t& time, edm::EventAuxiliary::ExperimentType&);
	virtual void produce(edm::Event & e);
	virtual void endJob() override;
	bool readHeader(void);
	bool readLines(void);
	runMap configuredRuns;
	std::string runEnergyMapFile;
	std::string inputPathFormat;
	std::string MWCInputPathFormat;
	std::vector<std::string> _fileNames;
	std::vector<std::string> _MWCFileNames;
	int mwcCounter;

	std::array< std::vector < unsigned int> , MAXLAYERS> m_lines;
	FILE* m_file;
	std::vector<MultiWireChambers> EventMultiWireChambers;
	unsigned int m_time;
	unsigned int m_event, m_run, m_spill, max_boards;
	bool _hasDanger;
	unsigned int NSpills;//Read while running how many spills we wish to run over
	int m_sourceId;

	std::map<int, int> eventsPerRun;
};
