// -------------------------------------------------------------------------
// Description: convert output of standalone simulator to TB 2016 test beam
// data format: SKIROC2DataFrames
// Created  April 2016 Harrison B. Prosper
// Updated: 04-23-2016 HBP use HGCCellMap to get mapping from (u, v) to
//                     (skiroc, channel id)
//          05-03-2016 HBP add noise code
// -------------------------------------------------------------------------
#include <iostream>
#include <algorithm>
#include "FWCore/Framework/interface/InputSourceMacros.h"
#include "HGCal/TBStandaloneSimulator/plugins/HGCSimDigiSource.h"
#include "HGCal/DataFormats/interface/HGCalTBDetId.h"
// -------------------------------------------------------------------------
using namespace std;

namespace
{
const size_t debugCount = 1;
};

HGCSimDigiSource::HGCSimDigiSource
(const edm::ParameterSet& pset, edm::InputSourceDescription const& desc)
	: edm::ProducerSourceFromFiles(pset, desc, true),
	  _run(pset.getUntrackedParameter<int>("runNumber", 101)),
	  _maxevents(pset.getUntrackedParameter<int>("maxEvents", -1)),
	  _minadccount(pset.getUntrackedParameter<int>("minADCCount", 1)),
	  _adcpermev(pset.getUntrackedParameter<double>("ADCperMeV", 200)),
	  _filenames(pset.getUntrackedParameter<vector<string> >
	             ("fileNames")),
	  _noisefilenames(pset.getUntrackedParameter<vector<string> >
	                  ("noiseFileNames")),
	  _chain(0),              // chain of files
	  _tree(0),               // sim tree
	  _entries(0),            // number of simulated events
	  _entry(0),              // entry number,
	  _cellmap(HGCCellMap()), // cell id to (u, v) map and (x, y) to (u, v)
	  _simhits(0),

	  _noisechain(0),         // chain of files
	  _noisetree(0),          // sim tree
	  _noiseentries(0),       // number of simulated events
	  _noiseentry(0),         // entry number
	  _random(TRandom3()),
	  _noise(vector<map<uint32_t, uint16_t> >())
{
	// collections to be saved
	produces<SKIROC2DigiCollection>();
	produces<HGCSSSimHitVec>();

	// ------------------------------------------------------
	// create a chain of files
	// ------------------------------------------------------
	_chain = new TChain("HGCSSTree");
	if ( !_chain )
		throw cms::Exception("ChainCreationFailed", "chain: HGCSSTree");

	for(size_t c = 0; c < _filenames.size(); c++)
		_chain->Add(_filenames[c].c_str());

	// determine number of events to read
	_entries = _chain->GetEntries();
	_entries = _maxevents < 0
	           ? _entries
	           : (_entries < (size_t)_maxevents ? _entries : _maxevents);
	cout << endl
	     << "==> Number of simulated events to read: "
	     << _entries
	     << endl;

	// map input tree to sim object pointers
	_tree = (TTree*)_chain;
	_tree->SetBranchAddress("HGCSSSimHitVec", &_simhits);

	// ------------------------------------------------------
	// read in noise mode
	// ------------------------------------------------------
	if ( _noisefilenames.size() > (size_t)0 ) {
		_noisechain = new TChain("Pedestals");
		if ( !_noisechain )
			throw cms::Exception("ChainCreationFailed",
			                     "chain: Pedestals");

		for(size_t c = 0; c < _noisefilenames.size(); c++)
			_noisechain->Add(_noisefilenames[c].c_str());

		// determine number of noise events to read
		_noiseentries = _noisechain->GetEntries();

		// map input tree to object pointers
		vector<uint32_t>* vkey = 0;
		vector<uint16_t>* vadc = 0;
		_noisetree = (TTree*)_noisechain;
		_noisetree->SetBranchAddress("vkey", &vkey);
		_noisetree->SetBranchAddress("vadc", &vadc);

		cout << "==> Creating noise model from file with "
		     << _noiseentries
		     << " pedestal events"
		     << endl;

		int N = 0;
		double a1 = 0.0;
		double a2 = 0.0;
		for(size_t entry = 0; entry < _noiseentries; entry++) {
			long localentry = _noisechain->LoadTree(entry);
			_noisetree->GetEntry(localentry);
			_noise.push_back(map<uint32_t, uint16_t>());
			if ( entry % 500 == 0 )
				cout << entry << "\t" << vkey->size() << endl;

			for(size_t c = 0; c < vkey->size(); c++) {
				long key  = (*vkey)[c];
				uint16_t adc = (*vadc)[c];
				_noise[entry][key] = adc;
				a1 += adc;
				a2 += adc * adc;
				N++;
			}
		}
		a1 /= N;
		a2 /= N;
		a2 = sqrt(a2 - a1 * a1);
		char rec[256];
		sprintf(rec, "==> pedestal = %8.1f +/-%-8.1f", a1, a2);
		cout << rec << endl;
		cout << "==> Done loading noise model"
		     << endl;
	} else {
		cout << "==> No noise will be added"
		     << _noiseentries
		     << endl;
	}
}

HGCSimDigiSource::~HGCSimDigiSource()
{
	if ( _chain ) delete _chain;
	if ( _noisechain ) delete _noisechain;
}

bool HGCSimDigiSource::
setRunAndEventInfo(edm::EventID& id,
                   edm::TimeValue_t& time,
                   edm::EventAuxiliary::ExperimentType&)
{
	if ( !_chain )
		throw cms::Exception("ChainNotFound")
		        << "sim file chain not open";

	if ( _entry >= _entries ) return false;

	// load sim objects into memory
	long localentry = _chain->LoadTree(_entry);
	_tree->GetEntry(localentry);
	_entry++;

	// construct event info
	id = edm::EventID(_run, 1, _entry);

	// FIXME: time hack
	edm::TimeValue_t present_time = presentTime();
	unsigned long time_between_events = timeBetweenEvents();
	time = present_time + time_between_events;

	return true;
}

void HGCSimDigiSource::produce(edm::Event& event)
{
	if ( _entry < debugCount ) {
		cout << endl
		     << "==> entry number: " << _entry
		     << endl;
	}

	// auto_ptr own objects they point to and are
	// automatically deleted when out of scope

	// add sim hits to event
	std::auto_ptr<HGCSSSimHitVec> simhits(new HGCSSSimHitVec());
	for(size_t c = 0; c < _simhits->size(); c++)
		simhits->push_back((*_simhits)[c]);
	event.put(simhits);

	// create skiroc digi objects and put in event
	std::auto_ptr<SKIROC2DigiCollection>
	digis(new SKIROC2DigiCollection(SKIROC::MAXSAMPLES));

	// sum energies of sim hits in each cell
	// and convert cell energies from MeV to ADC count
	vector<HGCSimDigiSource::Cell> channels;
	digitize(channels);

	// store digitized data
	if ( _entry < debugCount ) {
		cout << "    number of digi hits: " << channels.size() << endl;
	}

	for(size_t c = 0; c < channels.size(); c++) {
		HGCSimDigiSource::Cell& cell = channels[c];
		digis->addDataFrame(cell.detid);
		digis->backDataFrame().setSample(0,
		                                 cell.ADClow,
		                                 cell.ADChigh,
		                                 cell.TDC);
		if ( _entry < debugCount ) {
			char record[80];
			sprintf(record,
			        "%6d (ski,chan)=(%2d,%2d)"
			        " (u,v)=(%2d,%2d), (x,y)=(%6.2f,%6.2f), ADC=(%d)",
			        (int)(c + 1),
			        cell.skiroc, cell.channel,
			        cell.u, cell.v, cell.x, cell.y, cell.ADChigh);
			cout << record << endl;
			assert(cell.channel > -1);
		}
	}
	event.put(digis);
}

void HGCSimDigiSource::digitize(std::vector<HGCSimDigiSource::Cell>& channels)
{
	// Steps:
	// 1. for each cell, sum sim hit energies
	// 2. convert energies to adc counts
	// 3. randomly select a pedestal event and add to current event
	// 4. apply zero suppression
	// 5. sort cells

	if ( _entry < debugCount ) {
		cout << endl
		     << "    number of sim hits:  "
		     << _simhits->size() << endl;
	}

	// 1. for each cell, sum sim hit energies
	map<uint32_t, HGCSimDigiSource::Cell> hits;

	for(size_t c = 0; c < _simhits->size(); c++) {
		HGCSSSimHit& simhit = (*_simhits)[c];
		int cellid = simhit.cellid();
		int layer  = simhit.silayer() + 1;
		// FIXME: hard code for now
		int sensor_u = 0;
		int sensor_v = 0;
		pair<int, int> uv = _cellmap(cellid);
		int u = uv.first;
		int v = uv.second;
		assert(u > -1000);

		int celltype = _cellmap.celltype(layer,
		                                 sensor_u, sensor_v,
		                                 u, v);

		HGCalTBDetId detid(layer, sensor_u, sensor_v, u, v, celltype);
		uint32_t key = detid.rawId();

		if ( hits.find(key) == hits.end() ) {
			// this is a new cell, so initialize its data structure
			HGCSimDigiSource::Cell cell;
			cell.ADClow = 0;
			cell.ADChigh = 0;
			cell.TDC    = 0;
			cell.energy = 0.0;
			cell.layer  = layer;
			cell.sensor_u = sensor_u;
			cell.sensor_v = sensor_v;
			cell.u = u;
			cell.v = v;

			pair<double, double> xy = _cellmap.uv2xy(cell.u,
			                          cell.v);
			cell.x = xy.first;
			cell.y = xy.second;

			pair<int, int> eid = _cellmap.uv2eid(layer,
			                                     sensor_u,
			                                     sensor_v,
			                                     u, v);
			cell.skiroc   = eid.first;
			cell.channel  = eid.second;
			cell.celltype = celltype;
			cell.detid    = detid;

			hits[key] = cell;
		}

		// sum energy per cell
		hits[key].energy += simhit.energy();
	}

	// 2. convert cell energy to ADC counts
	for(map<uint32_t, HGCSimDigiSource::Cell>::iterator it = hits.begin();
	        it != hits.end(); it++) {
		uint32_t key  = it->first;
		double energy = hits[key].energy;
		uint16_t adc  = static_cast<uint16_t>(_adcpermev * energy);
		hits[key].ADChigh = adc;
	}

	// 3. add noise
	// randomly select a pedestal event and add to current
	// event, cell by cell
	if ( _noisechain ) {
		int index = _random.Integer(_noiseentries - 1);
		map<uint32_t, uint16_t>& noise = _noise[index];
		for(map<uint32_t, uint16_t>::iterator it = noise.begin();
		        it != noise.end(); it++) {
			uint32_t key = it->first; // this is the rawId
			uint16_t adc = it->second;
			if ( hits.find(key) != hits.end() )
				hits[key].ADChigh += adc;
			else {
				// this is a new cell with no signal count, only noise
				HGCalTBDetId detid(key);
				HGCSimDigiSource::Cell cell;
				cell.ADClow = 0;
				cell.ADChigh = adc;
				cell.TDC    = 0;
				cell.energy = 0.0;
				cell.layer    = detid.layer();
				cell.sensor_u = detid.sensorIU();
				cell.sensor_v = detid.sensorIV();
				cell.u        = detid.iu();
				cell.v        = detid.iv();
				cell.celltype = detid.cellType();
				cell.detid    = detid;

				pair<int, int> eid = _cellmap.uv2eid(cell.layer,
				                                     cell.sensor_u,
				                                     cell.sensor_v,
				                                     cell.u,
				                                     cell.v);
				cell.skiroc   = eid.first;
				cell.channel  = eid.second;

				pair<double, double> xy = _cellmap.uv2xy(cell.u, cell.v);
				cell.x = xy.first;
				cell.y = xy.second;

				hits[key] = cell;
			}
		}
	}

	// 4. apply zero suppression
	for(map<uint32_t, HGCSimDigiSource::Cell>::iterator it = hits.begin();
	        it != hits.end(); it++) {
		long key = it->first;
		HGCSimDigiSource::Cell& cell = hits[key];
		// apply "zero" suppression
		if ( cell.ADChigh < _minadccount ) continue;
		channels.push_back(cell);
	}

	// 5. sort so that SKIROC 2 comes before SKIROC 1
	// and channels increase monotonically
	sort(channels.begin(), channels.end());
}

void HGCSimDigiSource::fillDescriptions
(edm::ConfigurationDescriptions& descriptions)
{
	edm::ParameterSetDescription desc;
	desc.setComment("Test Beam 2016");
	desc.addUntracked<int>("runNumber", 101);
	desc.addUntracked<int>("maxEvents", -1);
	desc.addUntracked<int>("minADCCount", 1);
	desc.addUntracked<double>("ADCperMeV", 200);
	desc.addUntracked<std::vector<std::string> >("fileNames");
	desc.addUntracked<std::vector<std::string> >("noiseFileNames");
	descriptions.add("source", desc);
}

DEFINE_FWK_INPUT_SOURCE(HGCSimDigiSource);
