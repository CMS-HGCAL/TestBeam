# CMS HGCal Testbeam Analysis Framework

=============================================
* This new branch is to be used for analysing data from 2018 test beam. 


### Content

- [Download the code](#download-the-code)
- [Location of data files](#location-of-data-files)
- [Content of the code](#content-of-the-code)
- [New data types](#new-data-types)
- [DWC related plugins](#DWC-related-plugins)
- [Possibly helpful plugins](#possibly-helpful-plugins)
- [Exemplary analyzers](#exemplary-analyzers)
- [Running the code](#running-the-code)
- [More information](#more-information)


## Download the code
* CMSSW Version (!) 9.3.0

```bash
scram project ${RELEASE}
cd ${RELEASE}/src/
cmsenv
git cms-init
git clone git@github.com:CMS-HGCAL/TestBeam.git HGCal/
cd HGCal
git checkout CERN_TestBeam_2018_Ipbus_ntuples
git pull
scram b -j16
```

OUTDATED below here!!!


## Location of data files

#### Raw data files
* The RAW files can be found in **`/afs/cern.ch/work/b/barneyd/public/data_May_2017/disk2_2TB/eudaq_data`**.

#### Electronics map
* The EMAP for this run **`map_CERN_Hexaboard_OneLayers_May2017.txt`** can be found in CondObjects/data

## Content of the code
* The data are read using the EDProducer in **`RawToDigi/plugins/HGCalTBRawDataSource.cc(.h)`**. The data are rearranged to recreate the Skiroc2CMS raw data map (Skiroc2CMS information: [https://cms-docdb.cern.ch/cgi-bin/DocDB/ShowDocument?docid=13121](https://cms-docdb.cern.ch/cgi-bin/DocDB/ShowDocument?docid=13121)). The object storing those data is defined in **`DataFormats/interface/HGCalTBSkiroc2CMS.h`**. This producer is creating a collection of these objects. With the one layer data, the size of the collection should be 4 (4 Skiroc2CMS chips per layer).
* An example of EDAnalyzer, reading the collection of **`HGCalTBSkirocCMS`** is available in **`RawToDigi/plugins/RawDataPlotter.cc`**. It creates 2 .txt files, containing high (resp. low) gain pedestal mean and standard deviation for each channel, for each SCA. It also produces a root file output with basic histograms: 
    - high/low gain ADC distribution for all channels for all SCA
    - high/low gain as a function of the time for all channels
    - plots of mean and standard deviation of ADC counts using hexagonal geometry.
    - there is also an option to set to **`True`** to create event displays (One display per time sample). When the event display option is set, the code is slown down quite a lot, so my advice is to run only on few events.
* Another EDProducer **`RawToDigi/plugins/HGCalTBRawHitProducer.cc(.h)`** transforms the data format (HGCalTBSkiroc2CMS) into a digi format defined in **`DataFormats/interface/HGCalTBRawHit.h`** and stores collection these raw hits (the size of the collection should correspond to the number of channels = 4\*64 ). An option is also available to perform the pedestal subtraction before the raw hit creation. It is mandatory to have the pedestal subtraction here, since it is done for each SCA.
* The EDAnalyzer **`RawToDigi/plugins/RawHitPlotter.cc`** analyzes the collection of HGCalTBRawHit. It creates the same kind of plots as in RawDataPlotter (high/low gain ADC distribution for each time sample, versus time sample, event display). An option is available to run with common mode evaluation and subtraction.
* The EDAnalyzer **`RawToDigi/plugins/PulseShapePlotter`** reads the collection of HGCalTBRawHit, perform the common mode subtraction. It creates for each event, for each channel two histograms (high/low gain) of ADC counts as a function of the time. When signals are high enough, one should be able to see pulse shapes.
* To run the 2 previous analyzers (RawHitPlotter and PulseShapePlotter) with the common mode option, the pedestal subtraction option in the raw hit producer should be set.
* The rechit producer in **`Reco/plugins/HGCalTBRecHitProducer.cc(.h)`** transforms the raw hit into the HGCalTBRecHit format and save collection of rechits. It assumes the pedestal subtraction is already done and it performs the common mode subtraction (one common mode per chip). A first version of this producer use only the time sample 3 for the energy.
* The EDAnalyzer **`Reco/plugins/RecHitPlotter.cc`** reads the rechit collection and again produces basic plots.
## Running the code
* **`run2017_cfg.py`** is the script to run to start the analysis. This script is less complex (and less convenient) than the **`test_cfg.py`** in other branches. The chain sequence tool is absent, and will be re-implemented later. 
* All the EDProducer and EDAnalyzer (described above) are loaded with default parameters. Some paths might need modification depending on the user.
* 3 different cms.Path are pre-filled (2 are actually commented). One of them allows to run the RawDataPlotter plugin analysis, the 2nd runs the HGCalTBRawHitProducer and the RawHitPlotter analyzer and the last one runs the PulseShapePlotter analyzer. These 3 cms.Path creates 2 root files:
    - **`cmsswEvents.root`** which contains the cmssw objects (HGCalTBSkiroc2CMS,HGCalTBRawHit). Running the following command should return the list of cmssw collections in the events:```edmDumpEventContent cmsswEvents.root```
    - **`HexaOutput_${runnumber}.root`** contains the analysis output histograms

## New data types

### RunData
* Defined in `DataFormats/interface/HGCalTBRunData.h`. 
* An object that characterises a given event in terms of event numbers, trigger counters - independent of any event cuts if these variables are set at the unpacking stage -, run numbers, setup configurations ("1"=July 2017 TB - "4"=October 2017 with 20 layers), energy plus particle ID of the incident particles, and run types (pedestal, beam, sim.). 

```c++
    int run; //run number to which the event belongs. The edm's run object is always set to 1
    int trigger; //could be assigned to the trigger counter that is stored in the HGCal timing files
    int event; //is assigned to the actual event number starting at 1 for each new run
    int configuration; //1: July TB with 6 layers, 2: September TB with 17 layers, 4: October TB with 20 layers
    double energy; //beam energy
    int pdgID;	//particle ID of the beam particle
    RUNTYPES runType; //possible choices: HGCAL_TB_PEDESTAL, HGCAL_TB_BEAM, HGCAL_TB_SIM
```

* Arbitrary double and boolean values are assigned as so-called `UserRecords`. As an example, the time difference [ms] between a synchronised event to the matched event from the HGCal data stream could be stored as such a userrecord.

* The updated `HGCalTBRawDataSource` plugin is equipped with necessary configurable options to correctly set these quantities at the ORM file unpacking stage and to add them to the edm format.

### WireChamberData
* Defined in `DataFormats/interface/HGCalTBWireChamberData.h`. 
* A simple structure that contains reconstructed DWC impact positions (x,y,z). 
* Furthermore, lower level TDC data such as the number of received time stamps per DWC (<=4) as well as the averaged time stamp per DWC is stored.
* `WireChamberData` are instantiated in the `HGCalTBWireChamberSource` and `HGCalTBBeamWireChamberProducer` plugins.

### HGCalTBDWCTrack
* Defined in `DataFormats/interface/HGCalTBDWCTrack.h`. 
* An object that contains parameters from a straight-line fit using as many DWC measurements as possible in an event.
* Besides the four track parameters (slope and offset for each coordinate) and the quality information in terms of chi2's, the `int referenceType` indicates which DWCs have contributed to the track:

 * `referenceType += 1` for DWC E (CERN H2 TB area), or for the only DWC in H6A
 * `referenceType += 2` for DWC D (CERN H2 TB area)
 * `referenceType += 4` for DWC A (CERN H2 TB area)
 * `referenceType += 8` for DWC ext. (CERN H2 TB area)

* Under the provision of layer-distance files (see `CondObjects/data/layer_distances_CERN_September_7EELayers_10FHLayers_V0.txt`, values are given in mm, beginning EE part is at z=0mm.), the straight line trajectories are evaluated at each layer of the HGCal in the `DWCTrackProducer` plugin. Its extrapolations can be retrieved through `std::pair<double, double> DWCExtrapolation_XY(int layerIndex)` afterwards. Note that the layerIndex starts at "1".

## DWC related plugins

### HGCalTBTimingFileWriter
* Defined in `RawToDigi/plugins/HGCalTBTimingFileWriter.cc`. 
* A plugin that produces timing files similar to the ones created directly in the July 2017 beam test. These files are necessary to synchronise the HGCal data streams offline with e.g. the TDC and AHCal data streams.
* It runs on unpacked ORM files and prints the timing information of the skiroc indexed "0".

### HGCalTBWireChamberSource
* Defined in `RawToDigi/plugins/HGCalTBWireChamberSource.*`. This plugin runs on an arbitrary number of raw TDC files in which the hits from 4x4 TDC channels are stored in TTrees. The mapping to the DWCs is hard-coded. Similarly, the positions of the DWCs along the beam line (z-coordinate) are implemented with hard numbers both for H2 and H6.

* The position reconstruction is based on the time stamp difference and the multiplication by 0.2mm/ns. The TDC time stamp interval to ms is configurable in `TDCTriggerTimeStampConversionToMs`. 
* For the second half of the July2017 TB, all hits are stored. Consequently, more complex reconstruction algorithms making use of the entire collection of received hits (instead of just the first one) could be designed and implemented in a next iteration of this source plugin.

* Furthermore, running this plugin yields one printout per event to ensure the successful synchronisation with the HGCal data stream. Here, the timing file from the `HGCalTBTimingFileWriter` plugin must be provided for each run. Both event offsets in the HGCal and/or the TDC data stream at the beginning of each file are configurable and so is the time difference tolerance which sets a threshold for a successful and insufficient synchronisation of events. As the TDC in H6 has sometimes received random triggers, those are detected and skipped if the `allowForTDCEventSkipping` option is set true.

* Alignment parameters of the four DWCs as computed by using the millepede formalism (x/y translations and x-y rotations only) can be read and be applied.

* **In general, this plugin has been extensively tested and run on the data taken in the 2017 TB campaigns at CERN. The reconstructed data do not show any sign neither of a buggy reconstruction nor of a false synchronisation. Yet, its implementation is not general and its applicability to setups with other TDCs and/or more than 4 DWCs must be revised.**

### DWC_NTupelizer
* Defined in `Reco/plugins/DWC_NTupelizer`. 
* A simple plugin that writes out reconstructed impact positions on the DWCs into ROOT trees. 
* Residuals from simple tracking, position differences, and other information are written, too, if the `writeMinimal` option is set to false.

### HGCalTBBeamWireChamberProducer
* Defined in `RawToDigi/plugins/HGCalTBBeamWireChamberProducer`. 
* Reads from the reconstructed DWC tree (`DWC_NTupelizer`) and adds a `WireChamberData` collection to the edm file. 
* The matching of HGCal events to the DWC events makes use of the `RunData->event` parameter as this one is designed to be indepent to any event cuts that might be applied in between this merging of data streams and the initial synchronisation (i.e. comparison of time stamps) which is done right after the unpacking step.

* *It is suggested to run the merging of data streams which is done by this producer after the RecHit production.*

### DWCTrackProducer
* Defined in `Reco/plugins/DWCTrackProducer`. 
* `WireChamberData` are input to this producer. It computes `HGCalTBDWCTrack` objects. 
* Layer distance files must be provided to match beam line (=z-) positions to the individual layers where the beginning of the EE part is defined to be z=0mm.

### Suggested workflow for DWC reco, sync. and merging
 * **1.** Produce the timing files with `HGCalTBTimingFileWriter`. Determine event offsets.
 * **2.** Run the `DWC_NTupelizer` (`writeMinimal=True` suffices) with the `HGCalTBWireChamberSource` as the corresponding source to obtain reconstructed DWC files.
 * **3.** Run the `HGCalTBBeamWireChamberProducer` on edm files containing HGCal data (`skiroc2cms`, `rawhits` or `rechits`) and `RunData` - as its int event member is read for the event matching - under the provision of reconstructed DWC files. 
* **4.** The `DWCTrackProducer` computes tracks given the `WireChamberData` collection.

## Possibly helpful plugins

### HGCalTBGenSimSource
* Defined in `RawToDigi/plugins/HGCalTBGenSimSource`. 
* A source that reads from the simulated data files (=ROOT TTrees) and produces `HGCalTBRecHits`, `RunData` and `WireChamberData`.
* The hit IDs are converted into layer and cell numbers. The cell numbers can be uniquely mapped to their respective x,y positions on the sensor. This in turn yields the (u,v) tuple with which electronic IDs are instantiated. Hence, a simulated cell can be mapped to its physical counter part using skiroc and channel IDs. In particular, noisy channels may be removed using the same code and files as in the rawhit production.
* GeV to MIP conversion factors are configurable and must be retrieved from studies in simulation.
* Any smearing of the energy is not applied by default but configurable parameters to apply a gaussian smearing exist.
* DWC measurements are faked using the original beam positions as stored in the ROOT TTree.

* Ultimately, the simulated files are converted into files which may be input to other analyzer plugins that are designed to run on real reconstructed TB data. This should facilitate data to MC comparisons.

### NumpyConverter
* Defined in `Reco/plugins/NumpyConverter`. 
* A plugin that dumps the edm content into a specific numpy format. 
* Here, the data are stored as tensors in `.npz` files which may be loaded subsequently using standard `python`. 
* As an example, reconstructed energies from rechits are stored as NEvents x NLayers x 15 x 15 tensors. 
* The rechit energy of the cell with (u,v) = (-3, 4) in layer 4 and event 198 would be retrieved via 

```python
import numpy as np
data = np.load("<pathToDataFile>.npz") #just a pointer to the file, does not consume any RAM
rechits = data["rechits"] #actually loads the rechit data into the ram, ~300MB for 10k events with 17 modules
rechitOfInterest = rechits[198-1][4-1][-3+7][4+7]
```

* **(More detailed examples on how to retrieve the data are in preparation.)**

## Exemplary analyzers
### MIPFinder
* Defined in `Reco/plugins/MIPFinder`. 
* It is the analyzer with which the histograms for the MIP calibration with and without cuts on the extrapolated DWC tracks are made.
* It makes use of `RunData`, `WireChamberData`, `HGCalTBDWCTracks` and `HGCalTBRecHits` such that it should serve as a good reference for other persons intending to write own analyzers on any of these objects.

### Configuration files
* Step-by-step configration files are provided in the `runConfigurations2017`directory. 
* For each file, not more than two processes are chained to each other to keep the output minimal and maintain modularity of the targeted workflow.

## More Information
For general instructions on the framework, co-ordinate system, reconstruction sequence, etc please refer to
[http://cms-hgcal.github.io/TestBeam/](http://cms-hgcal.github.io/TestBeam/)
