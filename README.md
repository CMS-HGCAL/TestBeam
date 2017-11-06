CMS HGCal Testbeam Analysis Framework
=============================================
* This new branch is to be used for analysing data from 2017 test beam. 
* The code runs with the one layer data, and will need (major) updates when we will have more layers. 
* This readme assumes one layer data.
* The README.md in other branches contain useful information about 2016 test beam analysis which might be still useful.

### Content

- [Download the code](#download-the-code)
- [Location of data files](#location-of-data-files)
- [Content of the code](#content-of-the-code)
- [Running the code](#running-the-code)
- [More information](#more-information)


## Download the code
* CMSSW Version 8.0.1

```bash
scram project ${RELEASE}
cd ${RELEASE}/src/
cmsenv
git cms-init
git clone git@github.com:CMS-HGCAL/TestBeam.git HGCal/
git checkout CERN_TestBeam_2017
git pull
scram b -j16
```

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
    
## More Information
For general instructions on the framework, co-ordinate system, reconstruction sequence, etc please refer to
[http://cms-hgcal.github.io/TestBeam/](http://cms-hgcal.github.io/TestBeam/)
