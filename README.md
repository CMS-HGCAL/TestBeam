CMS HGCal Testbeam Analysis Framework
=============================================

### Content

- [Download the code](#download-the-code)
- [Location of data files](#location-of-data-files)
- [Running the code](#running-the-code)
- [positionResolution branch](#positionresolution-branch)
    - [TextInput Plugin](#textinput-plugin) 
    - [GenSim to RecHit Converter](#gensim-to-rechit-converter)
    - [HGCalTBRecHitProducer](#hgcaltbrechitproducer)
    - [PositionResolutionAnalyzer](#positionresolutionanalyzer)
    - [MillepedeBinaryWriter](#millepedebinarywriter)
    - [Helper Classes](#helper-classes) 
    - [Steering files](#steering-files) 
- [More information](#more-information)




## Download the code
* CMSMSW Version 8.0.1

```bash
>> scram project ${RELEASE}
>> cd ${RELEASE}/src/
>> cmsenv
>> git cms-init
>> git clone git@github.com:CMS-HGCAL/TestBeam.git HGCal/
>> git checkout FNAL_TB_4Layers
>> git pull
>> scram b -j16
```

## Location of data files

#### Raw data files
The RAW files can be found in **`/afs/cern.ch/user/r/rchatter/public/FNAL_TB_May/FNAL_TB_May_4Module_Runs`**.
This is structured by beam type and for the Electrons further by the energy of the electron beam.

#### Electronics map and pedestal files
* The EMAP for this run **`map_FNAL_Layer1234.txt`** can be found in CondObjects/data
The event average pedestal files for both low and high gain can be found in
**`Ped_LowGain_M1254.txt`** and **`Ped_HighGain_M1254.txt`** also in **`CondObjects/data`**.
* The information of the EMAP is needed by the framework when producing the digis from the 
RAW and is provided in the file **`RawToDigi/python/hgcaltbdigis_cfi.py`**.
It is also needed in analysis code if a conversion from Detector Id to Electronics Id is needed.
* May refer to the example **`RawToDigi/plugins/DigiPlotter.cc`** . In the **`EndJob()`** function of this same
code the format needed in which the pedestal files are to be written out can be found.
* The Pedestal information is needed in the step from Digis to Reco. And it is provided from
**`Reco/python/hgcaltbrechitproducer_cfi.py`**.

## Running the code
* In **`test_cfg.py`** the RAW data file to be run is to be provides in the line:

```python
process.source = 
cms.Source("HGCalTBTextSource",
  run=cms.untracked.int32(1),
  fileNames=cms.untracked.vstring("file:FNAL_TB_May_4Module_Runs/PED/HGC_Output_1600.txt")
)
```

* In **`cms.Path`**: 
**`process.hgcaltbdigis`** produces digis and **`process.hgcaltbrechits`** produces the rechits. The analyzer to be run can come after this. There are a few provided:
**`RawToDigi/plugins/DigiPlotter.cc`** to generate the pedestals
**`Reco/plugins/RecHitPlotter_HighGain_Correlation_CM.cc`** for the noise studies starting from pedestal subtracted digis.

* The output root file name is set in:

```python
process.TFileService = cms.Service("TFileService", fileName = cms.string("HGC_Output_1600_Reco.root") )
```

#### Event-based reconstruction (Outdated?)
```bash
sh scripts/rearrangeTxtFile.sh input.txt
>> CERN raw data rearranged in /eos/cms/store/group/upgrade/HGCAL/TestBeam/CERN/Sept2016/
cmsRun test_cfg_newEB.py chainSequence=6 pedestalsHighGain=CondObjects/data/Ped_HighGain_L8.txt pedestalsLowGain=CondObjects/data/Ped_LowGain_L8.txt runNumber=930 configuration=2 runType=HGCRun
>> chainSequence=<configuration>
>> test_cfg_newEB.py to load the reqarranged .txt:
-> mount eos
-> options.register('dataFolder','~/eos/cms/store/group/upgrade/HGCAL/TestBeam/CERN/Sept2016/',....
```

#### Overview of analysis configurations
* **`chainSequence=1`** Runs on Digis produces pedestal files in the path specifienfied under pedestalsHighGain and pedestalsLowGain
* **`chainSequence=4`** Runs event Display analyzer
* **`chainSequence=5`** Runs on Recos for each cell of the detector across layers. Use for pion, muon and pedestal run. Correlation across cells as implemented by Kai-Yu are evaluated.     
* **`chainSequence=6`** to run **`Layer_Sum_Analyzer`**: Evaluates for each layer as well as summed across layers: Max hit, Max hit + 6 nearest neighbours(7 cells), 19 cells and All cells --- **with each cell picked subject to a 2 MIP cut**. Results are available considering only the energy deposited in Silicon, and with dE/dx weights to recover the total energy deposited in the absorbers and the silicon sensors.

* **`configuration=-1`** ADCtoMIP for CERN (0 = ADCtoMIP for FNAL) 
* **`configuration=1`** (2) to select weights for 5X0 (25X0) 8 layers cern runs

## PositionResolution branch

### TextInput Plugin
The `HGCalTextSource` plugin has been modified to cope with the analysis of multiple runs / data files per execution. Except for the removed veto on `Danger=true` events, downward compatibility is granted. Application of the `BadSpillFilter` plugin should restore the compatibility.
The usage of the new features is optional. Besides the usual data readout, this plugin also adds a `RunData` collection:

```c++ 
produces<RunData>("RunData");
``` 
to the event, such that the information on beam energies, run numbers, ... is available for all subsequent plugins in the chain on an event-basis. Furthermore, MWC data are read and put to the event as `HGCalTBMultiWireChamberData`: 

```c++ 
produces<MultiWireChambers>("MultiWireChambers");
```
Hereby, it is assumed that the reconstructed coordinates are given in [mm].

 __**Options:**__

* **`fileNames`**: Array of file paths to be processed in that execution. Before, the plugin could only deal with one file at a time. If the paths to files are indicated, `runData` objects obtain dummy values.
* **`nSpills`**: Number of spills to be processed before the read-out is stopped. The spill number is read from the .txt files such that the according variable is reset for each file.
* **`runEnergyMapFile`**: Path to the configuration file that defines the mapping between the run numbers, H4DAQ run numbers, the appropriate beam momenta (labeled 'energy'), the run type (e.g. electron, muon, ...) and the setup configuration (e.g. 1, 2). An exemplary mapping file is added to `CondObjects/Data/runEnergiesPositionResolution.txt`. Missing H4 runs, i.e. such that do not have 6000 events, are marked by `xxx`. If  ```python  fileNames=cms.untracked.vstring(["file:DUMMY"])``` is set, the plugin reads the runs that are defined in this file. Otherwise, the **`fileNames`** has preference.
* **`inputPathFormat`**: Generic format of the path to the run data files. The setting of this parameter is obligatory if the read-out with the **`runEnergyMapFile`** is aimed. The placeholder for the run number is `<RUN>`. Up to five zeros are added to the run string, e.g. run=1 corresponds to 000001.
In summary, an exemplary setting could be: 

```python  
inputPathFormat=cms.untracked.string("file:%s/HGCRun_Output_<RUN>.txt"%options.dataFolder)
``` 

* **`MWCInputPathFormat`**: Similar to **`inputPathFormat`**. Zeros are not added to the run number string. Example: 

```python  
inputPathFormat=cms.untracked.string("file:%s/HGCRun_Output_<RUN>.txt"%options.dataFolder)
```

* **`readOnlyRuns`**: Array of run numbers to choose only a subset of runs for the processing. Those runs must be defined in the `runEnergyMapFile`.
* **`mwcRotation`**: Initial rotation of the multi-wire chamber coordinate system in degrees. The original coordinate information from the H4-files are rotated counter-clockwise along the z-axis after the readout. Set to 270 degrees in the position resolution analysis.
* **`mwc2DeltaX`**: Translation of the second multi-wire chamber with respect to the first one in x-direction AFTER the initial rotation. Set to zero in the position resolution analysis.
* **`mwc2DeltaY`**: Translation of the second multi-wire chamber with respect to the first one in y-direction AFTER the initial rotation. Set to zero in the position resolution analysis.



### GenSim to RecHit Converter
This plugin reads the simulated simulated data formats and converts the information into **`HGCalTBRecHit`** collections. Multi-wire chamber data are faked using the xBeam and yBeam branch of the input trees. Their distance w.r.t. the first module is hard-coded.
(x,y) position information of the cells are computed from the cell IDs as they are placed in the Geant4 simulation. So far, only the sensors with `133 cells` are implemented. Energies are converted from GeV to MIPs to ADC using the identical conversion factors in data.
Smearing effects in the fake MWC measurements and the noise have appropriate parameters.

 __**Options:**__
 
 * **`OutputCollectionName`**: Name of the HGCalTBRecHit collection that is put to the event. 
 * **`fileNames`**: See description of TextInput plugin.
 * **`e_mapFile_CERN`**: Electronics mapping file to perform the matching of the cell IDs to the hypothetical skiRocs on which MIP to ADC conversions depend.
 * **`runEnergyMapFile`**: See description of TextInput plugin. Since MWCs are faked and do not have their own files, according mapping files contain only four (instead of five) columns. An example is shown in `CondObjects/Data/runEnergiesPositionResolutionSimulation.txt`.
 
 * **`inputPathFormat`**: See description of the TextInput plugin. Paths to simulated files are assumed to have a placeholder for the simulated energy and a file index (=`run`). Therefore, an example is: 

```python 
inputPathFormat=cms.untracked.string("file:%s/<ENERGY>GeV/TBGenSim_<RUN>.root"%(options.dataFolder))
```
 
 * **`energyNoise`**: Energy noise that is added to the energy per cell. This parameter represents the mean value of a Gaussian in MIP unit from which the additional summands are diced. The parameter is set to zero in the position resolution analysis.
 * **`energyNoiseResolution`**: Energy noise that is added to the energy per cell. This parameter represents the RMS of a Gaussian in MIP unit from which the additional summands are diced. The parameter is set to zero in the position resolution analysis. Thus, no energy smearing is applied.
 * **`MWCSmearingResolution`**: Resolution of the multi-wire chamber fake measurements. The coordinates are smeared by a Gaussian centered at zero with width `MWCSmearingResolution` (in microns).


### HGCalTBRecHitProducer
* HGCalTBRecHits have the cell center coordinates as additional members `cellCenter_x`, `cellCenter_y`.
* The center's coordinates are set in this producer (lines 99 ff.) or in the `HGCalTBGenSimSource` plugin to grant a common interface for the position retrieval in the next steps.


### PositionResolutionAnalyzer
Plugin that reconstructs impact positions of the incident particles for each layer, extra- or interpolates them and computes the differences of both approaches. Both the reconstruction and the inter-/extrapolation schemes are configurable. Output of this plugin is a ROOT TTree (name: `deviations`) with information on the layerIndex, the impacts positions in x&y, the deviations to the reference, ... (see lines 351 ff in `Position_Resolution_Analyzer.cc`).
Z-positioning of the layers along the beam line in cm and radiation length is hard-coded. 

In principle, the procedure can be segmented into the following:

1. Loop through the `recHits`. Instantiate `sensor` objects, fill them with `recHits`. Alignment parameters are set if available.
2. Subtract the common mode noise for each layer (=`sensor`).
3. Reconstruction of the impact position for each layer.
4. Compute electron tracks using either the multi-wire chambers (extrapolation) or the impact positions from the stack ignoring the layer under investigation (interpolation). Implemented track models are straight lines (preferred, chosen model) and general broken lines (deprecated, since electron-induced shower axis' scattering angle cannot be computed trivially after certain raditiation lengths).
5. Calculation of deviations and writing to the output tree once per layer and event.

 __**Framework Options:**__

* **`RUNDATA`**: Input tag to the RunData collection. Name is fixed in the read-out plugins. E.g. 

```python 
RUNDATA = cms.InputTag("source","RunData","unpack" )
```

* **`MWCHAMBERS`**: Input tag to the multi-wire chambers collection. Name is fixed in the read-out plugins. E.g. 

```python 
MWCHAMBERS = cms.InputTag("source","MultiWireChambers","unpack" )
``` 

* **`HGCALTBRECHITS`**: Input tag to the recHits collection. E.g. 

```python 
HGCALTBRECHITS = cms.InputTag("hgcaltbrechits","","unpack" )
```

* **`HGCALTBCLUSTERS`**: Input tag to the cluster collection. E.g. 

```python 
HGCALTBCLUSTERS = cms.InputTag("hgcaltbclusters","","unpack" )
```

* **`HGCALTBCLUSTERS7`**: Input tag to the cluster (with seven closest cells) collection. E.g. 

```python 
HGCALTBCLUSTERS7 = cms.InputTag("hgcaltbclusters","7","unpack" )
```

* **`HGCALTBCLUSTERS19`**: Input tag to the cluster (with 19 closest cells) collection. E.g. 

```python 
HGCALTBCLUSTERS19 = cms.InputTag("hgcaltbclusters","19","unpack" )
```

* **`e_mapFile_CERN`**: Electronics mapping file to perform the matching of the cell IDs to the hypothetical skiRocs on which MIP to ADC conversions depend.

 __**Experimental Setup Options:**__

* **`layers_config`**: Index of the layer configuration (1: CERN <15 X0, 2: CERN <25X0). `Default: 1`.
* **`nLayers`**: Number of sensors/layers in the setup. `Default: 8`.
* **`SensorSize`**: Size of the sensors in terms of number of cells. `Supported: 133`.
* **`ADC_per_MIP`**: MIP to ADC conversion factors for each of the 2x8 skirocs. Example: 

```python 
ADC_per_MIP = cms.vdouble([17.31, 17.12, 16.37, 17.45, 17.31, 16.98, 16.45, 16.19, 17.55, 17.19, 16.99, 17.92, 15.95, 16.64, 16.79, 15.66])
```

* **`pedestalThreshold`**: Threshold in MIP unit to specify the boundary under which cell energy deposits are included in the common mode noise calculation and subsequent subtraction. `Default: 2`.
* **`alignmentParameterFiles`**: Array of paths to alignment files. The content of such files must follow the format as dictated by the `millepede` program. The parameter files are matched to each run. For this purpose, the path to the corresponding alignment parameter file must contain the substring `RUN_<theNumber>_`. Otherwise, i.e. if not available, values are assumed to be zero (=no alignment). The functionality of the readout is implemented in `Reco/src/PositionResolutionHelpers.cc`.

 __**Reconstruction Options:**__
 
 * **`considerationMethod`**: Specifies the set of cells that are considered in the impact position reconstruction for each layer. `all`, `closest7` (=one ring around the cell with highest energy deposit), `closest19`(=two rings) as well as the cells from the clustering algorithm (`clustersAll`, `clusters7`, `clusters19`) can be specified. `Default: closest19`.
 * **`weightingMethod`**: Specifies the weighting formula for the impact position reconstruction. Besides a logarithmic weightig scheme, linear and squared energy schemes are implemented. `Default: logWeighting_3.5_1.0`.
 
 __**Reference Options:**__
 
 * **`useMWCReference`**: Boolean flag indicating if multi-wire chamber extrapolation of HGCal stack interpolation is used as reference for the impact position.
 * **`fittingMethod`**: Corresponds to the track model of the electron shower's trajectory. Implemented are a straight line computed through an analytic chi2 minimisation (`lineAnalytical`), polynomial fits up to degree three using ROOT's minimisation (e.g. `pol3TGraphErrors`) and General Broken lines using the GBL core functions that are available in cmssw (`gblTrack`). The position resolution analysis makes exclusive use of the `lineAnalytical`. 
 * **`fitPointWeightingMethod`**: Rescaling scheme of errors on the reconstructed impact positions prior to the reference computation. This way, inaccurate positions, e.g. through uniform energy distribution across all cells in one layers, have less influence on the track fit. This procedure has been found to have no significant impact on the results. Therefore, the default configuration is `none`.


 

### MillepedeBinaryWriter
This plugin is very similar to the PositionResolutionAnalyzer as it also reconstructs impact positions and compares those references.
Yet, the output is not a voluminous ROOT TTree but a binary file that is subsequently input to the alignment procedure using `Millepede`.
Thus, the main difference is the track model derivation which is either done using the two MWCs in front (extrapolation) or all eight HGC stack layers only (for interpolation). Besides the residuals, their local derivatives w.r.t. the alignment parameters are computed and written to the binary file. [Link to the millepde manual.](http://www.desy.de/~blobel/Mptwo.pdf)

__**Options:**__

* **Same as for the PositionResolutionAnalyzer** 

__**Additional options: Quality cuts:**__

* **`totalEnergyThreshold`**: Energy threshold in units of MIPs for an event to be added to the binary file. If the energy sum across all cells in all layers after common mode subtraction is below this threshold, the entire event is discarded. `Default: 1000`.
* **`MWCQualityCut`**: Applies an event selection to the MWC measurement if it is set to `True`. The exact threshold of this cut is hard-coded for this analysis: 

```c++
double delta_x12 = mwcs->at(0).x - mwcs->at(1).x;
double delta_y12 = mwcs->at(0).y - mwcs->at(1).y;
//Todo: Hard coded numbers! Offsets correspond to angles
if (fabs(delta_x12-0.45) > 1.0 || fabs(delta_y12-0.076) > 1.0) return; 
```




### Helper Classes
#### HGCalTBRunData
An object that is created for each event in the data read-out with runEnergyMapFiles and passed to the central event objects to subsequent plugins. 

```c++
  int configuration; //index of the configuration
  int run; //number of the run
  int event; //event counter, counting occurs during the data readout prior to any cuts (e.g. before the bad-spill filtering)
  double energy; //more precisely: the set beam momentum
  std::string runType; //electrons, muons, pions
  bool hasValidMWCMeasurement; //an MWC measurement is valid if no coordinate shows -999
  bool hasDanger; //events that show Danger=true are not filtered but are flagged through this boolean. Therefore, the event filtering may occur later if necessary.
```
Note that certain variables, namely `hasValidMWCMeasurement`, `hasDanger`and `event` are properly set, even if runEnergyMap files are not indicated.

#### HGCalTBMultiWireChamberData
Simple structure to pass the read multi-wire chamber data to the event.

```c++
  explicit MultiWireChamberData(int _ID, double _x, double _y, double _z): ID(_ID), x(_x), y(_y), z(_z) {}; //default consructor

  int ID; //index of the MWC in the setup (e.g. 1, 2)
  double x; //reconstructed x-coordinate
  double y; //reconstructed y-coordinate
  double z; //position of the MWC in the setup, hard-coded values
};

typedef std::vector<MultiWireChamberData> MultiWireChambers; //two MWCs per event --> create a vector
```

#### Geometry
* **`HGCalWaferGeometry.cc`**: Conversion of cell IDs in simulation to x,y positions and cell types. Hard-coded mapping depending on the configuration in simulation.

* **`HGCalTBCellVertices::GetCellIUIVCoordinates`**: A function that returns the iu, iv coordinates of a cell on the 133 cell sensor with input of its x, y coordinates. It represents a useful function to retrieve iu, iv as those are input to the constructor of `HGCalTBDetId` which are required to instantiate `HGCalTBRecHits` like it is done in the GenSim to recHit conversion.


#### Sensors
Central class to perform all relevant steps on the sensors/layers for impact position reconstruction.

```c++
    SensorHitMap(int _label);
    //constructor, must indicate a label which could be the layer number (1-8)
    //most internal parameters are set to zero
    
    ~SensorHitMap();
    //default destructor
    
    /*************
    setters
    *************/     
    
    void setLabZ(double z_cm, double X0);
    //z-coordinate in cm and radiation lengths
    
    void setParticleEnergy(double e);
    //Sets the particle energy at the layer. For instance, one could approximate te electron energy using the PDG loss formula (implemented by `computeEnergyLoss` in Reco/interface/Tracks.h) . 
    //The energy is input to the kink angle RMW calculation for gbl tracks which are deprecated in the position resolution analysis.
    
    void setAlignmentParameters(double d_alpha, double d_beta, double d_gamma, 
    double d_x0, double d_y0, double d_z0);
    //Assigns up to six alignment parameters (three infinitesimal rotations and three translations) to each layer.
    
    void setResidualResolution(double r);
    //Sets the residual width to the layer. This information is required within the alignment binary writer.
    
    void setSensorSize(int s);
    //Sets the sensor size (default: 133). Used to boundary checks in position computations.
    
    void setPedestalThreshold(double t);
    //Sets the pedestal threshold under which all cells are input for the common mode noise subtraction.
    
    void setCenterHitPosition(double x, double y, double x_err, double y_err);
    //Manually fake the reconstructed impact position, e.g. if multi wire chambers are represented as instances of the sensor class.
    
    /*************
    computation
    *************/    
    void addHit(HGCalTBRecHit Rechit, double ADC_per_MIP);
    //Register a HGCalTBRecHit to the sensor. For the conversion of the energy to MIPs, the skiroc dependent ADC to MIP conversion factor has to be passed as an argument. Only cells of type '0' are incorporated.
    
    void registerClusterHit(HGCalTBDetId hit, int N_considered);
    //Registers a hit to be part of a cluster. 
    //The passed hit must have already been added via addHit(...) in order to be properly registered - otherwise it is ignored.
    //N_considered=-1 refers to the 'all' scheme, N_considered=7 for one ring and N_considered=19 for two rings around the cell with highest energy deposit.
    
    std::pair<int, double> subtractCM();  
    //Subtracts the common mode of each cell. Negative values are set to zero.
    //Returns the common mode noise and the number of contributing cells.
    
    void calculateCenterPosition(ConsiderationMethod considerationMethod, 
            WeightingMethod weightingMethod);
    //Applies the main impact position reconstruction for the layer.
    //Consideration- and WeightingMethod are dedicated enum data types.
    
    /*************
    getters
    *************/
    int label();
    //returns the initially set label
    
    double getTotalEnergy();
    //Returns the total energy summed over all cells after common mode subtraction in units of MIPs
    
    double getTotalClusterEnergy(int N_considered);
    //Returns the total energy summed over all cells which are entries in clusters before common mode subtraction in units of MIPs
    
    double getTotalWeight();
    //Returns the sum of all weights, equal to the denominator value in the weighting formula
    
    double getLabZ();
    //Returns the initially set z-coordinate in cm
    
    double getX0();
    //Returns the initially set z-coordinate in radiation lengths
    
    double getParticleEnergy();
    //Returns the initially energy at the layer in GeV
    
    double getIntrinsicHitZPosition();
    //Returns the z-coordinate of the main impact position in the sensor's reference frame. 
    //If no alignment is applied, it returns 0.
    
    double getResidualResolution();
    //Returns the initially set residual resolution or -1 if it has not been set.
    
    std::pair<double, double> getHitPosition(); 
    //Returns central main impact position in layer's own frame
    
    std::pair<double, double> getLabHitPosition();  
    //Returns central main impact position in the lab frame
    
    std::pair<double, double> getHitPositionError(); 
    //Returns the error of the reconstruction, calculated as the RMS of the weighting formula
    
    std::pair<double, double> getCenterOfClosestCell(std::pair<double, double> X_ref);
    //Returns the coordinates of the cell's center with minimal distance to X_ref. 
    //Function considers only those hits, that have been added to the sensor at this point. 
    //Has to be used with caution for running on simulation where not every cell on each sensor is available in the input data.
```


#### Tracks
Central class wrapping the functionality to make electron shower trajectory models in the position resolution analysis.

```c++
    ParticleTrack();
    //default constructor
    
    ~ParticleTrack();
    //default destructor
    
    void addFitPoint(SensorHitMap* sensor);
    //Adds a point to the model prior to any fit. 
    //Those points are the main impact position members of the sensor class. 
    //Thus, on the passed sensors, the reconstruction must be run in advance.
    //MWC measurements are added setting impact positions manually.
    
    void addReferenceSensor(SensorHitMap* reference);
     //Defines the reference sensor with which the trajectory model prediction is compared to the reconstruction.
    
    void weightFitPoints(FitPointWeightingMethod method);
    //Performs the energy-dependent rescaling of fit point errors.
    
    void fitTrack(TrackFittingMethod method);
    //Fits the trajectory. TrackFittingMethod is an enum type.
    
    std::pair<double, double> calculateReferenceXY();
    //Calculates the reference coordinates using the trajectory model on the reference sensor.
    
    std::pair<double, double> calculateReferenceErrorXY();
    //Calculates the reference coordinates uncertainties using the trajectory model on the reference sensor.
    
    std::pair<double, double> calculatePositionXY(double z, int layerLabel);
    //Calculates the reference coordinates using the trajectory model indicating any position on the z-axis (analytic models for the trajectory),
    //or through indication of a sensor label (e.g. for the general broken line trajectory).
    
    std::pair<double, double> calculatePositionErrorXY(double z, int layerLabel);
    //Calculates the reference coordinate uncertainties using the trajectory model indicating any position on the z-axis (analytic models for the trajectory),
    //or through indication of a sensor label (e.g. for the general broken line trajectory).
    
    double getSumOfEnergies();
    //Returns the sum over all energies of the fit point's total energies.
    
    /*****
    GBL specific members:
    *****/
    
    void addDummySensor(SensorHitMap* dummy);
     //Adds a sensor to the GBL formalism. 
     //Measurement is not specified but the absorbers/thin scatterers are placed.
     
    void gblTrackToMilleBinary(gbl::MilleBinary *milleBinary);
    //Practical, external interface to write out the necessary derivatives for alignment with millepede into the binary form.
```

### Steering files

**`test_cgf_newEB.py`**: Minimally adjusted old steering file to grant full compatibility with the existing code. Chains such as the LayerSumAnalyzer can be run as before. Only the arguments for the `HGCalTBTExtSource` plugin had to be extended.

**`position_resolution_cgf.py`**: Alternative steering file to run either the position resolution analyzer (producing the TTree) and/or the binary file creator for the alignment using Millepede.

__**Reading Options:**__

* **`repoFolder`**: Directory to the repository to correctly navigate to some file in CondObjects.

* **`dataFolder`**: Main directory containing the raw text input.

* **`pathToRunEnergyFile`**: (Absolute) path to the file that indicates the runs to run the analysis on.

* **`readOnlyRuns`**: Run the analysis only for the indicated runs.

* **`isData`**: Is the analysis run on real data (otherwise on simulated samples)?

* **`nSpills`**: Number of spills per run before read-out is stopped.

__**Pedestal Subtraction Options:**__

* **`pedestalsHighGain`**: Path to high gain pedestals file.

* **`pedestalsLowGain`**: Path to low gain pedestals file.

__**Analysis Execution Options:**__

* **`chainSequence`**: Here: only 8 (position resolution) and 9 (writing of millepede binary) is configured.

* **`reportEvery`**: Frequency of event count printouts on the console.

__**Experimental Setup Options:**__

* **`configuration`**: 1 if 8Layers with 5X0 sampling the center of the shower only; 2 if 8Layers with 25X0 sampling up to the tail of the shower.

* **`alignmentParameterFiles`**: Alignment parameter files as obtained from the (mille)pede framework.

__**Reference Trajectory Options:**__

* **`useMWCReference`**: Are the multi wire chamber information used for the exrapolation of referemce impact points?

* **`fittingMethod`**: Model of the electron tracks.

__**Reconstruction Options:**__

* **`considerationMethod`**: Possible arguments are: `all`, `closest7`, `closest19`, `clusters`.

* **`weightingMethod`**: Possible arguments are: `squaredWeighting`, `linearWeighting`, `logWeighting_3.5_1.0`, ... .

* **`pedestalThreshold`**: Threshold for the common mode noise subtraction. Unit is MIP.

__**Output Options:**__

* **`outputFolder`**: Main directory of the output files.

* **`outputPostfix`**: Postfix to the output file.

## More Information
For general instructions on the framework, co-ordinate system, reconstruction sequence, etc please refer to
[http://cms-hgcal.github.io/TestBeam/](http://cms-hgcal.github.io/TestBeam/)
