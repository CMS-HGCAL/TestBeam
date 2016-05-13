# TBStandaloneSimulator
A simple standalone simulator (a revamped version of the code developed
by Anne-Marie Magnan), whicn is intended for quick approximate studies. 
For sophisticated studies use the CMSSW test beam simulator. 

The code has been tested with CMSSW_8_0_6, 
slc6_amd64_gcc530, running within a CERNVM virtual machine on a mac. 
It should work on lxplus and cmslpc-sl6.

# Installation
```linux
  cd CMSSW_8_0_6/src
  git clone https://github.com/CMS-HGCAL/TestBeam.git HGCal
  cd HGCal/TBStandaloneSimulator
  cmsenv
  scram b
```
# Running the simulator
The simulator is called simulateTB. To simulate a couple of 32 GeV electron 
events, do
```linux
simulateTB geometry_1layer.py withvis.mac
```
which creates the files
```linux
PFcal.root
g4_00.wrl
g4_01.wrl
```
The Root file contains the results of the simulation, while the second and third
files contain graphical data that can be rendered using a VRML browser, 
such as freewrl.

To simulate 1000 32 GeV electrons events (without visualization) for a 1-layer
detector (like that investigated in the March 2016 test beam) do
```linux
simulateTB geometry_1layer.py gun.mac
```
which creates the file
```linux
PFcal.root
```

# Running producer (to create SKIROC2DataFrames)
```linux
  cd test
  cmsRun produceSKIROCCollection_cfg.py
```
This reads the sim file 
```linux
PFcal.root 
```
and creates and saves SKIROC dataframe objects to edm::Events. You 
should see the output file
```linux
HGC_Electrons_32GeV_2016_04_sim.root 
```
which can be analyzed in the same way as real test beam data.

If digitized events with noise are required, first create a noise file
from a pedestal run as in the following example
```linux
writePedestal.py HGC_Pedestals_2016_04_8272.root
```
This will create the files
```linux
HGC_Pedestals_2016_04_8272_Noise.root
HGC_Pedestals_2016_04_8272.txt
```
The first is a root file containing the noise to be added to the simulated 
digitized hits and the second is a text file containing the pedestals. In order
to activate the noise model during digitization, create the file noise_filelist,
which should contain the file name of the root file containing the noise data,
```linux
echo HGC_Pedestals_2016_04_8272_Noise.root > noise_filelist
```
then do
```linux
  cmsRun produceSKIROCCollection_cfg.py
```


