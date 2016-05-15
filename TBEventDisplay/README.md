# TBEventDisplay
A simple event display for the HGCal test beam experiments. The display has been tested with CMSSW_8_0_6, slc6_amd64_gcc530, running within a CERNVM virtual machine on a Mac. It should work on lxplus and cmslpc-sl6.

# Installation
```linux
  git clone https://github.com/CMS-HGCAL/TestBeam.git HGCal
  cd HGCal
  cmsenv
  scram b
```
# To Run
Once the cmsenv has been executed, it is possible to run TBEventDisplay.py from any directory as follows 
```linux
  TBEventDisplay.py <geometry-file> <root-file-name>
```
where the geometry file is of the same format as that used in TBStandaloneSimulator and the Root file is a file containing test beam digitized objects (SKIROC data frames).

