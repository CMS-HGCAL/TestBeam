# TBEventDisplay
A simple event dsisplay for the HGCal test beam experiments. The display has been tested with CMSSW_8_0_6, slc6_amd64_gcc530, running within a CERNVM virtual machine on a mac. It should work on lxplus and cmslpc-sl6.

# Installation
```linux
  git clone https://github.com/CMS-HGCAL/TestBeam.git HGcal
  cd HGCal
  cmsenv
  scram b

  cd TBEventDisplay
```
# To Run
Do 
```linux
  TBEventDisplay.py <geometry-file> <root-file-name>
```
where the geometry file is the same as that used in TBStandaloneSimulator
and the Root file is a file containing test beam digitized objects (SKIROC data frames)

