###############################################
scram project CMSSW_7_6_3_patch2
cd CMSSW_7_6_3_patch2/src/
cmsenv 
git cms-init
git clone git@github.com:CMS-HGCAL/TestBeam.git HGCal/
git checkout HGCal_SensorGeometry
scram b -j16
##################################################

To produce the fake Rechits(currently implemented for one sensor(128 cell) in one layer):
cd Reco
cmsRun test/produceRechitCollection_cfg.py

This produces a file test_RecHits_OneLayer_TB.root

To read this and plot the hits in the sensor:
cmsRun test/analyzeRechitCollection_cfg.py

This produces a file test_RecHitPlotter_OneLayer_TB.root which has the hits plotted in the cells of the hexagonal sensor.

The Geometry is in: HGCal/Geometry where the class HGCalTBCellVertices has the functions:
std::vector<std::pair<double,double>> GetCellCoordinates(int ix, int iv, int sensorsize) that returns the x,y vertices of a cell given ix,iv. Currently only full and half hexagons are implemented for 128 cell sensors.
The function std::pair<double,double> GetCellCentreCoordinates(int ix, int iv, int sensorsize) returns the x,y co-ordinates of a cell given ix,iv. This is useful in filling the full sensor histogram divided in bins shaped as cells for the plots. The former function is needed to set up this histogram.

In the same Geometry area the class HGCalTBTopology has a function bool ix_iv_valid(int ix, int iv, int sensorSize) that can say if the the values correspond to a cell with within a sensor(currently 128 cell)

Please Note that at the mouse-bitten corners the cells are represented as half hexagons and not pentagons yet.
