#include "HGCal/Reco/plugins/HGCalTBClustering.h"
#include "HGCal/Geometry/interface/HGCalTBTopology.h"
#include "HGCal/Geometry/interface/HGCalTBGeometryParameters.h"
#include "HGCal/Geometry/interface/HGCalTBCellVertices.h"

#include "algorithm"

void HGCalTBClustering::Run(HGCalTBRecHitCollection hitcol, std::vector<reco::HGCalTBCluster> &outClusterColl)
{
	HGCalTBCellVertices cellVertice;
	std::vector<HGCalTBDetId> temp;
	for( std::vector<HGCalTBRecHit>::iterator it=hitcol.begin(); it!=hitcol.end(); ++it){
		if( std::find(temp.begin(), temp.end(), (*it).id())!=temp.end() )
			continue;
		std::vector<HGCalTBDetId> clusterDetIDs;
		temp.push_back( (*it).id() );
		clusterDetIDs.push_back( (*it).id() );
		BuildCluster(hitcol, temp, clusterDetIDs);
		//std::cout << "clusterDetIDs.size() = " << clusterDetIDs.size() << std::endl;

		reco::HGCalTBCluster cluster;
		cluster.setLayer( clusterDetIDs.at(0).layer() );
		float energyHigh=0.;
		float energyLow=0.;
		float x,y,z; 
		x = y = z = 0.0;
		for( std::vector<HGCalTBDetId>::iterator jt=clusterDetIDs.begin(); jt!=clusterDetIDs.end(); ++jt){
			energyHigh+=(*hitcol.find(*jt)).energyHigh();
			energyLow+=(*hitcol.find(*jt)).energyLow();
			std::pair<double, double> CellCentreXY;
			CellCentreXY=cellVertice.GetCellCentreCoordinatesForPlots( ((*hitcol.find(*jt)).id()).layer(), 
					((*hitcol.find(*jt)).id()).sensorIU(), 
					((*hitcol.find(*jt)).id()).sensorIV(), 
					((*hitcol.find(*jt)).id()).iu(), 
					((*hitcol.find(*jt)).id()).iv(), 
					SENSORSIZE);
			x += CellCentreXY.first*(*hitcol.find(*jt)).energyHigh();
			y += CellCentreXY.second*(*hitcol.find(*jt)).energyHigh();
			//on verra plus tard pour z
		}
		cluster.setPosition( math::XYZPoint(x,y,z)/energyHigh );
		cluster.setEnergyLow(energyLow);
		cluster.setEnergyHigh(energyHigh);

		for( std::vector<HGCalTBDetId>::iterator jt=clusterDetIDs.begin(); jt!=clusterDetIDs.end(); ++jt)
			cluster.addHitAndFraction( (*jt), (*hitcol.find(*jt)).energyLow()/energyLow );

		outClusterColl.push_back(cluster);
	}
	//std::cout << "outClusterColl.size() " << outClusterColl.size() << std::endl;
}


void HGCalTBClustering::BuildCluster(HGCalTBRecHitCollection hitcol,
		std::vector<HGCalTBDetId> &temp,
		std::vector<HGCalTBDetId> &clusterDetIDs)
{ 
	HGCalTBTopology top;
	HGCalTBDetId detID=clusterDetIDs.back();
	std::set<HGCalTBDetId> neighbors=top.getNeighboringCellsDetID( detID, SENSORSIZE , settings.maxTransverse );
	for( std::set<HGCalTBDetId>::const_iterator it=neighbors.begin(); it!=neighbors.end(); ++it){
		if( settings.useDistanceInsteadCellID ) continue;
		if( std::find(temp.begin(), temp.end(), (*it))!=temp.end() || hitcol.find(*it)==hitcol.end() )
			continue;
		temp.push_back( (*it) );
		clusterDetIDs.push_back( (*it) );
		BuildCluster(hitcol, temp, clusterDetIDs);
	}
}

void HGCalTBClustering::RunSimple(HGCalTBRecHitCollection hitcol, reco::HGCalTBCluster &cluster)
{
	if( hitcol.size()==0 ) return;
	HGCalTBCellVertices cellVertice;
	HGCalTBTopology top;
	HGCalTBRecHit hitMax=(*hitcol.begin());
	for( std::vector<HGCalTBRecHit>::iterator it=hitcol.begin(); it!=hitcol.end(); ++it){
		if( (*it).energyHigh() > hitMax.energyHigh() )
			hitMax=(*it);
	}
	cluster.setLayer( hitMax.id().layer() );
	float energyHigh=hitMax.energyHigh();
	float energyLow=hitMax.energyLow();
	float x,y,z;
	x = y = z = 0.0;

	std::pair<double, double> CellCentreXY;
	CellCentreXY=cellVertice.GetCellCentreCoordinatesForPlots( hitMax.id().layer(), 
			hitMax.id().sensorIU(), 
			hitMax.id().sensorIV(), 
			hitMax.id().iu(), 
			hitMax.id().iv(), 
			SENSORSIZE);
	x += CellCentreXY.first*hitMax.energyHigh();
	y += CellCentreXY.second*hitMax.energyHigh();


	std::set<HGCalTBDetId> neighbors=top.getNeighboringCellsDetID( hitMax.id(), SENSORSIZE , settings.maxTransverse );
	for( std::set<HGCalTBDetId>::iterator jt=neighbors.begin(); jt!=neighbors.end(); ++jt){
		energyHigh+=(*hitcol.find(*jt)).energyHigh();
		energyLow+=(*hitcol.find(*jt)).energyLow();
		std::pair<double, double> CellCentreXY;
		CellCentreXY=cellVertice.GetCellCentreCoordinatesForPlots( ((*hitcol.find(*jt)).id()).layer(), 
				((*hitcol.find(*jt)).id()).sensorIU(), 
				((*hitcol.find(*jt)).id()).sensorIV(), 
				((*hitcol.find(*jt)).id()).iu(), 
				((*hitcol.find(*jt)).id()).iv(), 
				SENSORSIZE);
		x += CellCentreXY.first*(*hitcol.find(*jt)).energyHigh();
		y += CellCentreXY.second*(*hitcol.find(*jt)).energyHigh();
		//on verra plus tard pour z
	}
	cluster.setPosition( math::XYZPoint(x,y,z)/energyHigh );
	cluster.setEnergyLow(energyLow);
	cluster.setEnergyHigh(energyHigh);

	cluster.addHitAndFraction( hitMax.id(), hitMax.energyLow()/energyLow );
	for( std::set<HGCalTBDetId>::iterator jt=neighbors.begin(); jt!=neighbors.end(); ++jt)
		cluster.addHitAndFraction( (*jt), (*hitcol.find(*jt)).energyLow()/energyLow );

}



