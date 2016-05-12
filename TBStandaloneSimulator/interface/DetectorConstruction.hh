#ifndef DetectorConstruction_h
#define DetectorConstruction_h 1

#include "SamplingSection.hh"

#include "G4VUserDetectorConstruction.hh"
#include "globals.hh"

#include "HGCal/TBStandaloneSimulator/interface/TBGeometry.h"
#include "HGCal/TBStandaloneSimulator/interface/TBConfig.h"

#include <map>
#include <string>

class G4VSolid;
class G4Box;
class G4LogicalVolume;
class G4VPhysicalVolume;
class G4Material;
class G4UniformMagField;
class DetectorMessenger;
class G4Colour;

/**
   @class DetectorConstruction
   @short builds a simple detector
 */
class DetectorConstruction : public G4VUserDetectorConstruction
{
public:

	/**
	   @short CTOR
	 */
	DetectorConstruction(TBGeometry& geometry, TBConfig& config);

	/**
	   @short calorimeter structure (sampling sections)
	 */

	std::vector<SamplingSection>* getStructure()
	{
		return &m_caloStruct;
	}

	int getModel()
	{
		return m_geometry.model();
	}
	int getVersion()
	{
		return m_geometry.version();
	}

	const std::vector<G4LogicalVolume*>& getSiLogVol()
	{
		return m_logicSi;
	}
	const std::vector<G4LogicalVolume*>& getAbsLogVol()
	{
		return m_logicAbs;
	}

	/**
	   @short define the calorimeter materials
	 */
	void DefineMaterials();

	/**
	   @short set magnetic field
	 */
	void SetMagField(G4double fieldValue);

	/**
	   @short set detector model
	 */

	void SetDetModel(G4int model);

	/**
	   @short DTOR
	 */
	~DetectorConstruction();


	/**
	   @short getters
	 */
	G4double GetCalorSizeXY()
	{
		return m_CalorSizeXY;
	}
	G4double GetCalorSizeZ()
	{
		return m_CalorSizeZ;
	}
	G4double GetWorldSizeXY()
	{
		return m_WorldSizeXY;
	}
	G4double GetWorldSizeZ()
	{
		return m_WorldSizeZ;
	}
	G4double GetCellSize()
	{
		return m_CellSize;
	}

	/**
	   @short build the detector
	 */

	G4VPhysicalVolume* Construct();

	bool saveTracks()
	{
		return m_config.savetracks;
	}

private:

	/**
	   @short compute the calor dimensions
	 */
	void UpdateCalorSize();

	/**
	   @short build the calorimeter
	 */
	G4VPhysicalVolume* ConstructCalorimeter();

	void buildSectorStack();

	G4double getParameter(std::string baseName,
	                      TBGeometry::Element& e,
	                      std::string param);

	G4double getUnits(std::string baseName,
	                  TBGeometry::Element& e);

	G4VSolid* constructSolid(std::string baseName,
	                         TBGeometry::Element& e);

	TBGeometry m_geometry;
	TBConfig   m_config;

	std::vector<SamplingSection>  m_caloStruct;
	std::map<std::string, G4Material*> m_materials;
	std::map<std::string, G4double> m_dEdx;
	std::map<std::string, G4Colour> m_colours;
	G4UniformMagField* m_magField;      //pointer to the magnetic field

	std::vector<G4Material* > m_SensitiveMaterial;

	G4double           m_CalorSizeXY, m_CalorSizeZ, m_CellSize;
	G4double           m_WorldSizeXY, m_WorldSizeZ;
	G4VSolid*          m_solidWorld;    //pointer to the solid World
	G4LogicalVolume*   m_logicWorld;    //pointer to the logical World
	G4VPhysicalVolume* m_physWorld;     //pointer to the physical World

	std::vector<G4LogicalVolume*>   m_logicSi;//pointer to the logical Si volumes
	//pointer to the logical absorber volumes situated just before the si
	std::vector<G4LogicalVolume*>   m_logicAbs;

	DetectorMessenger* m_detectorMessenger;  //pointer to the Messenger
	int model_;
};


#endif

