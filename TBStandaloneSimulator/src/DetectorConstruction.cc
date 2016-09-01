// ----------------------------------------------------------------------------
// read geometry from a config file 04-28-2016 HBP
// ----------------------------------------------------------------------------
#include "HGCal/TBStandaloneSimulator/interface/DetectorConstruction.hh"
#include "HGCal/TBStandaloneSimulator/interface/DetectorMessenger.hh"
#include "HGCal/TBStandaloneSimulator/interface/HGCSSSimHit.hh"
#include "TMath.h"

#include "G4Material.hh"
#include "G4NistManager.hh"
#include "G4VSolid.hh"
#include "G4Box.hh"
#include "G4Tubs.hh"
#include "G4Polyhedra.hh"
#include "G4LogicalVolume.hh"
#include "G4PVPlacement.hh"
#include "G4PVReplica.hh"
#include "G4UniformMagField.hh"
#include "G4GeometryManager.hh"
#include "G4PhysicalVolumeStore.hh"
#include "G4LogicalVolumeStore.hh"
#include "G4SolidStore.hh"
#include "G4VisAttributes.hh"
#include "G4Colour.hh"
#include "G4PhysicalConstants.hh"
#include "G4SystemOfUnits.hh"
#include "G4FieldManager.hh"
#include "G4TransportationManager.hh"
#include "G4PhysicalConstants.hh"

using namespace std;

//
DetectorConstruction::DetectorConstruction(TBGeometry& geometry,
        TBConfig& config)
	: m_geometry(geometry),
	  m_config(config)
{
	//radiation lengths: cf. http://pdg.lbl.gov/2012/AtomicNuclearProperties/
	//W   3.504 mm
	//Pb  5.612 mm
	//Cu 14.36  mm

	G4cout << "BEGIN(DetectorConstruction)" << G4endl;

	vector<TBGeometry::Element> elements;
	for(size_t c = 0; c < m_geometry.size(); c++) {
		TBGeometry::Element e = m_geometry("geometry", c);
		if ( e.imap["first"] ) elements.clear();

		elements.push_back(e);

		if ( e.imap["last"] )
			m_caloStruct.push_back( SamplingSection(elements) );
	}

	DefineMaterials();
	SetMagField(0);
	UpdateCalorSize();
	m_detectorMessenger = new DetectorMessenger(this);

	G4cout << "END(DetectorConstruction)" << G4endl;
}

DetectorConstruction::~DetectorConstruction()
{
	delete m_detectorMessenger;
}

//
void DetectorConstruction::DefineMaterials()
{
	G4NistManager* nistManager = G4NistManager::Instance();
	m_materials["Abs"]
	    = nistManager->FindOrBuildMaterial(m_geometry.material("Abs"), false);
	m_materials["AbsHCAL"]
	    = nistManager->FindOrBuildMaterial(m_geometry.material("AbsHCAL"), false);
	m_materials["Al"] = nistManager->FindOrBuildMaterial("G4_Al", false);
	m_dEdx["Al"] = 0.4358;
	m_materials["W"] = nistManager->FindOrBuildMaterial("G4_W", false);
	m_dEdx["W"] = 2.210;
	m_materials["Pb"] = nistManager->FindOrBuildMaterial("G4_Pb", false);
	m_dEdx["Pb"] = 1.274;
	m_materials["Cu"] = nistManager->FindOrBuildMaterial("G4_Cu", false);
	m_dEdx["Cu"] = 1.257;
	m_materials["Si"] = nistManager->FindOrBuildMaterial("G4_Si", false);
	m_dEdx["Si"] = 0.3876;
	m_materials["Zn"] = nistManager->FindOrBuildMaterial("G4_Zn", false);
	m_dEdx["Zn"] = 1.007;
	m_materials["Air"] = nistManager->FindOrBuildMaterial("G4_AIR", false);
	m_dEdx["Air"] = 0;
	m_materials["Fe"] = nistManager->FindOrBuildMaterial("G4_Fe", false);
	m_dEdx["Fe"] = 1.143;
	m_materials["Mn"] = nistManager->FindOrBuildMaterial("G4_Mn", false);
	m_dEdx["Mn"] = 1.062 ;
	m_materials["C"] = nistManager->FindOrBuildMaterial("G4_C", false);
	m_dEdx["C"] = 0.3952;
	m_materials["H"] = nistManager->FindOrBuildMaterial("G4_H", false);
	m_dEdx["H"] =  0;
	m_materials["Cl"] = nistManager->FindOrBuildMaterial("G4_Cl", false);
	m_dEdx["Cl"] = 0;
	m_materials["Cr"] = nistManager->FindOrBuildMaterial("G4_Cr", false);
	m_dEdx["Cr"] = 1.046;
	m_materials["Ni"] = nistManager->FindOrBuildMaterial("G4_Ni", false);
	m_dEdx["Ni"] = 1.307;


	//from http://www.physi.uni-heidelberg.de/~adler/TRD/
	//TRDunterlagen/RadiatonLength/tgc2.htm

	G4Element* H  = new G4Element("Hydrogen", "H",  1., 1.01 * g / mole);
	G4Element* C  = new G4Element("Carbon",   "C",  6., 12.01 * g / mole);
	G4Element* O  = new G4Element("Oxygen",   "O" , 8., 16.00 * g / mole);
	G4Element* Si = new G4Element("Silicon",  "Si", 14., 28.09 * g / mole);

	m_materials["SiO2"] = new G4Material("SiO2", 2.200 * g / cm3, 2);
	m_materials["SiO2"]->AddElement(Si, 1);
	m_materials["SiO2"]->AddElement(O , 2);
	m_dEdx["SiO2"] = 5.938; // MeV/cm

	//Epoxy
	m_materials["Epoxy"] = new G4Material("Epoxy", 1.2 * g / cm3, 2);
	m_materials["Epoxy"]->AddElement(H, 2);
	m_materials["Epoxy"]->AddElement(C, 2);
	m_dEdx["Epoxy"] = 2.0; // JUST A GUESS

	//FR4
	m_materials["FR4"] = new G4Material("FR4", 1.86 * g / cm3, 2);
	m_materials["FR4"]->AddMaterial(m_materials["SiO2"], 0.528);
	m_materials["FR4"]->AddMaterial(m_materials["Epoxy"], 0.472);
	m_dEdx["FR4"] = 0.528 * m_dEdx["SiO2"] + 0.472 * m_dEdx["Epoxy"];

	//G10
	m_materials["G10"] = new G4Material("G10", 1.700 * g / cm3, 4);
	m_materials["G10"]->AddElement(nistManager->FindOrBuildElement(14), 1);
	m_materials["G10"]->AddElement(nistManager->FindOrBuildElement(8) , 2);
	m_materials["G10"]->AddElement(nistManager->FindOrBuildElement(6) , 3);
	m_materials["G10"]->AddElement(nistManager->FindOrBuildElement(1) , 3);
	m_dEdx["PCB"] = 3.179; // MeV/cm

	m_materials["PCB"] = m_materials["G10"];
	m_dEdx["PCB"] = m_dEdx["G10"];


	m_materials["Brass"] = new G4Material("Brass", 8.5 * g / cm3, 2);
	m_materials["Brass"]->AddMaterial(m_materials["Cu"]  , 70 * perCent);
	m_materials["Brass"]->AddMaterial(m_materials["Zn"]  , 30 * perCent);
	m_dEdx["Brass"] = 0.7 * m_dEdx["Cu"] + 0.3 * m_dEdx["Zn"];
	m_materials["Steel"] = new G4Material("Steel", 7.87 * g / cm3, 3);
	m_materials["Steel"]->AddMaterial(m_materials["Fe"]  , 0.9843);
	m_materials["Steel"]->AddMaterial(m_materials["Mn"], 0.014);
	m_materials["Steel"]->AddMaterial(m_materials["C"], 0.0017);
	m_dEdx["Steel"] = 0.9843 * m_dEdx["Fe"] + 0.014 * m_dEdx["Mn"] + 0.0017 * m_dEdx["C"];
	m_materials["SSteel"] = new G4Material("SSteel", 8.02 * g / cm3, 4);
	m_materials["SSteel"]->AddMaterial(m_materials["Fe"]  , 0.70);
	m_materials["SSteel"]->AddMaterial(m_materials["Mn"], 0.01);
	m_materials["SSteel"]->AddMaterial(m_materials["Cr"], 0.19);
	m_materials["SSteel"]->AddMaterial(m_materials["Ni"], 0.10);
	m_dEdx["SSteel"]
	    = 0.7 * m_dEdx["Fe"] + 0.01 * m_dEdx["Mn"] + 0.19 * m_dEdx["Cr"] + 0.1 * m_dEdx["Ni"];
	m_materials["Scintillator"]
	    = nistManager->FindOrBuildMaterial("G4_POLYSTYRENE", false);
	m_dEdx["Scintillator"] = m_dEdx["C"];
	G4cout << m_materials["Scintillator"] << G4endl;
	m_materials["Polystyrole"] = new G4Material("Polystyrole", 1.065 * g / cm3, 2);
	m_materials["Polystyrole"]->AddMaterial(m_materials["H"]  , 50 * perCent);
	m_materials["Polystyrole"]->AddMaterial(m_materials["C"]  , 50 * perCent);
	m_dEdx["Polystyrole"] = 0.5 * m_dEdx["C"];

	m_materials["PVC"] = new G4Material("PVC", 1.350 * g / cm3, 3);
	m_materials["PVC"]->AddMaterial(m_materials["H"]  , 50 * perCent);
	m_materials["PVC"]->AddMaterial(m_materials["C"]  , 33.33 * perCent);
	m_materials["PVC"]->AddMaterial(m_materials["Cl"]  , 16.67 * perCent);
	m_dEdx["PVC"] = 0.33 * m_dEdx["C"];

	m_materials["CFMix"] = new G4Material("CFMix", 0.120 * g / cm3, 3);
	m_materials["CFMix"]->AddMaterial(m_materials["Air"]  , 0.009);
	m_materials["CFMix"]->AddMaterial(m_materials["PVC"]  , 0.872);
	m_materials["CFMix"]->AddMaterial(m_materials["Polystyrole"]  , 0.119);
	m_dEdx["CFMix"] = 0;

	m_materials["Foam"] = new G4Material("Foam", 0.0999 * g / cm3, 2);
	m_materials["Foam"]->AddMaterial(m_materials["C"]  , 0.856);
	m_materials["Foam"]->AddMaterial(m_materials["H"]  , 0.144);
	m_dEdx["Foam"] = 1.749 * 0.856 * 0.0999 / 10.;

	m_materials["WCu"] = new G4Material("WCu", 14.979 * g / cm3, 2);
	m_materials["WCu"]->AddMaterial(m_materials["W"]  , 75 * perCent);
	m_materials["WCu"]->AddMaterial(m_materials["Cu"]  , 25 * perCent);
	m_dEdx["WCu"] = 0.75 * m_dEdx["W"] + 0.25 * m_dEdx["Cu"];

	m_materials["NeutMod"] = new G4Material("NeutMod", 0.950 * g / cm3, 2);
	m_materials["NeutMod"]->AddMaterial(m_materials["C"]  , 0.85628);
	m_materials["NeutMod"]->AddMaterial(m_materials["H"]  , 0.14372);
	m_dEdx["NeutMod"] = 1.749 * 0.86 * 0.950 / 10.;

}

//
void DetectorConstruction::UpdateCalorSize()
{
	size_t nsections =  m_caloStruct.size();

	// compute Z extent of detector
	TBGeometry::Element efirst = m_caloStruct[0].getElement(0);
	TBGeometry::Element elast  = m_caloStruct[nsections - 1].getElement(-1);

	G4double units = getUnits(m_caloStruct[0].ele_name[0], efirst);
	G4double center = (elast.dmap["z"] + efirst.dmap["z"]) * units / 2;
	m_CalorSizeZ
	    = (elast.dmap["z"] - efirst.dmap["z"]
	       + (elast.dmap["thickness"] + efirst.dmap["thickness"]) / 2) * units;

	// get XY extent of detector
	m_CellSize = 0;
	m_CalorSizeXY = 0;
	for(size_t i = 0; i < nsections; i++) {
		SamplingSection& section = m_caloStruct[i];
		size_t nEle = section.n_elements;
		for (unsigned ie(0); ie < nEle; ++ie) {
			TBGeometry::Element e = section.getElement(ie);



			// get the XY extent of the sensitive layers
			bool isSensitive = false;
			if ( e.imap.find("sensitive") != e.imap.end() )
				isSensitive = static_cast<bool>(e.imap["sensitive"]);

			G4double units = getUnits(section.ele_name[ie], e);
			if      ( e.smap["shape"] ==  "cylinder" ) {
				if ( isSensitive ) {
					G4double maxR = getParameter(section.ele_name[ie],
					                             e, "maxRadius") * units;
					if ( maxR > m_CalorSizeXY ) m_CalorSizeXY = maxR;
				}
			} else if ( e.smap["shape"] ==  "hexagon" ) {
				if ( isSensitive ) {
					m_CellSize = getParameter(section.ele_name[ie],
					                          e, "cellsize") * units;
					assert(m_CellSize > 0);

					// get center-to-corner distance
					G4double width = 2 * getParameter(section.ele_name[ie],
					                                  e, "side") * units;
					if ( width > m_CalorSizeXY ) m_CalorSizeXY = width;
				}
			} else if ( e.smap["shape"] ==  "square" ) {
				if ( isSensitive ) {
					G4double side = getParameter(section.ele_name[ie],
					                             e, "side") * units;
					if ( side > m_CalorSizeXY ) m_CalorSizeXY = side;
				}
			}
		}
	}

	std::string u = efirst.smap["units"];
	G4cout << "  XY  extent of detector: " << m_CalorSizeXY << u << G4endl;
	G4cout << "  Z   extent of detector: " << m_CalorSizeZ  << u << G4endl;
	G4cout << "  Z location of detector: " << center << u << G4endl;
	G4cout << "  cell size  of detector: " << m_CellSize  << u << G4endl;
}

//
G4VPhysicalVolume* DetectorConstruction::Construct()
{

	//clean old geometry
	G4GeometryManager::GetInstance()->OpenGeometry();
	G4PhysicalVolumeStore::GetInstance()->Clean();
	G4LogicalVolumeStore::GetInstance()->Clean();
	G4SolidStore::GetInstance()->Clean();

	//world
	TBGeometry::Element world = m_geometry("world");

	G4double units = getUnits("World", world);
	G4double expHall_x = world.dmap["xside"] * units;
	G4double expHall_y = world.dmap["yside"] * units;
	G4double expHall_z = world.dmap["zside"] * units;

	G4Box* experimentalHall_box = new G4Box("expHall_box",
	                                        expHall_x / 2,
	                                        expHall_y / 2,
	                                        expHall_z / 2);

	G4LogicalVolume* experimentalHall_log
	    = new G4LogicalVolume(experimentalHall_box, m_materials["Air"],
	                          "expHall_log");

	G4VPhysicalVolume* experimentalHall_phys
	    = new G4PVPlacement(0,                       // no rotation
	                        G4ThreeVector(0., 0., 0.), // translation position
	                        experimentalHall_log,    // its logical volume
	                        "expHall",               // its name
	                        0,                       // its mother volume
	                        false,                   // no boolean operations
	                        0);                      // its copy number

	// define detector's world volume, basically its
	// bounding box
	m_WorldSizeXY = 0.99 * expHall_x;
	m_WorldSizeZ  = 0.99 * expHall_z;
	m_solidWorld  = new G4Box("Wbox",
	                          m_WorldSizeXY / 2,
	                          m_WorldSizeXY / 2,
	                          m_WorldSizeZ / 2);

	m_logicWorld = new G4LogicalVolume(m_solidWorld,
	                                   m_materials["Air"],
	                                   "Wlog");
	m_logicWorld->SetVisAttributes(G4VisAttributes::Invisible);


	G4double xpos = 0.*mm;
	G4double ypos = 0.*mm;
	G4double zpos = 0.*mm;

	m_physWorld = new G4PVPlacement(0,
	                                G4ThreeVector(xpos, ypos, zpos),
	                                m_logicWorld,
	                                "Wphys",
	                                experimentalHall_log,
	                                false, 0);
	// build the solids that define the geometry
	buildSectorStack();

	return experimentalHall_phys;
}

void DetectorConstruction::buildSectorStack()
{
	char nameBuf[80];
	G4VSolid* solid = 0;
	G4double totalLengthX0 = 0;
	G4double totalLengthL0 = 0;

	// A section comprises an ordered sequence of elements.
	// An element is either an air gap, an absorber, or a sensor.

	G4cout << "== buildSectorStack == " << G4endl;
	G4cout << "number of sections: " << m_caloStruct.size()
	       << G4endl;

	for(size_t i = 0; i < m_caloStruct.size(); i++) {
		G4cout << "  ==> section: " << i << G4endl;
		const size_t nEle = m_caloStruct[i].n_elements;
		G4cout << "  number of elements/section: " << nEle << G4endl;


		//index for counting sensitive elements/section
		int idx = 0;

		// loop over elements of current section
		for (size_t ie(0); ie < nEle; ++ie) {
			// NOTE: must be a reference, NOT a copy!
			SamplingSection& section = m_caloStruct[i];

			TBGeometry::Element element = section.getElement(ie);
			G4double units = getUnits(section.ele_name[ie], element);

			// construct a name for the current element based on its
			// material
			std::string eleName = element.smap["material"];
			sprintf(nameBuf, "%s%d", eleName.c_str(), int(i + 1));
			// if the element is sensitive, give it a slightly different
			// name.
			bool isSensitive = static_cast<bool>(element.imap["sensitive"]);
			if ( isSensitive ) {
				sprintf(nameBuf, "%s%d_%d", eleName.c_str(), int(i + 1), idx);
				idx++;
			}
			std::string baseName(nameBuf);

			G4double thick = element.dmap["thickness"] * units;
			G4cout << "  \telement number: " << ie
			       << "\tname: " << nameBuf
			       << "\tthickness: " << thick <<  "mm" << G4endl;

			// build the solid corresponding to the current element
			solid = constructSolid(baseName, element);

			// build the logical volume that defines everything about
			// the element except its location
			G4LogicalVolume* logi = new G4LogicalVolume(solid,
			        m_materials[eleName],
			        baseName + "log");

			// cache various quantities for this element: radiation length,
			// dE/dx, and nuclear interaction length
			section.ele_X0[ie]   = m_materials[eleName]->GetRadlen();
			section.ele_dEdx[ie] = m_dEdx[eleName];
			section.ele_L0[ie]   = m_materials[eleName]->GetNuclearInterLength();
			totalLengthX0 += thick / section.ele_X0[ie];
			totalLengthL0 += thick / section.ele_L0[ie];

			G4cout << "  \t" << eleName
			       << " dEdx=" << section.ele_dEdx[ie]
			       << " X0=" << section.ele_X0[ie]
			       << G4endl
			       << "  \tL0=" << section.ele_L0[ie]
			       << G4endl
			       << "  \tTotX0=" << totalLengthX0
			       << " \tTotLambda=" << totalLengthL0 << G4endl;

			// keep track of volumes of sensitive elements
			if ( isSensitive ) m_logicSi.push_back(logi);

			// now place volume
			G4double xpos = getParameter(baseName, element, "x") * units;
			G4double ypos = getParameter(baseName, element, "y") * units;
			G4double zpos = getParameter(baseName, element, "z") * units;
			section.ele_vol[ie] =
			    new G4PVPlacement(0,
			                      G4ThreeVector(xpos, ypos, zpos),
			                      logi,
			                      baseName + "phys",
			                      m_logicWorld, false, 0);
			G4cout << "  \tpositioning element at (" << xpos
			       << ", " << ypos <<  ", " << zpos << ")"
			       << G4endl << G4endl;

			// make air invisible
			if ( element.smap["material"] ==  string("Air") )
				logi->SetVisAttributes(G4VisAttributes::Invisible);
			else {
				G4VisAttributes* simpleBoxVisAtt =
				    new G4VisAttributes(section.g4Colour(ie));
				simpleBoxVisAtt->SetVisibility(true);
				logi->SetVisAttributes(simpleBoxVisAtt);
			}


			//for sensitive volumes
			//add region to be able to set specific cuts for it
			//just for Si
			if ( isSensitive ) {
				unsigned nlogicsi = m_logicSi.size();
				G4Region* aRegion = new G4Region(baseName + "Reg");
				m_logicSi[nlogicsi - 1]->SetRegion(aRegion);
				aRegion->AddRootLogicalVolume(m_logicSi[nlogicsi - 1]);
			}
		}//loop over elements
	}//loop over sections
}//buildstack


void DetectorConstruction::SetMagField(G4double fieldValue)
{

	if(fieldValue <= 0) return;

	//apply a global uniform magnetic field along Z axis
	G4FieldManager* fieldMgr
	    = G4TransportationManager::GetTransportationManager()->GetFieldManager();
	if(m_magField) delete m_magField;   //delete the existing magn field
	m_magField = new G4UniformMagField(G4ThreeVector(0., 0., fieldValue));
	fieldMgr->SetDetectorField(m_magField);
	fieldMgr->CreateChordFinder(m_magField);
	fieldMgr->SetDetectorField(m_magField);
}

void DetectorConstruction::SetDetModel(G4int model)
{
	if (model <= 0) return;
	G4cout << " -- Setting detector model to " << model << G4endl;
	model_ = model;
}


double DetectorConstruction::getParameter(std::string name,
        TBGeometry::Element& e,
        std::string param)
{
	double x = 0;
	try {
		x = e.dmap[param];
	} catch (...) {
		G4cout << "** DetectorConstruction ** problem building "
		       << "detector element " << name << G4endl
		       << "** parameter " << param << " is missing. Check config file."
		       << G4endl;
		exit(0);
	}
	return x;
}


G4double DetectorConstruction::getUnits(std::string name,
                                        TBGeometry::Element& e)
{
	std::string units;
	try {
		units = e.smap["units"];
	} catch (...) {
		G4cout << "** DetectorConstruction ** no units given for"
		       << "detector element " << name << G4endl;
		exit(0);
	}

	G4double unit = mm;
	if      ( units == "m" )
		unit = m;
	else if ( units == "cm" )
		unit = cm;
	return unit;
}


G4VSolid* DetectorConstruction::constructSolid (std::string baseName,
        TBGeometry::Element& e)
{
	G4double units = getUnits(baseName, e);
	G4VSolid* solid = 0;
	if      ( e.smap["shape"] ==  "cylinder" ) {
		G4double minR     = getParameter(baseName, e, "minRadius") * units;
		G4double maxR     = getParameter(baseName, e, "maxRadius") * units;
		G4double thickness = getParameter(baseName, e, "thickness") * units;
		solid = new G4Tubs(baseName + "cyl",
		                   minR, maxR, thickness / 2, 0.*deg, 360.0 * deg);
	} else if ( e.smap["shape"] ==  "hexagon" ) {
		G4double thickness = getParameter(baseName, e, "thickness") * units;
		G4double zPlane[2] = { -thickness / 2, thickness / 2};
		// feed center-to-side distance to G4Polyhedra
		// Note: side is full length of hexagon side
		G4double side     = getParameter(baseName, e, "side") * units;
		double   minRadius = (side / 2) * sqrt(3.0) / 2;
		G4double rInner[2] = {0, 0};
		G4double rOuter[2] = {minRadius, minRadius};
		solid = new G4Polyhedra(baseName + "hexa",
		                        0,
		                        2 * TMath::Pi(),
		                        6,
		                        2, zPlane,
		                        rInner, rOuter);
	} else if ( e.smap["shape"] ==  "square" ) {
		G4double thickness = getParameter(baseName, e, "thickness") * units;
		G4double side     = getParameter(baseName, e, "side") * units;
		solid = new G4Box(baseName + "box", side / 2, side / 2, thickness / 2);
	} else {
		G4cout << "** DetectorConstruction ** problem building "
		       << "detector element " << baseName << G4endl
		       << "** shape " << e.smap["shape"] << " not yet implemented"
		       << G4endl;
		exit(0);
	}
	return solid;
}
