#ifndef PHASESPASEDETECTORCONSTRUCTION_HH
#define PHASESPASEDETECTORCONSTRUCTION_HH

#include "G4VUserDetectorConstruction.hh"

#include "G4SystemOfUnits.hh"

#include "G4Box.hh"
#include "globals.hh"
#include "G4VisAttributes.hh"
#include "G4LogicalVolume.hh"
#include "G4UnitsTable.hh"
#include "PhaseSpace_SD.hh"
//#include "HadrontherapyDetectorROGeometry.hh"

class G4VPhysicalVolume;
class G4LogicalVolume;
class G4PVPlacement;
//class HadrontherapyDetectorROGeometry;
//class HadrontherapyDetectorMessenger;
class PhaseSpaceDetector_SD;
//class HadrontherapyMatrix;

class PhaseSpaseDetectorConstruction: public G4VUserDetectorConstruction
{
public:

    PhaseSpaseDetectorConstruction(G4VPhysicalVolume*, G4double zPosition);

    ~PhaseSpaseDetectorConstruction();

    virtual G4VPhysicalVolume* Construct();
    virtual void ConstructSD();

    G4String sensitiveDetectorName = "PhaseSpaceLV";

public:
    static PhaseSpaseDetectorConstruction* GetInstance();
//    void InitializeDetectorROGeometry(HadrontherapyDetectorROGeometry*,
//                                      G4ThreeVector detectorToWorldPosition);
    G4VPhysicalVolume* motherPhys;
    PhaseSpaceDetector_SD*         detectorSD; // Pointer to sensitive detector

    //////////////////////////
    void VirtualLayer(G4bool Varbool);
    G4bool NewSource;
    void SetVirtualLayerPosition(G4ThreeVector);
    G4ThreeVector VirtualLayerPosition;

    //////////////////////////
private:

//    void ConstructPhSpDetector();
//    void ParametersCheck();
//    void CheckOverlaps();

public:

    void SetInWorldPosition(G4ThreeVector);
    void UpdateGeometry();

 //   G4LogicalVolume* GetDetectorLogicalVolume(){ return detectorLogicalVolume;}


private:
    static PhaseSpaseDetectorConstruction* instance;

    G4VisAttributes* skyBlue;
    G4VisAttributes* red;

 //   HadrontherapyMatrix*             matrix;

    G4Box *solid_phSpDet;
    G4LogicalVolume *logical_phSpDet;
    G4VPhysicalVolume *phys_phSpDet;

    G4Box* solid_voxel;
    G4LogicalVolume* logic_voxel;
    G4VPhysicalVolume*  phys_voxel;




    G4double voxel_sX = 1 *mm;
    G4double voxel_sY = 1 *mm;
    G4double voxel_sZ = 0.5 *mm;

    G4int Nx = 100;
    G4int Ny = 100;
    G4int Nz = 1;

    G4double phSp_posX = 0. *cm;
    G4double phSp_posY = 0. *cm;
    G4double phSp_posZ = 0. *cm;

    G4double voxel_HX = 0.5*voxel_sX;
    G4double voxel_HY = 0.5*voxel_sY;
    G4double voxel_HZ = 0.5*voxel_sZ;
    G4double gap = 0 *mm;

    G4double phSp_HX = 0.5 *mm;
    G4double phSp_HY = 0.5 *mm;
    G4double phSp_HZ = 0.5 *mm;

    G4Material* phSpDet_Material;

    G4LogicalVolume* phaseSpace_ScoringVolume;



};
#endif // PHASESPASEDETECTORCONSTRUCTION_HH
