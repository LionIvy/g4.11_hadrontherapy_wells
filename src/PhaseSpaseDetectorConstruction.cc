#include "PhaseSpaseDetectorConstruction.hh"


#include "G4GeometryManager.hh"
#include "G4PVPlacement.hh"
#include "G4NistManager.hh"
#include "PhaseSpace_SD.hh"


PhaseSpaseDetectorConstruction* PhaseSpaseDetectorConstruction::instance = 0;

PhaseSpaseDetectorConstruction::PhaseSpaseDetectorConstruction(G4VPhysicalVolume* physicalTreatmentRoom, G4double zPosition = 0)
: motherPhys(physicalTreatmentRoom),    // pointer to WORLD volume
  detectorSD(0)
{
    /* NOTE! that the PhaseSpaseDetectorConstruction class
     * does NOT inherit from G4VUserDetectorConstruction G4 class
     * So the Construct() mandatory virtual method is inside another geometric class
     * like the passiveProtonBeamLIne, ...
     */

    // Messenger to change parameters of the phantom/detector geometry
  //  detectorMessenger = new HadrontherapyDetectorMessenger(this);

    // Default detector voxels size
    // 200 slabs along the beam direction (X)
//    sizeOfVoxelAlongX = 200 *um;
//    sizeOfVoxelAlongY = 4 *cm;
//    sizeOfVoxelAlongZ = 4 *cm;

    // Define here the material of the water phantom and of the detector
   // SetPhantomMaterial("G4_WATER");
    // Construct geometry (messenger commands)
    // SetDetectorSize(4.*cm, 4.*cm, 4.*cm);
   // SetDetectorSize(4. *cm, 4. *cm, 4. *cm);
  //  SetPhantomSize(40. *cm, 40. *cm, 40. *cm);

  //  SetPhantomPosition(G4ThreeVector(20. *cm, 0. *cm, 0. *cm));
   // SetDetectorToPhantomPosition(G4ThreeVector(0. *cm, 0. *cm, 0. *cm));
  //  SetDetectorPosition();
    //GetDetectorToWorldPosition();

    // Write virtual parameters to the real ones and check for consistency
 //   UpdateGeometry();

    voxel_sX = 1 *mm;
    voxel_sY = 1 *mm;
    voxel_sZ = 0.5 *mm;

    Nx = 100;
    Ny = 100;
    Nz = 1;

    // смещения:
    // + 10*см
    // + 30*см
    // + 130*см
    phSp_posX = zPosition;
    phSp_posY = 0;
    phSp_posZ = 0;

    voxel_HX = 0.5*voxel_sX;
    voxel_HY = 0.5*voxel_sY;
    voxel_HZ = 0.5*voxel_sZ;
    gap = 0 *mm;


    phSp_HX = 0.5*((voxel_sX+gap)*Nx - gap);
    phSp_HY = 0.5*((voxel_sY+gap)*Ny - gap);
    phSp_HZ = voxel_HZ*Nz;


    phSpDet_Material = G4NistManager::Instance()->FindOrBuildMaterial("G4_AIR", false);

}


PhaseSpaseDetectorConstruction::~PhaseSpaseDetectorConstruction(){

}

PhaseSpaseDetectorConstruction* PhaseSpaseDetectorConstruction::GetInstance()
{
    return instance;
}


G4VPhysicalVolume* PhaseSpaseDetectorConstruction::Construct(){

    solid_phSpDet = new G4Box("solid_phSpDet",
                              phSp_HX,
                              phSp_HY,
                              phSp_HZ);

    // Definition of the logical volume of the solid_phSpDet
    logical_phSpDet = new G4LogicalVolume(solid_phSpDet,
                                          phSpDet_Material,
                                          "phSpDet_Log", 0, 0, 0);
    G4RotationMatrix* rMatrix = new G4RotationMatrix();
    rMatrix -> rotateY(-90.0*deg);
    rMatrix -> rotateZ(180.0*deg);

    // Definition of the physics volume of the Phantom
    phys_phSpDet = new G4PVPlacement(rMatrix,
                                     G4ThreeVector(phSp_posX,phSp_posY,phSp_posZ),
                                     "phSpDet_Phys",
                                     logical_phSpDet,
                                     motherPhys,
                                     false,
                                     0);

    // Visualisation attributes of the phantom
    G4VisAttributes *whiteBlock = new G4VisAttributes(G4Colour(1, 1 ,1));
    whiteBlock -> SetVisibility(true);
    //whiteBlock -> SetForceSolid(true);
    whiteBlock -> SetForceWireframe(true);
    logical_phSpDet -> SetVisAttributes(whiteBlock);

    solid_voxel = new G4Box("solidSD", voxel_HX, voxel_HY, voxel_HZ);

    sensitiveDetectorName = "PhaseSpaceLV";
    logic_voxel = new G4LogicalVolume(solid_voxel, phSpDet_Material, sensitiveDetectorName);


    G4bool check4overlaps = false;
    int sdNum = 0;

//    phys_voxel =new G4PVPlacement(0,
//                                  G4ThreeVector(-phSp_HX + voxel_HX + 1*(2*voxel_HX + gap),
//                                                -phSp_HY + voxel_HY + 1*(2*voxel_HY + gap),
//                                                0),
//                                  logic_voxel, "physSD",
//                                  logical_phSpDet,
//                                  true, 0 + 0 * Nx * Ny,
//                                  check4overlaps);

    for (int i=0; i<Nx; ++i){
        for (int j=0; j < Ny; ++j){
            ++sdNum;
            phys_voxel =new G4PVPlacement(0,
                                          G4ThreeVector(-phSp_HX + voxel_HX + i*(2*voxel_HX + gap),
                                                        -phSp_HY + voxel_HY + j*(2*voxel_HY + gap),
                                                        0),
                                          logic_voxel, "PhaseSpacePhysVol",
                                          logical_phSpDet,
                                          true, j + i * Ny,
                                          check4overlaps);
        }
    }
    G4VisAttributes *solid_gray = new G4VisAttributes(G4Colour(65, 65 ,65 ));
    solid_gray -> SetVisibility(true);
    solid_gray -> SetForceSolid(true);

    logic_voxel -> SetVisAttributes(solid_gray);
    phaseSpace_ScoringVolume = logic_voxel;

    return phys_voxel;

}

void PhaseSpaseDetectorConstruction::ConstructSD(){

    PhaseSpace_SD *sensDet = new PhaseSpace_SD(sensitiveDetectorName);
    phaseSpace_ScoringVolume ->SetSensitiveDetector(sensDet);

}
