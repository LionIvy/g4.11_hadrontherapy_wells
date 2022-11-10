#ifndef INRPassiveProtonBeamLine_H
#define INRPassiveProtonBeamLine_H 1

#include "globals.hh"
#include "G4VUserDetectorConstruction.hh"
#include "G4Box.hh"
#include "G4Tubs.hh"
#include "G4Cons.hh"
#include "G4SubtractionSolid.hh"
#include "G4VisAttributes.hh"
#include "G4LogicalVolume.hh"
#include <fstream>

#include "PhaseSpaceDataset.hh"

class G4VPhysicalVolume;
class HadrontherapyDetectorConstruction;
//class PhaseSpaceDetector;
//class PhaseSpaceDataset;
//class PhaseSpaceDetectorConstruction;



//class HadrontherapyCellDetectorConstruction;
//class HT_PetriDishDetectorConstruction;


class HadrontherapyModulator;
class HadrontherapyContouredScatterer;
class HadrontherapyRidgeFilter;

class INRPassiveProtonBeamLineMessenger;
class HadrontherapyDetectorROGeometry;

class G4ProductionCuts;
//class G4ProductionCut;



class INRPassiveProtonBeamLine : public G4VUserDetectorConstruction
{
public:
    G4int DetType;

    INRPassiveProtonBeamLine();
    ~INRPassiveProtonBeamLine();
  // static G4bool doCalculation;

    G4VPhysicalVolume* Construct();
    //void ConstructSDandField();
    //***************************** PW **************NON SERVE*************************

    static INRPassiveProtonBeamLine* GetInstance();

    //***************************** PW **************NON SERVE*************************

    void HadrontherapyAcceleratorTube();
    // задание ионопровода установки

    void HadrontherapyWedge();
    // Задание клина

    void BeamDump(bool dumpOn);

    void HadrontherapyGraphiteCollimator();

    void HadrontherapyFirstScatterer();
    // Задание первичного рассеивателя

    void HadrontherapyWall();
    // Задание параметров биологической защиты


    void Collimator(G4double CollimatorBoxX, G4double CollimatorBoxY, G4double CollimatorBoxZ, G4double CollimatorR, G4double CollimatorXPosition, G4Material* CollimatorMaterial);

   // void HadrontherapyPreCollimator();
    // Задание коллиматора перед ГФ

  //  void HadrontherapyRidgeFilter();
    // задание ГФ

  //  void HadrontherapyMainCollimator();
    // Задание коллиматора

    void HadrontherapyBolus();
    // Задание болюса


   void CellTestPhantom();

    void SetFirstScatteringFoilXSize(G4double);
    // This method allows to change the size of the first scattering foil
    // along the X axis

    
    G4Box *solid_phSpDet, *solid_voxel;
    G4LogicalVolume *logical_phSpDet, *logic_voxel;
    G4VPhysicalVolume *phys_phSpDet, *phys_voxel;
    void construct_PhaseSpace_detector();

protected:
  G4LogicalVolume* phaseSpace_ScoringVolume = nullptr;

private:

    std::ifstream File;
    static INRPassiveProtonBeamLine* instance;
    //passive proton line dimensions
    void SetDefaultDimensions();
    void ConstructINRProtonBeamLine();

    HadrontherapyModulator*          modulator; // Pointer to the modulator
    HadrontherapyContouredScatterer* ContouredScatterer;
    HadrontherapyRidgeFilter*        RidgeFilter;

    // geometry component
    INRPassiveProtonBeamLineMessenger* passiveMessenger;
    G4VPhysicalVolume* physicalTreatmentRoom;
    HadrontherapyDetectorConstruction* hadrontherapyDetectorConstruction;
//    PhaseSpaceDetector* phsp_det;
//    PhaseSpaceDetectorConstruction* PhaseSpaceDetector;
//    PhaseSpaceDataset* PhaseSpace_data_collection;
  //  HadrontherapyCellDetectorConstruction* hadrontherapyCellDetectorConstruction;
  //  HT_PetriDishDetectorConstruction* petriDishDetectorConstruction;


    G4double Area1XShift;


    G4double vacuumZoneXSize;
    G4double vacuumZoneRSize;
    G4double vacuumZoneXPosition;
    G4double kaptonWindowXSize;
   // G4double kaptonWindowXPosition;

    G4double GraphiteCollimatorXSize;
    G4double GraphiteCollimatorYSize;
    G4double GraphiteCollimatorZSize;
    G4double GraphiteCollimatorXPosition;
    G4double GraphiteCollimatorRSize;

    G4double firstScatteringFoilXSize;
    G4double firstScatteringFoilYSize;
    G4double firstScatteringFoilZSize;
    G4double firstScatteringFoilXPosition;

    G4double WallXSize;
    G4double WallYSize;
    G4double WallZSize;
    G4double WallAngle;
    G4double WallXPosition;
    G4double MetalTubeThickness;
    G4double PolyTubeThickness ;
    G4double AirWindowR;
    G4double WallWindowR;

    G4double ContScattBoxPosX;

    G4double RidgeFilterBoxPosX;


    G4double RFCollimatorBoxX;
    G4double RFCollimatorBoxY;
    G4double RFCollimatorBoxZ;
    G4double RFCollimatorBoxR;
    G4double RFCollimatorPosX;


    G4double FinalCollimatorBoxX;
    G4double FinalCollimatorBoxY;
    G4double FinalCollimatorBoxZ;
    G4double FinalCollimatorBoxR;
    G4double FinalCollimatorPosX;


    G4double PlasticWaterBoxX;
    G4double PlasticWaterBoxY;
    G4double PlasticWaterBoxZ;
    G4double PlasticWaterBoxPosX;

    G4int EpendorfNumber;
    G4double EpendorfRminBottom;
    G4double EpendorfRmaxBottom;
    G4double EpendorfRminTop;
    G4double EpendorfRmaxTop;
    G4double EpendorfHeight;
    G4double Perem_dx1;
    G4double Perem_dx2;
    G4double Perem_dy;
    G4double Perem_dz;


//    G4double innerRadiusStopper;
//    G4double heightStopper;
//    G4double startAngleStopper;
//    G4double spanningAngleStopper;
//    G4double stopperXPosition;
//    G4double stopperYPosition;
//    G4double stopperZPosition;
//    G4double outerRadiusStopper;

//    G4double secondScatteringFoilXSize;
//    G4double secondScatteringFoilYSize;
//    G4double secondScatteringFoilZSize;
//    G4double secondScatteringFoilXPosition;
//    G4double secondScatteringFoilYPosition;
//    G4double secondScatteringFoilZPosition;

//    G4double rangeShifterXSize;
//    G4double rangeShifterYSize;
//    G4double rangeShifterZSize;
//    G4double rangeShifterXPosition;
//    G4double rangeShifterYPosition;
//    G4double rangeShifterZPosition;




    //G4Box* GraphiteCollimatorBox;
    //G4Tubs* GraphiteAirGap;
    G4Box*  firstScatteringFoil;
    G4Box*  ConcreteWall;
    G4Tubs* ConcreteWallGAP;
    G4SubtractionSolid* ConcreteWallwGAP;
    G4Tubs* MetalTube;
    G4Tubs* PolyTube;
    G4Box*  ContouredScattererBox;

    G4LogicalVolume* logicAr1;

G4ProductionCuts* cutsA1;
    G4LogicalVolume* logicMetalTube;
    G4LogicalVolume* logicPolyTube;
    G4LogicalVolume* logicConcreteWall;
  //  G4LogicalVolume* logicContouredScattererBox;
  //  G4LogicalVolume* logicContScattCu;




    G4VPhysicalVolume* physiAr1;
    G4VPhysicalVolume* physiGraphiteCollimatorBox;
    G4VPhysicalVolume* physiGraphiteAirGap;
    G4VPhysicalVolume* physiFirstScatteringFoil;
    G4VPhysicalVolume* physiConcreteWall;
    G4VPhysicalVolume* physiMetalTube;
    G4VPhysicalVolume* physiPolyTube;
  //  G4VPhysicalVolume* physiContouredScattererBox;
   // G4VPhysicalVolume* physiContScattCu;
    G4VPhysicalVolume* physiCollimator;

    G4VPhysicalVolume* physPlasticWater;
    G4VPhysicalVolume* physPlasticWater2;
    G4VPhysicalVolume* physEpendorfHat;
    G4VPhysicalVolume* physEpendorfCone;
    G4VPhysicalVolume* physEpendorfBottom;
    G4VPhysicalVolume* physEpendorfInsideHat;
    G4VPhysicalVolume* physEpendorfInsideCone;
    G4VPhysicalVolume* physBTwinE;

   G4VPhysicalVolume* physiKaptonWindow;

   G4VPhysicalVolume* physCoreMainBox;
   G4VPhysicalVolume* physShellMainBox;
   G4VPhysicalVolume* physShellAirBox;
   G4VPhysicalVolume* physShellWallFront;
   G4VPhysicalVolume* physShellWallBack;










    G4VisAttributes* redWire;
    G4VisAttributes* blue;
    G4VisAttributes* gray;
    G4VisAttributes* white;
    G4VisAttributes* black;
    G4VisAttributes* red;
    G4VisAttributes* yellow;
    G4VisAttributes* green;
    G4VisAttributes* darkGreen;
    G4VisAttributes* darkOrange3;
    G4VisAttributes* skyBlue;
    G4VisAttributes* metalColor;
    G4VisAttributes* sandColor;
    G4VisAttributes* EmptyColor;



    G4Material* kapton;
    G4Material* vacuumZoneMaterial;
    G4Material* vacuumZoneWindowMaterial;
    G4Material* GraphiteCollimatorMaterial;
    G4Material* dumpShellMaterial;
    G4Material* dumpCoreMaterial;
    G4Material* firstScatteringFoilMaterial;
    G4Material* WallMaterial;
    G4Material* MetalTubeMaterial;
    G4Material* PolyTubeMaterial;
    G4Material* ContouredScattMaterial;
    G4Material* ContouredCompensatorMaterial;
    G4Material* Collimator1Material;
    G4Material* Collimator2Material;
    G4Material* PlasticWaterMaterial;
    G4Material* EpendorfMaterial;

    G4Material*  EpendorfInsideMaterial;

   // G4Material* rangeShifterMaterial;
   // G4Material* beamLineSupportMaterial;
    G4Material* kaptonWindowMaterial;
  //  G4Material* stopperMaterial;
  //  G4Material* secondScatteringFoilMaterial;
  //  G4Material* firstCollimatorMaterial;
  //  G4Material* holeFirstCollimatorMaterial;
  //  G4Material* modulatorBoxMaterial;
  //  G4Material* holeModulatorBoxMaterial;
//    G4Material* layer1MonitorChamberMaterial;
//    G4Material* layer2MonitorChamberMaterial;
//    G4Material* layer3MonitorChamberMaterial;
//    G4Material* layer4MonitorChamberMaterial;
//    G4Material* MOPIMotherVolumeMaterial;
//    G4Material* MOPIFirstKaptonLayerMaterial;
//    G4Material* MOPIFirstAluminumLayerMaterial;
//    G4Material* MOPIFirstAirGapMaterial;
//    G4Material* MOPICathodeMaterial;
//    G4Material* MOPISecondAirGapMaterial;
//    G4Material* MOPISecondAluminumLayerMaterial;
//    G4Material* MOPISecondKaptonLayerMaterial;
//    G4Material* nozzleSupportMaterial;
//    G4Material* holeNozzleSupportMaterial;

//    G4Material* brassTubeMaterial;
//    G4Material* brassTube2Material;
//    G4Material* brassTube3Material;
//    G4Material* finalCollimatorMaterial;


    HadrontherapyDetectorROGeometry* RO;


};
#endif
