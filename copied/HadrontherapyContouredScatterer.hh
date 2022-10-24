#ifndef HadrontherapyContouredScatterer_H
#define HadrontherapyContouredScatterer_H 1

#include "globals.hh"
#include <fstream>
#include "G4Material.hh"
#include "G4ThreeVector.hh"
#include "G4RotationMatrix.hh"
#include "HadrontherapyContouredScattererMessenger.hh"

// using namespace std;
class G4Tubs;
class G4Box;
class G4LogicalVolume;
class G4VPhysicalVolume;


class HadrontherapyContouredScatterer
{
public:

  HadrontherapyContouredScatterer();
  ~HadrontherapyContouredScatterer();

  void BuildScatterer(G4VPhysicalVolume*, G4double);
  void SetModulatorAngle(G4double);
  void SetModulatorMaterial(G4String);
  void SetCompensatorMaterial(G4String);
  void SetModulatorPosition(G4ThreeVector);
  void SetModulatorInnerRadius(G4double);
  void SetModulatorOuterRadius(G4double);
  void ModulatorDefaultProperties();
  void ModulatorPropertiesFromFile(G4String);
  void GetDataFromFile(G4String value);
  void GetStepInformation();
  void BuildCuLayers();

private:
  std::ifstream File;
   
   G4LogicalVolume * logicMotherMod ;
   G4VPhysicalVolume* physiMotherMod;


   G4LogicalVolume * logicScattMotherMod ;
   G4VPhysicalVolume* physiScattMotherMod;
   
  G4Material* Mod0Mater;
  G4Material* ModMater;
  G4Material* Mod1Mater;
  G4Material* Mod2Mater;

  G4Tubs*            ContScattCuLayer;
  G4LogicalVolume*   logicContScattCu;
  G4VPhysicalVolume* physiContScattCu;

  G4LogicalVolume*   logicPMMABoxBottomLayer;
  G4VPhysicalVolume* physiPMMABoxBottomLayer;
 // G4Box*             ContScattPMMABox;
  G4LogicalVolume*   logicContScattPMMA;
  G4VPhysicalVolume* physiContScattPMMA;
  G4LogicalVolume*   logicContScattAirTube;
  G4VPhysicalVolume* physiContScattAirTube;

  G4Tubs*            AirRing;
  G4LogicalVolume*   logicAirRing;
  G4VPhysicalVolume* physiAirRing;
  

  G4int CuStepNumber;
  G4double* CuHeight;
  G4double* CuR;
  G4double* PositionCu;

  G4double CuBoxWidth;
  G4double CuBoxHeigth;
  G4double CuBoxPosition;

  G4int PMMAStepNumber;
  G4double  PMMABoxWidth;
  G4double  PMMABoxBottomLayerHeight;
  G4double  PositionPMMABoxBottomLayerCenter;
  G4double  PMMABoxHeight;
  G4double  PositionPMMABoxCenter;
  G4double* AirHeight;
  G4double* AirR;
  G4double* PositionAirRing;


  G4RotationMatrix* rm;
  
  G4String FileName;
  HadrontherapyContouredScattererMessenger* ModulatorMessenger;

 

   
};
#endif
