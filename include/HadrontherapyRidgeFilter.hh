#ifndef HadrontherapyRidgeFilter_H
#define HadrontherapyRidgeFilter_H 1

#include "globals.hh"
#include <fstream>
#include "G4Material.hh"
#include "G4ThreeVector.hh"
#include "G4RotationMatrix.hh"
#include "HadrontherapyRidgeFilterMessenger.hh"

// using namespace std;
class G4Tubs;
class G4Box;
class G4LogicalVolume;
class G4VPhysicalVolume;
class G4VPVParameterisation;


class HadrontherapyRidgeFilter
{
public:

  HadrontherapyRidgeFilter();
  ~HadrontherapyRidgeFilter();

  void BuildFilter(G4VPhysicalVolume*, G4double);
  void SetFilterAngle(G4double);
  void SetFilterMaterial(G4String);
  void SetCompensatorMaterial(G4String);
  void SetFilterPosition(G4ThreeVector);

  void ModulatorDefaultProperties();
  void ModulatorReadPropertiesFromFile(G4String);
  void BuildModulatorFromFile(G4String);
  void GetDataFromFile(G4String value);
  void GetStepInformation();
  void BuildRFstairs();

private:
  std::ifstream File;

  G4Material* Mod0Mater;
  G4Material* ModMater;
  G4Material* Mod1Mater;
  G4Material* Mod2Mater;

   G4Box* solidMotherMod;
   G4LogicalVolume * logicMotherMod ;
   G4VPhysicalVolume* physiMotherMod;



   G4Box*             solidRFboxMotherMod;
   G4LogicalVolume*   logicRFboxMotherMod;
   G4VPhysicalVolume* physiRFboxMotherMod;
   



    G4Box*             RFstairBox;
    G4LogicalVolume*   logicRFstairBox;
    G4VPVParameterisation* chamberParam;
    G4VPhysicalVolume* physiRFstairBox;

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
  






  G4double FullXWidth;
  G4int NModules;
  G4int NStairs;
  G4int NElements;
  G4double* ElementFullHeight;
  G4double* ElementFullYWidth;
  G4double* ElementZPosition;
  G4double* ElementYPosition;
  G4double  ModulesYShift;
  G4double BoxFullWidth;
  G4double BoxFullHeigth;
  G4double BoxPosition;




//  G4int CuStepNumber;
//  G4double* CuHeight;
//  G4double* CuR;
//  G4double* PositionCu;
//  G4double CuBoxWidth;
//  G4double CuBoxHeigth;
//  G4double CuBoxPosition;
//  G4int PMMAStepNumber;
//  G4double  PMMABoxWidth;
//  G4double  PMMABoxBottomLayerHeight;
//  G4double  PositionPMMABoxBottomLayerCenter;
//  G4double  PMMABoxHeight;
//  G4double  PositionPMMABoxCenter;
//  G4double* AirHeight;
//  G4double* AirR;
//  G4double* PositionAirRing;


  G4RotationMatrix* rm;
  
  G4String FileName;
  HadrontherapyRidgeFilterMessenger* FilterMessenger;

 

   
};
#endif
