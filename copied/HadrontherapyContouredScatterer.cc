#include <fstream>

#include "globals.hh"
#include "G4SystemOfUnits.hh"
#include "G4Material.hh"
#include "G4Tubs.hh"
#include "G4Box.hh"
#include "G4SubtractionSolid.hh"
#include "TubesParameterisation.hh"
#include "G4LogicalVolume.hh"
#include "G4VPhysicalVolume.hh"
#include "G4ThreeVector.hh"
#include "G4PVPlacement.hh"
#include "G4PVParameterised.hh"
#include "G4Transform3D.hh"
#include "G4RotationMatrix.hh"
#include "G4VisAttributes.hh"
#include "G4Colour.hh"
#include "HadrontherapyContouredScatterer.hh"
#include "G4Transform3D.hh"
#include "G4ios.hh"
#include "G4RunManager.hh"
#include "G4NistManager.hh"
#include <iostream>

using namespace std;

HadrontherapyContouredScatterer::HadrontherapyContouredScatterer(): physiMotherMod(),
                         logicScattMotherMod(),physiScattMotherMod(),
                         ContScattCuLayer(),logicContScattCu(),physiContScattCu(),
                         logicContScattPMMA(),
                         AirRing(),logicAirRing(),physiAirRing(),
                         FileName("ContouredScatterers/Sc_01.txt")
{ 

   CuStepNumber=23;
   CuHeight=new G4double[CuStepNumber];
   CuR=new G4double[CuStepNumber];
   PositionCu=new G4double[CuStepNumber];

   CuBoxWidth=0;
   CuBoxHeigth=0;
   CuBoxPosition=0;


   PMMAStepNumber=10;
   PMMABoxWidth=4.5*cm;
   PMMABoxBottomLayerHeight=0.05*cm;
   PositionPMMABoxCenter =0*cm;
   AirHeight=new G4double[PMMAStepNumber];
   AirR=new G4double[PMMAStepNumber];
   PositionAirRing=new G4double[PMMAStepNumber];

     
   for (G4int i=0;i<CuStepNumber;++i)
  {
    CuHeight[i]=0;
    CuR[i]=0;
    PositionCu[i]=0;
  }
   for (G4int i=0;i<PMMAStepNumber;++i)
  {
    AirHeight[i]=0;
    AirR[i]=0;
    PositionAirRing[i]=0;
  }

	
  
  ModulatorMessenger = new  HadrontherapyContouredScattererMessenger(this);
  ModulatorDefaultProperties();

  rm = new G4RotationMatrix(); 
  G4double phi = 270. *deg;     
  rm -> rotateY(phi); 
}
/////////////////////////////////////////////////////////////////////////////
HadrontherapyContouredScatterer::~HadrontherapyContouredScatterer()
{
  delete rm;
  delete [] CuHeight;
  delete []	CuR;
  delete []	PositionCu;

  delete [] AirHeight;
  delete []	AirR;
  delete []	PositionAirRing;

  delete logicScattMotherMod;
  delete physiScattMotherMod;

  delete ContScattCuLayer;
  delete logicContScattCu;
  delete physiContScattCu;
  delete logicContScattPMMA;
  delete physiContScattPMMA;
  delete AirRing;
  delete logicAirRing;
  delete physiAirRing;






  delete ModulatorMessenger;
}

/////////////////////////////////////////////////////////////////////////////
void HadrontherapyContouredScatterer::ModulatorDefaultProperties()
{
/* Here we initialize the step properties of Modulator wheel, you can create your
specific modulator by changing the values in this class or writing them in an external
file and activate reading from file via a macrofile.	*/	

    CuR[0]  = 8.50  ;   CuHeight[0] =  0.1;
    CuR[1]  = 9.36  ;   CuHeight[1] =  0.1;
    CuR[2]  = 10.11 ;   CuHeight[2] =  0.1;
    CuR[3]  = 10.80 ;   CuHeight[3] =  0.1;
    CuR[4]  = 11.48 ;   CuHeight[4] =  0.1;
    CuR[5]  = 12.17 ;   CuHeight[5] =  0.1;
    CuR[6]  = 12.89 ;   CuHeight[6] =  0.1;
    CuR[7]  = 13.64 ;   CuHeight[7] =  0.1;
    CuR[8]  = 14.42 ;   CuHeight[8] =  0.1;
    CuR[9]  = 15.22 ;   CuHeight[9] =  0.1;
    CuR[10] = 16.01 ;   CuHeight[10] = 0.1;
    CuR[11] = 16.80 ;   CuHeight[11] = 0.1;
    CuR[12] = 17.55 ;   CuHeight[12] = 0.1;
    CuR[13] = 18.26 ;   CuHeight[13] = 0.1;
    CuR[14] = 18.93 ;   CuHeight[14] = 0.1;
    CuR[15] = 19.55 ;   CuHeight[15] = 0.1;
    CuR[16] = 20.14 ;   CuHeight[16] = 0.1;
    CuR[17] = 20.73 ;   CuHeight[17] = 0.1;
    CuR[18] = 21.37 ;   CuHeight[18] = 0.1;
    CuR[19] = 22.10 ;   CuHeight[19] = 0.1;
    CuR[20] = 23.01 ;   CuHeight[20] = 0.1;
    CuR[21] = 24.20 ;   CuHeight[21] = 0.1;
    CuR[22] = 25.81 ;   CuHeight[22] = 0.1;


      AirR[0]  = 2.69;   AirHeight[0]= 0.27;
      AirR[1]  = 5.38;   AirHeight[1]= 0.78;
      AirR[2]  = 8.07;   AirHeight[2]= 1.22;
      AirR[3]  = 10.76;  AirHeight[3]= 1.61;
      AirR[4]  = 13.45;  AirHeight[4]= 1.84;
      AirR[5]  = 16.14;  AirHeight[5]= 2.00;
      AirR[6]  = 18.83;  AirHeight[6]= 2.04;
      AirR[7]  = 21.52;  AirHeight[7]= 1.93;
      AirR[8]  = 24.21;  AirHeight[8]= 1.74;
      AirR[9]  = 26.90;  AirHeight[9]= 0.7;

 //Plex_H_set=[ .027 0.078 0.122 0.161 0.184 0.200 0.204 0.193 0.174 0.07]
//    Plex_R_set=[ .269  .538  .807 1.076 1.345 1.614 1.883 2.152 2.421 2.690];
 GetStepInformation();	
 
} 
/////////////////////////////////////////////////////////////////////////////
void HadrontherapyContouredScatterer:: ModulatorPropertiesFromFile(G4String Name)
{
  delete [] CuHeight;
  delete []	CuR;
  delete []	PositionCu;

  delete [] AirHeight;
  delete []	AirR;
  delete []	PositionAirRing;

  delete logicScattMotherMod ;
  delete physiScattMotherMod;

  delete  ContScattCuLayer;
  delete  logicContScattCu;

  delete logicContScattPMMA;
  delete physiContScattPMMA;
  delete AirRing;
  delete logicAirRing;




//    delete physiAirRing;
//  /ConrouredScatterers/Sc_01.txt
  File.open(Name,  std::ios::in);
  if(!File.is_open())
  {
  G4cout<<" WARNING: The File with name of "<<Name<<
 " doesn't exist to get modulator step properties. please modify it and try again"<<G4endl;
 
 G4Exception("HadrontherapyContouredScatterer::ContouredScattererPropertiesFromFile( )", "Hadrontherapy0009"
 , FatalException, "Error: No available external file for reading from");
    }	  

  G4String string;

  File >>string;
  File >>string>> CuStepNumber;
  File >>string>>string>>string;

  CuHeight=new G4double[CuStepNumber];
  CuR = new G4double[CuStepNumber];

  for(G4int i=0;i<CuStepNumber;++i)
   {
	 G4String stringX;
     File>>stringX>> CuR[i]>>CuHeight[i];
   }	
   
  File >>string;
  File >>string>> PMMAStepNumber;
  File >>string>> PMMABoxWidth;
  File >>string>> PMMABoxBottomLayerHeight;
  File >>string>>string>>string;

  AirHeight=new G4double[PMMAStepNumber];
  AirR=new G4double[PMMAStepNumber];

  for(G4int i=0;i<PMMAStepNumber;++i)
   {
     G4String stringX;
     File>>stringX>> AirR[i]>>AirHeight[i];
   }
File.close();
   GetStepInformation();
   BuildCuLayers();
   


}
////////////////////////////////////////////////////////////////////////////////
void HadrontherapyContouredScatterer::GetStepInformation()
{

G4double CuRmax=0;
G4double CuRingsTotalHeight=0;
for (G4int i=0;i<CuStepNumber;++i)
{
    CuRingsTotalHeight = CuRingsTotalHeight+ CuHeight[i];
    if(CuR[i]>CuRmax) {CuRmax=CuR[i];}
}

//G4double Starting_point=-2.5*cm;
G4double Starting_point=-CuRingsTotalHeight;

CuBoxWidth=CuRmax;
CuBoxHeigth=0.5*CuRingsTotalHeight;
//CuBoxPosition=Starting_point+0.5*CuBoxHeigth;
CuBoxPosition=Starting_point+0.5*CuRingsTotalHeight;


 PositionCu[0]= -CuBoxHeigth+ 0.5* CuHeight[0];
 for (G4int i=1;i<CuStepNumber;++i)
 {
        PositionCu[i]=  PositionCu[i-1]+0.5*CuHeight[i-1] + 0.5* CuHeight[i];
 }


  G4double AirRingsTotalHeight=0;
  for (G4int i=0;i<PMMAStepNumber;++i)
 {
        AirRingsTotalHeight = AirRingsTotalHeight+ AirHeight[i];
 }

  PositionPMMABoxBottomLayerCenter = CuBoxPosition+0.5*CuRingsTotalHeight+0.5*PMMABoxBottomLayerHeight;

  PMMABoxHeight=AirRingsTotalHeight;
  PositionPMMABoxCenter = CuBoxPosition+0.5*CuRingsTotalHeight+PMMABoxBottomLayerHeight+0.5*PMMABoxHeight;
// !!!
          //PositionCu[CuStepNumber]+0.5*CuHeight[CuStepNumber] + 0.5*PMMABoxHeight;
//  PositionAirRing[0]= -0.5*PMMABoxHeight+PMMABoxBottomLayerHeight+0.5* AirHeight[0];
  PositionAirRing[0]= -0.5*PMMABoxHeight+0.5* AirHeight[0];
  for (G4int i=1;i<PMMAStepNumber;++i)
 {
  PositionAirRing[i]=  PositionAirRing[i-1]+0.5*AirHeight[i-1] + 0.5* AirHeight[i];
  }
  //G4double test=1;
	
	
}
/////////////////////////////////////////////////////////////////////////////////
void HadrontherapyContouredScatterer::BuildScatterer(G4VPhysicalVolume* motherVolume, G4double BoxPosX)
{
  G4bool isotopes = false;
  G4Material* airNist =  G4NistManager::Instance()->FindOrBuildMaterial("G4_AIR", isotopes);
  G4Material* copperNistAsMaterial = G4NistManager::Instance()->FindOrBuildMaterial("G4_Cu", isotopes);
  G4Material* PMMANist = G4NistManager::Instance()->FindOrBuildMaterial("G4_PLEXIGLASS", isotopes);

   Mod0Mater = airNist;
   ModMater = airNist; // You have to change modulator material via a macrofile (default is air)
   Mod1Mater= copperNistAsMaterial;
   Mod2Mater= PMMANist;

  // Mother volume
  G4ThreeVector positionMotherMod = G4ThreeVector(BoxPosX, 0 *mm, 0 *mm);
 
  G4Box* solidMotherMod = new G4Box("ContScattMotherMod", 10.0 *cm, 10.0 *cm, 5.0 *cm);
 
  logicMotherMod = new G4LogicalVolume(solidMotherMod, Mod0Mater,"LogContScattMotherMod",0,0,0);

  physiMotherMod = new G4PVPlacement(rm,positionMotherMod,  "PhysContScattMotherMod",
				     logicMotherMod,    				  
				     motherVolume,      
				     false,           
                     0);


    BuildCuLayers();
				     
				     
				     
				 }            
 ///////////////////////////////////////////////////////////////////////////////////////
 void HadrontherapyContouredScatterer::BuildCuLayers()
 {
// Объем для параметаризации медных пластин
//==========================================================
     G4Box* solidScattMotherMod = new G4Box("SolidScattMotherMod", CuBoxWidth, CuBoxWidth, CuBoxHeigth);

     logicScattMotherMod = new G4LogicalVolume(solidScattMotherMod, Mod0Mater,"LogScattMotherMod",0,0,0);

     physiScattMotherMod = new G4PVPlacement(0 ,G4ThreeVector(0., 0.,CuBoxPosition),  "PhysScattMotherMod",
                        logicScattMotherMod,
                        physiMotherMod,
                        false,
                        0);


    G4double* CuRnull= new G4double[CuStepNumber];
    for (G4int i=0;i<CuStepNumber;++i)   {CuRnull[i]=0;}
     ContScattCuLayer
       = new G4Tubs("ContScatt",0, 1*mm, 1*m, 0.*deg, 360.*deg);
     logicContScattCu
       = new G4LogicalVolume(ContScattCuLayer,Mod1Mater,"ContScatt",0,0,0);
     G4VPVParameterisation* chamberParam
       = new TubesParameterisation(
                                     CuStepNumber,
                                     PositionCu,
                                     CuHeight,
                                     CuRnull,
                                     CuR);
     physiContScattCu
       =  new G4PVParameterised("ContScatt",       // their name
                                           logicContScattCu,   // their logical volume
                                           logicScattMotherMod ,       // Mother logical volume
                                           kZAxis,          // Are placed along this axis
                                           CuStepNumber,    // Number of chambers
                                           chamberParam,    // The parametrisation
                                           true); // checking overlaps


//==========================================================
   G4Box* PMMABoxBottomLayer
       = new G4Box("SolidPMMABoxBottomLayer",
                   PMMABoxWidth,PMMABoxWidth,0.5*PMMABoxBottomLayerHeight);
   logicPMMABoxBottomLayer
       = new G4LogicalVolume(PMMABoxBottomLayer,Mod2Mater,"LogPMMABoxBottomLayer",0,0,0);
   physiPMMABoxBottomLayer = new G4PVPlacement(0, G4ThreeVector(0., 0.,PositionPMMABoxBottomLayerCenter),
                                        "PhysPMMABoxBottomLayer", logicPMMABoxBottomLayer,physiMotherMod, false, 0);


//==========================================================
   // Общая коробка PPMA из которой вычитается цилиндр радусом максимального кольца
   G4Box* ContScattPMMABox
     = new G4Box("ContScatt",PMMABoxWidth,PMMABoxWidth,0.5*PMMABoxHeight);
//Находим максимальный радиус с запасом
   G4double AirRingRmax=0;
   for (G4int i=0;i<PMMAStepNumber;++i)  {  if(AirR[i]>AirRingRmax) {AirRingRmax=AirR[i];}  }
   AirRingRmax=AirRingRmax+1*mm;

   G4Tubs* ContScattPMMATube = new G4Tubs("Empty",0, AirRingRmax, 0.5*PMMABoxHeight+0.01*mm,0.,360*deg);
   G4SubtractionSolid* ContScattPMMA =  new G4SubtractionSolid("ContScattPMMA",
                                                                ContScattPMMABox, ContScattPMMATube,
                                                                0, G4ThreeVector(0,0,0));

   logicContScattPMMA
     = new G4LogicalVolume(ContScattPMMA,Mod2Mater,"ContScattPMMA",0,0,0);
   physiContScattPMMA = new G4PVPlacement(0, G4ThreeVector(0., 0.,PositionPMMABoxCenter),
                                         "ContScattPMMA", logicContScattPMMA,physiMotherMod, false, 0);


G4Tubs* ContScattAirTube = new G4Tubs("Empty",0, AirRingRmax, 0.5*PMMABoxHeight,0.,360*deg);
logicContScattAirTube = new G4LogicalVolume(ContScattAirTube,ModMater,"Empty",0,0,0);
physiContScattAirTube = new G4PVPlacement(0, G4ThreeVector(0., 0.,PositionPMMABoxCenter),
                                   "Empty", logicContScattAirTube,physiMotherMod, false, 0);



 G4double* PMMARmin= new G4double[PMMAStepNumber];
 for (G4int i=0;i<PMMAStepNumber;++i)   {PMMARmin[i]=AirRingRmax;}


   AirRing
      = new G4Tubs("AirRing",0, 0.1*mm, 0.1*mm, 0.*deg, 360.*deg);
    logicAirRing
      = new G4LogicalVolume(AirRing,Mod2Mater,"AirRing",0,0,0);
    G4VPVParameterisation* chamberParam2
      = new TubesParameterisation(
                                    PMMAStepNumber,
                                    PositionAirRing,
                                    AirHeight,
                                    AirR,
                                    PMMARmin);


    physiAirRing
      =  new G4PVParameterised("AirRing",       // their name
                               logicAirRing,   // their logical volume
                               logicContScattAirTube,//logicContScattPMMA,       // Mother logical volume
                               kZAxis,          // Are placed along this axis
                               PMMAStepNumber,    // Number of chambers
                               chamberParam2,    // The parametrisation
                               true); // checking overlaps




  // Inform the kernel about the new geometry
    G4RunManager::GetRunManager() -> GeometryHasBeenModified();
    G4RunManager::GetRunManager() -> PhysicsHasBeenModified();



  G4VisAttributes* darkOrange3= new G4VisAttributes( G4Colour(205/255. , 102/255. ,  000/255. ));
    darkOrange3 -> SetVisibility(true);
    darkOrange3 -> SetForceSolid(true);
  G4VisAttributes * red = new G4VisAttributes( G4Colour(1. ,0. ,0.));
  red-> SetVisibility(true);
  red-> SetForceSolid(true);

  G4VisAttributes* EmptyColor = new G4VisAttributes(G4Colour());
  EmptyColor -> SetVisibility(true);
  EmptyColor -> SetForceWireframe(true);

  G4VisAttributes* skyBlue = new G4VisAttributes( G4Colour(135/255. , 206/255. ,  235/255. ));
  skyBlue -> SetVisibility(true);
  skyBlue -> SetForceSolid(true);

  logicScattMotherMod -> SetVisAttributes (G4VisAttributes::Invisible);//-> SetVisAttributes(EmptyColor);
  logicMotherMod-> SetVisAttributes (G4VisAttributes::Invisible);// -> SetVisAttributes(EmptyColor);

  logicContScattCu -> SetVisAttributes(darkOrange3);//EmptyColor);//
  logicPMMABoxBottomLayer-> SetVisAttributes(skyBlue);//red);//EmptyColor);//
  logicContScattPMMA -> SetVisAttributes(skyBlue);//EmptyColor);//
  logicAirRing-> SetVisAttributes(skyBlue);//darkOrange3);//EmptyColor);//
  logicContScattAirTube-> SetVisAttributes (G4VisAttributes::Invisible);
 }





/////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////
// Messenger values
/////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////
void HadrontherapyContouredScatterer::SetModulatorAngle(G4double angle)
{
  G4double rotationAngle = angle;
  rm -> rotateZ(rotationAngle);
  physiMotherMod -> SetRotation(rm);  
  G4cout << "MODULATOR HAS BEEN ROTATED OF " << rotationAngle/deg 
	 << " deg" << G4endl;
  G4RunManager::GetRunManager() -> GeometryHasBeenModified(); 
}
/////////////////////////////////////////////////////////////////////////
// Change modulator material
void HadrontherapyContouredScatterer::SetModulatorMaterial(G4String Material)
{
    if (G4Material* NewMaterial = G4NistManager::Instance()->FindOrBuildMaterial(Material, false) )
    {
	if (NewMaterial) 
	{
        for(G4int i=1;i<CuStepNumber;++i)
	    {
        logicContScattCu -> SetMaterial(NewMaterial);
	  //  G4RunManager::GetRunManager() -> PhysicsHasBeenModified();
	    G4RunManager::GetRunManager() -> GeometryHasBeenModified();
	    
	  //  G4cout<<(logicMod[i]->GetMaterial()->GetName())<<G4endl;
	}
	G4cout << "The material of the Modulator wheel has been changed to " << Material << G4endl;
	}
    }
    else
    {
	G4cout << "WARNING: material \"" << Material << "\" doesn't exist in NIST elements/materials"
	    " table [located in $G4INSTALL/source/materials/src/G4NistMaterialBuilder.cc]" << G4endl; 
	G4cout << "Use command \"/parameter/nist\" to see full materials list!" << G4endl; 
	
	
    }
}
/////////////////////////////////////////////////////////////////////////
// Change compensator material
void HadrontherapyContouredScatterer::SetCompensatorMaterial(G4String Material)
{
    if (G4Material* NewMaterial = G4NistManager::Instance()->FindOrBuildMaterial(Material, false) )
    {
    if (NewMaterial)
    {
        for(G4int i=1;i<CuStepNumber;++i)
        {
        logicContScattPMMA -> SetMaterial(NewMaterial);
      //  G4RunManager::GetRunManager() -> PhysicsHasBeenModified();
        G4RunManager::GetRunManager() -> GeometryHasBeenModified();

      //  G4cout<<(logicMod[i]->GetMaterial()->GetName())<<G4endl;
    }
    G4cout << "The material of the Modulator wheel has been changed to " << Material << G4endl;
    }
    }
    else
    {
    G4cout << "WARNING: material \"" << Material << "\" doesn't exist in NIST elements/materials"
        " table [located in $G4INSTALL/source/materials/src/G4NistMaterialBuilder.cc]" << G4endl;
    G4cout << "Use command \"/parameter/nist\" to see full materials list!" << G4endl;


    }
}
////////////////////////////////////////////////////////////////////////////////
// Change modulator position in the beam line
void HadrontherapyContouredScatterer::SetModulatorPosition(G4ThreeVector Pos)
{
  G4ThreeVector NewModulatorPos=Pos; 
  physiMotherMod -> SetTranslation( NewModulatorPos); 
  G4RunManager::GetRunManager() -> GeometryHasBeenModified();
  G4cout << "The modulator wheel is translated to"<<  NewModulatorPos/mm <<"mm " <<G4endl;
	
}

/////////////////////////////////////////////////////////////////////////////////
void HadrontherapyContouredScatterer:: GetDataFromFile(G4String value)

{
G4String Name=value;
if(value=="default" )	
{
Name=FileName;
}
G4cout<<" Step properties of modulator will be get out from the external file "
 <<Name<<G4endl;
ModulatorPropertiesFromFile(Name);
}
