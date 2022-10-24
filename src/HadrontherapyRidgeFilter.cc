
#include <fstream>

#include "globals.hh"
#include "G4SystemOfUnits.hh"
#include "G4Material.hh"
#include "G4Tubs.hh"
#include "G4Box.hh"
#include "G4SubtractionSolid.hh"
#include "BoxesParameterisation.hh"
#include "G4LogicalVolume.hh"
#include "G4VPhysicalVolume.hh"
#include "G4ThreeVector.hh"
#include "G4PVPlacement.hh"
#include "G4PVParameterised.hh"
#include "G4Transform3D.hh"
#include "G4RotationMatrix.hh"
#include "G4VisAttributes.hh"
#include "G4Colour.hh"
#include "HadrontherapyRidgeFilter.hh"
#include "G4Transform3D.hh"
#include "G4ios.hh"
#include "G4RunManager.hh"
#include "G4NistManager.hh"
#include <iostream>

using namespace std;

HadrontherapyRidgeFilter::HadrontherapyRidgeFilter(): physiMotherMod(),
                         logicRFboxMotherMod(),physiRFboxMotherMod(),
                         RFstairBox(),logicRFstairBox(),chamberParam(),physiRFstairBox(),
                         FileName("RidgeFilters/RF_01.txt")
{ 

   FullXWidth=10*cm;                                    // Длина ступеней вдоль оси Х
   NModules=10;                                         // Количество одинаковых наборов ступеней
   NStairs=5;                                           // Количество ступеней в одном наборе

   NElements=NModules*NStairs;                          // Общее количество ступеней

   ElementFullHeight= new G4double[NElements];          // Полная высота каждой ступени
   ElementFullYWidth= new G4double[NElements];          // Полная ширина каждой ступени вдоль оси Y
   ElementZPosition= new G4double[NElements];           // Сдвиг каждой ступени по оси Z
   ElementYPosition= new G4double[NElements];           // Сдвиг каждой ступени по оси Y

   ModulesYShift=2.0*mm;                                // Расстояние между модулями гребенок



   BoxFullWidth=0;
   BoxFullHeigth=0;
   BoxPosition=0;


   for (G4int i=0;i<NElements;i++)
  {
    ElementFullHeight[i]=0;
    ElementFullYWidth[i]=0;
    ElementZPosition[i]=0;
    ElementYPosition[i]=0;
  }




  
  FilterMessenger = new  HadrontherapyRidgeFilterMessenger(this);

  G4String Name="RidgeFilters/default_RF_name.txt";
  File.open(Name,  std::ios::in);
  if(!File.is_open())
  {
  G4cout<<" WARNING: The File with name of "<<Name<<
" doesn't exist to get RidgeFilter properties. please modify it and try again"<<G4endl;

  ModulatorDefaultProperties();

  }else{
      G4String RFName;
      File >>  RFName;
      File.close();
      G4String RFName2="RidgeFilters/"+RFName;
      File.open(RFName2,  std::ios::in);
          if(!File.is_open())
          {
          G4cout<<" WARNING: The File with name of "<<Name<<
        " doesn't exist to get RidgeFilter properties. please modify it and try again"<<G4endl;
          ModulatorDefaultProperties();
          }else{
           File.close();
  ModulatorReadPropertiesFromFile(RFName2);
          }
  }


  rm = new G4RotationMatrix();
  G4double phi = 90. *deg;
  rm -> rotateY(phi);
}
/////////////////////////////////////////////////////////////////////////////
HadrontherapyRidgeFilter::~HadrontherapyRidgeFilter()
{
  delete rm;

  delete [] ElementFullHeight;
  delete [] ElementFullYWidth;
  delete [] ElementZPosition;
  delete [] ElementYPosition;

    delete solidRFboxMotherMod;
    delete logicRFboxMotherMod;
    delete physiRFboxMotherMod;

    delete  RFstairBox;
    delete  logicRFstairBox;
    delete  chamberParam;
    delete  physiRFstairBox;


  delete FilterMessenger;
}

/////////////////////////////////////////////////////////////////////////////
void HadrontherapyRidgeFilter::ModulatorDefaultProperties()
{
/* Here we initialize the step properties of Modulator wheel, you can create your
specific modulator by changing the values in this class or writing them in an external
file and activate reading from file via a macrofile.	*/	

    for (G4int i=0;i<NModules;i++)
   {

    ElementFullHeight[0+i*NStairs]=10*mm;ElementFullYWidth[0+i*NStairs]=1.6*mm;
    ElementFullHeight[1+i*NStairs]=20*mm;ElementFullYWidth[1+i*NStairs]=1.6*mm;
    ElementFullHeight[2+i*NStairs]=30*mm;ElementFullYWidth[2+i*NStairs]=1.6*mm;
    ElementFullHeight[3+i*NStairs]=40*mm;ElementFullYWidth[3+i*NStairs]=1.6*mm;
    ElementFullHeight[4+i*NStairs]=50*mm;ElementFullYWidth[4+i*NStairs]=1.6*mm;
    }

 GetStepInformation();
 
} 


/////////////////////////////////////////////////////////////////////////////
void HadrontherapyRidgeFilter::ModulatorReadPropertiesFromFile(G4String Name)
{
    delete [] ElementFullHeight;
    delete [] ElementFullYWidth;
    delete [] ElementYPosition;
    delete [] ElementZPosition;
          File.open(Name,  std::ios::in);
          if(!File.is_open())
          {
          G4cout<<" WARNING: The File with name of "<<Name<<
       " doesn't exist to get modulator step properties. please modify it and try again"<<G4endl;
        G4Exception("HadrontherapyRidgeFilter::RidgeFilterPropertiesFromFile( )", "Hadrontherapy0009"
         , FatalException, "Error: No available external file for reading from");
            }
                  G4String string;
                  File >> string>> FullXWidth;
                  File >> string>> NModules;
                  File >> string>> NStairs;
                  File >> string>> ModulesYShift;
                  File >>string;
                  File >>string>>string>>string;

                  NElements=NModules*NStairs;
                  ElementFullHeight=new G4double[NElements];
                  ElementFullYWidth=new G4double[NElements];
                  ElementYPosition=new G4double[NElements];
                  ElementZPosition=new G4double[NElements];
           //        G4double ElementFullYWidthI=0;
           //        G4double ElementFullHeightI=0;
                  for(G4int i=0;i<NStairs;++i)
                        {  G4String stringX;
                           File>>stringX>> ElementFullYWidth[i]>> ElementFullHeight[i];
           //                 ElementFullYWidthI=ElementFullYWidth[i];
           //                 ElementFullHeightI=ElementFullHeight[i];
                         }
          File.close();

                  for (G4int j=1;j<NModules;++j)
                 {    for(G4int i=0;i<NStairs;++i)
                         {  ElementFullHeight[i+j*NStairs]=ElementFullHeight[i];
                          ElementFullYWidth[i+j*NStairs]=ElementFullYWidth[i];
                         }
                 }
           GetStepInformation();
}
/////////////////////////////////////////////////////////////////////////////
void HadrontherapyRidgeFilter:: BuildModulatorFromFile(G4String Name)
{

//delete solidMotherMod;   нельзя удалять
//delete logicMotherMod;
//delete physiMotherMod;

delete solidRFboxMotherMod;
delete logicRFboxMotherMod;
delete physiRFboxMotherMod;

delete  RFstairBox;
delete  logicRFstairBox;
delete  chamberParam;
//delete  physiRFstairBox; //дает ошибку

    ModulatorReadPropertiesFromFile(Name);




       //BuildFilter(motherVolume, BoxPosX);
       // BoxFullWidth, BoxFullWidth, 0.5*BoxFullHeigth+1*mm
       solidMotherMod ->SetXHalfLength(BoxFullWidth);
       solidMotherMod ->SetYHalfLength(BoxFullWidth);
       solidMotherMod ->SetZHalfLength( 0.5*BoxFullHeigth+1*mm);
       BuildRFstairs();
}
////////////////////////////////////////////////////////////////////////////////
void HadrontherapyRidgeFilter::GetStepInformation()
{


 G4double StairMaxHeight=0;
 G4double ModuleYWidth=0*cm;
 for (G4int i=0;i<NStairs;++i)
 {
     ModuleYWidth=ModuleYWidth+ElementFullYWidth[i];
     if(ElementFullHeight[i]>StairMaxHeight) {StairMaxHeight=ElementFullHeight[i];}
 }
  BoxFullHeigth=StairMaxHeight;
  BoxFullWidth=(ModulesYShift+ModuleYWidth)*NModules;

  G4double Starting_point=-0.5*BoxFullWidth;
  G4int k=0;
  for (G4int i=0;i<NModules;++i)
  {
        ElementZPosition[k]= -0.5*(ElementFullHeight[0]-BoxFullHeigth); //-0.5*(BoxFullHeigth+ElementFullHeight[0]);           // Сдвиг каждой ступени по оси Z
        ElementYPosition[k]= Starting_point+ModulesYShift+(ModuleYWidth+ModulesYShift)*i+ 0.5*ElementFullYWidth[0];           // Сдвиг каждой ступени по оси Y
        k=k+1;
      for (G4int j=1;j<NStairs;++j)
      {

        //  ElementFullHeight= new G4double[NElements];          // Полная высота каждой ступени
        //  ElementFullYWidth= new G4double[NElements];          // Полная ширина каждой ступени вдоль оси Y
          ElementZPosition[k]= -0.5*(ElementFullHeight[j]-BoxFullHeigth);   //-0.5*(BoxFullHeigth+ElementFullHeight[j]);           // Сдвиг каждой ступени по оси Z
          ElementYPosition[k]= ElementYPosition[k-1]+0.5*ElementFullYWidth[j-1]+0.5*ElementFullYWidth[j];           // Сдвиг каждой ступени по оси Y
         k=k+1;
      }
  }

//  G4double NStairs0 = NStairs;
//  G4double NModules0 = NModules;
//  G4double BoxFullHeigth0 =BoxFullHeigth;
//  G4double BoxFullWidth0=BoxFullWidth;
//  G4double* ElementPosZ0=ElementZPosition;
//  G4double* ElementPosY0=ElementYPosition;
//  G4double* ElementFullYWidth0=ElementFullYWidth;
//  G4double* ElementFullHeight0=ElementFullHeight;
//  G4double emptyline=0*cm;





	
}
/////////////////////////////////////////////////////////////////////////////////
void HadrontherapyRidgeFilter::BuildFilter(G4VPhysicalVolume* motherVolume, G4double BoxPosX)
{
  G4bool isotopes = false;
  G4Material* airNist =  G4NistManager::Instance()->FindOrBuildMaterial("G4_AIR", isotopes);
  //G4Material* copperNistAsMaterial = G4NistManager::Instance()->FindOrBuildMaterial("G4_Cu", isotopes);
  G4Material* PMMANist = G4NistManager::Instance()->FindOrBuildMaterial("G4_PLEXIGLASS", isotopes);

   Mod0Mater = airNist;
   ModMater = airNist; // You have to change modulator material via a macrofile (default is air)
   Mod1Mater= PMMANist;

  // Mother volume
  G4ThreeVector positionMotherMod = G4ThreeVector(BoxPosX, 0 *mm, 0 *mm);

  G4double ZSz=BoxFullHeigth;
  //G4Box* solidMotherMod = new G4Box("MotherMod", 10.0 *cm, 10.0 *cm, 0.5*ZSz+1*mm);
  solidMotherMod = new G4Box("solidMotherMod", BoxFullWidth, BoxFullWidth, 0.5*ZSz+1*mm);



 
  logicMotherMod = new G4LogicalVolume(solidMotherMod, Mod0Mater,"logicMotherMod");//,0,0,0);

  physiMotherMod = new G4PVPlacement(rm,positionMotherMod,  "physiMotherMod",
				     logicMotherMod,    				  
				     motherVolume,      
				     false,           
                     0);


    BuildRFstairs();



				 }            
 ///////////////////////////////////////////////////////////////////////////////////////
 void HadrontherapyRidgeFilter::BuildRFstairs()
 {


// Объем для параметаризации ступеней
//==========================================================
     solidRFboxMotherMod = new G4Box("RFMotherMod", 0.5*BoxFullWidth, 0.5*BoxFullWidth, 0.5*BoxFullHeigth);

     logicRFboxMotherMod = new G4LogicalVolume(solidRFboxMotherMod, Mod0Mater,"logRFMotherMod");//,0,0,0);

     physiRFboxMotherMod = new G4PVPlacement(0 ,G4ThreeVector(0., 0.,BoxPosition),  "physRFMotherMod",
                        logicRFboxMotherMod,
                        physiMotherMod,
                        false,
                        0);




     RFstairBox
       = new G4Box("RFstairBox",1.*mm,1.*mm,1.*mm);
                   //,0, 1*mm, 1*m, 0.*deg, 360.*deg);
     logicRFstairBox
       = new G4LogicalVolume(RFstairBox,Mod1Mater,"RFstairBox");//,0,0,0);

     chamberParam
       = new BoxesParameterisation(  NElements,
                                     ElementYPosition,
                                     ElementZPosition,
                                     FullXWidth,
                                     ElementFullYWidth,
                                     ElementFullHeight);


     physiRFstairBox
       =  new G4PVParameterised("RFstairBox",       // their name
                                           logicRFstairBox,   // their logical volume
                                           logicRFboxMotherMod ,       // Mother logical volume
                                           kZAxis,          // Are placed along this axis
                                           NElements,    // Number of chambers
                                           chamberParam,    // The parametrisation
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

  logicRFboxMotherMod -> SetVisAttributes(EmptyColor);// -> SetVisAttributes (G4VisAttributes::Invisible);//
  logicMotherMod->SetVisAttributes(EmptyColor); // SetVisAttributes (G4VisAttributes::Invisible);// ->

  logicRFstairBox -> SetVisAttributes(skyBlue);//EmptyColor);//

//  logicPMMABoxBottomLayer-> SetVisAttributes(skyBlue);//red);//EmptyColor);//
//  logicContScattPMMA -> SetVisAttributes(skyBlue);//EmptyColor);//
//  logicAirRing-> SetVisAttributes(skyBlue);//darkOrange3);//EmptyColor);//
//  logicContScattAirTube-> SetVisAttributes (G4VisAttributes::Invisible);
 }





/////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////
// Messenger values
/////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////
void HadrontherapyRidgeFilter::SetFilterAngle(G4double angle)
{
  G4double rotationAngle = angle;
  rm -> rotateY(rotationAngle);
  physiMotherMod -> SetRotation(rm);  
  G4cout << "MODULATOR HAS BEEN ROTATED OF " << rotationAngle/deg 
	 << " deg" << G4endl;
  G4RunManager::GetRunManager() -> GeometryHasBeenModified(); 
  G4RunManager::GetRunManager() -> PhysicsHasBeenModified();
}
/////////////////////////////////////////////////////////////////////////
// Change modulator material
void HadrontherapyRidgeFilter::SetFilterMaterial(G4String Material)
{
    if (G4Material* NewMaterial = G4NistManager::Instance()->FindOrBuildMaterial(Material, false) )
    {
	if (NewMaterial) 
	{
        for(G4int i=0;i<NElements;++i)
	    {
        logicRFstairBox -> SetMaterial(NewMaterial);
	  //  G4RunManager::GetRunManager() -> PhysicsHasBeenModified();
	    G4RunManager::GetRunManager() -> GeometryHasBeenModified();
        G4RunManager::GetRunManager() -> PhysicsHasBeenModified();
	    
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

////////////////////////////////////////////////////////////////////////////////
// Change modulator position in the beam line
void HadrontherapyRidgeFilter::SetFilterPosition(G4ThreeVector Pos)
{
  G4ThreeVector NewModulatorPos=Pos; 
  physiMotherMod -> SetTranslation( NewModulatorPos); 
  G4RunManager::GetRunManager() -> GeometryHasBeenModified();
  G4RunManager::GetRunManager() -> PhysicsHasBeenModified();
  G4cout << "The RidgeFilter is translated to"<<  NewModulatorPos/mm <<"mm " <<G4endl;
	
}

/////////////////////////////////////////////////////////////////////////////////
void HadrontherapyRidgeFilter:: GetDataFromFile(G4String value)

{
G4String Name=value;
if(value=="default" )	
{
Name=FileName;
}
G4cout<<" Step properties of modulator will be get out from the external file "
 <<Name<<G4endl;
BuildModulatorFromFile(Name);
//G4RunManager::GetRunManager() -> GeometryHasBeenModified();
}
