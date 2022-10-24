//
// ********************************************************************
// * License and Disclaimer                                           *
// *                                                                  *
// * The  Geant4 software  is  copyright of the Copyright Holders  of *
// * the Geant4 Collaboration.  It is provided  under  the terms  and *
// * conditions of the Geant4 Software License,  included in the file *
// * LICENSE and available at  http://cern.ch/geant4/license .  These *
// * include a list of copyright holders.                             *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.  Please see the license in the file  LICENSE  and URL above *
// * for the full disclaimer and the limitation of liability.         *
// *                                                                  *
// * This  code  implementation is the result of  the  scientific and *
// * technical work of the GEANT4 collaboration.                      *
// * By using,  copying,  modifying or  distributing the software (or *
// * any work based  on the software)  you  agree  to acknowledge its *
// * use  in  resulting  scientific  publications,  and indicate your *
// * acceptance of all terms of the Geant4 Software license.          *
// ********************************************************************
//
// This is the *BASIC* version of Hadrontherapy, a Geant4-based application
// See more at: http://g4advancedexamples.lngs.infn.it/Examples/hadrontherapy
//
// Visit the Hadrontherapy web site (http://www.lns.infn.it/link/Hadrontherapy) to request 
// the *COMPLETE* version of this program, together with its documentation;
// Hadrontherapy (both basic and full version) are supported by the Italian INFN
// Institute in the framework of the MC-INFN Group
//


#include "HadrontherapyRidgeFilterMessenger.hh"
#include "HadrontherapyRidgeFilter.hh"
#include "G4UIdirectory.hh"
#include "G4UIcmdWithADoubleAndUnit.hh"
#include "G4UIcmdWithAString.hh"
#include "G4UIcmdWith3VectorAndUnit.hh"


    HadrontherapyRidgeFilterMessenger::HadrontherapyRidgeFilterMessenger(HadrontherapyRidgeFilter* Mod)
:RidgeFilter(Mod)

{
    filterDir = new G4UIdirectory("/RidgeFilter/");
    filterDir -> SetGuidance("Command to change the RidgeFilter proprties");

    filterMat1Cmd = new G4UIcmdWithAString("/RidgeFilter/ChangeFilterMaaterial",this);
    filterMat1Cmd -> SetGuidance("Set material of RidgeFilter");
    filterMat1Cmd -> SetParameterName("Material",false);
    filterMat1Cmd -> AvailableForStates(G4State_Idle);


    
    filterExternalFile=new G4UIcmdWithAString("/RidgeFilter/ReadData",this);
    filterExternalFile -> SetGuidance("set properties of RidgeFilter steps via a external file");
    filterExternalFile -> SetParameterName("FileName",true,false);
    filterExternalFile -> SetDefaultValue ("default");
    filterExternalFile -> AvailableForStates(G4State_Idle);

    filterPositionCmd = new G4UIcmdWith3VectorAndUnit("/RidgeFilter/position",this);
    filterPositionCmd -> SetGuidance("Set position of filter");
    filterPositionCmd -> SetParameterName("PositionAlongX",
						    "PositionAlongY", 
						    "PositionAlongZ",false);
    filterPositionCmd -> SetDefaultUnit("mm");
    filterPositionCmd -> SetUnitCandidates("mm cm m");
    filterPositionCmd -> AvailableForStates(G4State_Idle);


    filterAngleCmd = new G4UIcmdWithADoubleAndUnit("/RidgeFilter/angle",this);
    filterAngleCmd -> SetGuidance("Set filter Angle");
    filterAngleCmd -> SetParameterName("Size",false);
    filterAngleCmd -> SetRange("Size>=0.");
    filterAngleCmd -> SetUnitCategory("Angle");
    filterAngleCmd -> AvailableForStates(G4State_Idle);
}

 HadrontherapyRidgeFilterMessenger::~ HadrontherapyRidgeFilterMessenger()
{ 
    delete filterAngleCmd;
    delete filterMat1Cmd;
    delete filterPositionCmd;
    delete filterDir;
}




void  HadrontherapyRidgeFilterMessenger::SetNewValue(G4UIcommand* command,G4String newValue)
{ 
    if( command == filterAngleCmd )
    {RidgeFilter -> SetFilterAngle
     (filterAngleCmd -> GetNewDoubleValue(newValue));}

   else if( command == filterMat1Cmd )
   {RidgeFilter -> SetFilterMaterial(newValue);}
   
    else if (command== filterExternalFile)
     {RidgeFilter->GetDataFromFile(newValue);}

   else if( command == filterPositionCmd )
    { G4ThreeVector size = filterPositionCmd-> GetNew3VectorValue(newValue);
    RidgeFilter-> SetFilterPosition(size);}



}

