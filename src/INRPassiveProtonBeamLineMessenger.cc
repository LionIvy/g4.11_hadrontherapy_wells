

#include "INRPassiveProtonBeamLineMessenger.hh"
#include "INRPassiveProtonBeamLine.hh"
#include "G4UIdirectory.hh"
#include "G4UIcmdWithADoubleAndUnit.hh"
#include "G4UIcmdWithAString.hh"
#include "G4SystemOfUnits.hh"

    INRPassiveProtonBeamLineMessenger::INRPassiveProtonBeamLineMessenger(INRPassiveProtonBeamLine* beamLine)
:passiveProton(beamLine)

{
    changeTheBeamLineDir = new G4UIdirectory("/ChangeBeamLine/");
    changeTheBeamLineDir -> SetGuidance("Command to change the transport beam line");

    changeTheBeamLineNameCmd = new G4UIcmdWithAString("/ChangeBeamLine/beamLineName",this);
    changeTheBeamLineNameCmd -> SetGuidance("Insert the name of the beam line you want simulate");
    changeTheBeamLineNameCmd -> SetParameterName("List",false);
    changeTheBeamLineNameCmd -> AvailableForStates(G4State_PreInit);

    beamLineDir = new G4UIdirectory("/beamLine/");
    beamLineDir -> SetGuidance("set specification of range shifter");

    firstScatteringFoilDir = new G4UIdirectory("/beamLine/ScatteringFoil1/");
    firstScatteringFoilDir -> SetGuidance("set specification of first scattering foil");

    firstScatteringFoilXSizeCmd = new G4UIcmdWithADoubleAndUnit("/beamLine/ScatteringFoil1/thickness",this);
    firstScatteringFoilXSizeCmd -> SetGuidance("Set half thickness of first scattering foil");
    firstScatteringFoilXSizeCmd -> SetParameterName("Size",false);
    firstScatteringFoilXSizeCmd -> SetDefaultUnit("mm");
    firstScatteringFoilXSizeCmd -> SetUnitCandidates("mm cm m");
    firstScatteringFoilXSizeCmd -> AvailableForStates(G4State_Idle);

}

INRPassiveProtonBeamLineMessenger::~INRPassiveProtonBeamLineMessenger()
{

    delete firstScatteringFoilXSizeCmd;
    delete firstScatteringFoilDir;
    delete beamLineDir;

}




void INRPassiveProtonBeamLineMessenger::SetNewValue(G4UIcommand* command,G4String newValue)
{


if( command == firstScatteringFoilXSizeCmd )
    { passiveProton -> SetFirstScatteringFoilXSize
    (firstScatteringFoilXSizeCmd -> GetNewDoubleValue(newValue));}

}

