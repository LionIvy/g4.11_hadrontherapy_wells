

#ifndef INRPassiveProtonBeamLineMessenger_h
#define INRPassiveProtonBeamLineMessenger_h 1

#include "globals.hh"
#include "G4UImessenger.hh"


class INRPassiveProtonBeamLine;
class G4UIdirectory;
class G4UIcmdWithADoubleAndUnit;
class G4UIcmdWithAString;

class INRPassiveProtonBeamLineMessenger: public G4UImessenger
{
  public:
  INRPassiveProtonBeamLineMessenger(INRPassiveProtonBeamLine*);
  ~INRPassiveProtonBeamLineMessenger();

    void SetNewValue(G4UIcommand*, G4String);

private:

  // Pointer to the detector component
  INRPassiveProtonBeamLine* passiveProton;

  G4UIdirectory*       changeTheBeamLineDir;
  G4UIcmdWithAString*  changeTheBeamLineNameCmd; // Control the name of the beam line

  G4UIdirectory*       beamLineDir;  // Control of the beam line

  G4UIdirectory*       firstScatteringFoilDir;
  // Control of the first scattering foil component of the beam line

  G4UIcmdWithADoubleAndUnit* firstScatteringFoilXSizeCmd;
  // UI command to set half X size of the first scattering foil of
  // the beam line


};
#endif

