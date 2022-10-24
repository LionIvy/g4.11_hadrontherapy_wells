#ifndef PHASESPACE_SD_HH
#define PHASESPACE_SD_HH

#include "G4VSensitiveDetector.hh"
#include "G4SystemOfUnits.hh"
#include "INRPassiveProtonBeamLine.hh"

#include "G4AnalysisManager.hh"
#include "G4RunManager.hh"
#include "PhaseSpaceDataset.hh"

#include "HadrontherapyDetectorHit.hh"

typedef std::vector<G4double> G4double_vector;
typedef std::vector<G4int> G4int_2D_vector;
typedef std::vector<std::vector<G4int>> G4int_3D_vector;

class PhaseSpace_SD : public G4VSensitiveDetector
{
public:
    PhaseSpace_SD(G4String name);
    ~PhaseSpace_SD();
    void Initialize(G4HCofThisEvent*);
    virtual G4bool ProcessHits(G4Step *, G4TouchableHistory *);
private:
    HadrontherapyDetectorHitsCollection *HitsCollection;
    G4String sensitiveDetectorName;

    G4LogicalVolume* fScoringVolume = nullptr;

    G4int evnt;

    G4Track * track;
    G4StepPoint *preStepPoint;
    G4int X, Y;

    const G4VTouchable *touchable;
    G4int detectorCopyNo;

    G4ParticleDefinition *particleDef;
    G4int Z, A;
    G4double energy;
    G4ThreeVector momDirection;

    G4AnalysisManager* analysis;

    int Ny, Nz;
    G4int_2D_vector yz_matrix;

    int  Nenergy;
    int Nangle_yz;
    int Nangle_xz;





    //G4StepPoint *postStepPoint;
};

#endif // PHASESPACE_SD_HH
