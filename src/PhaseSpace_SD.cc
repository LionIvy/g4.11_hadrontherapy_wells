#include "PhaseSpace_SD.hh"
#include "G4AnalysisManager.hh"
#include "PhaseSpaceDataset.hh"


PhaseSpace_SD::PhaseSpace_SD(G4String name):
    G4VSensitiveDetector(name)
{
        analysis = G4AnalysisManager::Instance();
        G4String HCname;
        collectionName.insert(HCname="PhaseSpaceDetectorHitsCollection");
        HitsCollection = NULL;
        sensitiveDetectorName = name;
}

PhaseSpace_SD::~PhaseSpace_SD()
{

}


void PhaseSpace_SD::Initialize(G4HCofThisEvent*)
{
    HitsCollection = new HadrontherapyDetectorHitsCollection(sensitiveDetectorName,
                                                             collectionName[0]);

    Ny = 100;
    Nz = 100;
    G4cout << "PhaseSpace_SD init" << G4endl;
   // yz_matrix.

//    G4double min_energy = 150*MeV;
//    G4double max_energy = 160*MeV;
//    Nenergy = 80;

    //Nangle_yz = 100;
    //Nangle_xz = 100;


}

G4bool PhaseSpace_SD::ProcessHits(G4Step *aStep, G4TouchableHistory *ROhist)
{
    preStepPoint = aStep->GetPreStepPoint();
//==================================================================================
// Отсев по детектору
    if (preStepPoint -> GetPhysicalVolume() -> GetName() != "PhaseSpacePhysVol") return false;
//    if (preStepPoint -> GetPhysicalVolume() -> GetName() != "PhaseSpacePhysVol"){
//        track->SetTrackStatus(fStopAndKill);
//        return false;
//    }
//==================================================================================
// Отсев по частицам
    track = aStep  ->  GetTrack();
    particleDef = track -> GetDefinition();
    Z = particleDef-> GetAtomicNumber();
    A = particleDef-> GetAtomicMass();

    if (!(A == 1 && Z ==1)) return false; //только протоны
//==================================================================================
//
    evnt = G4RunManager::GetRunManager()->GetCurrentEvent()->GetEventID();
    touchable = preStepPoint->GetTouchable();
    detectorCopyNo = touchable->GetCopyNumber();

    G4cout << "PhaseSpace_SD progress hit" << G4endl;
//    PhaseSpaceDataset* phsp_data = PhaseSpaceDataset::getInstance();
//    phsp_data->fill_energy_dataset(detectorCopyNo, energy);
//    phsp_data->fill_fluence_dataset(detectorCopyNo);

    //int Nx = 100;
    Ny = 100;
    X = detectorCopyNo / Ny;
    Y = detectorCopyNo % Ny;

    momDirection =  preStepPoint->GetMomentumDirection();
    energy =  preStepPoint->GetKineticEnergy()/MeV;


//==================================================================================
//

    analysis->FillNtupleIColumn(0, evnt);
    analysis->FillNtupleIColumn(1, detectorCopyNo);
    analysis->FillNtupleIColumn(2, X);
    analysis->FillNtupleIColumn(3, Y);
    analysis->FillNtupleDColumn(4, energy);
    analysis->FillNtupleDColumn(5, momDirection[0]);
    analysis->FillNtupleDColumn(6, momDirection[1]);
    analysis->FillNtupleDColumn(7, momDirection[2]);


    analysis->FillH1(0, energy);

    analysis->AddNtupleRow(0);



    return true;
}
