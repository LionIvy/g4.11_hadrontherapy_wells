#include "G4RunManager.hh"
#include "G4UImanager.hh"
#include "G4PhysListFactory.hh"
#include "G4VModularPhysicsList.hh"
#include "HadrontherapyEventAction.hh"
#include "HadrontherapyPhysicsList.hh"
#include "HadrontherapyDetectorSD.hh"
#include "HadrontherapyPrimaryGeneratorAction.hh"
#include "HadrontherapyRunAction.hh"
#include "HadrontherapyMatrix.hh"
#include "Randomize.hh"

#include "G4UImessenger.hh"
#include "globals.hh"
#include "HadrontherapySteppingAction.hh"
#include "HadrontherapyGeometryController.hh"
#include "HadrontherapyGeometryMessenger.hh"
#include "HadrontherapyInteractionParameters.hh"
#include "HadrontherapyLet.hh"

#include "G4ScoringManager.hh"
#include "G4ParallelWorldPhysics.hh"
#include <time.h>
#include "G4Timer.hh"
#include "G4RunManagerFactory.hh"
#include "HadrontherapyActionInitialization.hh"

#include "G4VisExecutive.hh"
#include "G4UIExecutive.hh"

#include "PhaseSpaceDataset.hh"

//////////////////////////////////////////////////////////////////////////////////////////////
int main(int argc ,char ** argv)
{
        G4UIExecutive* ui = 0;
    if ( argc == 1 ) {
        ui = new G4UIExecutive(argc, argv);
    }
    
    //Instantiate the G4Timer object, to monitor the CPU time spent for
    //the entire execution
    G4Timer* theTimer = new G4Timer();
    //Start the benchmark
    theTimer->Start();
    
    // Set the Random engine
    // The following guarantees random generation also for different runs
    // in multithread
    CLHEP::RanluxEngine defaultEngine( 1234567, 4 );
    G4Random::setTheEngine( &defaultEngine );
    G4int seed = (G4int) time( NULL );
    G4Random::setTheSeed( seed );
 
 auto* runManager = G4RunManagerFactory::CreateRunManager();
 G4int nThreads = 4;
 runManager->SetNumberOfThreads(nThreads); 

    // Geometry controller is responsible for instantiating the
    // geometries. All geometry specific m tasks are now in class
    // HadrontherapyGeometryController.
    HadrontherapyGeometryController *geometryController = new HadrontherapyGeometryController();
    
    // Connect the geometry controller to the G4 user interface
    HadrontherapyGeometryMessenger *geometryMessenger = new HadrontherapyGeometryMessenger(geometryController);
    
    G4ScoringManager *scoringManager = G4ScoringManager::GetScoringManager();
    scoringManager->SetVerboseLevel(1);
    
    // Initialize the default Hadrontherapy geometry
    geometryController->SetGeometry("default");
    
    // Initialize the physics
    G4PhysListFactory factory;
    G4VModularPhysicsList* phys = 0;
    G4String physName = "";
    
    // Physics List name defined via environment variable
    char* path = std::getenv("PHYSLIST");
    if (path) { physName = G4String(path); }
    
    if(physName != "" && factory.IsReferencePhysList(physName))
    {
        phys = factory.GetReferencePhysList(physName);
    }
    if (phys)
    {
        G4cout << "Going to register G4ParallelWorldPhysics" << G4endl;
        phys->RegisterPhysics(new G4ParallelWorldPhysics("DetectorROGeometry"));
    }
    else
    {
        G4cout << "Using HadrontherapyPhysicsList()" << G4endl;
        phys = new HadrontherapyPhysicsList();
    }
    
    // Initialisations of physics
    runManager->SetUserInitialization(phys);
    
    // Initialisation of the Actions
    runManager->SetUserInitialization(new HadrontherapyActionInitialization);
    
    // Initialize command based scoring
    G4ScoringManager::GetScoringManager();
    
    // Interaction data: stopping powers
    HadrontherapyInteractionParameters* pInteraction = new HadrontherapyInteractionParameters(true);
    
    // Initialize analysis
    HadrontherapyAnalysis* analysis = HadrontherapyAnalysis::GetInstance();
    

    
// Initialise the Visualisation
    G4VisManager* visManager = new G4VisExecutive;
    visManager -> Initialize();
    
    //** Get the pointer to the User Interface manager
    G4UImanager* UImanager = G4UImanager::GetUIpointer();
    
    if ( !ui ) {
        // batch mode
        G4String command = "/control/execute ";
        G4String fileName = argv[1];
        UImanager->ApplyCommand(command+fileName);
        
    }
    
    else {

       // UImanager -> ApplyCommand("/control/execute macro/defaultMacro.mac");
        UImanager -> ApplyCommand("/control/execute macro/INR_default.mac");
        ui -> SessionStart();
        delete ui;
    }
    delete visManager;
 
    //Stop the benchmark here
    theTimer->Stop();
    
    G4cout << "The simulation took: " << theTimer->GetRealElapsed() << " s to run (real time)"
    << G4endl;
    
    // Job termination
    // Free the store: user actions, physics_list and detector_description are
    // owned and deleted by the run manager, so they should not be deleted
    // in the main() program !
    
    
        if ( HadrontherapyMatrix * pMatrix = HadrontherapyMatrix::GetInstance() )
    {
        // pMatrix -> TotalEnergyDeposit();
        pMatrix -> StoreDoseFluenceAscii();
        
    }
    
    if (HadrontherapyLet *let = HadrontherapyLet::GetInstance())
        if(let -> doCalculation)
        {
            let -> LetOutput(); 	// Calculate let
            let -> StoreLetAscii(); // Store it
        }

    G4cout << "main" << G4endl;
//    if(PhaseSpaceDataset* phsp_data = PhaseSpaceDataset::getInstance()){
//        //phsp_data
//        G4cout << "dataset found";
//    }
    
    delete geometryMessenger;
    delete geometryController;
    delete pInteraction;
    delete runManager;
    delete analysis;
    return 0;
    
}
