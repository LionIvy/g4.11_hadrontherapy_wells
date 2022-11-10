
#include "globals.hh"
#include "G4SystemOfUnits.hh"

#include "G4Box.hh"
#include "G4Tubs.hh"
#include "G4Trap.hh"
#include "G4Trd.hh"
#include "G4Cons.hh"
#include "G4SubtractionSolid.hh"
#include "G4VisAttributes.hh"
#include "G4Colour.hh"

#include "G4RunManager.hh"
#include "G4LogicalVolume.hh"
#include "G4PVPlacement.hh"
#include "G4RotationMatrix.hh"
#include "G4NistManager.hh"
#include "G4NistElementBuilder.hh"

#include "HadrontherapyDetectorConstruction.hh"


#include "INRPassiveProtonBeamLine.hh"
#include "INRPassiveProtonBeamLineMessenger.hh"
#include "TubesParameterisation.hh"
#include "G4PVParameterised.hh"
#include "HadrontherapyContouredScatterer.hh"
#include "HadrontherapyRidgeFilter.hh"

#include "G4Region.hh"
#include "G4RegionStore.hh"
#include "G4ProductionCuts.hh"

// под удаление
#include "HadrontherapyModulator.hh"

//#include "PhaseSpaceDetector.hh"
//#include "PhaseSpace_SD.hh"
//#include "PhaseSpaceDataset.hh"
//#include "PhaseSpaceDetectorConstruction.hh"



//G4bool PassiveProtonBeamLine::doCalculation = false;
/////////////////////////////////////////////////////////////////////////////
INRPassiveProtonBeamLine::INRPassiveProtonBeamLine():
    ContouredScatterer(),
    RidgeFilter(),
    physicalTreatmentRoom(),     // физ объем мира
    hadrontherapyDetectorConstruction(),      // конструкция детектора
 //   phsp_det(),
    firstScatteringFoil(),
    logicConcreteWall(),
    physiGraphiteCollimatorBox(),
    physiFirstScatteringFoil(),
    physiConcreteWall()

{

    // Messenger to change parameters of the passiveProtonBeamLine geometry
    passiveMessenger = new INRPassiveProtonBeamLineMessenger(this);

    //***************************** PW ***************************************
    static G4String ROGeometryName = "DetectorROGeometry";
    RO = new HadrontherapyDetectorROGeometry(ROGeometryName);

    G4cout << "Going to register Parallel world...";
    RegisterParallelWorld(RO);
    G4cout << "... done" << G4endl;

}
/////////////////////////////////////////////////////////////////////////////
INRPassiveProtonBeamLine::~INRPassiveProtonBeamLine()
{
    delete passiveMessenger;
    delete hadrontherapyDetectorConstruction;

}

/////////////////////////////////////////////////////////////////////////////
G4VPhysicalVolume* INRPassiveProtonBeamLine::Construct()
{
    // Sets default geometry and materials
    SetDefaultDimensions();

    // Construct the whole Passive Beam Line
    ConstructINRProtonBeamLine();

    //***************************** PW ***************************************
    if (!hadrontherapyDetectorConstruction){
    G4bool plasticWells = true;
    hadrontherapyDetectorConstruction = new HadrontherapyDetectorConstruction(physicalTreatmentRoom, plasticWells);
    hadrontherapyDetectorConstruction->InitializeDetectorROGeometry(RO,hadrontherapyDetectorConstruction->GetDetectorToWorldPosition());
    }

    // if (!phsp_det) phsp_det = new PhaseSpaceDetector(physicalTreatmentRoom);

    return physicalTreatmentRoom;
}

// In the following method the DEFAULTS used in the geometry of
// passive beam line are provided
// HERE THE USER CAN CHANGE THE GEOMETRY CHARACTERISTICS OF BEAM
// LINE ELEMENTS, ALTERNATIVELY HE/SHE CAN USE THE MACRO FILE (IF A
// MESSENGER IS PROVIDED)
//
// DEFAULT MATERIAL ARE ALSO PROVIDED
// and COLOURS ARE ALSO DEFINED
// ----------------------------------------------------------


//void INRPassiveProtonBeamLine::ConstructSDandField(){

//    PhaseSpace_SD *sensDet = new PhaseSpace_SD("PhaseSpaceLV");
//    phaseSpace_ScoringVolume ->SetSensitiveDetector(sensDet);


//}

/////////////////////////////////////////////////////////////////////////////
void INRPassiveProtonBeamLine::SetDefaultDimensions()
{
    // Set of coulors that can be used
    white = new G4VisAttributes( G4Colour());
    white -> SetVisibility(true);
    white -> SetForceSolid(true);

    black = new G4VisAttributes( G4Colour(1., 1., 1.));
    black -> SetVisibility(true);
    black -> SetForceSolid(true);


    blue = new G4VisAttributes(G4Colour(0. ,0. ,1.));
    blue -> SetVisibility(true);
    blue -> SetForceSolid(true);

    gray = new G4VisAttributes( G4Colour(0.5, 0.5, 0.5 ));
    gray-> SetVisibility(true);
    gray-> SetForceSolid(true);

    red = new G4VisAttributes(G4Colour(1. ,0. ,0.));
    red-> SetVisibility(true);
    red-> SetForceSolid(true);


    yellow = new G4VisAttributes(G4Colour(1., 1., 0. ));
    yellow-> SetVisibility(true);
    yellow-> SetForceSolid(true);

    green = new G4VisAttributes( G4Colour(25/255. , 255/255. ,  25/255. ));
    green -> SetVisibility(true);
    green -> SetForceSolid(true);

    darkGreen = new G4VisAttributes( G4Colour(0/255. , 100/255. ,  0/255. ));
    darkGreen -> SetVisibility(true);
    darkGreen -> SetForceSolid(true);

    darkOrange3 = new G4VisAttributes( G4Colour(205/255. , 102/255. ,  000/255. ));
    darkOrange3 -> SetVisibility(true);
    darkOrange3 -> SetForceSolid(true);

    skyBlue = new G4VisAttributes( G4Colour(135/255. , 206/255. ,  235/255. ));
    skyBlue -> SetVisibility(true);
    skyBlue -> SetForceSolid(true);

    metalColor = new G4VisAttributes( G4Colour(78/255. , 77/255. ,  74/255. ));
    metalColor -> SetVisibility(true);
    metalColor -> SetForceSolid(true);

    sandColor = new G4VisAttributes( G4Colour(194/255. , 178/255. ,  128/255. ));
    sandColor -> SetVisibility(true);
    sandColor -> SetForceWireframe(true);


    EmptyColor = new G4VisAttributes(G4Colour(1.,1.,1.));
    EmptyColor -> SetVisibility(true);
    EmptyColor -> SetForceWireframe(true);









    //_____________________________________________________________
    // Стена

    G4double defaultWallXSize = 120.0*cm;
    WallXSize = defaultWallXSize;
    G4double defaultWallYSize = 520.0*cm;
    WallYSize = defaultWallYSize;
    G4double defaultWallZSize = 520.0*cm;
    WallZSize = defaultWallZSize;
    //G4double defaultWallXPosition = GraphiteCollimatorXPosition + GraphiteCollimatorXSize + WallXSize+ 50.0 *cm;

    G4double defaultWallAngle=0.3764*rad;
    WallAngle=defaultWallAngle;
    G4double defaultWallXPosition = -0.5*WallXSize/cos(WallAngle);

    WallXPosition = defaultWallXPosition;


    G4double defaultMetalTubeThickness = 3.0*mm;
    MetalTubeThickness = defaultMetalTubeThickness;

    G4double defaultPolyTubeThickness = 57.0*mm;
    PolyTubeThickness = defaultPolyTubeThickness;

    G4double defaultAirWindowR = 5.0*cm;
    AirWindowR = defaultAirWindowR;

    G4double defaultWallWindowR = defaultAirWindowR + 2*defaultMetalTubeThickness + defaultPolyTubeThickness;
    WallWindowR = defaultWallWindowR;

    //_____________________________________________________________
    //              ОБЪЕКТЫ ЗА СТЕНОЙ:
    //_____________________________________________________________

    //_____________________________________________________________
    // Параметры выходного канала

    G4double defaultVacuumZoneRSize = 3.0 *cm;
    vacuumZoneRSize = defaultVacuumZoneRSize;

    G4double defaultVacuumZoneXSize = 50.0 *mm;
    vacuumZoneXSize = defaultVacuumZoneXSize;


    // KAPTON WINDOW: it prmits the passage of the beam from vacuum to air
    G4double defaultKaptonWindowXSize = 0.1045 *cm;//0.005225*mm;//0.5*0.1045 *cm= 0.05225 *cm
    kaptonWindowXSize = defaultKaptonWindowXSize;

    //G4double defaultVacuumZoneXPosition = deltaXPosition-5950.0 *mm;
    G4double defaultVacuumZoneXPosition = -WallXSize-288.50*cm-0.5*defaultVacuumZoneXSize;
    vacuumZoneXPosition = defaultVacuumZoneXPosition;

    //      G4double defaultKaptonWindowXPosition = vacuumZoneXPosition+0.5*vacuumZoneXSize+0.5*defaultKaptonWindowXSize;
    //      kaptonWindowXPosition = defaultKaptonWindowXPosition;



    //_____________________________________________________________
    // Параметры графитового коллиматора
    G4double defaultGraphiteCollimatorXSize = 27.5 *cm;
    GraphiteCollimatorXSize = defaultGraphiteCollimatorXSize;

    G4double defaultGraphiteCollimatorYSize = 15.0 *cm;
    GraphiteCollimatorYSize = defaultGraphiteCollimatorYSize;

    G4double defaultGraphiteCollimatorZSize = 15.0 *cm;
    GraphiteCollimatorZSize = defaultGraphiteCollimatorZSize;

    G4double defaultGraphiteCollimatorXPosition = vacuumZoneXPosition+0.5*vacuumZoneXSize + 55.0*cm +  0.5*GraphiteCollimatorXSize ;
    GraphiteCollimatorXPosition = defaultGraphiteCollimatorXPosition;

    G4double defaultGraphiteCollimatorRSize = 10.0*mm;
    GraphiteCollimatorRSize = defaultGraphiteCollimatorRSize;
    //_____________________________________________________________
    // Параметры первичного рассеивателя

    // FIRST SCATTERING FOIL: a thin foil performing a first scattering
    // of the original beam
    G4double defaultFirstScatteringFoilXSize = 0.10 *mm; // значение полной толщины пластины
    firstScatteringFoilXSize = defaultFirstScatteringFoilXSize;

    G4double defaultFirstScatteringFoilYSize = 52.5   *mm;
    firstScatteringFoilYSize = defaultFirstScatteringFoilYSize;

    G4double defaultFirstScatteringFoilZSize = 52.5   *mm;
    firstScatteringFoilZSize = defaultFirstScatteringFoilZSize;

    G4double defaultFirstScatteringFoilXPosition = GraphiteCollimatorXPosition + 0.5*GraphiteCollimatorXSize + 0.5*firstScatteringFoilXSize + 1.5*cm;
    firstScatteringFoilXPosition = defaultFirstScatteringFoilXPosition;


    //_____________________________________________________________
    //              ОБЪЕКТЫ В ПРОЦЕДУРНОЙ:
    //_____________________________________________________________


    //_____________________________________________________________
    // Фигурный рассеиватель


    G4double defaultContScattBoxPosX = 21.3*cm;//+ 30*cm;
    ContScattBoxPosX = defaultContScattBoxPosX;


    //_____________________________________________________________
    // Гребенчатый фильтр ГФ

    G4double defaultRigdeFilterPosX = defaultContScattBoxPosX+100.5*cm;//+ 30*cm;
    RidgeFilterBoxPosX = defaultRigdeFilterPosX;

    //-------------------------------------------------------------
    // Коллиматор перед ГФ
    G4double defaultRFCollimatorBoxX=7.5*cm;
    RFCollimatorBoxX=defaultRFCollimatorBoxX;

    G4double defaultRFCollimatorBoxY=12.0*cm;
    RFCollimatorBoxY=defaultRFCollimatorBoxY;

    G4double defaultRFCollimatorBoxZ=12.0*cm;;
    RFCollimatorBoxZ=defaultRFCollimatorBoxZ;

    G4double defaultRFCollimatorBoxR=2.6*cm;
    RFCollimatorBoxR=defaultRFCollimatorBoxR;

    G4double defaultRFCollimatorPosX= defaultContScattBoxPosX + 75.0*cm;//RidgeFilterBoxPosX - 20*cm;
    RFCollimatorPosX = defaultRFCollimatorPosX;

    //-------------------------------------------------------------
    // Финальный коллиматор

    G4double defaultFinalCollimatorBoxX=7.0*cm;
    FinalCollimatorBoxX=defaultFinalCollimatorBoxX;

    G4double defaultFinalCollimatorBoxY=11.0*cm;
    FinalCollimatorBoxY=defaultFinalCollimatorBoxY;

    G4double defaultFinalCollimatorBoxZ=11.0*cm;
    FinalCollimatorBoxZ=defaultFinalCollimatorBoxZ;

    G4double defaultFinalCollimatorBoxR=1.5*cm;
    FinalCollimatorBoxR=defaultFinalCollimatorBoxR;

    G4double defaultFinalCollimatorPosX= defaultContScattBoxPosX + 140.0*cm;//RidgeFilterBoxPosX - 20*cm;
    FinalCollimatorPosX = defaultFinalCollimatorPosX;


    //-------------------------------------------------------------
    // Фантом для клеток

    G4double defaultPlasticWaterBoxX=10.0*cm;
    PlasticWaterBoxX=defaultPlasticWaterBoxX;

    G4double defaultPlasticWaterBoxY=5.0*cm;
    PlasticWaterBoxY=defaultPlasticWaterBoxY;

    G4double defaultPlasticWaterBoxZ=2.5*cm;
    PlasticWaterBoxZ=defaultPlasticWaterBoxZ;

    //G4double defaultFinalCollimatorBoxR=1.5*cm;
    //FinalCollimatorBoxR=defaultFinalCollimatorBoxR;

    G4double defaultPlasticWaterBoxPosX= defaultContScattBoxPosX + 170.0*cm;//RidgeFilterBoxPosX - 20*cm;
    PlasticWaterBoxPosX = defaultPlasticWaterBoxPosX;

    G4int defaultEpendorfNumber=1;
    G4double defaultEpendorfRminBottom=3.0*mm;
    G4double defaultEpendorfRmaxBottom=3.5*mm;
    G4double defaultEpendorfRminTop=4.0*mm;
    G4double defaultEpendorfRmaxTop=4.5*mm;
    G4double defaultEpendorfHeight=1.0*cm;

    EpendorfNumber=defaultEpendorfNumber;
    EpendorfRminBottom=defaultEpendorfRminBottom;
    EpendorfRmaxBottom=defaultEpendorfRmaxBottom;
    EpendorfRminTop=defaultEpendorfRminTop;
    EpendorfRmaxTop=defaultEpendorfRmaxTop;
    EpendorfHeight=defaultEpendorfHeight;


    G4double defaultPerem_dx1=1.0*mm;
    // G4double defaultPerem_dx2=defaultPerem_dx1+2.0*(EpendorfRmaxTop-EpendorfRmaxBottom);
    G4double defaultPerem_dy=1.5*mm;
    // G4double defaultPerem_dz=EpendorfHeight;

    Perem_dx1=defaultPerem_dx1;
    Perem_dx2=Perem_dx1+2.0*(EpendorfRmaxTop-EpendorfRmaxBottom);
    Perem_dy=defaultPerem_dy;
    Perem_dz=EpendorfHeight;


    // DEFAULT DEFINITION OF THE MATERIALS
    // All elements and compound definition follows the NIST database

    // ELEMENTS
    G4double A;  // atomic mass
    G4double Z;  // atomic number
    G4bool isotopes = false;
    G4Material* aluminumNist = G4NistManager::Instance()->FindOrBuildMaterial("G4_Al", isotopes);
    //G4Material* tantalumNist = G4NistManager::Instance()->FindOrBuildMaterial("G4_Ta", isotopes);
    G4Material* copperNistAsMaterial = G4NistManager::Instance()->FindOrBuildMaterial("G4_Cu", isotopes);
    G4Element* zincNist = G4NistManager::Instance()->FindOrBuildElement("Zn");
    G4Element* copperNist = G4NistManager::Instance()->FindOrBuildElement("Cu");

    A = 16.00*g/mole;
    G4Element* elO =     new G4Element("Oxygen",   "O",  Z = 8.,A);
    A = 22.99*g/mole;
    G4Element* elNa =    new G4Element("Sodium",   "Na", Z = 11.,A);
    A=26.98*g/mole;
    G4Element* elAl =    new G4Element("Aluminum", "Al", Z = 13.,A);
    A = 28.09*g/mole;
    G4Element* elSi  =   new G4Element("Silicon",  "Si", Z = 14.,A);
    A = 40.08*g/mole;
    G4Element* elCa =    new G4Element("Calcium",  "Ca", Z = 20.,A);
    A = 52.00*g/mole;
    G4Element* elCr  =   new G4Element("Chromium", "Cr", Z = 24.,A);
    A  =  54.94*g/mole;
    G4Element* elMn   =  new G4Element("Manganese","Mn", Z = 25.,A);
    A = 55.85*g/mole;
    G4Element* elFe  =   new G4Element("Iron",     "Fe", Z = 26.,A);
    A = 58.70*g/mole;
    G4Element* elNi  =   new G4Element("Nickel",   "Ni", Z = 28.,A);

    // COMPOUND
    //G4Material* airNist =  G4NistManager::Instance()->FindOrBuildMaterial("G4_AIR", isotopes);
    //G4Material* kaptonNist = G4NistManager::Instance()->FindOrBuildMaterial("G4_KAPTON", isotopes);
    G4Material* galacticNist = G4NistManager::Instance()->FindOrBuildMaterial("G4_Galactic", isotopes);
    G4Material* PMMANist = G4NistManager::Instance()->FindOrBuildMaterial("G4_PLEXIGLASS", isotopes);
    //G4Material* mylarNist = G4NistManager::Instance()->FindOrBuildMaterial("G4_MYLAR", isotopes);
    G4Material* polyethilene = G4NistManager::Instance()->FindOrBuildMaterial("G4_POLYETHYLENE");
    G4Material* lead = G4NistManager::Instance()->FindOrBuildMaterial("G4_LEAD_OXIDE");
    G4Material* water = G4NistManager::Instance()->FindOrBuildMaterial("G4_WATER", false);

    G4double d; // Density
    G4int nComponents;// Number of components
    G4double fractionmass; // Fraction in mass of an element in a material

    d = 8.40*g/cm3;
    nComponents = 2;
    G4Material* brass = new G4Material("Brass", d, nComponents);
    brass -> AddElement(zincNist, fractionmass = 30 *perCent);
    brass -> AddElement(copperNist, fractionmass = 70 *perCent);


    // Stainless steel (Medical Physics, Vol 25, No 10, Oct 1998)
    d = 8.02*g/cm3 ;
    G4Material* matsteel = new G4Material("Stainless steel",d,5);
    matsteel->AddElement(elMn, 0.02);
    matsteel->AddElement(elSi, 0.01);
    matsteel->AddElement(elCr, 0.19);
    matsteel->AddElement(elNi, 0.10);
    matsteel->AddElement(elFe, 0.68);

    // graphite
    G4Isotope* C12 = new G4Isotope("C12", 6, 12);
    G4Element* elC = new G4Element("TS_C_of_Graphite","C", 1);
    elC->AddIsotope(C12, 100.*perCent);
    G4Material* graphite =
            new G4Material("graphite", 2.27*g/cm3, 1,
                           kStateSolid, 293*kelvin, 1*atmosphere);
    graphite->AddElement(elC, 1);

    // Concrete
    d = 2.5*g/cm3;
    G4Material* Concrete = new G4Material("Concrete",d,6);
    Concrete->AddElement(elO, 0.52);
    Concrete->AddElement(elSi, 0.325);
    Concrete->AddElement(elCa, 0.06);
    Concrete->AddElement(elNa, 0.015);
    Concrete->AddElement(elFe, 0.04);
    Concrete->AddElement(elAl, 0.04);
    //***************************** PW ***************************************

    // DetectorROGeometry Material
    new G4Material("dummyMat", 1., 1.*g/mole, 1.*g/cm3);

    //***************************** PW ***************************************



    // MATERIAL ASSIGNMENT

    // Vacuum pipe
    vacuumZoneMaterial = galacticNist;

    // Vacuum pipe
    vacuumZoneWindowMaterial = aluminumNist;

    // Graphite collimator
    GraphiteCollimatorMaterial = graphite;

    // dumpShell
    dumpShellMaterial = lead;
    dumpCoreMaterial = graphite;


    // Material of the fisrt scattering foil
    firstScatteringFoilMaterial = copperNistAsMaterial; //brass;
    //holeFirstCollimatorMaterial = airNist;

    //Wall
    WallMaterial = Concrete;
    MetalTubeMaterial = matsteel;
    PolyTubeMaterial = polyethilene;

    //Фигурный рассеиватель
    //holeModulatorBoxMaterial = airNist;
    ContouredScattMaterial = copperNistAsMaterial;
    ContouredCompensatorMaterial = PMMANist;


    // Коллиматор перед ГФ
    Collimator1Material = lead;

    // Последний коллиматор
    Collimator2Material = lead;

    //Фантом для клеток
    PlasticWaterMaterial=PMMANist;
    EpendorfMaterial=PMMANist;

    EpendorfInsideMaterial=water;

}

/////////////////////////////////////////////////////////////////////////////
void INRPassiveProtonBeamLine::ConstructINRProtonBeamLine()
{
    // -----------------------------
    // Treatment room - World volume
    //------------------------------
    // Treatment room sizes
    const G4double worldX = 700.0 *cm;
    const G4double worldY = 400.0 *cm;
    const G4double worldZ = 400.0 *cm;
    G4bool isotopes = false;

    G4Material* airNist =  G4NistManager::Instance()->FindOrBuildMaterial("G4_AIR", isotopes);
    G4Box* treatmentRoom = new G4Box("TreatmentRoom",worldX,worldY,worldZ);
    G4LogicalVolume* logicTreatmentRoom = new G4LogicalVolume(treatmentRoom,
                                                              airNist,
                                                              "logicTreatmentRoom",
                                                              0,0,0);
    physicalTreatmentRoom = new G4PVPlacement(0,
                                              G4ThreeVector(),
                                              "physicalTreatmentRoom",
                                              logicTreatmentRoom,
                                              0,false,0);


    // The treatment room is invisible in the Visualisation
    logicTreatmentRoom -> SetVisAttributes(EmptyColor);// (G4VisAttributes::Invisible);


    // Components of the Passive Proton Beam Line
    G4String Name="macro/INR_setup_list.mac";
    bool WedgeOn = 0;
    bool GraphiteCollimatorOn = 0;//
    bool FirstScattererOn = 0;//0;//
    bool WallOn = 1;//0;//
    bool ContouredScattererOn = 1;//0;//
    bool Collimator1On = 0;//0;//
    bool RidgeFilterOn = 1;//0;//
    bool Collimator2On = 0;//0;//
    bool BolusOn = 0;
    bool BeamDumpSet = 1;
    bool BeamDumpOn = 1;
    //Name='w';
    File.open(Name,  std::ios::in);
    if(!File.is_open())
    {      }
    else{
        G4String string;
        File >> WedgeOn              >> string>> string;
        File >> GraphiteCollimatorOn >> string>> string;
        File >> FirstScattererOn     >> string>> string;
        File >> WallOn               >> string>> string;
        File >> ContouredScattererOn >> string>> string;
        File >> Collimator1On        >> string>> string;
        File >> RidgeFilterOn        >> string>> string;
        File >> Collimator2On        >> string>> string;
        File >> BolusOn              >> string>> string;
        File >> BeamDumpSet          >> string>> string;
        File >> BeamDumpOn           >> string>> string;

        File.close();
    }





    //   Изменить трапецию на БОКС ислючающий водный фантом,
    //   расстояние до ВФ ~10 см, как и ограничение по пробегу - 10 см
    G4double WaterPhantomStartX=237.0*cm;
    G4double AirGap=30.0 *cm;
    G4double FormSysStartPosX=vacuumZoneXPosition-0.5*vacuumZoneXSize;//
    G4double pX= 0.5*( WaterPhantomStartX-AirGap-FormSysStartPosX);
    G4double pY= worldY-0.5*cm;
    G4double pZ= worldZ-0.5*cm;
    Area1XShift=-(FormSysStartPosX+pX);

    G4Box* Ar1 = new G4Box("Area1Solid",pX,pY,pZ);
    logicAr1 = new G4LogicalVolume(Ar1, airNist, "logicLowQualityArea");
    physiAr1 = new G4PVPlacement(0, G4ThreeVector(-Area1XShift, 0., 0.),
                                 "VacuumZone", logicAr1,
                                 physicalTreatmentRoom, false, 0);

    // logicAr1 -> SetVisAttributes (G4VisAttributes::Invisible);
    red -> SetForceWireframe(true);
   // red-> SetVisibility(false);
    logicAr1 -> SetVisAttributes (red);
    logicTreatmentRoom  -> SetVisAttributes (red);
    //logicAr1=logicTreatmentRoom;
    //physiAr1=physicalTreatmentRoom;
    //Area1XShift=0;

    // Create a region //
    G4Region* RegionOfLowQuality = new G4Region("RegionOfLowQuality");
    // Attach a logical volume to the region
    RegionOfLowQuality->AddRootLogicalVolume(logicAr1);// [...]



    HadrontherapyAcceleratorTube();
    //if (WedgeOn) {HadrontherapyWedge();}
    if (GraphiteCollimatorOn) { HadrontherapyGraphiteCollimator();}
    if (!GraphiteCollimatorOn && BeamDumpSet) { BeamDump(BeamDumpOn);}

    if (FirstScattererOn    ) { HadrontherapyFirstScatterer()    ;}
    if (WallOn              ) { HadrontherapyWall();
        //   RegionOfLowQuality->AddRootLogicalVolume(logicConcreteWall);// [...]
        //   RegionOfLowQuality->AddRootLogicalVolume(logicMetalTube);// [...]
    }

    if (ContouredScattererOn)
    {
        ContouredScatterer = new HadrontherapyContouredScatterer();
        //ContouredScatterer -> BuildScatterer(physicalTreatmentRoom,ContScattBoxPosX);
        ContouredScatterer -> BuildScatterer(physiAr1,ContScattBoxPosX+Area1XShift);
    }

    if (RidgeFilterOn)
    {
        RidgeFilter = new HadrontherapyRidgeFilter();
        // RidgeFilter ->   BuildFilter(physicalTreatmentRoom,RidgeFilterBoxPosX);
        RidgeFilter ->   BuildFilter(physiAr1,RidgeFilterBoxPosX+Area1XShift);
    }


    if (Collimator1On)
    {Collimator(   RFCollimatorBoxX,   RFCollimatorBoxY,   RFCollimatorBoxZ,   RFCollimatorBoxR,   RFCollimatorPosX, Collimator1Material);}
    if (Collimator2On)
    {Collimator(FinalCollimatorBoxX,FinalCollimatorBoxY,FinalCollimatorBoxZ,FinalCollimatorBoxR,FinalCollimatorPosX, Collimator2Material );}


    //construct_PhaseSpace_detector();
//    if (!PhaseSpaceDetector) PhaseSpaceDetector = new PhaseSpaceDetectorConstruction(physiAr1, Area1XShift + 30*cm);



    // Retrieve the region by its name
    G4Region* regionA1    = G4RegionStore::GetInstance()->GetRegion("RegionOfLowQuality");
    //    // Create production cuts
    cutsA1 = new G4ProductionCuts;
    cutsA1->SetProductionCut(20.0*cm, G4ProductionCuts::GetIndex("gamma"));
    cutsA1->SetProductionCut(30.0*cm, G4ProductionCuts::GetIndex("e-"));
    cutsA1->SetProductionCut(30.0*cm, G4ProductionCuts::GetIndex("e+"));
// cutsA1->SetProductionCut(200.0*m, G4ProductionCuts::GetIndex("neutron"));

    // Attach cuts to the region
    regionA1->SetProductionCuts(cutsA1);


    //HadrontherapyContouredScatterer();
    // HadrontherapyPreCollimator();
    // HadrontherapyRidgeFilter();
    // HadrontherapyMainCollimator();
    // HadrontherapyBolus();

//    G4Box* nCtop = new G4Box("nCtop", 45.0*cm,0.5*coreWidth,0.5*coreWidth);
//    G4LogicalVolume* logicNCtop = new G4LogicalVolume(nCtop, airNist, "nCtop");
//    G4VPhysicalVolume* physCoreMainBox =  new G4PVPlacement(0,G4ThreeVector(shellPosX,shellPosY,shellPosZ-switchOn*0.5*(airGap+coreWidth)),
//                                          "nCtop",logicNCtop,
//                                          physiAr1,
//                                          false,0);


}



////////////////////////////////////////////////////////////////////////////
void INRPassiveProtonBeamLine::HadrontherapyAcceleratorTube()
{

    G4Tubs* vacuumZone = new G4Tubs("VacuumZone",
                                    0, vacuumZoneRSize,
                                    0.5*vacuumZoneXSize,
                                    0.0*deg,360.0*deg);
    G4RotationMatrix* rMatrix = new G4RotationMatrix();
    rMatrix -> rotateY(-90*deg);
    G4LogicalVolume* logicVacuumZone = new G4LogicalVolume(vacuumZone, vacuumZoneMaterial, "VacuumZone");//
    G4VPhysicalVolume* physiVacuumZone = new G4PVPlacement(rMatrix, G4ThreeVector(Area1XShift+vacuumZoneXPosition, 0., 0.),
                                                           "VacuumZone", logicVacuumZone,
                                                           physiAr1, true, 0);

    // -------------------//
    // THE KAPTON WINDOWS //
    //--------------------//
    //It prmits the passage of the beam from vacuum to air

    G4Tubs* solidKaptonWindow = new G4Tubs("KaptonWindow",
                                           0,vacuumZoneRSize,
                                           0.5*kaptonWindowXSize,
                                           0.0*deg,360.0*deg);

    G4LogicalVolume* logicKaptonWindow = new G4LogicalVolume(solidKaptonWindow,vacuumZoneWindowMaterial,
                                                             "KaptonWindow");

    physiKaptonWindow = new G4PVPlacement(0, G4ThreeVector(0., 0., 0.5*vacuumZoneXSize-0.5*kaptonWindowXSize),
                                          "KaptonWindow", logicKaptonWindow,
                                          physiVacuumZone, false,	0);

    logicKaptonWindow -> SetVisAttributes(darkOrange3);




    logicVacuumZone -> SetVisAttributes(skyBlue);
    //  logicVacuumZone -> SetVisAttributes(EmptyColor);


}

/////////////////////////////////////////////////////////////////////////////
void INRPassiveProtonBeamLine::HadrontherapyGraphiteCollimator()
{

    // Задаем форму ящика графитового коллиматора
    G4Box* GraphiteCollimatorBox = new G4Box("GraphiteCollimatorBox",
                                             0.5*GraphiteCollimatorXSize,
                                             0.5*GraphiteCollimatorYSize,
                                             0.5*GraphiteCollimatorZSize);


    G4Tubs* GraphiteAirGap = new G4Tubs("GraphiteAirGap",
                                        0, GraphiteCollimatorRSize,
                                        0.5*GraphiteCollimatorXSize+0.01*mm,
                                        0.0*deg,360.0*deg);

    G4RotationMatrix* rMatrix = new G4RotationMatrix();
    rMatrix -> rotateY(-90*deg);
    G4SubtractionSolid* GraphiteCollimatorBody = new G4SubtractionSolid("GraphiteCollimatorBody",
                                                                        GraphiteCollimatorBox, GraphiteAirGap,
                                                                        rMatrix,G4ThreeVector(0,0,0));



    G4LogicalVolume* logicGraphiteCollimatorBox = new G4LogicalVolume(GraphiteCollimatorBody,
                                                                      GraphiteCollimatorMaterial,
                                                                      "GraphiteCollimatorBox");




    physiGraphiteCollimatorBox = new G4PVPlacement(0, G4ThreeVector(Area1XShift+GraphiteCollimatorXPosition, 0.,0.),
                                                   "GraphiteCollimatorBox", logicGraphiteCollimatorBox,
                                                   physiAr1,
                                                   false, 0);




    //G4LogicalVolume* logicGraphiteAirGap = new G4LogicalVolume(GraphiteAirGap, holeFirstCollimatorMaterial, "GraphiteAirGap");

    //physiGraphiteAirGap = new G4PVPlacement(rMatrix, G4ThreeVector(0., 0.,0.),
    //                                        logicGraphiteAirGap,"GraphiteAirGap",
    //                                        logicGraphiteCollimatorBox,
    //                                        false,0);

    //physiGraphiteAirGap = new G4PVPlacement(rMatrix, G4ThreeVector(0., 0.,0.),
    //                                        "GraphiteAirGap",logicGraphiteAirGap,
    //                                        physiGraphiteCollimatorBox,
    //                                        false,0);

    //G4VPhysicalVolume** physiGraphiteAirGap = new G4PVPlacement(rMatrix, 0,
    //                                                   "GraphiteAirGap", logicGraphiteAirGap, physiGraphiteCollimatorBox, false, 0);

    logicGraphiteCollimatorBox -> SetVisAttributes(gray);
    //logicGraphiteAirGap -> SetVisAttributes(EmptyColor);
}
/////////////////////////////////////////////////////////////////////////////
void INRPassiveProtonBeamLine::HadrontherapyFirstScatterer()
{
    firstScatteringFoil = new G4Box("FirstScatteringFoil",
                                    0.5*firstScatteringFoilXSize,
                                    0.5*firstScatteringFoilYSize,
                                    0.5*firstScatteringFoilZSize);

    G4LogicalVolume* logicFirstScatteringFoil = new G4LogicalVolume(firstScatteringFoil,
                                                                    firstScatteringFoilMaterial,
                                                                    "FirstScatteringFoil");

    physiFirstScatteringFoil = new G4PVPlacement(0, G4ThreeVector(Area1XShift+firstScatteringFoilXPosition, 0.,0.),
                                                 "FirstScatteringFoil", logicFirstScatteringFoil,
                                                 physiAr1,
                                                 false, 0);
    logicFirstScatteringFoil -> SetVisAttributes(darkOrange3);

}

/////////////////////////////////////////////////////////////////////////////
void INRPassiveProtonBeamLine::HadrontherapyWall()
{
    ConcreteWall = new G4Box("ConcreteWall",
                             0.5*WallZSize,
                             0.5*WallYSize,
                             0.5*WallXSize);
    ConcreteWallGAP = new G4Tubs("ConcreteWall",
                                 0, WallWindowR,
                                 WallXSize,//WallXSize,//
                                 0.0*deg,360.0*deg);




    G4RotationMatrix* rMatrix = new G4RotationMatrix();
    //
    //rMatrix -> rotateY(-90*deg);//-90*deg
    rMatrix -> rotateY(WallAngle);
    ConcreteWallwGAP =  new G4SubtractionSolid("ConcreteWall",
                                               ConcreteWall, ConcreteWallGAP,
                                               rMatrix, G4ThreeVector(0,0,0));





    //G4LogicalVolume* logicConcreteWall
    logicConcreteWall = new G4LogicalVolume(ConcreteWallwGAP,WallMaterial,"ConcreteWall");

    G4RotationMatrix* rMatrix2 = new G4RotationMatrix();
    rMatrix2 -> rotateY(-90*deg-WallAngle);
    physiConcreteWall  = new G4PVPlacement(rMatrix2, G4ThreeVector(WallXPosition+Area1XShift, 0.,0.),
                                           "ConcreteWall", logicConcreteWall, physiAr1,
                                           false, 0);






    G4RotationMatrix* rMatrix3 = new G4RotationMatrix();
    //rMatrix -> rotateX(WallAngle);
    rMatrix3 -> rotateY(-90*deg);//-90*deg
    MetalTube = new G4Tubs("MetalTube",
                           AirWindowR, WallWindowR,
                           0.5*(8.3*cm+WallXSize/cos(WallAngle)),
                           0.0*deg,360.0*deg);
    logicMetalTube = new G4LogicalVolume(MetalTube,MetalTubeMaterial, "MetalTube");

    physiMetalTube = new G4PVPlacement(rMatrix3, G4ThreeVector(WallXPosition+Area1XShift, 0.,0.),
                                       "MetalTube", logicMetalTube, physiAr1,
                                       false, 0);


    PolyTube = new G4Tubs("PolyTube",
                          AirWindowR+MetalTubeThickness, WallWindowR-MetalTubeThickness,
                          0.5*WallXSize+2*cm-MetalTubeThickness,
                          0.0*deg,360.0*deg);
    logicPolyTube = new G4LogicalVolume(PolyTube,PolyTubeMaterial, "PolyTube");

    physiPolyTube = new G4PVPlacement(0, G4ThreeVector(0., 0.,0.),
                                      "PolyTube", logicPolyTube, physiMetalTube,
                                      false, 0);



    logicConcreteWall -> SetVisAttributes(sandColor);//EmptyColor);//
    logicMetalTube -> SetVisAttributes(metalColor);
    //logicMetalTube -> SetVisAttributes(darkGreen);
}
/////////////////////////////////////////////////////////////////////////////

void INRPassiveProtonBeamLine::Collimator(G4double CollimatorBoxX,G4double CollimatorBoxY,G4double CollimatorBoxZ,
                                          G4double CollimatorBoxR,
                                          G4double CollimatorXPosition,
                                          G4Material* CollimatorMaterial)
{
    G4Box* CollimatorBox = new G4Box("CollimatorBox", 0.5*CollimatorBoxX,0.5*CollimatorBoxY,0.5*CollimatorBoxZ);
    G4Tubs* CollimatorTunnel = new G4Tubs("CollimatorTunnel", 0, CollimatorBoxR, 0.5*CollimatorBoxX+0.01*mm, 0.0*deg, 360.0*deg);
    G4RotationMatrix* rMatrix = new G4RotationMatrix();
    rMatrix -> rotateY(-90.0*deg);
    G4SubtractionSolid* CollimatorBody = new G4SubtractionSolid("CollimatorBody",
                                                                CollimatorBox, CollimatorTunnel,
                                                                rMatrix,G4ThreeVector(0,0,0));
    G4LogicalVolume* logicCollimator = new G4LogicalVolume(CollimatorBody, CollimatorMaterial, "logCollimator");


    //     CollimatorMaterial,
    //     "CollimatorLogic");
    physiCollimator =  new G4PVPlacement(0,G4ThreeVector(CollimatorXPosition+Area1XShift,0,0),
                                         "physiCollimator",logicCollimator,
                                         physiAr1,
                                         false,0);
    logicCollimator -> SetVisAttributes(gray);

}

/////////////////////////////////////////////////////////////////////////////
void INRPassiveProtonBeamLine::HadrontherapyBolus()
{

}

/////////////////////////////////////////////////////////////////////////////
void INRPassiveProtonBeamLine::CellTestPhantom()
{
    // Пластиковый фантом
    G4Box* PlasticWaterBox = new G4Box("PlasticWaterBox", 0.5*PlasticWaterBoxX,0.5*PlasticWaterBoxY,0.5*PlasticWaterBoxZ);
    G4LogicalVolume* logicPlasticWater = new G4LogicalVolume(PlasticWaterBox, PlasticWaterMaterial, "logPlasticWater");

    physPlasticWater =  new G4PVPlacement(0,G4ThreeVector(PlasticWaterBoxPosX,-0.5*PlasticWaterBoxY,0),
                                          "physiPlasticWater",logicPlasticWater,
                                          physiAr1,
                                          false,0);
    G4Box* PlasticWaterBox2 = new G4Box("PlasticWaterBox2", 0.5*PlasticWaterBoxX,0.5*PlasticWaterBoxY,0.5*PlasticWaterBoxZ);
    G4LogicalVolume* logicPlasticWater2 = new G4LogicalVolume(PlasticWaterBox2, PlasticWaterMaterial, "logPlasticWater2");

    physPlasticWater2 =  new G4PVPlacement(0,G4ThreeVector(PlasticWaterBoxPosX,+0.5*PlasticWaterBoxY,0),
                                           "physiPlasticWater2",logicPlasticWater2,
                                           physiAr1,
                                           false,0);

    //Дополнительные переменные
    G4double EpendorfHatHeight=0.5*mm;
    G4double EpendorfHatPosY=0.5*(PlasticWaterBoxY-EpendorfHatHeight);

    G4double EpendorfConePosY=0.5*(PlasticWaterBoxY-EpendorfHeight)-EpendorfHatHeight;

    G4double EpendorfBottomHeight=1.0*mm;

    G4double EpendorfBottomRminTop = (EpendorfBottomHeight/EpendorfHeight)*(EpendorfRminTop-EpendorfRminBottom)+EpendorfRminBottom;
    G4RotationMatrix* rMatrix = new G4RotationMatrix();
    rMatrix -> rotateX(90*deg);

    // "шляпка эпендорфа"
    G4Tubs* EpendorfHatSolid=new G4Tubs("EpndrfHatSolid",
                                        0.0,EpendorfRmaxTop+0.5*mm,
                                        0.5*EpendorfHatHeight,
                                        0.0*deg,360.0*deg);
    G4LogicalVolume* logicEpendorfHat = new G4LogicalVolume(EpendorfHatSolid, EpendorfMaterial, "logEpendorfHat");
    physEpendorfHat =  new G4PVPlacement(rMatrix,G4ThreeVector(0,EpendorfHatPosY,0),
                                         "physEpendorfHat",logicEpendorfHat,
                                         physPlasticWater,
                                         false,0);

    G4Tubs* EpendorfInsideHatSolid=new G4Tubs("EpndrfInsideHatSolid",
                                              0.0,EpendorfRminTop,
                                              0.5*EpendorfHatHeight,
                                              0.0*deg,360.0*deg);
    G4LogicalVolume* logicEpendorfInsideHat = new G4LogicalVolume(EpendorfInsideHatSolid, EpendorfInsideMaterial, "logEpendorfInsideHat");
    physEpendorfInsideHat =  new G4PVPlacement(0,G4ThreeVector(0,0,0),
                                               "physEpendorfInsideHat",logicEpendorfInsideHat,
                                               physEpendorfHat,
                                               false,0);


    // "тело эпендорфа"
    G4Cons* EpendorfConeSolid = new G4Cons("EpndrfCone",
                                           0.0,EpendorfRmaxBottom,
                                           0.0,EpendorfRmaxTop,
                                           0.5*EpendorfHeight,
                                           0.0*deg,360.0*deg);
    G4LogicalVolume* logicEpendorfCone = new G4LogicalVolume(EpendorfConeSolid, EpendorfMaterial, "logEpendorfCone");
    physEpendorfCone =  new G4PVPlacement(rMatrix,G4ThreeVector(0.0,EpendorfConePosY,0.0),
                                          "physEpendorfCone",logicEpendorfCone,
                                          physPlasticWater,
                                          false,0);

    G4Cons* EpendorfInsideConeSolid = new G4Cons("EpndrfInsideCone",
                                                 0.0,EpendorfBottomRminTop,
                                                 0.0,EpendorfRminTop,
                                                 0.5*(EpendorfHeight-EpendorfBottomHeight),
                                                 0.0*deg,360.0*deg);
    G4LogicalVolume* logicInsideEpendorfCone = new G4LogicalVolume(EpendorfInsideConeSolid, EpendorfInsideMaterial, "logEpendorfCone");
    physEpendorfInsideCone =  new G4PVPlacement(0,G4ThreeVector(0,0,0.5*EpendorfBottomHeight),
                                                "physEpendorfCone",logicInsideEpendorfCone,
                                                physEpendorfCone,
                                                false,0);

    // перемычка
    // /vis/viewer/set/targetPoint 93 2.5 0 cm
    //   G4double Perem_dx1=1.0*mm;
    //   G4double Perem_dx2=Perem_dx1+2.0*(EpendorfRmaxTop-EpendorfRmaxBottom);
    //   G4double Perem_dy=1.5*mm;
    //   G4double Perem_dz=EpendorfHeight;
    G4Trd* BTwinEsolid = new G4Trd("BTwinEsolid",
                                   0.5*Perem_dx2,0.5*Perem_dx1,       //hX
                                   0.5*Perem_dy ,0.5*Perem_dy, //hY
                                   0.5*Perem_dz);         //hZ
    G4LogicalVolume* logicBTwinE = new G4LogicalVolume(BTwinEsolid, PlasticWaterMaterial, "logBTwin");

    physBTwinE =  new G4PVPlacement(rMatrix ,G4ThreeVector(EpendorfRmaxTop+0.5*Perem_dx1,EpendorfConePosY,0.0),
                                    "physBTwinE",logicBTwinE,
                                    physPlasticWater,
                                    true,0);



    gray -> SetForceWireframe(true);
    logicPlasticWater -> SetVisAttributes(gray);
    //yellow -> SetForceWireframe(true);
    logicEpendorfHat -> SetVisAttributes(yellow);
    logicEpendorfInsideHat-> SetVisAttributes(blue);
    //logicEpendorfInsideHat -> SetVisAttributes(blue);

    logicEpendorfCone ->SetVisAttributes(yellow);//SetVisAttributes(gray);//
    logicInsideEpendorfCone -> SetVisAttributes(blue);

    logicBTwinE  -> SetVisAttributes(yellow);

}
/////////////////////////////////////////////////////////////////////////////

void INRPassiveProtonBeamLine::BeamDump(bool dumpOn = false)
{
    //G4Material* air =  G4NistManager::Instance()->FindOrBuildMaterial("G4_AIR", false);

    G4double coreLength = 30.0 *cm;
    G4double coreWidth  = 25.0 *cm;

    G4double airGap  = 5.0 *cm;
    G4double shellThickness = 5.0*cm;
    G4double shellWindowRadius = 3.5*cm;

    G4double shellWidth = coreWidth*2 + airGap*3+2*shellThickness;
    G4double shellLength = coreLength + 2*airGap;
    G4double shellHeight = coreWidth + 2*airGap+2*shellThickness;

    G4double shellPosZ  =  0.5*(coreWidth+airGap);//coreWidth*1.5 + airGap*2+shellThickness;//0.5*(shellWidth+airGap);//
    //G4double a = coreWidth*1.5 + airGap*2+shellThickness;
    G4double shellPosY = 0*cm;
    G4double shellPosX =  Area1XShift+GraphiteCollimatorXPosition;

    //dumpShellMaterial = lead;
    //dumpCoreMaterial = graphite;

    // Внешняя оболочка - свинец
    //Тело с воздушной полостью
    //Тело
    G4Box* shellBox = new G4Box("Box", 0.5*shellLength,0.5*shellHeight,0.5*shellWidth);
    G4Box* shellAirBox = new G4Box("AirBox", 0.5*shellLength+1*mm,0.5*shellHeight-shellThickness,0.5*shellWidth-shellThickness);

    G4SubtractionSolid* shellMainBox = new G4SubtractionSolid("MainBox",
                                                             shellBox, shellAirBox,
                                                             0,G4ThreeVector(0,0,0));


    G4LogicalVolume* logicShellMainBox = new G4LogicalVolume(shellMainBox, dumpShellMaterial, "logShellMainBox");

    physShellMainBox =  new G4PVPlacement(0,G4ThreeVector(shellPosX,shellPosY,shellPosZ),
                                          "physiShellMainBox",logicShellMainBox,
                                          physiAr1,
                                          false,0);

    logicShellMainBox -> SetVisAttributes(gray);


    //Передняя и задняя стенки с отверстием
    G4Box* shWall = new G4Box("shWall", 0.5*shellThickness,0.5*shellHeight,0.5*shellWidth);
    G4Tubs* shellWindow = new G4Tubs("shWin",
                                    0, shellWindowRadius,
                                    0.5*shellThickness+0.01*mm,
                                    0.0*deg,360.0*deg);


    G4RotationMatrix* rMatrix = new G4RotationMatrix();
    rMatrix -> rotateY(-90*deg);


    G4SubtractionSolid* shellWall = new G4SubtractionSolid("shellWall",
                                                             shWall, shellWindow,
                                                             rMatrix,G4ThreeVector(0,0,-shellPosZ));


    G4LogicalVolume* logicShellWall = new G4LogicalVolume(shellWall, dumpShellMaterial, "logShellWall");

    physShellWallFront =  new G4PVPlacement(0,G4ThreeVector(shellPosX-0.5*shellLength-0.5*shellThickness,shellPosY,shellPosZ),
                                          "physiShellMainBox",logicShellWall,
                                          physiAr1,
                                          false,0);
    physShellWallBack =  new G4PVPlacement(0,G4ThreeVector(shellPosX+0.5*shellLength+0.5*shellThickness,shellPosY,shellPosZ),
                                          "physiShellMainBox",logicShellWall,
                                          physiAr1,
                                          false,0);

   // logicShellWall -> SetVisAttributes(gray);
    logicShellWall -> SetVisAttributes(EmptyColor);


    //Графитовая мишень
    int switchOn = -1;
    if(dumpOn) switchOn = 1;

    G4Box* coreBox = new G4Box("Box", 0.5*coreLength,0.5*coreWidth,0.5*coreWidth);
    G4LogicalVolume* logicCoreBox = new G4LogicalVolume(coreBox, dumpCoreMaterial, "logCoreMainBox");
    physCoreMainBox =  new G4PVPlacement(0,G4ThreeVector(shellPosX,shellPosY,shellPosZ-switchOn*0.5*(airGap+coreWidth)),
                                          "physiCoreMainBox",logicCoreBox,
                                          physiAr1,
                                          false,0);

    logicCoreBox -> SetVisAttributes(gray);

//    G4Box* nSDtop = new G4Box("Box", 0.5*coreLength,0.5*coreWidth,0.5*coreWidth);
//    logicSDtop = new G4LogicalVolume(nSDtop, air, "logicSDtop");
//    physSDtop =  new G4PVPlacement(0,G4ThreeVector(shellPosX,shellPosY,shellPosZ-switchOn*0.5*(airGap+coreWidth)),
//                                          "physSDtop",logicSDtop,
//                                          physiAr1,
//                                          false,0);



}
/////////////////////////////////////////////////////////////////////////////
void INRPassiveProtonBeamLine::construct_PhaseSpace_detector(){

    G4double voxel_sX = 1 *mm;
    G4double voxel_sY = 1 *mm;
    G4double voxel_sZ = 0.5 *mm;

    G4int Nx = 100;
    G4int Ny = 100;
    G4int Nz = 1;

    // смещения:
    // + 10*см
    // + 30*см
    // + 130*см
    G4double phSp_posX = Area1XShift + 130. *cm;
    G4double phSp_posY = 0;
    G4double phSp_posZ = 0;

    G4double voxel_HX = 0.5*voxel_sX;
    G4double voxel_HY = 0.5*voxel_sY;
    G4double voxel_HZ = 0.5*voxel_sZ;
    G4double gap = 0 *mm;


    G4double phSp_HX = 0.5*((voxel_sX+gap)*Nx - gap);
    G4double phSp_HY = 0.5*((voxel_sY+gap)*Ny - gap);
    G4double phSp_HZ = voxel_HZ*Nz;


    G4Material* phSpDet_Material = G4NistManager::Instance()->FindOrBuildMaterial("G4_AIR", false);
    //physiAr1

    solid_phSpDet = new G4Box("solid_phSpDet",
                              phSp_HX,
                              phSp_HY,
                              phSp_HZ);

    // Definition of the logical volume of the solid_phSpDet
    logical_phSpDet = new G4LogicalVolume(solid_phSpDet,
                                          phSpDet_Material,
                                          "phSpDet_Log", 0, 0, 0);
    G4RotationMatrix* rMatrix = new G4RotationMatrix();
    rMatrix -> rotateY(-90.0*deg);
    rMatrix -> rotateZ(180.0*deg);

    // Definition of the physics volume of the Phantom
    phys_phSpDet = new G4PVPlacement(rMatrix,
                                     G4ThreeVector(phSp_posX,phSp_posY,phSp_posZ),
                                     "phSpDet_Phys",
                                     logical_phSpDet,
                                     physiAr1,
                                     false,
                                     0);

    // Visualisation attributes of the phantom
    G4VisAttributes *whiteBlock = new G4VisAttributes(G4Colour(1, 1 ,1));
    whiteBlock -> SetVisibility(true);
    //whiteBlock -> SetForceSolid(true);
    whiteBlock -> SetForceWireframe(true);
    logical_phSpDet -> SetVisAttributes(whiteBlock);

    solid_voxel = new G4Box("solidSD", voxel_HX, voxel_HY, voxel_HZ);
    logic_voxel = new G4LogicalVolume(solid_voxel, phSpDet_Material, "PhaseSpaceLV");




    G4bool check4overlaps = false;
    int sdNum = 0;

//    phys_voxel =new G4PVPlacement(0,
//                                  G4ThreeVector(-phSp_HX + voxel_HX + 1*(2*voxel_HX + gap),
//                                                -phSp_HY + voxel_HY + 1*(2*voxel_HY + gap),
//                                                0),
//                                  logic_voxel, "physSD",
//                                  logical_phSpDet,
//                                  true, 0 + 0 * Nx * Ny,
//                                  check4overlaps);

    for (int i=0; i<Nx; ++i){
        for (int j=0; j < Ny; ++j){
            ++sdNum;
            phys_voxel =new G4PVPlacement(0,
                                          G4ThreeVector(-phSp_HX + voxel_HX + i*(2*voxel_HX + gap),
                                                        -phSp_HY + voxel_HY + j*(2*voxel_HY + gap),
                                                        0),
                                          logic_voxel, "PhaseSpacePhysVol",
                                          logical_phSpDet,
                                          true, j + i * Ny,
                                          check4overlaps);
        }
    }
    G4VisAttributes *solid_gray = new G4VisAttributes(G4Colour(65, 65 ,65 ));
    solid_gray -> SetVisibility(true);
    solid_gray -> SetForceSolid(true);

    logic_voxel -> SetVisAttributes(solid_gray);
    phaseSpace_ScoringVolume = logic_voxel;

//    G4cout << "INR_beamline" << G4endl;
//    PhaseSpace_data_collection = PhaseSpaceDataset::getInstance(Nx, -phSp_HX, phSp_HX,
//                                                                Ny, -phSp_HY, phSp_HZ,
//                                                                100, 150 *MeV, 160* MeV);


}

/////////////////////////////////////////////////////////////////////////////

/////////////////////////////////////////////////////////////////////////////

/////////////////////////////////////////////////////////////////////////////

/////////////////////////////////////////////////////////////////////////////

/////////////////////////////////////////////////////////////////////////////

/////////////////////////////////////////////////////////////////////////////








/////////////////////////// MESSENGER ///////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////


/////////////////////////////////////////////////////////////////////////////
void INRPassiveProtonBeamLine::SetFirstScatteringFoilXSize(G4double value)
{
    firstScatteringFoil -> SetXHalfLength(value);
    G4RunManager::GetRunManager() -> GeometryHasBeenModified();
    G4cout <<"The X size of the first scattering foil is (mm):"<<
             ((firstScatteringFoil -> GetXHalfLength())*2.)/mm
          << G4endl;
}


