
#include "DetectorConstruction.hh"

#include "G4Material.hh"
#include "G4Element.hh"
#include "G4LogicalBorderSurface.hh"
#include "G4LogicalSkinSurface.hh"
#include "G4OpticalSurface.hh"
#include "G4Polycone.hh"
#include "G4Box.hh"
#include "G4Trd.hh"
#include "G4Cons.hh"
#include "G4Tubs.hh"
#include "G4LogicalVolume.hh"
#include "G4ThreeVector.hh"
#include "G4PVPlacement.hh"
#include "G4SystemOfUnits.hh"
#include "G4NistManager.hh"
#include "G4VisAttributes.hh"
#include "G4SubtractionSolid.hh"

G4int SiPM_n = 1;
G4int SiPM_m = 1;
G4double SiPM_pitch = 3; //SiPM pixel pitch/2
G4double SiPM_sep = 0.05;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

DetectorConstruction::DetectorConstruction()
 : G4VUserDetectorConstruction()
{
  fExpHall_x = fExpHall_y = fExpHall_z = 10*m;//set expermiental hall size
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

DetectorConstruction::~DetectorConstruction(){;}//destructor is empty

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4VPhysicalVolume* DetectorConstruction::Construct()
{

// ------------- Materials -------------

  G4double a, z, density; //these are the labels of the stuff
  G4int nelements;
  G4NistManager* man = G4NistManager::Instance();//make a new manager, works well for modular

// Elements
//
  G4Element* N = new G4Element("Nitrogen", "N", z=7 , a=14.01*g/mole);
  G4Element* O = new G4Element("Oxygen"  , "O", z=8 , a=16.00*g/mole);
  G4Element* C = new G4Element("Carbon"  , "C",  z=6 , a=12.011*g/mole);
  G4Element* H = new G4Element("Hydrogen" , "H", z=1, a=1.008*g/mole);
  G4Element* Cs = new G4Element("Caesium",  "Cs", z=55, a=132.905*g/mole);
  G4Element* Li = new G4Element("Lithium",  "Li", z=3,  a=6.94*g/mole);
  G4Element* Y  = new G4Element("Yttrium",  "Y",  z=39, a=88.905*g/mole);
  G4Element* Cl = new G4Element("Chlorine", "Cl", z=17, a=35.45*g/mole);
  G4Element* Ce = new G4Element("Cerium",   "Ce", z=58, a=140.116*g/mole);
  G4Element* Li_6 = new G4Element("Lithium", "Li", z=3, a=6.*g/mole);
  G4Element* Gd = new G4Element("Gadolinium", "Gd", z=64, a=157.25*g/mole);
  G4Element* Al_elt = new G4Element("Aluminium", "Al", z=13, a=26.98*g/mole);
  G4Element* Ga = new G4Element("Gallium", "Ga", z=31 , a=69.723*g/mole);
  //Air
  //
  G4Material* air = new G4Material("Air", density=1.29*mg/cm3, nelements=2);
  air->AddElement(N, 70.*perCent);
  air->AddElement(O, 30.*perCent);

// Water
//
  //G4Element* H = new G4Element("Hydrogen", "H", z=1 , a=1.01*g/mole);

  G4Material* water = new G4Material("Water", density= 1.0*g/cm3, nelements=2);
  water->AddElement(H, 2);
  water->AddElement(O, 1);
// perspex
//
  //G4Element* H = new G4Element("Hydrogen", "H", z=1 , a=1.01*g/mole);

  G4Material* perspex = new G4Material("Perspex", density= 1.18*g/cm3, nelements=3);
  perspex->AddElement(H, 8);
  perspex->AddElement(O, 2);
  perspex->AddElement(C, 5);

    G4Material* SiPM_mat = new G4Material("SiPM_mat", density= 1.18*g/cm3, nelements=3);
  SiPM_mat->AddElement(H, 8);
  SiPM_mat->AddElement(O, 2);
  SiPM_mat->AddElement(C, 5);
//wrapping, density about perspex for now
  G4Material* wrapping = new G4Material("wrapping", density= 1.18*g/cm3, nelements=3);
  wrapping->AddElement(H, 8);
  wrapping->AddElement(O, 2);
  wrapping->AddElement(C, 5);

//epoxy
  G4Material* Epoxy = new G4Material("Epoxy", density= 1.69*g/cm3, nelements=4);
  Epoxy->AddElement(H, 25);
  Epoxy->AddElement(O, 5);
  Epoxy->AddElement(C, 21);
  Epoxy->AddElement(Cl, 1);

// CLYC
//

  G4Material* CLYC = new G4Material("CLYC", density= 3.31*g/cm3, nelements=5);
  CLYC->AddElement(Cs, 2);
  CLYC->AddElement(Li_6, 1);
  CLYC->AddElement(Y, 1);
  CLYC->AddElement(Cl,6);
  CLYC->AddElement(Ce,1);
  //Currently wrong proportion of Ce, but does not seem to cause much interferance
  
  //GAGG
  G4Material* GAGG = new G4Material("GAGG", density = 6.63*g/cm3,nelements=5);
  GAGG->AddElement(Gd,3);
  GAGG->AddElement(Al_elt,2);
  GAGG->AddElement(Ga,3);
  GAGG->AddElement(O,12);
  GAGG->AddElement(Ce,1);
  //again wrong Ce, quicker to fix than comment but CBA
  
  G4Material* quartz = man->FindOrBuildMaterial("G4_SILICON_DIOXIDE");
  G4Material* Al = man->FindOrBuildMaterial("G4_Al");
  G4Material* metal = man->FindOrBuildMaterial("G4_STAINLESS-STEEL"); 
//
// ------------ Generate & Add Material Properties Table ------------

// Photon energies (to define refractive index at a given energy)
 
G4double photonEnergy[] =
            { 2.034*eV, 2.068*eV, 2.103*eV, 2.139*eV,
              2.177*eV, 2.216*eV, 2.256*eV, 2.298*eV,
              2.341*eV, 2.386*eV, 2.433*eV, 2.481*eV,
              2.532*eV, 2.585*eV, 2.640*eV, 2.697*eV,
              2.757*eV, 2.820*eV, 2.885*eV, 2.954*eV,
              3.026*eV, 3.102*eV, 3.181*eV, 3.265*eV,
              3.353*eV, 3.446*eV, 3.545*eV, 3.649*eV,
              3.760*eV, 3.877*eV, 4.002*eV, 4.136*eV };
 
  const G4int nEntries = sizeof(photonEnergy)/sizeof(G4double);
 
//
// Air
 
  G4double refractiveIndex2[] =
            { 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00,
              1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00,
              1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00,
              1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00,
              1.00, 1.00, 1.00, 1.00 };

  G4MaterialPropertiesTable* Air = new G4MaterialPropertiesTable();
  Air->AddProperty("RINDEX", photonEnergy, refractiveIndex2, nEntries)
    ->SetSpline(true);
  //only adding R-index for air, its all it deservs
  G4cout << "Air G4MaterialPropertiesTable" << G4endl;
  Air->DumpTable();
  //can dump the tables, this writes it to the running file
  air->SetMaterialPropertiesTable(Air);

//
// Perspex
//
  G4double refractiveIndex3[] =
            { 1.486, 1.486, 1.486, 1.486, 1.486, 1.486, 1.486,
              1.486, 1.486, 1.486, 1.486, 1.486, 1.486, 1.486,
              1.486, 1.486, 1.486, 1.486, 1.486, 1.486, 1.486,
              1.486, 1.486, 1.486, 1.486, 1.486, 1.486, 1.486,
              1.486, 1.486, 1.486, 1.486 };

  G4MaterialPropertiesTable* Perspex = new G4MaterialPropertiesTable();
  Perspex->AddProperty("RINDEX", photonEnergy, refractiveIndex3, nEntries)->SetSpline(true);

  G4cout << "PERRSPEX G4MaterialPropertiesTable" << G4endl;
  Perspex->DumpTable();

  perspex->SetMaterialPropertiesTable(Perspex);
  SiPM_mat->SetMaterialPropertiesTable(Perspex);
  //about right... I know, lazy right

// GAGG
//
    G4MaterialPropertiesTable* myMPTGAGG = new G4MaterialPropertiesTable();
  G4double refractiveIndexGAGG[] =
            { 1.9, 1.9,  1.9, 1.9,  1.9,
              1.9,  1.9, 1.9,  1.9, 1.9,
              1.9, 1.9, 1.9,   1.9, 1.9,
              1.9, 1.9, 1.9, 1.9, 1.9,
              1.9, 1.9,  1.9, 1.9,  1.9,
              1.9, 1.9,  1.9, 1.9,  1.9,
              1.9,   1.9};

  assert(sizeof(refractiveIndexGAGG) == sizeof(photonEnergy));


  G4double scintilFastGAGG[] =
            { 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00,
              1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00,
              1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00,
              1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00,
              1.00, 1.00, 1.00, 1.00 };

  assert(sizeof(scintilFastGAGG) == sizeof(photonEnergy));

 
  G4double scintilSlowGAGG[] =
            { 1., 1., 1., 1., 1., 1., 1.,
              1., 1., 1., 1., 1., 1., 1.,
              1., 1., 1., 1., 1., 1., 1.,
              1., 1., 1., 1., 1., 1., 1.,
              1., 1., 1., 1. };
  assert(sizeof(scintilSlowGAGG) == sizeof(photonEnergy));

  myMPTGAGG->AddProperty("RINDEX",       photonEnergy, refractiveIndexGAGG,nEntries)
        ->SetSpline(true);
  myMPTGAGG->AddProperty("FASTCOMPONENT",photonEnergy, scintilFastGAGG,     nEntries)
        ->SetSpline(true);
  myMPTGAGG->AddProperty("SLOWCOMPONENT",photonEnergy, scintilSlowGAGG,     nEntries)
        ->SetSpline(true);

  myMPTGAGG->AddConstProperty("SCINTILLATIONYIELD",40000./MeV);
  myMPTGAGG->AddConstProperty("RESOLUTIONSCALE",1.0);
  myMPTGAGG->AddConstProperty("FASTTIMECONSTANT", 50.*ns);
  myMPTGAGG->AddConstProperty("SLOWTIMECONSTANT",150.*ns);
  myMPTGAGG->AddConstProperty("YIELDRATIO",0.91);
  GAGG->SetMaterialPropertiesTable(myMPTGAGG);
  //I assume here 40k photons/MeV and that the fast slow ration is 91%, don't know
  //if this is true, but for now, close enough
  //quartz
  G4double refractiveIndex5[] =
            { 1.46, 1.46,  1.46, 1.46,  1.46,
              1.46,  1.46, 1.46,  1.46, 1.46,
              1.46, 1.46, 1.46,   1.46, 1.46,
              1.46, 1.46, 1.46, 1.46, 1.46,
              1.46, 1.46,  1.46, 1.46,  1.46,
              1.46, 1.46,  1.46, 1.46,  1.46,
              1.46,   1.46};

  assert(sizeof(refractiveIndex4) == sizeof(photonEnergy));
  G4MaterialPropertiesTable* myMPT5 = new G4MaterialPropertiesTable();
  myMPT5->AddProperty("RINDEX",photonEnergy, refractiveIndex5,nEntries)->SetSpline(true);
  quartz->SetMaterialPropertiesTable(myMPT5);
//
// ------------- Volumes --------------

// The experimental Hall
//
  G4Box* expHall_box = new G4Box("World",fExpHall_x,fExpHall_y,fExpHall_z);

  G4LogicalVolume* expHall_log
    = new G4LogicalVolume(expHall_box,air,"World",0,0,0);

  //G4VPhysicalVolume* expHall_phys
  expHall_phys
  = new G4PVPlacement(0,G4ThreeVector(),expHall_log,"World",0,false,0);
    expHall_log->SetVisAttributes(G4VisAttributes::Invisible);
    
// The Scintillator
//
    G4double r1[6] = {0*mm,5*mm,10*mm,15*mm,20*mm,24.01*mm};
    G4double r2[6] = {0,0,0,0,0,0};
    G4double r3[6] = {28.25*mm,26.65*mm,24.779*mm ,22.525*mm,19.488*mm ,13.34*mm};
    G4double r4[6] = {28.35*mm,26.75*mm,24.879*mm,22.625*mm,19.588*mm,13.44*mm};
    G4double r5[6] = {28.55*mm,26.95*mm,25.079*mm,22.825*mm,19.788*mm,13.64*mm};
    G4double zpos[2] = {0,10*mm};
    G4double innerd[2] = {0*mm,0*mm};
    G4double outerd[2] = {1.5*mm, 0.15*mm};
    G4double wrapin[2] = {1.51*mm, 0.16*mm};
    G4double wrapout[2] = {1.56*mm, 0.21*mm};

    // Cylinder
    // G4Cons* Guide_box = new G4Cons("Guide",0*mm,25.5*mm,0*mm,25.4/2*mm,30*mm,0*deg,360*deg);
    // Cuboid - 12x12x30mm for Beth's experiment
    G4Box* Guide_box=new G4Box("Guide",6,6,15);
    // Cone
    // G4Polycone* Guide_box = new G4Polycone("Guide",0,360*deg,2,zpos,innerd,outerd);
    // G4Polycone* Eminem = new G4Polycone("Tupac",0,360*deg,2,zpos,wrapin,wrapout);

  // Setting Materials for the scintillator and the scintillator cover (Guide_log = scintillator, Snooop = Al coating) 
  G4LogicalVolume* Guide_log
    = new G4LogicalVolume(Guide_box,GAGG,"Guide",0,0,0);
 
 
  //G4LogicalVolume* Snoop 
   // = new G4LogicalVolume(Eminem,Al,"Tupac",0,0,0);

   // Setting scintillator colour to orange (just for funsies really)
   G4VisAttributes* boxVisAtt = new G4VisAttributes(G4Colour(1.0,.43,0.0));
   Guide_log->SetVisAttributes(boxVisAtt);

   fScint = new G4PVPlacement(0,G4ThreeVector(0,0,0.3*mm),Guide_log,"Guide",
                        expHall_log,false,0);

//G4VPhysicalVolume * Nelly = 
     // new G4PVPlacement(0,G4ThreeVector(0,0,0.3*mm),Snoop,"Tupac", expHall_log,false,0);

 
// Makes a smaller volume inside the larger volume so you can wrap the scintillator, or whatever.
 
G4Box* CookieCutter = new G4Box("Cutter", 6.016,6.016,15.016);
//G4Tubs* EpoxyFill = new G4Tubs("Dummy", 0. *mm, 0.912*mm, 1.75*mm, 0*deg, 360*deg);
G4SubtractionSolid* Turtle = new G4SubtractionSolid("Turtle", CookieCutter, Guide_box, 0, G4ThreeVector(0,0,0.017));
G4LogicalVolume * Turtle_log = new G4LogicalVolume(Donatello, Al, "Turtle", 0,0,0);
G4VPhysicalVolume * Turtle_phys = new G4PVPlacement(0,G4ThreeVector(0,0,15.3*mm), Turtle_log, "Turtle", expHall_log, false, 0);



//-----------making the SiPM array----------
//
//SiPM object, only the optical part for now, rest plays havok with array
 G4Box* SiPM_box = new G4Box("SiPM", SiPM_pitch,SiPM_pitch, 0.3*mm);
 SiPM_log = new G4LogicalVolume(SiPM_box, SiPM_mat, "SiPM",0,0,0);
 //now to make an array, the nxm matrix values should be at the top of this as SiPM_n and m
 if(SiPM_n==1&&SiPM_m==1){
fSiPM = new G4PVPlacement(0,G4ThreeVector(0,0,0), SiPM_log,"SiPM", expHall_log,false,0);
 }
 else{
 G4double iter_x;
 G4double iter_y;
 int cpy_n; 
 for(int i = 0; i<SiPM_n;i++){
   for(int j = 0; j<SiPM_n;j++){
     iter_x = 2*i*(SiPM_pitch+SiPM_sep)*mm - (SiPM_n-1)*(SiPM_pitch+SiPM_sep);
     iter_y = 2*j*(SiPM_pitch+SiPM_sep)*mm - (SiPM_m-1)*(SiPM_pitch+SiPM_sep);
     cpy_n = SiPM_n*i+j;
     fSiPM = new G4PVPlacement(0,G4ThreeVector(iter_x,iter_y,-0.3*mm),SiPM_log,"SiPM",expHall_log,false,cpy_n);
   
   }
 }
 }
 
// ------------- Surfaces --------------
//
//Wrapping
//
  G4OpticalSurface* opWrapSurface = new G4OpticalSurface("WrapSurface");
  opWrapSurface->SetType(dielectric_metal);
  opWrapSurface->SetFinish(polishedfrontpainted);
  opWrapSurface->SetModel(glisur);

  G4double pp[] = {2.0*eV, 3.5*eV};
  const G4int num2 = sizeof(pp)/sizeof(G4double);
  G4double reflectivity3[] = {1.,1.};
  assert(sizeof(reflectivity2) == sizeof(pp));
  G4double efficiency2[] = {0.0,0.0};
  assert(sizeof(effieciency2) == sizeof(pp));
  G4MaterialPropertiesTable* ScintWrap_MPT = new G4MaterialPropertiesTable();
  ScintWrap_MPT->AddProperty("Reflectivity",pp,reflectivity3,num2);
  ScintWrap_MPT->AddProperty("Efficiency",pp,efficiency2,num2);
  opWrapSurface->SetMaterialPropertiesTable(ScintWrap_MPT);
  new G4LogicalSkinSurface("WrapMirror",Snoop,opWrapSurface);
//new G4LogicalBorderSurface("WrapSurface",
  //                             Scint_Wrap_phys,expHall_phys,opWrapSurface);
  G4OpticalSurface* opWrapSurface2 = new G4OpticalSurface("WrapSurface2");
  opWrapSurface2->SetType(dielectric_metal);
  opWrapSurface2->SetFinish(polishedfrontpainted);
  opWrapSurface2->SetModel(glisur);
   opWrapSurface2->SetMaterialPropertiesTable(ScintWrap_MPT);
   // new G4LogicalBorderSurface("WrapSurface",
   //                            Scint_Wrap_phys,expHall_phys,opWrapSurface);
//Guide_Wrapping
//
  G4OpticalSurface* opGuideWrapSurface = new G4OpticalSurface("GuideWrapSurface");
  opGuideWrapSurface->SetType(dielectric_metal);
  opGuideWrapSurface->SetFinish(polishedfrontpainted);
  opGuideWrapSurface->SetModel(glisur);
  opGuideWrapSurface->SetMaterialPropertiesTable(ScintWrap_MPT);
  
  //new G4LogicalBorderSurface("GuideWrapSurface",
  //                             Guide_Wrap_phys,expHall_phys,opGuideWrapSurface);

// Light Guide
//
  G4OpticalSurface* opGuideSurface = new G4OpticalSurface("GuideSurface");
  opGuideSurface->SetType(dielectric_dielectric);
  opGuideSurface->SetFinish(polishedbackpainted);
  opGuideSurface->SetModel(unified);

  //new G4LogicalBorderSurface("GuideSurface",
  //                             Guide_phys,expHall_phys,opGuideSurface);

//Window
//
  G4OpticalSurface* opWindowSurface = new G4OpticalSurface("WindowSurface");
  opWindowSurface->SetType(dielectric_dielectric);
  opWindowSurface->SetFinish(polishedbackpainted);
  opWindowSurface->SetModel(unified);

  //new G4LogicalBorderSurface("WindowSurface",
  //                             Window_Pys,expHall_phys,opWindowSurface);
  
//SiPM
//
  G4OpticalSurface* opSiPMSurface = new G4OpticalSurface("SiPMSurface");
  opSiPMSurface->SetType(dielectric_dielectric);
  opSiPMSurface->SetFinish(polishedbackpainted);
  opSiPMSurface->SetModel(unified);

  //new G4LogicalBorderSurface("SiPMSurface", fSiPM,Guide_phys,opSiPMSurface);
// Scintillator
//
  G4OpticalSurface* opScintSurface = new G4OpticalSurface("ScintSurface");
  opScintSurface->SetType(dielectric_dielectric);
  opScintSurface->SetFinish(polished);
  opScintSurface->SetModel(unified);

  //G4LogicalSkinSurface* airSurface =
  //      new G4LogicalSkinSurface("ScintSurface", Scint_log, opScintSurface);

  //G4OpticalSurface* opticalSurface = dynamic_cast <G4OpticalSurface*>
  //    (airSurface->GetSurface(Scint_log)->GetSurfaceProperty());

  //if (opticalSurface) opticalSurface->DumpInfo();

//
// Generate & Add Material Properties Table attached to the optical surfaces
//
  const G4int num = 2;
  G4double ephoton[num] = {2.034*eV, 4.136*eV};

  //OpticalGuideSurface
  G4double refractiveIndex[num] = {1.48, 1.48};
  G4double specularLobe[num]    = {0.3, 0.3};
  G4double specularSpike[num]   = {0.2, 0.2};
  G4double backScatter[num]     = {0.2, 0.2};

  G4MaterialPropertiesTable* myST1 = new G4MaterialPropertiesTable();

  myST1->AddProperty("RINDEX",                ephoton, refractiveIndex, num)->SetSpline(true);
  myST1->AddProperty("SPECULARLOBECONSTANT",  ephoton, specularLobe,    num);
  myST1->AddProperty("SPECULARSPIKECONSTANT", ephoton, specularSpike,   num);
  myST1->AddProperty("BACKSCATTERCONSTANT",   ephoton, backScatter,     num);

  G4cout << "SiPM Surface G4MaterialPropertiesTable" << G4endl;
  myST1->DumpTable();

  opGuideSurface->SetMaterialPropertiesTable(myST1);
  opWindowSurface->SetMaterialPropertiesTable(myST1);
  opSiPMSurface->SetMaterialPropertiesTable(myST1);

  //OpticalScintSurface
  G4double reflectivity[num] = {0.3, 0.5};
  G4double efficiency[num]   = {0.8,1.0};//{0.8, 1.0};
  G4double reflectivity2[num] = {1.0, 1.0};

  G4MaterialPropertiesTable *myST2 = new G4MaterialPropertiesTable();

  myST2->AddProperty("REFLECTIVITY", ephoton, reflectivity, num)->SetSpline(true);
  myST2->AddProperty("EFFICIENCY",   ephoton, efficiency,   num)->SetSpline(true);

  G4cout << "Air Surface G4MaterialPropertiesTable" << G4endl;
  myST2->DumpTable();

  opScintSurface->SetMaterialPropertiesTable(myST2);
  G4MaterialPropertiesTable *myST3 = new G4MaterialPropertiesTable();

  myST3->AddProperty("REFLECTIVITY", ephoton, reflectivity2, num)->SetSpline(true);
  myST3->AddProperty("EFFICIENCY",   ephoton, efficiency,   num)->SetSpline(true);

  G4cout << "Air Surface G4MaterialPropertiesTable" << G4endl;
  myST3->DumpTable();

  opWrapSurface->SetMaterialPropertiesTable(myST3);
  opWrapSurface2->SetMaterialPropertiesTable(myST3);
  G4MaterialPropertiesTable *myST4 = new G4MaterialPropertiesTable();

  myST4->AddProperty("REFLECTIVITY", ephoton, reflectivity2, num)->SetSpline(true);
  myST4->AddProperty("EFFICIENCY",   ephoton, efficiency,   num)->SetSpline(true);

  G4cout << "Air Surface G4MaterialPropertiesTable" << G4endl;
  myST4->DumpTable();

  opGuideWrapSurface->SetMaterialPropertiesTable(myST4);



//always return the physical World
  return expHall_phys;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
