#include "BACDetectorConstruction_167_2.hh"
#include "AeroSD.hh"
#include "MPPCSD.hh"



#include "G4Colour.hh"
#include "G4NistManager.hh"
#include "G4Box.hh"
#include "G4ExtrudedSolid.hh"
#include "G4LogicalVolume.hh"
#include "G4PVPlacement.hh"
#include "G4SystemOfUnits.hh"
#include "G4Element.hh"
#include "G4Material.hh"
#include "G4UnionSolid.hh"
#include "G4Polycone.hh"
#include "G4SubtractionSolid.hh"
#include "G4Trap.hh"
#include "G4SDManager.hh"
#include "G4VSensitiveDetector.hh"
#include "G4RunManager.hh"
#include "G4VPhysicalVolume.hh"
#include "G4VisAttributes.hh"
#include "G4ios.hh"
#include "TString.h"
#include "TMath.h"
#include "G4Tubs.hh"
#include "G4Trd.hh"
#include "G4OpticalSurface.hh"
#include "G4LogicalSkinSurface.hh"
#include "G4RotationMatrix.hh"



BACDetectorConstruction_167_2::BACDetectorConstruction_167_2()
  : G4VUserDetectorConstruction()
{
}

BACDetectorConstruction_167_2::~BACDetectorConstruction_167_2()
{
  for (auto visAttributes: fVisAttributes){
    delete visAttributes;
  }
}


G4VPhysicalVolume* BACDetectorConstruction_167_2::Construct()
{


  G4NistManager *nist = G4NistManager::Instance();
  G4bool checkOverlaps = true;

  //Material---------------------------------------------------------------

  G4Material* world_mat = nist-> FindOrBuildMaterial("G4_AIR");
  G4Material* Al = nist-> FindOrBuildMaterial("G4_Al");

  G4Element* C = new G4Element("Carbon","C",6.,12.011*g/mole);
  G4Element* H = new G4Element("Hydrogen","H",1.,1.00794*g/mole);
  G4Element* O = new G4Element("Oxygen","O",15.,15.9994*g/mole);
  //G4Element* Al = new G4Element("Aluminum","Al",13.,26.981538*g/mole);
  G4Element* B = new G4Element("Boron","B",5.,10.811*g/mole);
  G4Element* Na = new G4Element("Na","Na",11.,23.0*g/mole);
  G4Element* Si = new G4Element("Silicon","Si",14.,28.0844*g/mole);
  G4Element* Cl = new G4Element("Chlorine","Cl",17.,35.453*g/mole);
  G4Element* K = new G4Element("Potassium","K",19.,39.093*g/mole);
  G4Element* F = new G4Element("Fluorine","F",9.,18.998*g/mole);
  


  G4Material *Aerogel1 = new G4Material("Aerogel1",0.2000*g/cm3,2);
  Aerogel1->AddElement(Si,1);
  Aerogel1->AddElement(O,2);

  G4Material *Aerogel2 = new G4Material("Aerogel2",0.2000*g/cm3,2);
  Aerogel2->AddElement(Si,1);
  Aerogel2->AddElement(O,2);

  G4Material *Aerogel3 = new G4Material("Aerogel3",0.2000*g/cm3,2);
  Aerogel3->AddElement(Si,1);
  Aerogel3->AddElement(O,2);

  G4Material* blacksheet = new G4Material("blacksheet",0.95*g/cm3,2);
  blacksheet->AddElement(C,1);
  blacksheet->AddElement(H,2);

  G4Material *Mylar = new G4Material("Mylar", 1.39*g/cm3, 3);
  Mylar->AddElement(C, 5);
  Mylar->AddElement(H, 4);
  Mylar->AddElement(O, 2);

  G4Material* Tyvek = new G4Material("Tyvek",0.38*g/cm3,2);
  Tyvek->AddElement(C,1);
  Tyvek->AddElement(H,2);

  G4Material* Epoxi = new G4Material("Epoxi",1.1*g/cm3,4);
  Epoxi->AddElement(C, 21);
  Epoxi->AddElement(H, 25);
  Epoxi->AddElement(O,  5);
  Epoxi->AddElement(Cl, 1);


  G4Material* Acrylic = new G4Material("Acrylic", 1.19*g/cm3, 3);
  Acrylic->AddElement(C, 5);
  Acrylic->AddElement(H, 8);
  Acrylic->AddElement(O, 2);

  G4Material* Teflon = new G4Material("Teflon",2.2*g/cm3,2);
  Teflon->AddElement(C,2);
  Teflon->AddElement(F,4);


  //Property--------------------------------------------------------------

  //air property---------------------------------------------------------

   
  G4MaterialPropertiesTable* prop_air = new G4MaterialPropertiesTable();
  G4double air_ep[] = {1.3*eV,7.*eV};
  G4double air_rindex[] = {1.0,1.0};
  prop_air->AddProperty("RINDEX",air_ep,air_rindex,2)->SetSpline(true);
  world_mat->SetMaterialPropertiesTable(prop_air);



  //Aerogel1 Property------------------------------------------------------
  G4MaterialPropertiesTable* prop_aerogel1 = new G4MaterialPropertiesTable();


  G4double factor = 8;

  G4double aerogel1_ep[] = {1.3*eV,7.*eV};
  G4double aerogel_ep[] = {1.3*eV,1.56*eV,1.68*eV,1.84*eV,2.06*eV,2.26*eV,2.54*eV,2.90*eV,3.10*eV,3.28*eV,3.94*eV,4.94*eV,7.0*eV};
  
  G4double aerogel1_abs[] = {500*mm,128*mm,120*mm,97*mm,77*mm,59*mm,41*mm,26*mm,20*mm,17*mm,8*mm,4*mm,1*mm};
  for(int i=0;i<13;i++)aerogel1_abs[i]*=factor;
  G4double aerogel1_rindex[]={1.167,1.167};

  //G4double aerogel_ray[] = {6.16*pow(10,10),6.16*pow(10,10)};

  assert(sizeof(aerogel1_ep_abs)==sizeof(aerogel1_abs));
  
  prop_aerogel1->AddProperty("RINDEX",aerogel1_ep,aerogel1_rindex,2)->SetSpline(true);
  prop_aerogel1->AddProperty("ABSLENGTH",aerogel_ep,aerogel1_abs,13)->SetSpline(true);
  prop_aerogel1->AddConstProperty("RESOLUTIONSCALE",1.0);
  prop_aerogel1->AddConstProperty("FASTTIMECONSTANT", 1.*ns);
  prop_aerogel1->AddConstProperty("SLOWTIMECONSTANT",2.8*ns);
  prop_aerogel1->AddConstProperty("YIELDRATIO",0.8);
  Aerogel1->SetMaterialPropertiesTable(prop_aerogel1);

  //Aerogel2 Property------------------------------------------------------
  G4MaterialPropertiesTable* prop_aerogel2 = new G4MaterialPropertiesTable();

  

  G4double aerogel2_ep[] = {1.3*eV,7.*eV};
  G4double aerogel2_abs[] = {500*mm,128*mm,120*mm,97*mm,77*mm,59*mm,41*mm,26*mm,20*mm,17*mm,8*mm,4*mm,1*mm};
  for(int i=0;i<13;i++)aerogel2_abs[i]*=factor;

  G4double aerogel2_rindex[]={1.167,1.167};

  //G4double aerogel_ray[] = {6.16*pow(10,10),6.16*pow(10,10)};

  assert(sizeof(aerogel2_ep_abs)==sizeof(aerogel2_abs));
  
  prop_aerogel2->AddProperty("RINDEX",aerogel2_ep,aerogel2_rindex,2)->SetSpline(true);
  prop_aerogel2->AddProperty("ABSLENGTH",aerogel_ep,aerogel2_abs,13)->SetSpline(true);
  prop_aerogel2->AddConstProperty("RESOLUTIONSCALE",1.0);
  prop_aerogel2->AddConstProperty("FASTTIMECONSTANT", 1.*ns);
  prop_aerogel2->AddConstProperty("SLOWTIMECONSTANT",2.8*ns);
  prop_aerogel2->AddConstProperty("YIELDRATIO",0.8);
  Aerogel2->SetMaterialPropertiesTable(prop_aerogel2);

  //Aerogel3 Property------------------------------------------------------
  G4MaterialPropertiesTable* prop_aerogel3 = new G4MaterialPropertiesTable();

  

  G4double aerogel3_ep[] = {1.3*eV,7.*eV};
  G4double aerogel3_abs[] = {500*mm,128*mm,120*mm,97*mm,77*mm,59*mm,41*mm,26*mm,20*mm,17*mm,8*mm,4*mm,1*mm};
  for(int i=0;i<13;i++)aerogel3_abs[i]*=factor;
  G4double aerogel3_rindex[]={1.167,1.167};

  //G4double aerogel_ray[] = {6.16*pow(10,10),6.16*pow(10,10)};

  assert(sizeof(aerogel3_ep_abs)==sizeof(aerogel3_abs));
  
  prop_aerogel3->AddProperty("RINDEX",aerogel3_ep,aerogel3_rindex,2)->SetSpline(true);
  prop_aerogel3->AddProperty("ABSLENGTH",aerogel_ep,aerogel3_abs,13)->SetSpline(true);
  prop_aerogel3->AddConstProperty("RESOLUTIONSCALE",1.0);
  prop_aerogel3->AddConstProperty("FASTTIMECONSTANT", 1.*ns);
  prop_aerogel3->AddConstProperty("SLOWTIMECONSTANT",2.8*ns);
  prop_aerogel3->AddConstProperty("YIELDRATIO",0.8);
  Aerogel3->SetMaterialPropertiesTable(prop_aerogel3);


  
  

  //Black sheet property-------------------------------------------------
  G4MaterialPropertiesTable* prop_bs = new G4MaterialPropertiesTable();

  G4double bs_ep[] = {1.3*eV,7.*eV};
  G4double bs_abs[] = {1.0e-9*cm,1.0e-9*cm};
  G4double bs_rindex[]={1.6,1.6};
  
  prop_bs->AddProperty("RINDEX",bs_ep,bs_rindex,2)->SetSpline(true);
  prop_bs->AddProperty("ABSLENGTH",bs_ep,bs_abs,2)->SetSpline(true);
  blacksheet->SetMaterialPropertiesTable(prop_bs);


  //Mylar property
  G4double mylar_ep[]={1.4*eV,1.48*eV,1.52*eV,1.56*eV,1.6*eV,1.7*eV,1.8*eV,1.9*eV,2*eV,2.2*eV,2.4*eV,2.6*eV,2.8*eV,3*eV,3.4*eV,3.8*eV,4*eV,5*eV,6*eV,7*eV};

  G4double mylar_real[]={2.2802,2.6945,2.7668,2.7675,2.6154,2.1606,1.8301,1.5724,1.366,1.0728,0.8734,0.7278,0.6079,0.52135,0.39877,0.31474,0.28003,0.18137,0.12677,0.094236};
  G4double mylar_ima[]={8.1134,8.1878,8.2573,8.3866,8.4914,8.3565,8.0601,7.7354,7.4052,6.7839,6.2418,5.7781,5.3676,5.0008,4.3957,3.9165,3.7081,2.9029,2.3563,1.9519};

  G4double mylar_ep1[] = {1.3*eV,7.*eV};

  assert (sizeof(mylar_ep) == sizeof(mylar_real));
  assert (sizeof(mylar_ep) == sizeof(mylar_ima));
  
  const G4int numentries_mylar = sizeof(mylar_ep)/sizeof(G4double);
  G4double mylar_abs[] = {1.0e-9*cm,1.0e-9*cm};
  G4MaterialPropertiesTable* prop_mylar = new G4MaterialPropertiesTable();
  prop_mylar->AddProperty("REALRINDEX",mylar_ep,mylar_real,numentries_mylar)->SetSpline(true);
  prop_mylar->AddProperty("IMAGINARYRINDEX",mylar_ep,mylar_ima,numentries_mylar)->SetSpline(true);
  prop_mylar->AddProperty("ABSLENGTH",mylar_ep1,mylar_abs,2)->SetSpline(true);
  Mylar->SetMaterialPropertiesTable(prop_mylar);



  //Teflon property
  G4MaterialPropertiesTable* prop_teflon = new G4MaterialPropertiesTable();
  
  G4double teflon_rindex[2]={1.35,1.35};
  G4double teflon_ep[] = {1.6*eV,7.*eV};
  G4double teflon_abs[] = {1.0e-9*cm,1.0e-9*cm};

  prop_teflon->AddProperty("RINDEX",teflon_ep,teflon_rindex,2)->SetSpline(true);
  prop_teflon->AddProperty("ABSLENGTH",teflon_ep,teflon_abs,2)->SetSpline(true);
  Teflon->SetMaterialPropertiesTable(prop_teflon);

  
  //MPPC property
  G4MaterialPropertiesTable* prop_mppc = new G4MaterialPropertiesTable();
  
  G4double mppc_rindex[2]={1.,1.};
  G4double mppc_ep[] = {1.6*eV,7.*eV};
  G4double mppc_abs[] = {1.0*cm,1.0*cm};


  prop_mppc->AddProperty("RINDEX",mppc_ep,mppc_rindex,2)->SetSpline(true);
  prop_mppc->AddProperty("ABSLENGTH",mppc_ep,mppc_abs,2)->SetSpline(true);
  Al->SetMaterialPropertiesTable(prop_mppc);

  //MPPC surface (Epoxi) property
  G4MaterialPropertiesTable* prop_epoxi = new G4MaterialPropertiesTable();
 
  G4double epoxi_rindex[2]={1.5,1.5};
  G4double epoxi_ep[] = {1.3*eV,7.*eV};
  G4double epoxi_abs[] = {1.0*cm,1.0*cm};


  prop_epoxi->AddProperty("RINDEX",epoxi_ep,epoxi_rindex,2)->SetSpline(true);
  prop_epoxi->AddProperty("ABSLENGTH",epoxi_ep,epoxi_abs,2)->SetSpline(true);
  Epoxi->SetMaterialPropertiesTable(prop_epoxi);



  //Geometry---------------------------------------------------------------

  //World------------------------------------------------------------------

  
  G4double world_size = 1*m;
  G4Box* solidWorld = new G4Box("World",world_size, world_size, world_size); 
  G4LogicalVolume* logicWorld = new G4LogicalVolume(solidWorld, world_mat, "World");
  physWorld = new G4PVPlacement(0,G4ThreeVector(), logicWorld, "World",0,false,0,checkOverlaps);


  G4double Aero_space = 0.25*mm;
  
  G4double Aerox_real = 58 *mm;
  G4double Aeroy_real = 105*mm;
  G4double Aeroz_real = 10.0*mm;
  G4double Aerox = Aerox_real+Aero_space;
  G4double Aeroy = Aeroy_real+Aero_space;
  G4double Aeroz = Aeroz_real+Aero_space;
  
  
  G4double Ref_space = 25*mm;
  G4double Ref_spacez = 5*mm;


  G4double mppc_thick = 1*mm;
  G4double re_thick = 1*mm;

  G4Box* Reflect = new G4Box("Reflect",(Aerox*2+Ref_space)/2,(Aeroy+Ref_space)/2,(Aeroz*2+Ref_spacez)/2);
  ReflectLW = new G4LogicalVolume(Reflect,Teflon,"Reflect");
  new G4PVPlacement(0,G4ThreeVector(0*mm,0*mm,0*mm),ReflectLW,"Reflect",logicWorld,false,0,checkOverlaps);


  //Aerogel
  G4Box* Aero_air = new G4Box("Aero_air",Aerox/2,Aeroy/2,Aeroz/2);
  Aero_airLW = new G4LogicalVolume(Aero_air,world_mat,"Aero_air");
  new G4PVPlacement(0,G4ThreeVector(Aerox/2,0*mm,Aeroz/2),Aero_airLW,"Aero_air",ReflectLW,false,10,checkOverlaps);
  new G4PVPlacement(0,G4ThreeVector(-Aerox/2,0*mm,Aeroz/2),Aero_airLW,"Aero_air",ReflectLW,false,20,checkOverlaps);
  new G4PVPlacement(0,G4ThreeVector(Aerox/2+5*mm,0*mm,-Aeroz/2),Aero_airLW,"Aero_air",ReflectLW,false,30,checkOverlaps);
  new G4PVPlacement(0,G4ThreeVector(-Aerox/2+5*mm,0*mm,-Aeroz/2),Aero_airLW,"Aero_air",ReflectLW,false,40,checkOverlaps);
  
  G4Box* Aero = new G4Box("Aero",Aerox_real/2,Aeroy_real/2,Aeroz_real/2);
  AeroLW = new G4LogicalVolume(Aero,Aerogel1,"Aero");
  new G4PVPlacement(0,G4ThreeVector(0*mm,0*mm,0*mm),AeroLW,"Aero",Aero_airLW,false,123,checkOverlaps);
  




   
  //MPPC---------------------------------------------------------------------------

  //G4double sensor_pos = 36*mm;
  G4double sensor_pos = 50*mm;
  G4Box* MPPC = new G4Box("MPPC",12*mm,12*mm,mppc_thick/2);
  MPPCLW = new G4LogicalVolume(MPPC,Epoxi,"MPPC");
  G4RotationMatrix *rotM = new G4RotationMatrix();
  rotM->rotateX(90*degree);
  new G4PVPlacement(rotM,G4ThreeVector(-sensor_pos/2,Aeroy/2+mppc_thick/2,0),MPPCLW,"MPPC",ReflectLW,false,1,checkOverlaps);
  new G4PVPlacement(rotM,G4ThreeVector(sensor_pos/2,Aeroy/2+mppc_thick/2,0),MPPCLW,"MPPC",ReflectLW,false,2,checkOverlaps);
  //new G4PVPlacement(rotM,G4ThreeVector(40*mm,Aeroy_real/2+mppc_thick,0),MPPCLW,"MPPC",logicWorld,false,3,checkOverlaps);
  new G4PVPlacement(rotM,G4ThreeVector(-sensor_pos/2,-Aeroy/2-mppc_thick/2,0),MPPCLW,"MPPC",ReflectLW,false,4,checkOverlaps);
  new G4PVPlacement(rotM,G4ThreeVector(sensor_pos/2,-Aeroy/2-mppc_thick/2,0),MPPCLW,"MPPC",ReflectLW,false,5,checkOverlaps);
  //new G4PVPlacement(rotM,G4ThreeVector(40*mm,-Aeroy_real/2-mppc_thick,0),MPPCLW,"MPPC",logicWorld,false,6,checkOverlaps);
 

 

  //visattributes------------------------------------------------
  auto visAttributes = new G4VisAttributes(G4Colour(1.0, 1.0, 1.0));
  visAttributes -> SetVisibility(false);
  logicWorld->SetVisAttributes(visAttributes);
  fVisAttributes.push_back(visAttributes);

  
  visAttributes = new G4VisAttributes(G4Color::Blue()); 
  AeroLW->SetVisAttributes(visAttributes);
  fVisAttributes.push_back(visAttributes);



  

  visAttributes = new G4VisAttributes(G4Colour(0.5,0.5,0.0));
  fVisAttributes.push_back(visAttributes);





  //Surface--------------------------------------------------------------
  //black sheet surface------------------
  
  G4OpticalSurface* surface_bs = new G4OpticalSurface("surface_bs");
  surface_bs->SetType(dielectric_dielectric);
  surface_bs->SetModel(unified);
  surface_bs->SetFinish(polishedtyvekair);

  G4MaterialPropertiesTable* sp_bs = new G4MaterialPropertiesTable();

  G4double bs_effi[] = {0.0,0.0};
  G4double bs_reflec[]= {0.05,0.05};
  G4double bs_specularLobe[] = {0.3,0.3};
  G4double bs_specularSpike[] = {0.2,0.2};
  G4double bs_backScatter[] = {0,0};

  sp_bs->AddProperty("EFFICIENCY",bs_ep,bs_effi,2)->SetSpline(true);
  sp_bs->AddProperty("REFLECTIVITY",bs_ep,bs_reflec,2)->SetSpline(true);
  sp_bs->AddProperty("SPECULARLOBECONSTANT",bs_ep,bs_specularLobe,2)->SetSpline(true);
  sp_bs->AddProperty("SPECULARSPIKECONSTANT",bs_ep,bs_specularSpike,2)->SetSpline(true);
  sp_bs->AddProperty("BACKSCATTERCONSTANT",bs_ep,bs_backScatter,2)->SetSpline(true);
  surface_bs->SetMaterialPropertiesTable(sp_bs);


  
  
  

  //mylar_al surface-------------------
  G4OpticalSurface* surface_mylar = new G4OpticalSurface("surface_mylar");
  surface_mylar->SetType(dielectric_metal);
  surface_mylar->SetFinish(polished);    
  surface_mylar->SetModel(unified);
  
  G4MaterialPropertiesTable* sp_mylar = new G4MaterialPropertiesTable();
  G4double mylar_reflec[] = {0.80,0.80};  //for metal, reflectivity is calculated using rindex, they use polarization, angle, energy
  G4double mylar_effi[] = {0.0,0.0};
  G4double mylar_specularLobe[] = {0.85,0.85};
  G4double mylar_specularSpike[]={0.87,0.87};
  G4double mylar_backScatter[] = {0,0};
  sp_mylar->AddProperty("EFFICIENCY",mylar_ep1,mylar_effi,2)->SetSpline(true);
  sp_mylar->AddProperty("REFLECTIVITY",air_ep,mylar_reflec,2)->SetSpline(true);
  sp_mylar->AddProperty("SPECULARLOBECONSTANT",mylar_ep1,mylar_specularLobe,2)->SetSpline(true);
  sp_mylar->AddProperty("SPECULARSPIKECONSTANT",mylar_ep1,mylar_specularSpike,2)->SetSpline(true);
  sp_mylar->AddProperty("BACKSCATTERCONSTANT",mylar_ep1,mylar_backScatter,2)->SetSpline(true);
  surface_mylar->SetMaterialPropertiesTable(sp_mylar);




  
  //tyvek surface--------------------------------------------------
  G4OpticalSurface* surface_tyvek = new G4OpticalSurface("surface_tyvek");
  surface_tyvek->SetType(dielectric_dielectric);
  surface_tyvek->SetModel(unified);
  surface_tyvek->SetFinish(groundfrontpainted);

  //teflon surface--------------------------------------------------
  G4OpticalSurface* surface_teflon = new G4OpticalSurface("surface_teflon");
  /*
  surface_teflon->SetType(dielectric_dielectric);
  surface_teflon->SetModel(unified);
  surface_teflon->SetFinish(groundfrontpainted);
  */
  surface_teflon->SetModel(DAVIS);
  surface_teflon->SetType(dielectric_LUTDAVIS);
  surface_teflon->SetFinish(RoughTeflon_LUT);

  G4MaterialPropertiesTable* sp_teflon = new G4MaterialPropertiesTable();

  G4double teflon_ep1[] = {1.3*eV,1.56*eV,1.61*eV,1.74*eV,1.90*eV,2.05*eV,2.22*eV,2.34*eV,5.42*eV,7*eV};
  G4double teflon_reflec[] = {0.85,0.91,0.93,0.95,0.97,0.98,1,1,1,1};

  G4double teflon_effi[] = {0.0,0.0};

  //sp_teflon->AddProperty("EFFICIENCY",air_ep,teflon_effi,2)->SetSpline(true);
  sp_teflon->AddProperty("REFLECTIVITY",teflon_ep1,teflon_reflec,10)->SetSpline(true);

  surface_teflon->SetMaterialPropertiesTable(sp_teflon);
  

  new G4LogicalSkinSurface("teflon_surface",ReflectLW,surface_teflon);
  //new G4LogicalSkinSurface("teflon_surface",Reflect1LW,surface_teflon);
  //new G4LogicalSkinSurface("teflon_surface",FrameLW,surface_teflon);

  

 


  return physWorld;
}		    


void BACDetectorConstruction_167_2::ConstructSDandField()
{
  
  auto aeroSD = new AeroSD("aeroSD");
  G4SDManager::GetSDMpointer()->AddNewDetector(aeroSD);
  AeroLW->SetSensitiveDetector(aeroSD);


  auto mppcSD = new MPPCSD("mppcSD");
  G4SDManager::GetSDMpointer()->AddNewDetector(mppcSD);
  MPPCLW->SetSensitiveDetector(mppcSD);


  





}



  
