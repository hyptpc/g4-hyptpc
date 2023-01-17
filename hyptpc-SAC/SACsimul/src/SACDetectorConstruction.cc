#include "SACDetectorConstruction.hh"
#include "AeroSD.hh"
#include "PMT15SD.hh"



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


SACDetectorConstruction::SACDetectorConstruction()
  : G4VUserDetectorConstruction()
{
}


SACDetectorConstruction::~SACDetectorConstruction()
{
  for (auto visAttributes: fVisAttributes){
    delete visAttributes;
  }
}


G4VPhysicalVolume* SACDetectorConstruction::Construct()
{

  G4NistManager *nist = G4NistManager::Instance();
  G4bool checkOverlaps = true;

  //Material---------------------------------------------------------------

  G4Material* world_mat = nist-> FindOrBuildMaterial("G4_AIR");

  G4Element* C = new G4Element("Carbon","C",6.,12.011*g/mole);
  G4Element* H = new G4Element("Hydrogen","H",1.,1.00794*g/mole);
  G4Element* O = new G4Element("Oxygen","O",15.,15.9994*g/mole);
  G4Element* B = new G4Element("Boron","B",5.,10.811*g/mole);
  G4Element* Na = new G4Element("Na","Na",11.,23.0*g/mole);
  G4Element* Si = new G4Element("Silicon","Si",14.,28.0844*g/mole);
  G4Element* Cl = new G4Element("Chlorine","Cl",17.,35.453*g/mole);
  G4Element* K = new G4Element("Potassium","K",19.,39.093*g/mole);
  G4Element* F = new G4Element("Fluorine","F",9.,18.998*g/mole);
  

  G4Material *Aerogel = new G4Material("Aerogel",0.2000*g/cm3,2);
  Aerogel->AddElement(Si,1);
  Aerogel->AddElement(O,2);

  G4Material* Teflon = new G4Material("Teflon",2.2*g/cm3,2);
  Teflon->AddElement(C,2);
  Teflon->AddElement(F,4);

  G4Material *Boro = new G4Material("Boro",1.19*g/cm3,3);
  Boro->AddElement(C,5);
  Boro->AddElement(H,8);
  Boro->AddElement(O,2);

  G4Material *UVT = new G4Material("UVT",1.19*g/cm3,3);
  UVT->AddElement(C,5);
  UVT->AddElement(H,8);
  UVT->AddElement(O,2);
  



  //Property--------------------------------------------------------------

  //air property---------------------------------------------------------

   
  G4MaterialPropertiesTable* prop_air = new G4MaterialPropertiesTable();
  G4double air_ep[] = {1.3*eV,7.*eV};
  G4double air_rindex[] = {1.0,1.0};
  prop_air->AddProperty("RINDEX",air_ep,air_rindex,2)->SetSpline(true);
  world_mat->SetMaterialPropertiesTable(prop_air);


  //Aerogel Property------------------------------------------------------
  G4MaterialPropertiesTable* prop_aerogel = new G4MaterialPropertiesTable();
  G4double aerogel_ep[] = {1.3*eV,7.*eV};
  G4double aerogel_ep1[] = {1.3*eV,1.56*eV,1.68*eV,1.84*eV,2.06*eV,2.26*eV,2.54*eV,2.90*eV,3.10*eV,3.28*eV,3.94*eV,4.94*eV,7.0*eV};
  G4double aerogel_abs[] = {500*mm,128*mm,120*mm,97*mm,77*mm,59*mm,41*mm,26*mm,20*mm,17*mm,8*mm,4*mm,1*mm};
  G4double factor = 8;
  for(int i=0;i<13;i++)aerogel_abs[i]*=factor;
  G4double aerogel_rindex[]={1.05,1.05};
  assert(sizeof(aerogel_ep_abs)==sizeof(aerogel_abs));
  
  prop_aerogel->AddProperty("RINDEX",aerogel_ep,aerogel_rindex,2)->SetSpline(true);
  prop_aerogel->AddProperty("ABSLENGTH",aerogel_ep1,aerogel_abs,13)->SetSpline(true);
  Aerogel->SetMaterialPropertiesTable(prop_aerogel);


  //Teflon property
  G4MaterialPropertiesTable* prop_teflon = new G4MaterialPropertiesTable();
  
  G4double teflon_rindex[2]={1.35,1.35};
  G4double teflon_ep[] = {1.3*eV,7.*eV};
  G4double teflon_abs[] = {1.0e-9*cm,1.0e-9*cm};

  prop_teflon->AddProperty("RINDEX",teflon_ep,teflon_rindex,2)->SetSpline(true);
  prop_teflon->AddProperty("ABSLENGTH",teflon_ep,teflon_abs,2)->SetSpline(true);
  Teflon->SetMaterialPropertiesTable(prop_teflon);

  
  
  //PMT window (glass) property -> UV enhenced
  const int n_uvt_abs = 23;
  G4double refractiveIndex_ch2[2] = { 1.49, 1.49};
  double uvt_pe[n_uvt_abs] = { 3.28216*eV, 3.38243*eV, 3.46465*eV,
			       3.57658*eV, 3.68828*eV, 3.80135*eV, 3.91291*eV, 4.03121*eV,
			       4.13058*eV, 4.22631*eV, 4.3463*eV, 4.43332*eV, 4.5239*eV,
			       4.61653*eV, 4.7041*eV, 4.78211*eV, 4.83997*eV, 4.89925*eV,
			       4.94812*eV, 4.98392*eV, 5.03451*eV, 5.05916*eV, 5.08405*eV };
  double uvt_abs[n_uvt_abs] = { 14.715*mm, 14.6688*mm, 14.6229*mm,
				11.4544*mm, 9.42164*mm, 6.96188*mm, 5.25449*mm, 3.77347*mm,
				2.86657*mm, 2.36169*mm, 1.86815*mm, 1.5681*mm, 1.31777*mm,
				1.11303*mm, 0.971912*mm, 0.835886*mm, 0.748513*mm, 0.661041*mm,
				0.595097*mm, 0.538428*mm, 0.476456*mm, 0.415598*mm, 0.357541*mm  };

  G4MaterialPropertiesTable* mpt_uvt = new G4MaterialPropertiesTable();
  mpt_uvt->AddProperty("RINDEX",teflon_ep,refractiveIndex_ch2,2)->SetSpline(true);
  mpt_uvt->AddProperty("ABSLENGTH",uvt_pe,uvt_abs,n_uvt_abs)->SetSpline(true);
  UVT->SetMaterialPropertiesTable(mpt_uvt);
  

  // PMT window (glass) property -> borosilicate (normal)
  const int n_bsg_abs = 27;
  double bsg_pe[n_bsg_abs] = {2.52649*eV, 2.56198*eV, 
			      2.60559*eV, 2.65071*eV, 2.70801*eV, 2.77157*eV, 2.83105*eV, 
			      2.89652*eV, 2.9566*eV, 3.02365*eV, 3.07922*eV, 3.1329*eV, 
			      3.18766*eV, 3.23422*eV, 3.27609*eV, 3.32529*eV, 3.37048*eV, 
			      3.41034*eV, 3.4454*eV, 3.47534*eV, 3.49986*eV, 3.52473*eV, 
			      3.54894*eV, 3.56835*eV, 3.58693*eV, 3.59942*eV, 3.61833*eV  };
  double bsg_abs[n_bsg_abs] = {15.3173*mm, 17.8103*mm, 
			       21.2352*mm, 21.1411*mm, 17.6112*mm, 13.2111*mm, 10.5909*mm, 
			       8.82763*mm, 6.63541*mm, 4.8361*mm, 3.80526*mm, 3.04797*mm, 
			       2.42786*mm, 2.01705*mm, 1.75316*mm, 1.44645*mm, 1.21704*mm, 
			       1.05022*mm, 0.931761*mm, 0.831044*mm, 0.749793*mm, 0.683042*mm, 
			       0.616324*mm, 0.552945*mm, 0.503814*mm, 0.45659*mm, 0.355462*mm};
  
  G4MaterialPropertiesTable* mpt_bsg = new G4MaterialPropertiesTable();
  mpt_bsg->AddProperty("RINDEX",teflon_ep,refractiveIndex_ch2,2)->SetSpline(true);  // temporary
  mpt_bsg->AddProperty("ABSLENGTH",bsg_pe,bsg_abs,n_bsg_abs)->SetSpline(true); // 
  Boro->SetMaterialPropertiesTable(mpt_bsg);
  

  



  //Geometry---------------------------------------------------------------

  //World------------------------------------------------------------------
  //Total World Size
  G4double world_size = 1*m;
  G4Box* solidWorld = new G4Box("World",world_size, world_size, world_size); 
  G4LogicalVolume* logicWorld = new G4LogicalVolume(solidWorld, world_mat, "World");
  physWorld = new G4PVPlacement(0,G4ThreeVector(), logicWorld, "World",0,false,0,checkOverlaps);


  //Various size information
  //Real Aerogel Size 
  G4double Aerox_real = 145.0 *mm;    
  G4double Aeroy_real = 113.0 *mm;
  G4double Aeroz_real = 33*mm;

  //Since there is air gap between Aerogel blocks and the Teflon sheet, I also considered that effect.
  G4double air_gap = 1*mm;
  G4double Aerox = Aerox_real + air_gap*2; //Air gap is at both side
  G4double Aeroy = Aeroy_real + air_gap*2;
  G4double Aeroz = Aeroz_real + air_gap*2;

  //Teflon size
  G4double teflon_thick = 0.5*mm;

  //Size related to PMT
  G4double pmt_dia = 12.5*mm;
  G4double pmt_thick = 1*mm;
  G4double pmt_position = 41.5*mm;
  

  //G4int numRZ = 30;
  G4int numRZ = 150;
  G4double win_thick = 1*mm;
  
  G4double re_thick = 1*mm;
  G4double mppc_place = 70*mm;
    



  //Aerogel 
  G4Box* Aero = new G4Box("Aero",Aerox_real/2,Aeroy_real/2,Aeroz_real/2);
  AeroLW = new G4LogicalVolume(Aero,Aerogel,"Aero");
  new G4PVPlacement(0,G4ThreeVector(0*mm,0*mm,0*mm),AeroLW,"Aero",logicWorld,false,9999,checkOverlaps);


  //Reflector - Teflon
  G4Box* Reflect = new G4Box("Reflect",Aerox/2,Aeroy/2,teflon_thick/2);
  ReflectLW = new G4LogicalVolume(Reflect,Teflon,"Reflect");

  G4Box* Reflect1 = new G4Box("Reflect1",Aeroz/2,Aeroy/2,teflon_thick/2);
  Reflect1LW = new G4LogicalVolume(Reflect1,Teflon,"Reflect1");

  G4RotationMatrix *rotR = new G4RotationMatrix();
  rotR->rotateY(90*degree);
    
  new G4PVPlacement(0,G4ThreeVector(0*mm,0*mm,Aeroz/2+teflon_thick/2),ReflectLW,"Reflect",logicWorld,false,0,checkOverlaps);
  new G4PVPlacement(0,G4ThreeVector(0*mm,0*mm,-Aeroz/2-teflon_thick/2),ReflectLW,"Reflect",logicWorld,false,0,checkOverlaps);

  new G4PVPlacement(rotR,G4ThreeVector(Aerox/2+teflon_thick/2,0*mm,0*mm),Reflect1LW,"Reflect1",logicWorld,false,0,checkOverlaps);
  new G4PVPlacement(rotR,G4ThreeVector(-Aerox/2-teflon_thick/2,0*mm,0*mm),Reflect1LW,"Reflect1",logicWorld,false,0,checkOverlaps);

  
  

  
  //PMT
  G4Tubs* PMT15 = new G4Tubs("PMT15",0*mm,pmt_dia,pmt_thick/2,0,2*TMath::Pi());
  PMT15LW = new G4LogicalVolume(PMT15,Boro,"PMT15");
  G4RotationMatrix *rotM = new G4RotationMatrix();
  rotM->rotateX(90*degree);
  
  new G4PVPlacement(rotM,G4ThreeVector(-pmt_position,Aeroy/2+pmt_thick/2,0),PMT15LW,"PMT15",logicWorld,false,1,checkOverlaps);
  new G4PVPlacement(rotM,G4ThreeVector(0*mm,Aeroy/2+pmt_thick/2,0),PMT15LW,"PMT15",logicWorld,false,2,checkOverlaps);
  new G4PVPlacement(rotM,G4ThreeVector(pmt_position,Aeroy/2+pmt_thick/2,0),PMT15LW,"PMT15",logicWorld,false,3,checkOverlaps);
  new G4PVPlacement(rotM,G4ThreeVector(-pmt_position,-Aeroy/2-pmt_thick/2,0),PMT15LW,"PMT15",logicWorld,false,4,checkOverlaps);
  new G4PVPlacement(rotM,G4ThreeVector(0*mm,-Aeroy/2-pmt_thick/2,0),PMT15LW,"PMT15",logicWorld,false,5,checkOverlaps);
  new G4PVPlacement(rotM,G4ThreeVector(pmt_position,-Aeroy/2-pmt_thick/2,0),PMT15LW,"PMT15",logicWorld,false,6,checkOverlaps);


  
  //Frame around PMT - Since Teflon sheet is also attached at the frame around PMT, I put the matieral as teflon, you can change it!
  G4Box* Frame_be = new G4Box("Frame_be",Aerox/2,pmt_thick/2,Aeroz/2);
  G4Tubs* PMThole = new G4Tubs("PMThole",0*mm,pmt_dia,pmt_thick/2,0,2*TMath::Pi());
    
  G4SubtractionSolid* Frame_1 = new G4SubtractionSolid("Frame_1",Frame_be,PMThole,rotM,G4ThreeVector(-pmt_position,0,0));
  G4SubtractionSolid* Frame_2 = new G4SubtractionSolid("Frame_2",Frame_1,PMThole,rotM,G4ThreeVector(0,0,0));
  G4SubtractionSolid* Frame = new G4SubtractionSolid("Frame",Frame_2,PMThole,rotM,G4ThreeVector(pmt_position,0,0));

  FrameLW = new G4LogicalVolume(Frame,Teflon,"Frame");
  new G4PVPlacement(0,G4ThreeVector(0,Aeroy/2+pmt_thick/2,0),FrameLW,"Frame",logicWorld,false,0,checkOverlaps);
  new G4PVPlacement(0,G4ThreeVector(0,-Aeroy/2-pmt_thick/2,0),FrameLW,"Frame",logicWorld,false,0,checkOverlaps);  
  



  //visattributes------------------------------------------------
  auto visAttributes = new G4VisAttributes(G4Colour(1.0, 1.0, 1.0));
  visAttributes -> SetVisibility(false);
  logicWorld->SetVisAttributes(visAttributes);
  fVisAttributes.push_back(visAttributes);

  
  visAttributes = new G4VisAttributes(G4Color::Blue()); 
  AeroLW->SetVisAttributes(visAttributes);
  PMT15LW->SetVisAttributes(visAttributes);
  fVisAttributes.push_back(visAttributes);
  

  visAttributes = new G4VisAttributes(G4Colour(0.5,0.5,0.0));
  fVisAttributes.push_back(visAttributes);





  //Surface--------------------------------------------------------------


  //teflon surface--------------------------------------------------
  G4OpticalSurface* surface_teflon = new G4OpticalSurface("surface_teflon");

  surface_teflon->SetModel(DAVIS);
  surface_teflon->SetType(dielectric_LUTDAVIS);
  surface_teflon->SetFinish(RoughTeflon_LUT);

  G4MaterialPropertiesTable* sp_teflon = new G4MaterialPropertiesTable();

  G4double teflon_ep1[] = {1.3*eV,1.56*eV,1.61*eV,1.74*eV,1.90*eV,2.05*eV,2.22*eV,2.34*eV,5.42*eV,7*eV};
  G4double teflon_reflec[] = {0.85,0.91,0.93,0.95,0.97,0.98,1,1,1,1};

  sp_teflon->AddProperty("REFLECTIVITY",teflon_ep1,teflon_reflec,10)->SetSpline(true);

  surface_teflon->SetMaterialPropertiesTable(sp_teflon);
  

  new G4LogicalSkinSurface("teflon_surface",ReflectLW,surface_teflon);
  new G4LogicalSkinSurface("teflon_surface",Reflect1LW,surface_teflon);
  new G4LogicalSkinSurface("teflon_surface",FrameLW,surface_teflon);

  


 


  return physWorld;
}		    


void SACDetectorConstruction::ConstructSDandField()
{
  
  auto aeroSD = new AeroSD("aeroSD");
  G4SDManager::GetSDMpointer()->AddNewDetector(aeroSD);
  AeroLW->SetSensitiveDetector(aeroSD);

  auto pmt15SD = new PMT15SD("pmt15SD");
  G4SDManager::GetSDMpointer()->AddNewDetector(pmt15SD);
  PMT15LW->SetSensitiveDetector(pmt15SD);


}



  
