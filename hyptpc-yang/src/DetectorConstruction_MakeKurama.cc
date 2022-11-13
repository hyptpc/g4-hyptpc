/*
  DetectorConstruction_MakeKurama.cc
  make KURAMA magnet, ch1, dc1, dc2, dc3, ftof, fac, dummy_detectors
  2018/8  Yang
*/

#include "DetectorConstruction.hh"
#include "MaterialList.hh"
#include "SCField.hh"
#include "UNIField.hh"

#include "G4FieldManager.hh"
#include "G4ChordFinder.hh"
#include "G4TransportationManager.hh"

#include "G4Box.hh"
#include "G4Trd.hh"
#include "G4Tubs.hh"
#include "G4Sphere.hh"
#include "G4Paraboloid.hh"
#include "G4Torus.hh"
#include "G4NistManager.hh"
#include "G4SubtractionSolid.hh"
#include "G4UnionSolid.hh"
#include "G4IntersectionSolid.hh"
#include "G4BooleanSolid.hh"
#include "G4ThreeVector.hh"
#include "G4RotationMatrix.hh"
#include "G4LogicalVolume.hh"
#include "G4Polyhedra.hh"
#include "G4PVPlacement.hh"
#include "G4UniformMagField.hh"
#include "G4Material.hh"
#include "G4PhysicalConstants.hh"
#include "G4SystemOfUnits.hh"
#include "G4UserLimits.hh"

#include "ConfMan.hh"
//#include "DCGeomMan.hh"

//#include "DetectorSize.hh"
//#include "Tpc.hh"
#include "G4SDManager.hh"
#include "G4VSensitiveDetector.hh"
#include "TargetSD.hh"
#include "TofSD.hh"
#include "TpcSD.hh"
#include "FACSD.hh"
#include "CH1SD.hh"
#include "DCSD.hh"
#include "FTOFSD.hh"
#include "ScintSD.hh"

#include "G4VisAttributes.hh"
#include "G4Colour.hh"

#include "TString.h"


/// ///
///// KURAMA Magnet E07/////

void DetectorConstruction::MakeKuramaMagnetE07(G4VPhysicalVolume *pMother, G4ThreeVector &pos, G4RotationMatrix &rot)
{

  /// KURAMA Magnet Physics World ///
  /// Physics World ///
  // param //
  const G4double X_PW = 2200.0/2.0*mm;
  const G4double Y_PW=  1600.0/2.0*mm;
  const G4double Z_PW  = 1640.0/2.0*mm; 
  std::string fullNamePW = "KuramaMagnet_PhysicsWorld";
  // //
  
  G4Box* solidDetectorPW = new G4Box(G4String(fullNamePW.c_str()),
				     X_PW,
				     Y_PW,
				     Z_PW
				     );

  G4LogicalVolume* logicDetectorPW = new G4LogicalVolume(solidDetectorPW,
							 mat_air,
							 G4String(fullNamePW.c_str()));
  //logicDetectorPW->SetVisAttributes(new G4VisAttributes(G4Color::Blue()));
  
  G4RotationMatrix rot_PW = rot;
  G4ThreeVector TVPW = rot*pos;
  G4VPhysicalVolume *physMother = new G4PVPlacement(G4Transform3D(rot_PW, TVPW),
						    G4String(fullNamePW.c_str()),
						    logicDetectorPW,
						    pMother,
						    false,
						    0
						    );


  // Magnet Parameters // 
  const double kurama_gap=400.0;
  int XCOORD = 0;
  int YCOORD = 1;
  int ZCOORD = 2;
  G4double size_MFIELD[3];
  size_MFIELD[XCOORD] = 500.0*mm;
  size_MFIELD[YCOORD] = kurama_gap*mm;
  size_MFIELD[ZCOORD] = 400.0*mm;
  //size_MFIELD[XCOORD] = 0.0*mm;
  //size_MFIELD[YCOORD] = 0.0*mm;
  //size_MFIELD[ZCOORD] = 0.0*mm;

  G4double size_COIL1[3];
  size_COIL1[XCOORD] = 900./2*mm;
  size_COIL1[YCOORD] = 193./2*mm;
  size_COIL1[ZCOORD] = 280./2*mm;

  G4double size_COIL2[3];
  size_COIL2[XCOORD] = 193./2.*mm;
  size_COIL2[YCOORD] = 58.5*mm;
  size_COIL2[ZCOORD] = 280.0/2.*mm;

  G4double size_COIL3[3];
  size_COIL3[XCOORD] = 193.0/2.*mm;
  size_COIL3[YCOORD] = 68.5*mm;
  size_COIL3[ZCOORD] = 280.0/2.*mm;

  G4double size_COIL4[3];
  size_COIL4[XCOORD] = 193.0/2.*mm;
  size_COIL4[YCOORD] = 280.0/2.*mm;
  size_COIL4[ZCOORD] = 740.0/2.*mm;

  G4double size_COIL5[3];
  size_COIL5[XCOORD] = 900.0/2.*mm;
  size_COIL5[YCOORD] = 214./2.*mm;
  size_COIL5[ZCOORD] = 280./2.*mm;



  G4double size_YOKE_UD[3];
  size_YOKE_UD[XCOORD] = 2200.0/2.*mm;
  size_YOKE_UD[YCOORD] = 370.0/2.*mm;
  size_YOKE_UD[ZCOORD] = size_MFIELD[ZCOORD]*mm;

  G4double size_YOKE_LR[3];
  size_YOKE_LR[XCOORD] = 200.0*mm;
  size_YOKE_LR[YCOORD] = size_MFIELD[YCOORD]*mm;
  size_YOKE_LR[ZCOORD] = size_MFIELD[ZCOORD]*mm;

  G4double size_YOKE_LR_GapSpace[3]; //E42
  size_YOKE_LR_GapSpace[ZCOORD] = 400.0*mm;
  size_YOKE_LR_GapSpace[XCOORD] = 220.0/2.*mm;
  size_YOKE_LR_GapSpace[YCOORD] = 120.0/2.*mm;

  G4double size_UGUARD_UD[3];
  size_UGUARD_UD[XCOORD] = 950.0*mm;
  size_UGUARD_UD[YCOORD] = 310.0*mm;
  size_UGUARD_UD[ZCOORD] = 50.0*mm;

  G4double size_UGUARD_LR[3];
  size_UGUARD_LR[XCOORD] = 400.0*mm;
  size_UGUARD_LR[YCOORD] = 150.0*mm;
  size_UGUARD_LR[ZCOORD] = 50.0*mm;

  G4double size_DGUARD_UD[3];
  size_DGUARD_UD[XCOORD] = 800.0*mm;
  size_DGUARD_UD[YCOORD] = 210.0*mm;
  size_DGUARD_UD[ZCOORD] = 50.0*mm;

  G4double size_DGUARD_LR[3];
  size_DGUARD_LR[XCOORD] = 125.0*mm;
  size_DGUARD_LR[YCOORD] = size_MFIELD[YCOORD]*mm+150.0*mm;
  size_DGUARD_LR[ZCOORD] = 50.0*mm;

  G4double size_UGUARD_R1[3];
  size_UGUARD_R1[ZCOORD] =  50.0*mm;
  size_UGUARD_R1[XCOORD] = 590.0/2*mm;
  size_UGUARD_R1[YCOORD] = 30.0*mm;

  G4double size_UGUARD_R2[3];
  size_UGUARD_R2[ZCOORD] =  50.0*mm;
  size_UGUARD_R2[XCOORD] = 325.0*mm;
  size_UGUARD_R2[YCOORD] = 60.0*mm;
  
  G4double size_CH[3];
  size_CH[XCOORD] = 11.5*mm;
  size_CH[YCOORD] = 400.0*mm;
  size_CH[ZCOORD] = 2.0*mm;

  G4double size_FTOF[3];
  size_FTOF[XCOORD] = 40.0*mm;
  size_FTOF[YCOORD] = 900.0*mm;
  size_FTOF[ZCOORD] = 15.0*mm;


  G4double k_move_x=5.;
  G4double k_pos_z=1500.;

  G4double pos_MFIELD[3];
  pos_MFIELD[XCOORD] = 0.0*mm;
  pos_MFIELD[YCOORD] = 0.0*mm;
  pos_MFIELD[ZCOORD] = 0.0*mm;
  //pos_MFIELD[XCOORD] = k_move_x*mm;
  //pos_MFIELD[YCOORD] = 0.0*mm;
  //pos_MFIELD[ZCOORD] = k_pos_z*mm;
  //  pos_MFIELD[ZCOORD] = (1105.+400.)*mm;
  //  pos_MFIELD[ZCOORD] = 0.*mm;


  G4double pos_COIL1[3];
  pos_COIL1[XCOORD] = pos_MFIELD[XCOORD];
  pos_COIL1[YCOORD] = size_MFIELD[YCOORD]+size_YOKE_UD[YCOORD]*2. -(size_COIL1[YCOORD]+20.*mm);
  pos_COIL1[ZCOORD] = pos_MFIELD[ZCOORD]-size_MFIELD[ZCOORD]-(size_COIL1[ZCOORD]+20.*mm);

  G4double pos_COIL4L[3];
  pos_COIL4L[XCOORD] = pos_MFIELD[XCOORD]+size_MFIELD[XCOORD]+size_COIL4[XCOORD];
  pos_COIL4L[YCOORD] = size_MFIELD[YCOORD]/2;
  pos_COIL4L[ZCOORD] = pos_MFIELD[ZCOORD];


  G4double pos_COIL4R[3];
  pos_COIL4R[XCOORD] = pos_MFIELD[XCOORD]-size_MFIELD[XCOORD]-size_COIL4[XCOORD];
  pos_COIL4R[YCOORD] = size_MFIELD[YCOORD]/2;
  pos_COIL4R[ZCOORD] = pos_MFIELD[ZCOORD];


  G4double pos_COIL2L[3];
  pos_COIL2L[XCOORD] = pos_MFIELD[XCOORD]+size_MFIELD[XCOORD]+size_COIL4[XCOORD];
  pos_COIL2L[YCOORD] = (pos_COIL4L[YCOORD]+pos_COIL1[YCOORD])/2+(size_COIL4[YCOORD]-size_COIL1[YCOORD])/2;
  pos_COIL2L[ZCOORD] = pos_MFIELD[ZCOORD]-size_MFIELD[ZCOORD]-(size_COIL1[ZCOORD]+20.*mm);

  G4double pos_COIL2R[3];
  pos_COIL2R[XCOORD] = pos_MFIELD[XCOORD]-size_MFIELD[XCOORD]-size_COIL4[XCOORD];
  pos_COIL2R[YCOORD] = (pos_COIL4R[YCOORD]+pos_COIL1[YCOORD])/2+(size_COIL4[YCOORD]-size_COIL1[YCOORD])/2;
  pos_COIL2R[ZCOORD] = pos_MFIELD[ZCOORD]-size_MFIELD[ZCOORD]-(size_COIL1[ZCOORD]+20.*mm);


  G4double pos_COIL5[3];
  pos_COIL5[XCOORD] = pos_MFIELD[XCOORD];
  pos_COIL5[YCOORD] = size_MFIELD[YCOORD]+size_YOKE_UD[YCOORD]*2. -(size_COIL5[YCOORD]);
  pos_COIL5[ZCOORD] = pos_MFIELD[ZCOORD] + size_MFIELD[ZCOORD]+(size_COIL5[ZCOORD]+21.*mm);


  G4double pos_COIL3L[3];
  pos_COIL3L[XCOORD] = pos_MFIELD[XCOORD]+size_MFIELD[XCOORD]+size_COIL4[XCOORD];
  pos_COIL3L[YCOORD] = (pos_COIL4L[YCOORD]+pos_COIL5[YCOORD])/2+(size_COIL4[YCOORD]-size_COIL5[YCOORD])/2;
  pos_COIL3L[ZCOORD] = pos_MFIELD[ZCOORD]+size_MFIELD[ZCOORD]+(size_COIL5[ZCOORD]+21.*mm);

  G4double pos_COIL3R[3];
  pos_COIL3R[XCOORD] = pos_MFIELD[XCOORD]-size_MFIELD[XCOORD]-size_COIL4[XCOORD];
  pos_COIL3R[YCOORD] = (pos_COIL4R[YCOORD]+pos_COIL5[YCOORD])/2+(size_COIL4[YCOORD]-size_COIL5[YCOORD])/2;
  pos_COIL3R[ZCOORD] = pos_MFIELD[ZCOORD]+size_MFIELD[ZCOORD]+(size_COIL5[ZCOORD]+21.*mm);



  G4double pos_YOKE_U[3];
  pos_YOKE_U[XCOORD] = pos_MFIELD[XCOORD];
  pos_YOKE_U[YCOORD] = size_MFIELD[YCOORD]+size_YOKE_UD[YCOORD];
  pos_YOKE_U[ZCOORD] = pos_MFIELD[ZCOORD];

  G4double pos_YOKE_D[3];
  pos_YOKE_D[XCOORD] = pos_MFIELD[XCOORD];
  pos_YOKE_D[YCOORD] = -(size_MFIELD[YCOORD]+size_YOKE_UD[YCOORD]);
  pos_YOKE_D[ZCOORD] = pos_MFIELD[ZCOORD];

  G4double pos_YOKE_L[3];
  pos_YOKE_L[XCOORD] = pos_MFIELD[XCOORD] + size_MFIELD[XCOORD]+size_YOKE_LR[XCOORD]+200.*mm;
  pos_YOKE_L[YCOORD] = 0.0*mm;
  pos_YOKE_L[ZCOORD] = pos_MFIELD[ZCOORD];

  G4double pos_YOKE_R[3];
  pos_YOKE_R[XCOORD] = pos_MFIELD[XCOORD] - (size_MFIELD[XCOORD]+size_YOKE_LR[XCOORD])-200.*mm;
  pos_YOKE_R[YCOORD] = 0.0*mm;
  pos_YOKE_R[ZCOORD] = pos_MFIELD[ZCOORD];

  G4double pos_YOKE_L_GapSpace[3];
  pos_YOKE_L_GapSpace[ZCOORD] = pos_YOKE_L[ZCOORD];
  pos_YOKE_L_GapSpace[XCOORD] = pos_YOKE_L[XCOORD]-size_YOKE_LR[XCOORD]-size_YOKE_LR_GapSpace[XCOORD];
  pos_YOKE_L_GapSpace[YCOORD] = 0.0*mm;

  G4double pos_YOKE_R_GapSpace[3];
  pos_YOKE_R_GapSpace[ZCOORD] = pos_YOKE_R[ZCOORD];
  pos_YOKE_R_GapSpace[XCOORD] = pos_YOKE_R[XCOORD]+size_YOKE_LR[XCOORD]+size_YOKE_LR_GapSpace[XCOORD];
  pos_YOKE_R_GapSpace[YCOORD] = 0.0*mm;


  ///up guard
  G4double pos_UGUARD_U[3];
  pos_UGUARD_U[XCOORD] = 0.0*mm;
  pos_UGUARD_U[YCOORD] = size_UGUARD_LR[YCOORD]+size_UGUARD_UD[YCOORD];
  pos_UGUARD_U[ZCOORD] = pos_MFIELD[ZCOORD] - 820.0*mm + size_UGUARD_UD[ZCOORD];

  G4double pos_UGUARD_D[3];
  pos_UGUARD_D[XCOORD] = 0.0*mm;
  pos_UGUARD_D[YCOORD] = -(size_UGUARD_LR[YCOORD]+size_UGUARD_UD[YCOORD]);
  pos_UGUARD_D[ZCOORD] = pos_MFIELD[ZCOORD]- 820.0*mm + size_UGUARD_UD[ZCOORD];

  G4double pos_UGUARD_L[3];
  pos_UGUARD_L[XCOORD] = 0.0*mm + (300.0*mm + size_UGUARD_LR[XCOORD]); // 150 --> gap
  pos_UGUARD_L[YCOORD] = 0.0*mm;
  pos_UGUARD_L[ZCOORD] = pos_MFIELD[ZCOORD]- 820.0*mm + size_UGUARD_LR[ZCOORD];

  G4double pos_UGUARD_R[3];
  pos_UGUARD_R[XCOORD] = 0.0*mm - (300.0*mm + size_UGUARD_LR[XCOORD]);
  pos_UGUARD_R[YCOORD] = 0.0*mm;
  pos_UGUARD_R[ZCOORD] = pos_MFIELD[ZCOORD]- 820.0*mm + size_UGUARD_LR[ZCOORD];

  G4double pos_UGUARD_R1[3];
  pos_UGUARD_R1[ZCOORD] = pos_MFIELD[ZCOORD]- 820.0*mm + size_UGUARD_LR[ZCOORD];
  pos_UGUARD_R1[XCOORD] = pos_UGUARD_U[XCOORD] -(size_UGUARD_UD[XCOORD] - size_UGUARD_R1[XCOORD]);
  pos_UGUARD_R1[YCOORD] = 0.0*mm;

  G4double pos_UGUARD_R2U[3];
  pos_UGUARD_R2U[ZCOORD] = pos_MFIELD[ZCOORD]- 820.0*mm + size_UGUARD_LR[ZCOORD];
  pos_UGUARD_R2U[XCOORD] = pos_UGUARD_U[XCOORD] -(size_UGUARD_UD[XCOORD] - size_UGUARD_R2[XCOORD]);
  pos_UGUARD_R2U[YCOORD] = size_UGUARD_R1[YCOORD]+size_UGUARD_R2[YCOORD];

  G4double pos_UGUARD_R2D[3];
  pos_UGUARD_R2D[ZCOORD] = pos_MFIELD[ZCOORD]- 820.0*mm + size_UGUARD_LR[ZCOORD];
  pos_UGUARD_R2D[XCOORD] = pos_UGUARD_U[XCOORD] -(size_UGUARD_UD[XCOORD] - size_UGUARD_R2[XCOORD]);
  pos_UGUARD_R2D[YCOORD] = -(size_UGUARD_R1[YCOORD]+size_UGUARD_R2[YCOORD]);

  G4double pos_DGUARD_U[3];
  pos_DGUARD_U[XCOORD] = pos_MFIELD[XCOORD];
  pos_DGUARD_U[YCOORD] = size_DGUARD_LR[YCOORD]+size_DGUARD_UD[YCOORD];
  //pos_DGUARD_U[YCOORD] = size_DGUARD_LR[YCOORD]+size_DGUARD_UD[YCOORD] -200.; //Hwang
  pos_DGUARD_U[ZCOORD] = pos_MFIELD[ZCOORD] + 820.0*mm - size_DGUARD_UD[ZCOORD];

  G4double pos_DGUARD_D[3];
  pos_DGUARD_D[XCOORD] = pos_MFIELD[XCOORD];
  pos_DGUARD_D[YCOORD] = -(size_DGUARD_LR[YCOORD]+size_DGUARD_UD[YCOORD]);
  //pos_DGUARD_D[YCOORD] = -(size_DGUARD_LR[YCOORD]+size_DGUARD_UD[YCOORD]) + 200.; //Hwang
  pos_DGUARD_D[ZCOORD] = pos_MFIELD[ZCOORD]+ 820.0*mm - size_DGUARD_UD[ZCOORD];

  G4double pos_DGUARD_L[3];
  pos_DGUARD_L[XCOORD] = pos_MFIELD[XCOORD] + (size_MFIELD[XCOORD] + 50*mm + size_DGUARD_LR[XCOORD]);
  //pos_DGUARD_L[XCOORD] = pos_MFIELD[XCOORD] + (size_MFIELD[XCOORD] + 300*mm + size_DGUARD_LR[XCOORD]); //Hwang
  pos_DGUARD_L[YCOORD] = 0.0*mm;
  pos_DGUARD_L[ZCOORD] = pos_MFIELD[ZCOORD]+ 820.0*mm - size_DGUARD_LR[ZCOORD];

  G4double pos_DGUARD_R[3];
  pos_DGUARD_R[XCOORD] = pos_MFIELD[XCOORD] - (size_MFIELD[XCOORD] + 50*mm + size_DGUARD_LR[XCOORD]);
  //pos_DGUARD_R[XCOORD] = pos_MFIELD[XCOORD] - (size_MFIELD[XCOORD] + 300*mm + size_DGUARD_LR[XCOORD]); //Hwang
  pos_DGUARD_R[YCOORD] = 0.0*mm;
  pos_DGUARD_R[ZCOORD] = pos_MFIELD[ZCOORD]+ 820.0*mm - size_DGUARD_LR[ZCOORD];


  //G4double pos_DGUARD_U[3];
  //pos_DGUARD_U[XCOORD] = pos_MFIELD[XCOORD];
  //pos_DGUARD_U[YCOORD] = size_DGUARD_LR[YCOORD]+size_DGUARD_UD[YCOORD]/2.-200.;
  //pos_DGUARD_U[ZCOORD] = pos_MFIELD[ZCOORD] + 820.0*mm - size_DGUARD_UD[ZCOORD];

  //G4double pos_DGUARD_D[3];
  //pos_DGUARD_D[XCOORD] = pos_MFIELD[XCOORD];
  //pos_DGUARD_D[YCOORD] = -(size_DGUARD_LR[YCOORD]+size_DGUARD_UD[YCOORD]/2.)+200.;
  //pos_DGUARD_D[ZCOORD] = pos_MFIELD[ZCOORD]+ 820.0*mm - size_DGUARD_UD[ZCOORD];

  //G4double pos_DGUARD_L[3];
  //pos_DGUARD_L[XCOORD] = pos_MFIELD[XCOORD] + (size_MFIELD[XCOORD]+300.0*mm + size_DGUARD_LR[XCOORD]);
  //pos_DGUARD_L[YCOORD] = 0.0*mm;
  //pos_DGUARD_L[ZCOORD] = pos_MFIELD[ZCOORD]+ 820.0*mm - size_DGUARD_LR[ZCOORD];

  //G4double pos_DGUARD_R[3];
  //pos_DGUARD_R[XCOORD] = pos_MFIELD[XCOORD] - (size_MFIELD[XCOORD]+300.0*mm + size_DGUARD_LR[XCOORD]);
  //pos_DGUARD_R[YCOORD] = 0.0*mm;
  //pos_DGUARD_R[ZCOORD] = pos_MFIELD[ZCOORD]+ 820.0*mm - size_DGUARD_LR[ZCOORD];

  char nam[128];
  int  m_model=0;


  G4double maxStep;

  //G4RotationMatrix *rotForwardSp = new G4RotationMatrix();
  G4RotationMatrix rotForwardSp;

  //G4String spectrometer_ang = getenv("Spectrometer_ang");
  //G4double fSpectrometerAngle=atof( spectrometer_ang.c_str())*deg;//out side less 100 mm. 10+5*x < 100 mm is pad_in_num
  G4double fSpectrometerAngle = 0.0*deg;
  //  G4double fSpectrometerAngle=10*deg;
  //rotForwardSp->rotateY(-fSpectrometerAngle);
  rotForwardSp.rotateY(-fSpectrometerAngle);


  /////////////
  // Construct KURAMA Magnet
  //////////////coil1U
  G4Box* Coil1_box = new G4Box("Coil1_box",
      size_COIL1[XCOORD],size_COIL1[YCOORD],size_COIL1[ZCOORD]);
  G4LogicalVolume*  Coil1_log = new G4LogicalVolume(Coil1_box, NistMan->FindOrBuildMaterial("G4_Cu"), "Coil1_log",0,0,0);
  Coil1_log->SetVisAttributes(new G4VisAttributes(G4Color::Red()));

  G4VPhysicalVolume* Coil1U_phys = new G4PVPlacement(G4Transform3D(rotForwardSp, G4ThreeVector(pos_COIL1[XCOORD],pos_COIL1[YCOORD],pos_COIL1[ZCOORD]).rotateY(fSpectrometerAngle)),
						     "Coil1U_phys",
						     Coil1_log,
						     physMother,
						     false,
						     0);
  G4VPhysicalVolume* Coil1D_phys = new G4PVPlacement(G4Transform3D(rotForwardSp, G4ThreeVector(pos_COIL1[XCOORD],-pos_COIL1[YCOORD],pos_COIL1[ZCOORD]).rotateY(fSpectrometerAngle)),
      "Coil1D_phys",
      Coil1_log,
      physMother,
      false,
      0);

  //////////////coil4RLUD
  G4Box* Coil4_box = new G4Box("Coil4_box",
      size_COIL4[XCOORD],size_COIL4[YCOORD],size_COIL4[ZCOORD]);
  G4LogicalVolume*  Coil4_log = new G4LogicalVolume(Coil4_box, NistMan->FindOrBuildMaterial("G4_Cu"), "Coil4_log",0,0,0);
  Coil4_log->SetVisAttributes(new G4VisAttributes(G4Color::Red()));


  G4VPhysicalVolume* Coil4UR_phys = new G4PVPlacement(G4Transform3D(rotForwardSp, G4ThreeVector(pos_COIL4L[XCOORD],pos_COIL4L[YCOORD],pos_COIL4L[ZCOORD]).rotateY(fSpectrometerAngle)),
      "Coil4UR_phys",
      Coil4_log,
      physMother,
      false,
      0);
  G4VPhysicalVolume* Coil4UL_phys = new G4PVPlacement(G4Transform3D(rotForwardSp, G4ThreeVector(pos_COIL4R[XCOORD],pos_COIL4R[YCOORD],pos_COIL4R[ZCOORD]).rotateY(fSpectrometerAngle)),
      "Coil4UL_phys",
      Coil4_log,
      physMother,
      false,
      0);
  G4VPhysicalVolume* Coil4DR_phys = new G4PVPlacement(G4Transform3D(rotForwardSp, G4ThreeVector(pos_COIL4L[XCOORD],-pos_COIL4L[YCOORD],pos_COIL4L[ZCOORD]).rotateY(fSpectrometerAngle)),
      "Coil4DR_phys",
      Coil4_log,
      physMother,
      false,
      0);
  G4VPhysicalVolume* Coil4DL_phys = new G4PVPlacement(G4Transform3D(rotForwardSp, G4ThreeVector(pos_COIL4R[XCOORD],-pos_COIL4R[YCOORD],pos_COIL4R[ZCOORD]).rotateY(fSpectrometerAngle)),
      "Coil4DL_phys",
      Coil4_log,
      physMother,
      false,
      0);



  //////////////coil5UD
  G4Box* Coil5_box = new G4Box("Coil5_box",
      size_COIL5[XCOORD],size_COIL5[YCOORD],size_COIL5[ZCOORD]);
  G4LogicalVolume*  Coil5_log = new G4LogicalVolume(Coil5_box, NistMan->FindOrBuildMaterial("G4_Cu"), "Coil5_log",0,0,0);
  Coil5_log->SetVisAttributes(new G4VisAttributes(G4Color::Red()));


  G4VPhysicalVolume* Coil5U_phys = new G4PVPlacement(G4Transform3D(rotForwardSp, G4ThreeVector(pos_COIL5[XCOORD],pos_COIL5[YCOORD],pos_COIL5[ZCOORD]).rotateY(fSpectrometerAngle)),
      "Coil5U_phys",
      Coil5_log,
      physMother,
      false,
      0);
  G4VPhysicalVolume* Coil5D_phys = new G4PVPlacement(G4Transform3D(rotForwardSp, G4ThreeVector(pos_COIL5[XCOORD],-pos_COIL5[YCOORD],pos_COIL5[ZCOORD]).rotateY(fSpectrometerAngle)),
      "Coil5D_phys",
      Coil5_log,
      physMother,
      false,
      0);

  //    padSolid[i] = new G4Tubs("E07 pad", pad_in[i]*mm, pad_out[i]*mm,
  //           250.*mm, 0., angle[i]*deg);

  //////////////coil6RLUD
  G4double size_COIL6[4];
  //0:in
  //1:out
  //2:z
  //3:angle
  size_COIL6[0] = 50.0*mm;
  size_COIL6[1] = size_COIL1[1]*2+size_COIL6[0];
  //  size_COIL6[1] = 330.*mm;
  size_COIL6[2] = 280./2*mm;
  size_COIL6[3] = 90.*deg;

  G4double pos_COIL6LU[3];
  G4double pos_COIL6RU[3];
  G4double pos_COIL6LD[3];
  G4double pos_COIL6RD[3];
  //LU
  pos_COIL6LU[XCOORD] = pos_MFIELD[XCOORD] +size_MFIELD[XCOORD]-size_COIL6[0];
  pos_COIL6LU[YCOORD] = pos_COIL1[YCOORD]  -(size_COIL6[0]+size_COIL1[1]);
  pos_COIL6LU[ZCOORD] = pos_MFIELD[ZCOORD] - size_MFIELD[ZCOORD]-(size_COIL6[ZCOORD]+21.*mm);
  //RU
  pos_COIL6RU[XCOORD] = pos_MFIELD[XCOORD] -size_MFIELD[XCOORD]+size_COIL6[0];
  pos_COIL6RU[YCOORD] = pos_COIL1[YCOORD]  -(size_COIL6[0]+size_COIL1[1]);
  pos_COIL6RU[ZCOORD] = pos_MFIELD[ZCOORD] - size_MFIELD[ZCOORD]-(size_COIL6[ZCOORD]+21.*mm);
  //LD
  pos_COIL6LD[XCOORD] = pos_MFIELD[XCOORD]  +size_MFIELD[XCOORD]-size_COIL6[0];
  pos_COIL6LD[YCOORD] = -pos_COIL1[YCOORD] +(size_COIL6[0]+size_COIL1[1]);
  pos_COIL6LD[ZCOORD] = pos_MFIELD[ZCOORD]  - size_MFIELD[ZCOORD]-(size_COIL6[ZCOORD]+21.*mm);
  //RD
  pos_COIL6RD[XCOORD] = pos_MFIELD[XCOORD]  -size_MFIELD[XCOORD]+size_COIL6[0];
  pos_COIL6RD[YCOORD] = -pos_COIL1[YCOORD] +(size_COIL6[0]+size_COIL1[1]);
  pos_COIL6RD[ZCOORD] = pos_MFIELD[ZCOORD]  - size_MFIELD[ZCOORD]-(size_COIL6[ZCOORD]+21.*mm);


  G4Tubs* Coil6_tub = new G4Tubs("Coil6_tubs",
      size_COIL6[0],size_COIL6[1],size_COIL6[2],0.,size_COIL6[3]);
  G4LogicalVolume*  Coil6_log = new G4LogicalVolume(Coil6_tub, NistMan->FindOrBuildMaterial("G4_Cu"), "Coil6_log",0,0,0);
  Coil6_log->SetVisAttributes(new G4VisAttributes(G4Color::Red()));
  maxStep=0.00001*mm;
  //upGuard_L_log->SetUserLimits(new G4UserLimits(maxStep));

  //G4RotationMatrix *rotcoil6lu = new G4RotationMatrix();
  //G4RotationMatrix *rotcoil6ru = new G4RotationMatrix();
  //G4RotationMatrix *rotcoil6ld = new G4RotationMatrix();
  //G4RotationMatrix *rotcoil6rd = new G4RotationMatrix();
  //rotcoil6lu->rotateY(-fSpectrometerAngle);
  //rotcoil6ru->rotateY(-fSpectrometerAngle);
  //rotcoil6ld->rotateY(-fSpectrometerAngle);
  //rotcoil6rd->rotateY(-fSpectrometerAngle);


  G4RotationMatrix rotcoil6lu;
  G4RotationMatrix rotcoil6ru;
  G4RotationMatrix rotcoil6ld;
  G4RotationMatrix rotcoil6rd;
  rotcoil6lu.rotateY(-fSpectrometerAngle);
  rotcoil6ru.rotateY(-fSpectrometerAngle);
  rotcoil6ld.rotateY(-fSpectrometerAngle);
  rotcoil6rd.rotateY(-fSpectrometerAngle);


  rotcoil6lu.rotateZ(0.*deg);
  G4VPhysicalVolume* Coil6LU_phys
    = new G4PVPlacement(G4Transform3D(rotcoil6lu,
				      G4ThreeVector(pos_COIL6LU[XCOORD],
						    pos_COIL6LU[YCOORD],
						    pos_COIL6LU[ZCOORD]).rotateY(fSpectrometerAngle)),
			"Coil6LU_phys",
			Coil6_log,
			physMother,
			false,
			0);
  //rotcoil6ru.rotateZ(-90.*deg);
  rotcoil6ru.rotateZ(90.*deg);
  G4VPhysicalVolume* Coil6RU_phys
    = new G4PVPlacement(G4Transform3D(rotcoil6ru,
				      G4ThreeVector(pos_COIL6RU[XCOORD],
						    pos_COIL6RU[YCOORD],
						    pos_COIL6RU[ZCOORD]).rotateY(fSpectrometerAngle)),
			"Coil6RU_phys",
			Coil6_log,
			physMother,
			false,
			0);

  //rotcoil6ld.rotateZ(-180.*deg);
  rotcoil6ld.rotateZ(180.*deg);
  G4VPhysicalVolume* Coil6RD_phys
    = new G4PVPlacement(G4Transform3D(rotcoil6ld,
				      G4ThreeVector(pos_COIL6RD[XCOORD],
						    pos_COIL6RD[YCOORD],
						    pos_COIL6RD[ZCOORD]).rotateY(fSpectrometerAngle)),
			"Coil6RD_phys",
			Coil6_log,
			physMother,
			false,
			0);
  //rotcoil6rd.rotateZ(-270.*deg);
  rotcoil6rd.rotateZ(270.*deg);
  G4VPhysicalVolume* Coil6LD_phys
    = new G4PVPlacement(G4Transform3D(rotcoil6rd,
				      G4ThreeVector(pos_COIL6LD[XCOORD],
						    pos_COIL6LD[YCOORD],
						    pos_COIL6LD[ZCOORD]).rotateY(fSpectrometerAngle)),
			"Coil6LD_phys",
			Coil6_log,
			physMother,
			false,
			0);
  //  rotForwardSp->rotateZ(90.*deg);
  //  rotForwardSp->rotateZ(90.*deg);



  //////////////coil8RLUD
  G4double size_COIL8[4];
  //0:in
  //1:out
  //2:z
  //3:angle
  size_COIL8[0] = 50.0*mm;
  size_COIL8[1] = size_COIL5[1]*2.+size_COIL8[0];
  size_COIL8[2] = 280./2*mm;
  size_COIL8[3] = 90.*deg;

  G4double pos_COIL8LU[3];
  G4double pos_COIL8RU[3];
  G4double pos_COIL8LD[3];
  G4double pos_COIL8RD[3];
  //LU
  pos_COIL8LU[XCOORD] = pos_MFIELD[XCOORD] +size_MFIELD[XCOORD]-size_COIL8[0];
  pos_COIL8LU[YCOORD] = pos_COIL5[YCOORD]  -(size_COIL8[0]+size_COIL5[1]);
  pos_COIL8LU[ZCOORD] = pos_MFIELD[ZCOORD] + size_MFIELD[ZCOORD]+(size_COIL8[ZCOORD]+21.*mm);
  //RU
  pos_COIL8RU[XCOORD] = pos_MFIELD[XCOORD] -size_MFIELD[XCOORD]+size_COIL8[0];
  pos_COIL8RU[YCOORD] = pos_COIL5[YCOORD]  -(size_COIL8[0]+size_COIL5[1]);
  pos_COIL8RU[ZCOORD] = pos_MFIELD[ZCOORD] + size_MFIELD[ZCOORD]+(size_COIL8[ZCOORD]+21.*mm);
  //LD
  pos_COIL8LD[XCOORD] = pos_MFIELD[XCOORD]  +size_MFIELD[XCOORD]-size_COIL8[0];
  pos_COIL8LD[YCOORD] = -pos_COIL5[YCOORD] +(size_COIL8[0]+size_COIL5[1]);
  pos_COIL8LD[ZCOORD] = pos_MFIELD[ZCOORD]  + size_MFIELD[ZCOORD]+(size_COIL8[ZCOORD]+21.*mm);
  //RD
  pos_COIL8RD[XCOORD] = pos_MFIELD[XCOORD]  -size_MFIELD[XCOORD]+size_COIL8[0];
  pos_COIL8RD[YCOORD] = -pos_COIL5[YCOORD] +(size_COIL8[0]+size_COIL5[1]);
  pos_COIL8RD[ZCOORD] = pos_MFIELD[ZCOORD]  + size_MFIELD[ZCOORD]+(size_COIL8[ZCOORD]+21.*mm);


  G4Tubs* Coil8_tub = new G4Tubs("Coil8_tubs",
      size_COIL8[0],size_COIL8[1],size_COIL8[2],0.,size_COIL8[3]);
  G4LogicalVolume*  Coil8_log = new G4LogicalVolume(Coil8_tub, NistMan->FindOrBuildMaterial("G4_Cu"), "Coil8_log",0,0,0);
  Coil8_log->SetVisAttributes(new G4VisAttributes(G4Color::Red()));
  maxStep=0.00001*mm;
  //upGuard_L_log->SetUserLimits(new G4UserLimits(maxStep));

  G4RotationMatrix rotcoil8lu;
  G4RotationMatrix rotcoil8ru;
  G4RotationMatrix rotcoil8ld;
  G4RotationMatrix rotcoil8rd;
  rotcoil8lu.rotateY(-fSpectrometerAngle);
  rotcoil8ru.rotateY(-fSpectrometerAngle);
  rotcoil8ld.rotateY(-fSpectrometerAngle);
  rotcoil8rd.rotateY(-fSpectrometerAngle);

  rotcoil8lu.rotateZ(0.*deg);
  G4VPhysicalVolume* Coil8LU_phys
    = new G4PVPlacement(G4Transform3D(rotcoil8lu,
				      G4ThreeVector(pos_COIL8LU[XCOORD],
						    pos_COIL8LU[YCOORD],
						    pos_COIL8LU[ZCOORD]).rotateY(fSpectrometerAngle)),
			"Coil8LU_phys",
			Coil8_log,
			physMother,
			false,
			0);
  //rotcoil8ru.rotateZ(-90.*deg);
  rotcoil8ru.rotateZ(90.*deg);
  G4VPhysicalVolume* Coil8RU_phys
    = new G4PVPlacement(G4Transform3D(rotcoil8ru,
				      G4ThreeVector(pos_COIL8RU[XCOORD],
						    pos_COIL8RU[YCOORD],
						    pos_COIL8RU[ZCOORD]).rotateY(fSpectrometerAngle)),
			"Coil8RU_phys",
			Coil8_log,
			physMother,
			false,
			0);

  //rotcoil8ld.rotateZ(-180.*deg);
  //rotcoil8rd.rotateZ(270.*deg);
  rotcoil8rd.rotateZ(180.*deg);
  G4VPhysicalVolume* Coil8RD_phys
    = new G4PVPlacement(G4Transform3D(rotcoil8rd,
				      G4ThreeVector(pos_COIL8RD[XCOORD],
						    pos_COIL8RD[YCOORD],
						    pos_COIL8RD[ZCOORD]).rotateY(fSpectrometerAngle)),
			"Coil8RD_phys",
			Coil8_log,
			physMother,
			false,
			0);
  //rotcoil8rd.rotateZ(-270.*deg);
  //rotcoil8ld.rotateZ(180.*deg);
  rotcoil8ld.rotateZ(270.*deg);
  G4VPhysicalVolume* Coil8LD_phys
    = new G4PVPlacement(G4Transform3D(rotcoil8ld,
				      G4ThreeVector(pos_COIL8LD[XCOORD],
						    pos_COIL8LD[YCOORD],
						    pos_COIL8LD[ZCOORD]).rotateY(fSpectrometerAngle)),
			"Coil8LD_phys",
			Coil8_log,
			
			physMother,
			false,
			0);
  //  rotForwardSp->rotateZ(90.*deg);
  //  rotForwardSp->rotateZ(90.*deg);
  


  //////////////coil7RLUD
  G4double size_COIL7[4];
  //0:in
  //1:out
  //2:z
  //3:angle
  size_COIL7[0] = 50.0*mm;
  size_COIL7[1] = size_COIL4[1]*2.+size_COIL7[0];
  size_COIL7[2] = size_COIL4[0];
  size_COIL7[3] = 90.*deg;

  G4double pos_COIL7ULU[3];
  G4double pos_COIL7URU[3];
  G4double pos_COIL7ULD[3];
  G4double pos_COIL7URD[3];
  G4double pos_COIL7DLU[3];
  G4double pos_COIL7DRU[3];
  G4double pos_COIL7DLD[3];
  G4double pos_COIL7DRD[3];
  //ULU
  //  pos_COIL7ULU[XCOORD] = pos_MFIELD[XCOORD] +size_MFIELD[XCOORD]+size_COIL7[0];
  pos_COIL7ULU[XCOORD] = pos_COIL4L[XCOORD];
  pos_COIL7ULU[YCOORD] = pos_COIL4L[YCOORD] +(size_COIL7[0]+size_COIL4[1]);
  pos_COIL7ULU[ZCOORD] = pos_MFIELD[ZCOORD] - size_COIL4[ZCOORD];
  //URU
  pos_COIL7URU[XCOORD] = pos_COIL4R[XCOORD];
  pos_COIL7URU[YCOORD] = pos_COIL4R[YCOORD] +(size_COIL7[0]+size_COIL4[1]);
  pos_COIL7URU[ZCOORD] = pos_MFIELD[ZCOORD] - size_COIL4[ZCOORD];
  //ULD
  pos_COIL7ULD[XCOORD] = pos_COIL4L[XCOORD];
  pos_COIL7ULD[YCOORD] = -pos_COIL4L[YCOORD] -(size_COIL7[0]+size_COIL4[1]);
  pos_COIL7ULD[ZCOORD] = pos_MFIELD[ZCOORD]  - size_COIL4[ZCOORD];
  //URD
  pos_COIL7URD[XCOORD] = pos_COIL4R[XCOORD];
  pos_COIL7URD[YCOORD] = -pos_COIL4R[YCOORD] -(size_COIL7[0]+size_COIL4[1]);
  pos_COIL7URD[ZCOORD] = pos_MFIELD[ZCOORD] - size_COIL4[ZCOORD];


  //DLU
  //  pos_COIL7ULU[XCOORD] = pos_MFIELD[XCOORD] +size_MFIELD[XCOORD]+size_COIL7[0];
  pos_COIL7DLU[XCOORD] = pos_COIL4L[XCOORD];
  pos_COIL7DLU[YCOORD] = pos_COIL4L[YCOORD] +(size_COIL7[0]+size_COIL4[1]);
  pos_COIL7DLU[ZCOORD] = pos_MFIELD[ZCOORD] + size_COIL4[ZCOORD];
  //DRU
  pos_COIL7DRU[XCOORD] = pos_COIL4R[XCOORD];
  pos_COIL7DRU[YCOORD] = pos_COIL4R[YCOORD] +(size_COIL7[0]+size_COIL4[1]);
  pos_COIL7DRU[ZCOORD] = pos_MFIELD[ZCOORD] + size_COIL4[ZCOORD];
  //DLD
  pos_COIL7DLD[XCOORD] = pos_COIL4L[XCOORD];
  pos_COIL7DLD[YCOORD] = -pos_COIL4L[YCOORD] -(size_COIL7[0]+size_COIL4[1]);
  pos_COIL7DLD[ZCOORD] = pos_MFIELD[ZCOORD]  + size_COIL4[ZCOORD];
  //DRD
  pos_COIL7DRD[XCOORD] = pos_COIL4R[XCOORD];
  pos_COIL7DRD[YCOORD] = -pos_COIL4R[YCOORD] -(size_COIL7[0]+size_COIL4[1]);
  pos_COIL7DRD[ZCOORD] = pos_MFIELD[ZCOORD] + size_COIL4[ZCOORD];


  G4Tubs* Coil7_tub = new G4Tubs("Coil7_tubs",
      size_COIL7[0],size_COIL7[1],size_COIL7[2],0.,size_COIL7[3]);
  G4LogicalVolume*  Coil7_log = new G4LogicalVolume(Coil7_tub, NistMan->FindOrBuildMaterial("G4_Cu"), "Coil7_log",0,0,0);
  Coil7_log->SetVisAttributes(new G4VisAttributes(G4Color::Red()));
  maxStep=0.00001*mm;
  //upGuard_L_log->SetUserLimits(new G4UserLimits(maxStep));

  G4RotationMatrix rotcoil7ulu;
  G4RotationMatrix rotcoil7uru;
  G4RotationMatrix rotcoil7uld;
  G4RotationMatrix rotcoil7urd;
  G4RotationMatrix rotcoil7dlu;
  G4RotationMatrix rotcoil7dru;
  G4RotationMatrix rotcoil7dld;
  G4RotationMatrix rotcoil7drd;
  rotcoil7ulu.rotateY(-fSpectrometerAngle);
  rotcoil7uru.rotateY(-fSpectrometerAngle);
  rotcoil7uld.rotateY(-fSpectrometerAngle);
  rotcoil7urd.rotateY(-fSpectrometerAngle);
  rotcoil7dlu.rotateY(-fSpectrometerAngle);
  rotcoil7dru.rotateY(-fSpectrometerAngle);
  rotcoil7dld.rotateY(-fSpectrometerAngle);
  rotcoil7drd.rotateY(-fSpectrometerAngle);

  
  rotcoil7ulu.rotateZ(0.*deg);
  rotcoil7ulu.rotateX(180.*deg);
  rotcoil7ulu.rotateY(90.*deg);
  G4VPhysicalVolume* Coil7ULU_phys
    = new G4PVPlacement(G4Transform3D(rotcoil7ulu,
				      G4ThreeVector(pos_COIL7ULU[XCOORD],
						    pos_COIL7ULU[YCOORD],
						    pos_COIL7ULU[ZCOORD]).rotateY(fSpectrometerAngle)),
			"Coil7ULU_phys",
			Coil7_log,
			physMother,
			false,
			0);


  rotcoil7uru.rotateZ(0.*deg);
  rotcoil7uru.rotateX(180.*deg);
  rotcoil7uru.rotateY(90.*deg);
  G4VPhysicalVolume* Coil7URU_phys
    = new G4PVPlacement(G4Transform3D(rotcoil7uru,
				      G4ThreeVector(pos_COIL7URU[XCOORD],
						    pos_COIL7URU[YCOORD],
						    pos_COIL7URU[ZCOORD]).rotateY(fSpectrometerAngle)),
			"Coil7URU_phys",
			Coil7_log,
			physMother,
			false,
			0);

  rotcoil7urd.rotateZ(0.*deg);
  //rotcoil7urd.rotateX(90.*deg);
  rotcoil7urd.rotateX(0.*deg);
  rotcoil7urd.rotateY(90.*deg);
  G4VPhysicalVolume* Coil7URD_phys
    = new G4PVPlacement(G4Transform3D(rotcoil7urd,
				      G4ThreeVector(pos_COIL7URD[XCOORD],
						    pos_COIL7URD[YCOORD],
						    pos_COIL7URD[ZCOORD]).rotateY(fSpectrometerAngle)),
			"Coil7URD_phys",
			Coil7_log,
			physMother,
			false,
			0);

  rotcoil7uld.rotateZ(0.*deg);
  //rotcoil7uld.rotateX(90.*deg);
  rotcoil7uld.rotateX(0.*deg);
  rotcoil7uld.rotateY(90.*deg);
  G4VPhysicalVolume* Coil7ULD_phys
    = new G4PVPlacement(G4Transform3D(rotcoil7uld,
				      G4ThreeVector(pos_COIL7ULD[XCOORD],
						    pos_COIL7ULD[YCOORD],
						    pos_COIL7ULD[ZCOORD]).rotateY(fSpectrometerAngle)),
			"Coil7ULD_phys",
			Coil7_log,
			physMother,
			false,
			0);


  ///down
  rotcoil7dlu.rotateX(-90.*deg);
  rotcoil7dlu.rotateY(180.*deg);
  rotcoil7dlu.rotateZ(90.*deg);
  G4VPhysicalVolume* Coil7DLU_phys
    = new G4PVPlacement(G4Transform3D(rotcoil7dlu,
				      G4ThreeVector(pos_COIL7DLU[XCOORD],
						    pos_COIL7DLU[YCOORD],
						    pos_COIL7DLU[ZCOORD]).rotateY(fSpectrometerAngle)),
			"Coil7DLU_phys",
			Coil7_log,
			physMother,
			false,
			0);



  rotcoil7dru.rotateX(-90.*deg);
  rotcoil7dru.rotateY(180.*deg);
  rotcoil7dru.rotateZ(90.*deg);
  G4VPhysicalVolume* Coil7DRU_phys
    = new G4PVPlacement(G4Transform3D(rotcoil7dru,
				      G4ThreeVector(pos_COIL7DRU[XCOORD],
						    pos_COIL7DRU[YCOORD],
						    pos_COIL7DRU[ZCOORD]).rotateY(fSpectrometerAngle)),
			"Coil7DRU_phys",
			Coil7_log,
			physMother,
			false,
			0);

  rotcoil7drd.rotateZ(0.*deg);
  rotcoil7drd.rotateX(0.*deg);
  rotcoil7drd.rotateY(-90.*deg);
  G4VPhysicalVolume* Coil7DRD_phys
    = new G4PVPlacement(G4Transform3D(rotcoil7drd,
				      G4ThreeVector(pos_COIL7DRD[XCOORD],
						    pos_COIL7DRD[YCOORD],
						    pos_COIL7DRD[ZCOORD]).rotateY(fSpectrometerAngle)),
			"Coil7DRD_phys",
			Coil7_log,
			physMother,
			false,
			0);

  rotcoil7dld.rotateZ(0.*deg);
  rotcoil7dld.rotateX(0.*deg);
  rotcoil7dld.rotateY(-90.*deg);
  G4VPhysicalVolume* Coil7DLD_phys
    = new G4PVPlacement(G4Transform3D(rotcoil7dld,
				      G4ThreeVector(pos_COIL7DLD[XCOORD],
						    pos_COIL7DLD[YCOORD],
						    pos_COIL7DLD[ZCOORD]).rotateY(fSpectrometerAngle)),
			"Coil7DLD_phys",
			Coil7_log,
			physMother,
			false,
			0);
  
  ///coil2
  //////////////coil2
  G4Box* Coil2_box = new G4Box("Coil2_box",
      size_COIL2[XCOORD],size_COIL2[YCOORD],size_COIL2[ZCOORD]);
  G4LogicalVolume*  Coil2_log = new G4LogicalVolume(Coil2_box, NistMan->FindOrBuildMaterial("G4_Cu"), "Coil2_log",0,0,0);
  Coil2_log->SetVisAttributes(new G4VisAttributes(G4Color::Red()));
  maxStep=0.00001*mm;
  //upGuard_L_log->SetUserLimits(new G4UserLimits(maxStep));

  G4VPhysicalVolume* Coil2UL_phys = new G4PVPlacement(G4Transform3D(rotForwardSp,
								    G4ThreeVector(pos_COIL2L[XCOORD],
										  pos_COIL2L[YCOORD],
										  pos_COIL2L[ZCOORD]).rotateY(fSpectrometerAngle)),
						      "Coil2UL_phys",
						      Coil2_log,
						      physMother,
						      false,
						      0);
  G4VPhysicalVolume* Coil2UR_phys = new G4PVPlacement(G4Transform3D(rotForwardSp,
								    G4ThreeVector(pos_COIL2R[XCOORD],
										  pos_COIL2R[YCOORD],
										  pos_COIL2R[ZCOORD]).rotateY(fSpectrometerAngle)),
						      "Coil2UR_phys",
						      Coil2_log,
						      physMother,
						      false,
						      0);


  G4VPhysicalVolume* Coil2DL_phys = new G4PVPlacement(G4Transform3D(rotForwardSp,
								    G4ThreeVector(pos_COIL2L[XCOORD],
										  -pos_COIL2L[YCOORD],
										  pos_COIL2L[ZCOORD]).rotateY(fSpectrometerAngle)),
						      "Coil2DL_phys",
						      Coil2_log,
						      physMother,
						      false,
						      0);
  G4VPhysicalVolume* Coil2DR_phys = new G4PVPlacement(G4Transform3D(rotForwardSp,
								    G4ThreeVector(pos_COIL2R[XCOORD],
										  -pos_COIL2R[YCOORD],
										  pos_COIL2R[ZCOORD]).rotateY(fSpectrometerAngle)),
						      "Coil2DR_phys",
						      Coil2_log,
						      physMother,
						      false,
						      0);


  ///coil3
  //////////////coil3
  G4Box* Coil3_box = new G4Box("Coil3_box",
      size_COIL3[XCOORD],size_COIL3[YCOORD],size_COIL3[ZCOORD]);
  G4LogicalVolume*  Coil3_log = new G4LogicalVolume(Coil3_box, NistMan->FindOrBuildMaterial("G4_Cu"), "Coil3_log",0,0,0);
  Coil3_log->SetVisAttributes(new G4VisAttributes(G4Color::Red()));
  maxStep=0.00001*mm;
  //upGuard_L_log->SetUserLimits(new G4UserLimits(maxStep));

  G4VPhysicalVolume* Coil3UL_phys = new G4PVPlacement(G4Transform3D(rotForwardSp,
								    G4ThreeVector(pos_COIL3L[XCOORD],
										  pos_COIL3L[YCOORD],
										  pos_COIL3L[ZCOORD]).rotateY(fSpectrometerAngle)),
						      "Coil3UL_phys",
						      Coil3_log,
						      physMother,
						      false,
						      0);
  G4VPhysicalVolume* Coil3UR_phys = new G4PVPlacement(G4Transform3D(rotForwardSp,
								    G4ThreeVector(pos_COIL3R[XCOORD],
										  pos_COIL3R[YCOORD],
										  pos_COIL3R[ZCOORD]).rotateY(fSpectrometerAngle)),
						      "Coil3UR_phys",
						      Coil3_log,
						      physMother,
						      false,
						      0);
  
  
  G4VPhysicalVolume* Coil3DL_phys = new G4PVPlacement(G4Transform3D(rotForwardSp,
								    G4ThreeVector(pos_COIL3L[XCOORD],
										  -pos_COIL3L[YCOORD],
										  pos_COIL3L[ZCOORD]).rotateY(fSpectrometerAngle)),
						      "Coil3DL_phys",
						      Coil3_log,
						      physMother,
						      false,
						      0);

  G4VPhysicalVolume* Coil3DR_phys = new G4PVPlacement(G4Transform3D(rotForwardSp,
								    G4ThreeVector(pos_COIL3R[XCOORD],
										  -pos_COIL3R[YCOORD],
										  pos_COIL3R[ZCOORD]).rotateY(fSpectrometerAngle)),
						      "Coil3DR_phys",
						      Coil3_log,
						      physMother,
						      false,
						      0);




  //-------------------- Upstream End Guard
#if 1      
    G4Box* upGuard_UD_box = new G4Box("upGuard_UD_box",
        size_UGUARD_UD[XCOORD],size_UGUARD_UD[YCOORD],size_UGUARD_UD[ZCOORD]);
    G4LogicalVolume*  upGuard_U_log = new G4LogicalVolume(upGuard_UD_box, NistMan->FindOrBuildMaterial("G4_Fe"), "upGuard_U_log",0,0,0);
    upGuard_U_log->SetVisAttributes(new G4VisAttributes(G4Color::Blue()));
    maxStep=0.00001*mm;
    //upGuard_U_log->SetUserLimits(new G4UserLimits(maxStep));
    
    G4VPhysicalVolume*  upGuard_U_phys = new G4PVPlacement(G4Transform3D(rotForwardSp,
									 G4ThreeVector(pos_UGUARD_U[XCOORD],
										       pos_UGUARD_U[YCOORD],
										       pos_UGUARD_U[ZCOORD]).rotateY(fSpectrometerAngle)),
							   "upGuard_U_phys",
							   upGuard_U_log,
							   physMother,
							   false,
							   0);
    
    G4LogicalVolume*  upGuard_D_log = new G4LogicalVolume(upGuard_UD_box, NistMan->FindOrBuildMaterial("G4_Fe"), "upGuard_D_log",0,0,0);
    upGuard_D_log->SetVisAttributes(new G4VisAttributes(G4Color::Blue()));
    maxStep=0.00001*mm;
    //upGuard_D_log->SetUserLimits(new G4UserLimits(maxStep));
    
    G4VPhysicalVolume*  upGuard_D_phys = new G4PVPlacement(G4Transform3D(rotForwardSp,
									 G4ThreeVector(pos_UGUARD_D[XCOORD],
										       pos_UGUARD_D[YCOORD],
										       pos_UGUARD_D[ZCOORD]).rotateY(fSpectrometerAngle)),
							   "upGuard_D_phys",
							   upGuard_D_log,
							   physMother,
							   false,
							   0);
    
    G4Box* upGuard_LR_box = new G4Box("upGuard_LR_box",
        size_UGUARD_LR[XCOORD],size_UGUARD_LR[YCOORD],size_UGUARD_LR[ZCOORD]);
    G4LogicalVolume*  upGuard_L_log = new G4LogicalVolume(upGuard_LR_box, NistMan->FindOrBuildMaterial("G4_Fe"), "upGuard_L_log",0,0,0);
    upGuard_L_log->SetVisAttributes(new G4VisAttributes(G4Color::Blue()));
    maxStep=0.00001*mm;
    //upGuard_L_log->SetUserLimits(new G4UserLimits(maxStep));
    G4VPhysicalVolume* upGuard_L_phys = new G4PVPlacement(G4Transform3D(rotForwardSp,
									G4ThreeVector(pos_UGUARD_L[XCOORD],
										      pos_UGUARD_L[YCOORD],
										      pos_UGUARD_L[ZCOORD]).rotateY(fSpectrometerAngle)),
							  "upGuard_L_phys",
							  upGuard_L_log,
							  physMother,
							  false,
							  0);
    /*
    G4LogicalVolume*  upGuard_R_log = new G4LogicalVolume(upGuard_LR_box, NistMan->FindOrBuildMaterial("G4_Fe"), "upGuard_R_log",0,0,0);
    upGuard_R_log->SetVisAttributes(new G4VisAttributes(G4Color::Blue()));
    maxStep=0.00001*mm;
    //upGuard_R_log->SetUserLimits(new G4UserLimits(maxStep));
    G4VPhysicalVolume* upGuard_R_phys = new G4PVPlacement(G4Transform3D(rotForwardSp,
									G4ThreeVector(pos_UGUARD_R[XCOORD],
										      pos_UGUARD_R[YCOORD],
										      pos_UGUARD_R[ZCOORD]).rotateY(fSpectrometerAngle)),
							  "upGuard_R_phys",							  
							  upGuard_R_log,
							  physMother,
							  false,
							  0);
    */
    G4Box* upGuard_R1_box =
      new G4Box("upGuard_R1_box",
		size_UGUARD_R1[XCOORD],size_UGUARD_R1[YCOORD],size_UGUARD_R1[ZCOORD]);
    
    G4Box* upGuard_R2_box =
      new G4Box("upGuard_R2_box",
		size_UGUARD_R2[XCOORD],size_UGUARD_R2[YCOORD],size_UGUARD_R2[ZCOORD]);
    
    
    G4LogicalVolume* upGuard_R1_log = new G4LogicalVolume(upGuard_R1_box,
							  NistMan->FindOrBuildMaterial("G4_Fe"),
							  "upGuard_R1_log",0,0,0);
    upGuard_R1_log->SetVisAttributes(new G4VisAttributes(G4Color::Blue()));
                           
    G4VPhysicalVolume* upGuard_R1_phys =
      new G4PVPlacement(G4Transform3D(rotForwardSp,
				      G4ThreeVector(pos_UGUARD_R1[XCOORD],
						    pos_UGUARD_R1[YCOORD],
						    pos_UGUARD_R1[ZCOORD]).rotateY(fSpectrometerAngle)),
			"upGuard_R1_phys",
			upGuard_R1_log,
			physMother,
			false,
			0);
    
    G4LogicalVolume* upGuard_R2U_log = new G4LogicalVolume(upGuard_R2_box,
							  NistMan->FindOrBuildMaterial("G4_Fe"),
							   "upGuard_R2U_log",0,0,0);
    upGuard_R2U_log->SetVisAttributes(new G4VisAttributes(G4Color::Blue()));
    G4VPhysicalVolume* upGuard_R2U_phys =
      new G4PVPlacement(0,
			G4ThreeVector(pos_UGUARD_R2U[XCOORD],pos_UGUARD_R2U[YCOORD],pos_UGUARD_R2U[ZCOORD]),
			"upGuard_R2U_phys",
			upGuard_R2U_log,
			physMother,
			false,
			0);
    
    G4LogicalVolume* upGuard_R2D_log = new G4LogicalVolume(upGuard_R2_box,
							   NistMan->FindOrBuildMaterial("G4_Fe"),
							   "upGuard_R2D_log",0,0,0);
    upGuard_R2D_log->SetVisAttributes(new G4VisAttributes(G4Color::Blue()));
    maxStep=0.00001*mm;
    G4VPhysicalVolume* upGuard_R2D_phys =
      new G4PVPlacement(0,
			G4ThreeVector(pos_UGUARD_R2D[XCOORD],pos_UGUARD_R2D[YCOORD],pos_UGUARD_R2D[ZCOORD]),
			"upGuard_R2D_phys",
			upGuard_R2D_log,
			physMother,
			false,
			0);
    


#endif
  
   //temp

  //-------------------- Yoke
  
  G4Box* Yoke_UD_box = new G4Box("Yoke_UD_box",
      size_YOKE_UD[XCOORD],size_YOKE_UD[YCOORD],size_YOKE_UD[ZCOORD]);
  G4LogicalVolume*  Yoke_U_log = new G4LogicalVolume(Yoke_UD_box, NistMan->FindOrBuildMaterial("G4_Fe"), "Yoke_U_log",0,0,0);

  Yoke_U_log->SetVisAttributes(new G4VisAttributes(G4Color::Blue()));
  maxStep=0.00001*mm;
  //Yoke_U_log->SetUserLimits(new G4UserLimits(maxStep));
  G4VPhysicalVolume* Yoke_U_phys = new G4PVPlacement(G4Transform3D(rotForwardSp,
								   G4ThreeVector(pos_YOKE_U[XCOORD],
										 pos_YOKE_U[YCOORD],
										 pos_YOKE_U[ZCOORD]).rotateY(fSpectrometerAngle)),
						     "Yoke_U_phys",
						     Yoke_U_log,
						     physMother,
						     false,
						     0);

  G4LogicalVolume*  Yoke_D_log = new G4LogicalVolume(Yoke_UD_box, NistMan->FindOrBuildMaterial("G4_Fe"), "Yoke_D_log",0,0,0);
  Yoke_D_log->SetVisAttributes(new G4VisAttributes(G4Color::Blue()));
  maxStep=0.00001*mm;
  //Yoke_D_log->SetUserLimits(new G4UserLimits(maxStep));
  G4VPhysicalVolume* Yoke_D_phys = new G4PVPlacement(G4Transform3D(rotForwardSp,
								   G4ThreeVector(pos_YOKE_D[XCOORD],
										 pos_YOKE_D[YCOORD],
										 pos_YOKE_D[ZCOORD]).rotateY(fSpectrometerAngle)),
						     "Yoke_D_phys",
						     Yoke_D_log,
						     physMother,
						     false,
						     0);

  G4Box* Yoke_LR_box = new G4Box("Yoke_LR_box",
      size_YOKE_LR[XCOORD],size_YOKE_LR[YCOORD],size_YOKE_LR[ZCOORD]);
  G4LogicalVolume*  Yoke_L_log = new G4LogicalVolume(Yoke_LR_box, NistMan->FindOrBuildMaterial("G4_Fe"), "Yoke_L_log",0,0,0);
  Yoke_L_log->SetVisAttributes(new G4VisAttributes(G4Color::Blue()));
  maxStep=0.00001*mm;
  //Yoke_L_log->SetUserLimits(new G4UserLimits(maxStep));
  G4VPhysicalVolume* Yoke_L_phys = new G4PVPlacement(G4Transform3D(rotForwardSp,
								   G4ThreeVector(pos_YOKE_L[XCOORD],
										 pos_YOKE_L[YCOORD],
										 pos_YOKE_L[ZCOORD]).rotateY(fSpectrometerAngle)),
						     "Yoke_L_phys",
						     Yoke_L_log,
						     physMother,
						     false,
						     0);

  G4LogicalVolume*  Yoke_R_log = new G4LogicalVolume(Yoke_LR_box, NistMan->FindOrBuildMaterial("G4_Fe"), "Yoke_R_log",0,0,0);
  Yoke_R_log->SetVisAttributes(new G4VisAttributes(G4Color::Blue()));
  maxStep=0.00001*mm;
  //Yoke_R_log->SetUserLimits(new G4UserLimits(maxStep));
  G4VPhysicalVolume* Yoke_R_phys = new G4PVPlacement(G4Transform3D(rotForwardSp,
								   G4ThreeVector(pos_YOKE_R[XCOORD],
										 pos_YOKE_R[YCOORD],
										 pos_YOKE_R[ZCOORD]).rotateY(fSpectrometerAngle)),
						     "Yoke_R_phys",
						     Yoke_R_log,
						     physMother,
						     false,
						     0);
  
  //E42 Yoke spacer addition
  #if 1
  G4Box* Yoke_LR_GapSpace_box = new G4Box("Yoke_LR_GapSpace_box",
                                          size_YOKE_LR_GapSpace[XCOORD],
                                          size_YOKE_LR_GapSpace[YCOORD],
                                          size_YOKE_LR_GapSpace[ZCOORD]);
  G4LogicalVolume* Yoke_L_GapSpace_log = new G4LogicalVolume(Yoke_LR_GapSpace_box,
                                                             NistMan->FindOrBuildMaterial("G4_Fe"),
                                                             "Yoke_L_GapSpace_log",0,0,0);

  Yoke_L_GapSpace_log->SetVisAttributes(new G4VisAttributes(G4Color::Blue()));
  G4VPhysicalVolume* Yoke_L_GapSpace_phys =
    new G4PVPlacement(0,
                      G4ThreeVector(pos_YOKE_L_GapSpace[XCOORD],pos_YOKE_L_GapSpace[YCOORD],pos_YOKE_L_GapSpace[ZCOORD]),
                      "Yoke_L_GapSpace_phys",
                      Yoke_L_GapSpace_log,
                      physMother,
                      false,
                      0);
  
  G4LogicalVolume* Yoke_R_GapSpace_log = new G4LogicalVolume(Yoke_LR_GapSpace_box,
                                                             NistMan->FindOrBuildMaterial("G4_Fe"),
                                                             "Yoke_R_GapSpace_log",0,0,0);

  Yoke_R_GapSpace_log->SetVisAttributes(new G4VisAttributes(G4Color::Blue()));
  G4VPhysicalVolume* Yoke_R_GapSpace_phys =
    new G4PVPlacement(0,
                      G4ThreeVector(pos_YOKE_R_GapSpace[XCOORD],pos_YOKE_R_GapSpace[YCOORD],pos_YOKE_R_GapSpace[ZCOORD]),
                      "Yoke_R_GapSpace_phys",
                      Yoke_R_GapSpace_log,
                      physMother,
                      false,
                      0);
  #endif
  //-------------------- Downstream End Guard
  #if 1
  G4Box* downGuard_UD_box = new G4Box("downGuard_UD_box",
      size_DGUARD_UD[XCOORD],size_DGUARD_UD[YCOORD],size_DGUARD_UD[ZCOORD]);
  G4LogicalVolume*  downGuard_U_log = new G4LogicalVolume(downGuard_UD_box, NistMan->FindOrBuildMaterial("G4_Fe"), "downGuard_U_log",0,0,0);
  downGuard_U_log->SetVisAttributes(new G4VisAttributes(G4Color::Blue()));
  maxStep=0.00001*mm;
  //downGuard_U_log->SetUserLimits(new G4UserLimits(maxStep));
  G4VPhysicalVolume* downGuard_U_phys = new G4PVPlacement(G4Transform3D(rotForwardSp,
									G4ThreeVector(pos_DGUARD_U[XCOORD],
										      pos_DGUARD_U[YCOORD],
										      pos_DGUARD_U[ZCOORD]).rotateY(fSpectrometerAngle)),
							  "downGuard_U_phys",
							  downGuard_U_log,
							  physMother,
							  false,
							  0);

  G4LogicalVolume*  downGuard_D_log = new G4LogicalVolume(downGuard_UD_box, NistMan->FindOrBuildMaterial("G4_Fe"), "downGuard_D_log",0,0,0);
  downGuard_D_log->SetVisAttributes(new G4VisAttributes(G4Color::Blue()));
  maxStep=0.00001*mm;
  //downGuard_D_log->SetUserLimits(new G4UserLimits(maxStep));
  G4VPhysicalVolume* downGuard_D_phys = new G4PVPlacement(G4Transform3D(rotForwardSp,
									G4ThreeVector(pos_DGUARD_D[XCOORD],
										      pos_DGUARD_D[YCOORD],
										      pos_DGUARD_D[ZCOORD]).rotateY(fSpectrometerAngle)),
							  "downGuard_D_phys",
							  downGuard_D_log,
							  physMother,
							  false,
							  0);

  G4Box* downGuard_LR_box = new G4Box("downGuard_LR_box",
      size_DGUARD_LR[XCOORD],size_DGUARD_LR[YCOORD],size_DGUARD_LR[ZCOORD]);
  G4LogicalVolume*  downGuard_L_log = new G4LogicalVolume(downGuard_LR_box, NistMan->FindOrBuildMaterial("G4_Fe"), "downGuard_L_log",0,0,0);
  downGuard_L_log->SetVisAttributes(new G4VisAttributes(G4Color::Blue()));
  maxStep=0.00001*mm;
  //downGuard_L_log->SetUserLimits(new G4UserLimits(maxStep));
  G4VPhysicalVolume* downGuard_L_phys = new G4PVPlacement(G4Transform3D(rotForwardSp,
									G4ThreeVector(pos_DGUARD_L[XCOORD],
										      pos_DGUARD_L[YCOORD],
										      pos_DGUARD_L[ZCOORD]).rotateY(fSpectrometerAngle)),
      "downGuard_L_phys",
      downGuard_L_log,
      physMother,
      false,
      0);

  G4LogicalVolume*  downGuard_R_log = new G4LogicalVolume(downGuard_LR_box, NistMan->FindOrBuildMaterial("G4_Fe"), "downGuard_R_log",0,0,0);
  downGuard_R_log->SetVisAttributes(new G4VisAttributes(G4Color::Blue()));
  maxStep=0.00001*mm;
  //downGuard_R_log->SetUserLimits(new G4UserLimits(maxStep));
  G4VPhysicalVolume*  downGuard_R_phys = new G4PVPlacement(G4Transform3D(rotForwardSp,
									 G4ThreeVector(pos_DGUARD_R[XCOORD],
										       pos_DGUARD_R[YCOORD],
										       pos_DGUARD_R[ZCOORD]).rotateY(fSpectrometerAngle)),
      "downGuard_R_phys",
      downGuard_R_log,
      physMother,
      false,
      0);
  
#endif
  
  //Detectors inside KURAMA//
  G4ThreeVector DC1Pos ( -pos.x(), -pos.y(), - 635.4675*mm);
  G4RotationMatrix DC1Rot;

  G4ThreeVector CH1Pos ( 37.- pos.x(), 1. - pos.y(), - 387.*mm);
  G4RotationMatrix CH1Rot;

  G4ThreeVector DummyDetectorPos2 ( 0., 0., -820.0 + 100.0 + 0.01/2.0);// End HC 
  G4RotationMatrix DummyDetectorRot2;

  MakeDummyDetector2( physMother, DummyDetectorPos2, DummyDetectorRot2 );
  MakeDC1( physMother, DC1Pos, DC1Rot);
  MakeCH1( physMother, CH1Pos, CH1Rot);

  G4ThreeVector DummyDetectorPos4 ( 0., 0., - 635.4675*mm + 22.6*mm + 0.01/2.0*mm);// End DC1
  G4RotationMatrix DummyDetectorRot4;
  //MakeDummyDetector4( physMother, DummyDetectorPos4, DummyDetectorRot4 );


}


void DetectorConstruction::MakeDC1(G4VPhysicalVolume *pMother, G4ThreeVector &pos, G4RotationMatrix &rot)
{
  /// DC1 Physics World ///
  /// Physics World ///
  // param //
  const G4double X_PW = 289.0*mm;
  const G4double Y_PW=  215.0*mm;
  const G4double Z_PW  = 146.0/2.0*mm; 
  std::string fullNamePW = "DC1_PhysicsWorld";
  // //
  
  G4Box* solidDetectorPW = new G4Box(G4String(fullNamePW.c_str()),
				     X_PW,
				     Y_PW,
				     Z_PW
				     );

  G4LogicalVolume* logicDetectorPW = new G4LogicalVolume(solidDetectorPW,
							 mat_air,
							 G4String(fullNamePW.c_str()));
  logicDetectorPW->SetVisAttributes(new G4VisAttributes(G4Color::Blue()));
  
  G4RotationMatrix rot_PW = rot;
  G4ThreeVector TVPW = rot*pos;
  G4VPhysicalVolume *physMother = new G4PVPlacement(G4Transform3D(rot_PW, TVPW),
						    G4String(fullNamePW.c_str()),
						    logicDetectorPW,
						    pMother,
						    false,
						    0
						    );

  int XCOORD = 0;
  int YCOORD = 1;
  int ZCOORD = 2;

  G4double size_DC1[3];
  size_DC1[XCOORD] = 289.0*2.*mm;
  size_DC1[YCOORD] = 215.0*2.*mm;
  size_DC1[ZCOORD] = 146.0*mm;

  G4double size_DC1Plane[3];
  size_DC1Plane[XCOORD] = 6.0*67.0*0.5*mm;
  size_DC1Plane[YCOORD] = 6.0*42.0*0.5*mm;
  size_DC1Plane[ZCOORD] = 0.0001*mm;

  // DC1 //                                                                                        
  G4Box* DC1_box = new G4Box("DC1_box",size_DC1[XCOORD]/2,size_DC1[YCOORD]/2,size_DC1[ZCOORD]/2);
  G4LogicalVolume*  DC1_log = new G4LogicalVolume(DC1_box, mat_ArGas, "DC1_log",0,0,0);
  DC1_log->SetVisAttributes(G4Color::Green());
  G4RotationMatrix dummy_rot;
  G4VPhysicalVolume * DC1_phys = new G4PVPlacement(G4Transform3D(dummy_rot, G4ThreeVector(0., 0., 0.)),
					      "DC1_phys",
					      DC1_log,
					      physMother,
					      false,
					      0);

  // DC1 Planes //                                                                                      
  G4Box* DC1Plane_box = new G4Box("DC1Plane_box",
				  size_DC1Plane[XCOORD],size_DC1Plane[YCOORD],size_DC1Plane[ZCOORD]);
  G4LogicalVolume*    DC1Plane_log[6];
  G4VPhysicalVolume*   DC1Plane_phys[6];
  G4double pos_DC1Plane[3]={0.};
  char name1[30], name2[30];
  for (int i=0; i<6; i++) {
    switch (i) {
    case 0:
      sprintf(name1, "DC1U_log");
      sprintf(name2, "DC1U_phys");
      pos_DC1Plane[XCOORD] = 0.0*mm;
      pos_DC1Plane[YCOORD] = 0.0*mm;
      pos_DC1Plane[ZCOORD] = -22.5985*mm;
      break;
    case 1:
      sprintf(name1, "DC1Up_log");
      sprintf(name2, "DC1Up_phys");
      pos_DC1Plane[XCOORD] = 0.0*mm;
      pos_DC1Plane[YCOORD] = 0.0*mm;
      pos_DC1Plane[ZCOORD] = -17.4015*mm;
      break;
    case 2:
      sprintf(name1, "DC1X_log");
      sprintf(name2, "DC1X_phys");
      pos_DC1Plane[XCOORD] = 0.0*mm;
      pos_DC1Plane[YCOORD] = 0.0*mm;
      pos_DC1Plane[ZCOORD] = -2.5985*mm;
      break;
    case 3:
      sprintf(name1, "DC1Xp_log");
      sprintf(name2, "DC1Xp_phys");
      pos_DC1Plane[XCOORD] = 0.0*mm;
      pos_DC1Plane[YCOORD] = 0.0*mm;
      pos_DC1Plane[ZCOORD] = 2.5985*mm;
      break;
    case 4:
      sprintf(name1, "DC1V_log");
      sprintf(name2, "DC1V_phys");
      pos_DC1Plane[XCOORD] = 0.0*mm;
      pos_DC1Plane[YCOORD] = 0.0*mm;
      pos_DC1Plane[ZCOORD] = 17.4015*mm;
      break;
    case 5:
      sprintf(name1, "DC1Vp_log");
      sprintf(name2, "DC1Vp_phys");
      pos_DC1Plane[XCOORD] = 0.0*mm;
      pos_DC1Plane[YCOORD] = 0.0*mm;
      pos_DC1Plane[ZCOORD] = 22.5985*mm;
      break;
    }
    DC1Plane_log[i] = new G4LogicalVolume(DC1Plane_box, mat_ArGas, name1,0,0,0);
    DC1Plane_phys[i] = new G4PVPlacement(G4Transform3D(dummy_rot,
						       G4ThreeVector(pos_DC1Plane[XCOORD],pos_DC1Plane[YCOORD],pos_DC1Plane[ZCOORD])),
					 name2,
					 DC1Plane_log[i],
					 DC1_phys,
					 false,
					 i);
  }


  // sensitive detector //
  G4SDManager *SDMan = G4SDManager::GetSDMpointer();
  G4VSensitiveDetector *dcSD = SDMan -> FindSensitiveDetector("/DC");

  if(!dcSD)
    {
      DCSD* dcSD = new DCSD("/DC");
      SDMan -> AddNewDetector( dcSD );
    }

  G4VSensitiveDetector *dcSD1 = SDMan -> FindSensitiveDetector("/DC");
  for(int i=0; i<6; i++)
    {
      DC1Plane_log[i]->SetSensitiveDetector( dcSD1 );
    }

}


void DetectorConstruction::MakeCH1(G4VPhysicalVolume *pMother, G4ThreeVector &pos, G4RotationMatrix &rot)
{
  /// CH1 Physics World ///
  /// Physics World ///
  // param //
  const G4double X_PW = (11.5*64.-63.0)/2.0*mm;
  const G4double Y_PW=  400.0/2.0*mm;
  const G4double Z_PW  = 4.0/2.0*mm; 
  std::string fullNamePW = "CH1_PhysicsWorld";
  // //
  
  G4Box* solidDetectorPW = new G4Box(G4String(fullNamePW.c_str()),
				     X_PW,
				     Y_PW,
				     Z_PW
				     );

  G4LogicalVolume* logicDetectorPW = new G4LogicalVolume(solidDetectorPW,
							 mat_air,
							 G4String(fullNamePW.c_str()));
  logicDetectorPW->SetVisAttributes(new G4VisAttributes(G4Color::Blue()));
  
  G4RotationMatrix rot_PW = rot;
  G4ThreeVector TVPW = rot*pos;
  G4VPhysicalVolume *physMother = new G4PVPlacement(G4Transform3D(rot_PW, TVPW),
						    G4String(fullNamePW.c_str()),
						    logicDetectorPW,
						    pMother,
						    false,
						    0
						    );
  // CH1 //
  int CHMAX = 64;
  int XCOORD = 0;
  int YCOORD = 1;
  int ZCOORD = 2;

  G4double size_CH[3];
  size_CH[XCOORD] = 11.5*mm;
  size_CH[YCOORD] = 400.0*mm;
  size_CH[ZCOORD] = 2.0*mm;

  //--------------CH                                                                                         
  G4Box* CH_box = new G4Box("CH_box",size_CH[XCOORD]/2,size_CH[YCOORD]/2,size_CH[ZCOORD]/2);
  G4double CH_Overlap = 1.*mm;
  G4LogicalVolume* CH_log[CHMAX];
  G4VPhysicalVolume* CH_phys[CHMAX];
  double pos_CH_orgX = 0.0;
  G4Colour aqua(0.247, 0.8, 1.0);
  char name1[30], name2[30];
  G4double pos_CH[3] = {0.};
  G4RotationMatrix dummy_rot;

  for (int i=0; i<CHMAX; i++) 
    {
      sprintf(name1, "CH%d_log", i);
      CH_log[i] = new G4LogicalVolume(CH_box, mat_Scin, name1,0,0,0);
      CH_log[i]->SetVisAttributes(aqua);
      pos_CH[XCOORD] = pos_CH_orgX - CHMAX/2*(size_CH[XCOORD]-CH_Overlap)- CH_Overlap/2 + size_CH[XCOORD]/2 + (size_CH[XCOORD]-CH_Overlap)*i;
      sprintf(name1, "CH%d", i);
      if(i%2==0)
	{
	  CH_phys[i] = new G4PVPlacement(G4Transform3D(dummy_rot,
					 G4ThreeVector(pos_CH[XCOORD],pos_CH[YCOORD],pos_CH[ZCOORD]-1.*mm)),
					 name1,
					 CH_log[i],
					 physMother,
					 false,
					 i);
	}
      else if(i%2==1)
	{
	  CH_phys[i] = new G4PVPlacement(G4Transform3D(dummy_rot,
					 G4ThreeVector(pos_CH[XCOORD],pos_CH[YCOORD],pos_CH[ZCOORD]+1.*mm)),
					 name1,
					 CH_log[i],
					 physMother,
					 false,
					 i);
	}
    }

  // sensitive detector //
  G4SDManager *SDMan = G4SDManager::GetSDMpointer();
  CH1SD* ch1SD = new CH1SD("/CH1");
  SDMan -> AddNewDetector( ch1SD );
  for(int i=0; i<CHMAX; i++)
    {
      CH_log[i]->SetSensitiveDetector( ch1SD );
    }
 
}



void DetectorConstruction::MakeDC2(G4VPhysicalVolume *pMother, G4ThreeVector &pos, G4RotationMatrix &rot)
{
  /// DC2 Physics World ///
  /// Physics World ///
  // param //
  const G4double X_PW = 1186.5/2.0*mm;
  const G4double Y_PW=  1186.5/2.0*mm;
  const G4double Z_PW  = 100.0/2.0*mm; 
  std::string fullNamePW = "DC2_PhysicsWorld";
  // //
  
  G4Box* solidDetectorPW = new G4Box(G4String(fullNamePW.c_str()),
				     X_PW,
				     Y_PW,
				     Z_PW
				     );

  G4LogicalVolume* logicDetectorPW = new G4LogicalVolume(solidDetectorPW,
							 mat_air,
							 G4String(fullNamePW.c_str()));
  logicDetectorPW->SetVisAttributes(new G4VisAttributes(G4Color::Blue()));
  
  G4RotationMatrix rot_PW = rot;
  G4ThreeVector TVPW = rot*pos;
  G4VPhysicalVolume *physMother = new G4PVPlacement(G4Transform3D(rot_PW, TVPW),
						    G4String(fullNamePW.c_str()),
						    logicDetectorPW,
						    pMother,
						    false,
						    0
						    );

  int XCOORD = 0;
  int YCOORD = 1;
  int ZCOORD = 2;

  G4double size_DC2[3];
  size_DC2[XCOORD] = 1186.5*mm;//# of wires are 128                                                          
  size_DC2[YCOORD] = 1186.5*mm;//# of wires are 128                                                          
  //size_DC2[ZCOORD] = 45.0*mm;                                                                              
  size_DC2[ZCOORD] = 100.0*mm;

  G4double size_DC2Plane[3];
  size_DC2Plane[XCOORD] = 9.0*128.0*0.5*mm;
  size_DC2Plane[YCOORD] = 9.0*128.0*0.5*mm;
  size_DC2Plane[ZCOORD] = 0.0001*mm;

  // DC2 //                                                                                        
  G4Box* DC2_box = new G4Box("DC2_box",size_DC2[XCOORD]/2,size_DC2[YCOORD]/2,size_DC2[ZCOORD]/2);
  G4LogicalVolume*  DC2_log = new G4LogicalVolume(DC2_box, mat_ArGas, "DC2_log",0,0,0);
  DC2_log->SetVisAttributes(G4Color::Green());
  G4RotationMatrix dummy_rot;
  G4VPhysicalVolume * DC2_phys = new G4PVPlacement(G4Transform3D(dummy_rot, G4ThreeVector(0.0, 0.0, 0.0)),
					      "DC2_phys",
					      DC2_log,
					      physMother,
					      false,
					      0);
  // DC2 Planes //                                                                                      
  G4Box* DC2Plane_box = new G4Box("DC2Plane_box",
				  size_DC2Plane[XCOORD],size_DC2Plane[YCOORD],size_DC2Plane[ZCOORD]);

  G4LogicalVolume* DC2Plane_log[4];
  G4VPhysicalVolume* DC2Plane_phys[4];
  G4double pos_DC2Plane[3]={0.};
  char name1[30], name2[30];
  for (int i=0; i<4; i++) {
    switch (i) {
    case 0:
      sprintf(name1, "DC2X_log");
      sprintf(name2, "DC2X_phys");
      pos_DC2Plane[XCOORD] = 0.0*mm;
      pos_DC2Plane[YCOORD] = 0.0*mm;
      pos_DC2Plane[ZCOORD] = -16*mm;
      break;
    case 1:
      sprintf(name1, "DC2Xp_log");
      sprintf(name2, "DC2Xp_phys");
      pos_DC2Plane[XCOORD] = 0.0*mm;
      pos_DC2Plane[YCOORD] = 0.0*mm;
      pos_DC2Plane[ZCOORD] = -8.206*mm;
      break;
    case 2:
      sprintf(name1, "DC2Y_log");
      sprintf(name2, "DC2Y_phys");
      pos_DC2Plane[XCOORD] = 0.0*mm;
      pos_DC2Plane[YCOORD] = 0.0*mm;
      pos_DC2Plane[ZCOORD] = 8.206*mm;
      break;
    case 3:
      sprintf(name1, "DC2Yp_log");
      sprintf(name2, "DC2Yp_phys");
      pos_DC2Plane[XCOORD] = 0.0*mm;
      pos_DC2Plane[YCOORD] = 0.0*mm;
      pos_DC2Plane[ZCOORD] = 16*mm;
      break;
    }
    DC2Plane_log[i] = new G4LogicalVolume(DC2Plane_box, mat_ArGas, name1,0,0,0);
    DC2Plane_phys[i] = new G4PVPlacement(G4Transform3D(dummy_rot, G4ThreeVector(pos_DC2Plane[XCOORD],pos_DC2Plane[YCOORD],pos_DC2Plane[ZCOORD])),
					 name2,
					 DC2Plane_log[i],
					 DC2_phys,
					 false,
					 i+10);
  }

  // sensitive detector //
  G4SDManager *SDMan = G4SDManager::GetSDMpointer();
  G4VSensitiveDetector *dcSD = SDMan -> FindSensitiveDetector("/DC");
  if(!dcSD)
    {
      DCSD* dcSD = new DCSD("/DC");
      SDMan -> AddNewDetector( dcSD );
    }

  G4VSensitiveDetector *dcSD1 = SDMan -> FindSensitiveDetector("/DC");
  for(int i=0; i<4; i++)
    {
      DC2Plane_log[i]->SetSensitiveDetector( dcSD1 );
    }
  
}


void DetectorConstruction::MakeDC3(G4VPhysicalVolume *pMother, G4ThreeVector &pos, G4RotationMatrix &rot)
{
  /// DC3 Physics World ///
  /// Physics World ///
  // param //
  const G4double X_PW = 1980.0/2.0*mm;
  const G4double Y_PW=  1300.0/2.0*mm;
  const G4double Z_PW  = 150.0/2.0*mm; 
  std::string fullNamePW = "DC3_PhysicsWorld";
  // //
  
  G4Box* solidDetectorPW = new G4Box(G4String(fullNamePW.c_str()),
				     X_PW,
				     Y_PW,
				     Z_PW
				     );

  G4LogicalVolume* logicDetectorPW = new G4LogicalVolume(solidDetectorPW,
							 mat_air,
							 G4String(fullNamePW.c_str()));
  logicDetectorPW->SetVisAttributes(new G4VisAttributes(G4Color::Blue()));
  
  G4RotationMatrix rot_PW = rot;
  G4ThreeVector TVPW = rot*pos;
  G4VPhysicalVolume *physMother = new G4PVPlacement(G4Transform3D(rot_PW, TVPW),
						    G4String(fullNamePW.c_str()),
						    logicDetectorPW,
						    pMother,
						    false,
						    0
						    );

  int XCOORD = 0;
  int YCOORD = 1;
  int ZCOORD = 2;
  G4double size_DC3[3];
  size_DC3[XCOORD] = 1980.0*mm;
  size_DC3[YCOORD] = 1300.0*mm;
  size_DC3[ZCOORD] = 150.*mm;

  G4double size_DC3Plane[3];
  size_DC3Plane[XCOORD] = 20.0*96.0*0.5*mm;
  size_DC3Plane[YCOORD] = 20.0*64.0*0.5*mm;
  size_DC3Plane[ZCOORD] = 0.0001*mm;

  //--------------DC3                                                                                        
  G4Box* DC3_box = new G4Box("DC3_box",size_DC3[XCOORD]/2,size_DC3[YCOORD]/2,size_DC3[ZCOORD]/2);
  G4LogicalVolume*  DC3_log = new G4LogicalVolume(DC3_box, mat_ArGas, "DC3_log",0,0,0);
  DC3_log->SetVisAttributes(G4Color::Green());
  //  pos_DC3[XCOORD] = par_cham->get_DCPlaneCenter(DC3X, XCOORD)*mm;                                        
  //  pos_DC3[YCOORD] = par_cham->get_DCPlaneCenter(DC3X, YCOORD)*mm;                                        
  //  pos_DC3[ZCOORD] = (par_cham->get_DCPlaneCenter(DC3X, ZCOORD)                                           
  //         +par_cham->get_DCPlaneCenter(DC3Yp, ZCOORD))*0.5*mm;                                            
  G4RotationMatrix dummy_rot;
  G4VPhysicalVolume* DC3_phys = new G4PVPlacement(G4Transform3D(dummy_rot,
					      G4ThreeVector(0., 0., 0.)),
					      "DC3_phys",
					      DC3_log,
					      physMother,
					      false,
					      0);

  
  //---------DC3 Planes                                                                                      
  G4Box* DC3Plane_box = new G4Box("DC3Plane_box",
				  size_DC3Plane[XCOORD],size_DC3Plane[YCOORD],size_DC3Plane[ZCOORD]);
  G4LogicalVolume* DC3Plane_log[4];
  G4VPhysicalVolume* DC3Plane_phys[4];
  G4double pos_DC3Plane[3]={0.};
  char name1[30], name2[30];

  for (int i=0; i<4; i++) {
    switch (i) {
    case 0:
      sprintf(name1, "DC3X_log");
      sprintf(name2, "DC3X_phys");
      pos_DC3Plane[XCOORD] = 0.0*mm;
      pos_DC3Plane[YCOORD] = 0.0*mm;
      pos_DC3Plane[ZCOORD] = -34.868*mm;
      break;
    case 1:
      sprintf(name1, "DC3Xp_log");
      sprintf(name2, "DC3Xp_phys");
      pos_DC3Plane[XCOORD] = 0.0*mm;
      pos_DC3Plane[YCOORD] = 0.0*mm;
      pos_DC3Plane[ZCOORD] = -17.547*mm;
      break;
    case 2:
      sprintf(name1, "DC3Y_log");
      sprintf(name2, "DC3Y_phys");
      pos_DC3Plane[XCOORD] = 0.0*mm;
      pos_DC3Plane[YCOORD] = 0.0*mm;
      pos_DC3Plane[ZCOORD] = 17.547*mm;
      break;
    case 3:
      sprintf(name1, "DC3Yp_log");
      sprintf(name2, "DC3Yp_phys");
      pos_DC3Plane[XCOORD] = 0.0*mm;
      pos_DC3Plane[YCOORD] = 0.0*mm;
      pos_DC3Plane[ZCOORD] = 34.868*mm;
      break;
    }
    DC3Plane_log[i] = new G4LogicalVolume(DC3Plane_box, mat_ArGas, name1,0,0,0);
    DC3Plane_phys[i] = new G4PVPlacement(G4Transform3D(dummy_rot,
						       G4ThreeVector(pos_DC3Plane[XCOORD],pos_DC3Plane[YCOORD],pos_DC3Plane[ZCOORD])),
					 name2,
					 DC3Plane_log[i],
					 DC3_phys,
					 false,
					 i+20);
  }

  // sensitive detector //
  G4SDManager *SDMan = G4SDManager::GetSDMpointer();
  G4VSensitiveDetector *dcSD = SDMan -> FindSensitiveDetector("/DC");
  if(!dcSD)
    {
      DCSD* dcSD = new DCSD("/DC");
      SDMan -> AddNewDetector( dcSD );
    }
  G4VSensitiveDetector *dcSD1 = SDMan -> FindSensitiveDetector("/DC");
  for(int i=0; i<4; i++)
    {
      DC3Plane_log[i]->SetSensitiveDetector( dcSD1 );
    }
}


void DetectorConstruction::MakeFTOF(G4VPhysicalVolume *pMother, G4ThreeVector &pos, G4RotationMatrix &rot)
{

  int XCOORD = 0;
  int YCOORD = 1;
  int ZCOORD = 2;
  int FTOFMAX = 24;
  
  G4double size_FTOF[3];
  size_FTOF[XCOORD] = 40.0*mm;
  size_FTOF[YCOORD] = 900.0*mm;
  size_FTOF[ZCOORD] = 15.0*mm;
  
  G4RotationMatrix rot_PW = rot;
  G4ThreeVector TVPW = rot*pos;
  G4RotationMatrix rot_in_PW;
  rot_in_PW.rotateY(12*deg); //E07
  rot_PW = rot_PW*rot_in_PW;
  
  //--------------FTOF                                                                                       
  G4Box* FTOF_mo_box = new G4Box("FTOF_mo_box",size_FTOF[XCOORD]*24+50.*mm,size_FTOF[YCOORD]+50.*mm,size_FTOF[ZCOORD]+50.*mm);
  G4LogicalVolume* FTOF_mo_log = new G4LogicalVolume(FTOF_mo_box, mat_air, "FTOF_mo_log",0,0,0);
  G4VPhysicalVolume* FTOF_mo_phys = new G4PVPlacement(G4Transform3D(rot_PW, TVPW),
						      "FTOF_mo_phys",
						      FTOF_mo_log,
						      pMother,
						      false,
						      0);
  
  G4VisAttributes* FTOF_mo_VisAtt= new G4VisAttributes(false, G4Colour(1.,0.,0.));
  FTOF_mo_log-> SetVisAttributes(FTOF_mo_VisAtt);
  
  ////end mother volume of FTOF                                                                              
  G4Box* FTOF_box = new G4Box("FTOF_box",size_FTOF[XCOORD],size_FTOF[YCOORD],size_FTOF[ZCOORD]);
  G4VPhysicalVolume* FTOF_phys[FTOFMAX];
  G4LogicalVolume* FTOF_log[FTOFMAX];
  
  G4double FTOF_Overlap=5.*mm;
  G4double pos_FTOF_bar[3];
  G4Colour aqua(0.247, 0.8, 1.0);

  G4RotationMatrix dummy_rot;
  char name1[30], name2[30];

  for (int i=0; i<FTOFMAX; i++) 
    {
      sprintf(name1, "FTOF%d_log", i);
      FTOF_log[i] = new G4LogicalVolume(FTOF_box, mat_Scin, name1,0,0,0);
      FTOF_log[i]->SetVisAttributes(aqua);
      pos_FTOF_bar[XCOORD]=-FTOFMAX/2*(size_FTOF[XCOORD]*2-FTOF_Overlap) - FTOF_Overlap/2 + size_FTOF[XCOORD] + (size_FTOF[XCOORD]*2-FTOF_Overlap)*i;
      sprintf(name1, "FTOF%d", i);
      if(i%2==0)
	{
	  FTOF_phys[i] = new G4PVPlacement(G4Transform3D(dummy_rot,
							 G4ThreeVector(pos_FTOF_bar[XCOORD],0.,size_FTOF[ZCOORD]+5.0*mm)),
					   name1,
					   FTOF_log[i],
					   FTOF_mo_phys,
					   false,
					   i);
	}
      else if(i%2==1)
	{
	  FTOF_phys[i] = new G4PVPlacement(G4Transform3D(dummy_rot,
							 G4ThreeVector(pos_FTOF_bar[XCOORD],0.,-size_FTOF[ZCOORD]-5.0*mm)),
					   name1,
					   FTOF_log[i],
					   FTOF_mo_phys,
					   false,
					   i);
	}
    }


  // sensitive detector //
  G4SDManager *SDMan = G4SDManager::GetSDMpointer();
  FTOFSD* ftofSD = new FTOFSD("/FTOF");
  SDMan -> AddNewDetector( ftofSD );
  for(int i=0; i<24; i++)
    {
      FTOF_log[i]->SetSensitiveDetector( ftofSD );
    }

 
}

// Make Dummy Detector //
void DetectorConstruction::MakeDummyDetector1(G4VPhysicalVolume *pMother, G4ThreeVector &pos, G4RotationMatrix &rot)
{
  /// Dummy Detector Physics World ///
  /// Physics World ///
  // param //
  const G4double X_PW = 1550.0/2.0*mm;
  const G4double Y_PW=  950.0/2.0*mm;
  const G4double Z_PW  = 0.01/2.0*mm; 
  std::string fullNamePW = "Dummy_Detector_PhysicsWorld";
  // //
  
  G4Box* solidDetectorPW = new G4Box(G4String(fullNamePW.c_str()),
				     X_PW,
				     Y_PW,
				     Z_PW
				     );

  G4LogicalVolume* logicDetectorPW = new G4LogicalVolume(solidDetectorPW,
							 mat_air,
							 G4String(fullNamePW.c_str()));

  logicDetectorPW->SetVisAttributes(new G4VisAttributes(G4Color::Yellow()));
  
  G4RotationMatrix rot_PW = rot;
  G4ThreeVector TVPW = rot*pos;
  G4VPhysicalVolume *physMother = new G4PVPlacement(G4Transform3D(rot_PW, TVPW),
						    G4String(fullNamePW.c_str()),
						    logicDetectorPW,
						    pMother,
						    false,
						    0
						    );

  // Dummy Detector //                                                                                        
  G4Box* dummy_box = new G4Box("dummy_box", X_PW, Y_PW, Z_PW);
  G4LogicalVolume*  dummy_log = new G4LogicalVolume(dummy_box, mat_air, "dummy_log", 0, 0, 0);
  dummy_log->SetVisAttributes(G4Color::Yellow());

  G4RotationMatrix dummy_rot;
  G4VPhysicalVolume* dummy_phys = new G4PVPlacement(G4Transform3D(dummy_rot,
								  G4ThreeVector(0., 0., 0.)),
						    "dummy_phys",
						    dummy_log,
						    physMother,
						    false,
						    100);

  // sensitive detector //
  G4SDManager *SDMan = G4SDManager::GetSDMpointer();
  G4VSensitiveDetector *dcSD = SDMan -> FindSensitiveDetector("/DC");

  if(!dcSD)
    {
      DCSD* dcSD = new DCSD("/DC");
      SDMan -> AddNewDetector( dcSD );
    }

  G4VSensitiveDetector *dcSD1 = SDMan -> FindSensitiveDetector("/DC");
  dummy_log->SetSensitiveDetector( dcSD1 );

}



void DetectorConstruction::MakeDummyDetector2(G4VPhysicalVolume *pMother, G4ThreeVector &pos, G4RotationMatrix &rot)
{
  /// Dummy Detector Physics World ///
  /// Physics World ///
  // param //
  const G4double X_PW = 2200.0/2.0*mm;
  const G4double Y_PW=  1500.0/2.0*mm;
  const G4double Z_PW  = 0.01/2.0*mm; 
  std::string fullNamePW = "Dummy_Detector_PhysicsWorld";
  // //
  
  G4Box* solidDetectorPW = new G4Box(G4String(fullNamePW.c_str()),
				     X_PW,
				     Y_PW,
				     Z_PW
				     );

  G4LogicalVolume* logicDetectorPW = new G4LogicalVolume(solidDetectorPW,
							 mat_air,
							 G4String(fullNamePW.c_str()));

  logicDetectorPW->SetVisAttributes(new G4VisAttributes(G4Color::Yellow()));
  
  G4RotationMatrix rot_PW = rot;
  G4ThreeVector TVPW = rot*pos;
  G4VPhysicalVolume *physMother = new G4PVPlacement(G4Transform3D(rot_PW, TVPW),
						    G4String(fullNamePW.c_str()),
						    logicDetectorPW,
						    pMother,
						    false,
						    0
						    );

  // Dummy Detector //                                                                                        
  G4Box* dummy_box = new G4Box("dummy_box", X_PW, Y_PW, Z_PW);
  G4LogicalVolume*  dummy_log = new G4LogicalVolume(dummy_box, mat_air, "dummy_log", 0, 0, 0);
  dummy_log->SetVisAttributes(G4Color::Yellow());

  G4RotationMatrix dummy_rot;
  G4VPhysicalVolume* dummy_phys = new G4PVPlacement(G4Transform3D(dummy_rot,
								  G4ThreeVector(0., 0., 0.)),
						    "dummy_phys",
						    dummy_log,
						    physMother,
						    false,
						    101);

  // sensitive detector //
  G4SDManager *SDMan = G4SDManager::GetSDMpointer();
  G4VSensitiveDetector *dcSD = SDMan -> FindSensitiveDetector("/DC");

  if(!dcSD)
    {
      DCSD* dcSD = new DCSD("/DC");
      SDMan -> AddNewDetector( dcSD );
    }

  G4VSensitiveDetector *dcSD1 = SDMan -> FindSensitiveDetector("/DC");
  dummy_log->SetSensitiveDetector( dcSD1 );

}

void DetectorConstruction::MakeDummyDetector3(G4VPhysicalVolume *pMother, G4ThreeVector &pos, G4RotationMatrix &rot)
{
  /// Dummy Detector Physics World ///
  /// Physics World ///
  // param //
  const G4double X_PW = 2200.0/2.0*mm;
  const G4double Y_PW=  1500.0/2.0*mm;
  const G4double Z_PW  = 0.01/2.0*mm; 
  std::string fullNamePW = "Dummy_Detector_PhysicsWorld";
  // //
  
  G4Box* solidDetectorPW = new G4Box(G4String(fullNamePW.c_str()),
				     X_PW,
				     Y_PW,
				     Z_PW
				     );

  G4LogicalVolume* logicDetectorPW = new G4LogicalVolume(solidDetectorPW,
							 mat_air,
							 G4String(fullNamePW.c_str()));

  logicDetectorPW->SetVisAttributes(new G4VisAttributes(G4Color::Yellow()));
  
  G4RotationMatrix rot_PW = rot;
  G4ThreeVector TVPW = rot*pos;
  G4VPhysicalVolume *physMother = new G4PVPlacement(G4Transform3D(rot_PW, TVPW),
						    G4String(fullNamePW.c_str()),
						    logicDetectorPW,
						    pMother,
						    false,
						    0
						    );

  // Dummy Detector //                                                                                        
  G4Box* dummy_box = new G4Box("dummy_box", X_PW, Y_PW, Z_PW);
  G4LogicalVolume*  dummy_log = new G4LogicalVolume(dummy_box, mat_air, "dummy_log", 0, 0, 0);
  dummy_log->SetVisAttributes(G4Color::Yellow());

  G4RotationMatrix dummy_rot;
  G4VPhysicalVolume* dummy_phys = new G4PVPlacement(G4Transform3D(dummy_rot,
								  G4ThreeVector(0., 0., 0.)),
						    "dummy_phys",
						    dummy_log,
						    physMother,
						    false,
						    102);

  // sensitive detector //
  G4SDManager *SDMan = G4SDManager::GetSDMpointer();
  G4VSensitiveDetector *dcSD = SDMan -> FindSensitiveDetector("/DC");

  if(!dcSD)
    {
      DCSD* dcSD = new DCSD("/DC");
      SDMan -> AddNewDetector( dcSD );
    }

  G4VSensitiveDetector *dcSD1 = SDMan -> FindSensitiveDetector("/DC");
  dummy_log->SetSensitiveDetector( dcSD1 );

}


// KURAMA Magnet E42 Setup //
void DetectorConstruction::MakeKuramaMagnetE42(G4VPhysicalVolume *pMother, G4ThreeVector &pos, G4RotationMatrix &rot)
{

  /// KURAMA Magnet Physics World ///
  /// Physics World ///
  // param //
  const G4double X_PW = 2200.0/2.0*mm;
  const G4double Y_PW=  1500.0/2.0*mm;
  const G4double Z_PW  = 1640.0/2.0*mm; 
  std::string fullNamePW = "KuramaMagnet_PhysicsWorld";
  // //
  
  G4Box* solidDetectorPW = new G4Box(G4String(fullNamePW.c_str()),
				     X_PW,
				     Y_PW,
				     Z_PW
				     );

  G4LogicalVolume* logicDetectorPW = new G4LogicalVolume(solidDetectorPW,
							 mat_air,
							 G4String(fullNamePW.c_str()));
  logicDetectorPW->SetVisAttributes(new G4VisAttributes(G4Color::Blue()));
  
  G4RotationMatrix rot_PW = rot;
  G4ThreeVector TVPW = rot*pos;
  G4VPhysicalVolume *physMother = new G4PVPlacement(G4Transform3D(rot_PW, TVPW),
						    G4String(fullNamePW.c_str()),
						    logicDetectorPW,
						    pMother,
						    false,
						    0
						    );


  // Magnet Parameters // 
  const double kurama_gap=400.0;
  int XCOORD = 0;
  int YCOORD = 1;
  int ZCOORD = 2;
  G4double size_MFIELD[3];
  size_MFIELD[XCOORD] = 500.0*mm;
  size_MFIELD[YCOORD] = kurama_gap*mm;
  size_MFIELD[ZCOORD] = 400.0*mm;
  //size_MFIELD[XCOORD] = 0.0*mm;
  //size_MFIELD[YCOORD] = 0.0*mm;
  //size_MFIELD[ZCOORD] = 0.0*mm;

  G4double size_COIL1[3];
  size_COIL1[XCOORD] = 900./2*mm;
  size_COIL1[YCOORD] = 193./2*mm;
  size_COIL1[ZCOORD] = 280./2*mm;

  G4double size_COIL2[3];
  size_COIL2[XCOORD] = 193./2.*mm;
  size_COIL2[YCOORD] = 58.5*mm;
  size_COIL2[ZCOORD] = 280.0/2.*mm;

  G4double size_COIL3[3];
  size_COIL3[XCOORD] = 193.0/2.*mm;
  size_COIL3[YCOORD] = 68.5*mm;
  size_COIL3[ZCOORD] = 280.0/2.*mm;

  G4double size_COIL4[3];
  size_COIL4[XCOORD] = 193.0/2.*mm;
  size_COIL4[YCOORD] = 280.0/2.*mm;
  size_COIL4[ZCOORD] = 740.0/2.*mm;

  G4double size_COIL5[3];
  size_COIL5[XCOORD] = 900.0/2.*mm;
  size_COIL5[YCOORD] = 214./2.*mm;
  size_COIL5[ZCOORD] = 280./2.*mm;



  G4double size_YOKE_UD[3];
  size_YOKE_UD[XCOORD] = 2200.0/2.*mm;
  size_YOKE_UD[YCOORD] = 370.0/2.*mm;
  size_YOKE_UD[ZCOORD] = size_MFIELD[ZCOORD]*mm;

  G4double size_YOKE_LR[3];
  size_YOKE_LR[XCOORD] = 200.0*mm;
  size_YOKE_LR[YCOORD] = size_MFIELD[YCOORD]*mm;
  size_YOKE_LR[ZCOORD] = size_MFIELD[ZCOORD]*mm;

  G4double size_YOKE_LR_GapSpace[3]; //E42
  size_YOKE_LR_GapSpace[ZCOORD] = 400.0*mm;
  size_YOKE_LR_GapSpace[XCOORD] = 220.0/2.*mm;
  size_YOKE_LR_GapSpace[YCOORD] = 120.0/2.*mm;

  G4double size_UGUARD_UD[3];
  size_UGUARD_UD[XCOORD] = 950.0*mm;
  size_UGUARD_UD[YCOORD] = 310.0*mm;
  size_UGUARD_UD[ZCOORD] = 50.0*mm;

  G4double size_UGUARD_LR[3];
  //size_UGUARD_LR[XCOORD] = 400.0*mm;
  size_UGUARD_LR[XCOORD] = 325.0*mm;
  size_UGUARD_LR[YCOORD] = 150.0*mm;
  size_UGUARD_LR[ZCOORD] = 50.0*mm;

  G4double size_UGUARD_R1[3];
  size_UGUARD_R1[ZCOORD] =  50.0*mm;
  size_UGUARD_R1[XCOORD] = 590.0/2*mm;
  size_UGUARD_R1[YCOORD] = 30.0*mm;

  G4double size_UGUARD_R2[3];
  size_UGUARD_R2[ZCOORD] =  50.0*mm;
  size_UGUARD_R2[XCOORD] = 325.0*mm;
  size_UGUARD_R2[YCOORD] = 60.0*mm;

  G4double size_DGUARD_UD[3];
  size_DGUARD_UD[XCOORD] = 800.0*mm;
  size_DGUARD_UD[YCOORD] = 210.0*mm;
  size_DGUARD_UD[ZCOORD] = 50.0*mm;

  G4double size_DGUARD_LR[3];
  size_DGUARD_LR[XCOORD] = 125.0*mm;
  size_DGUARD_LR[YCOORD] = size_MFIELD[YCOORD]*mm+150.0*mm;
  size_DGUARD_LR[ZCOORD] = 50.0*mm;

  
  G4double size_CH[3];
  size_CH[XCOORD] = 11.5*mm;
  size_CH[YCOORD] = 400.0*mm;
  size_CH[ZCOORD] = 2.0*mm;

  G4double size_FTOF[3];
  size_FTOF[XCOORD] = 40.0*mm;
  size_FTOF[YCOORD] = 900.0*mm;
  size_FTOF[ZCOORD] = 15.0*mm;


  G4double k_move_x=5.;
  G4double k_pos_z=1500.;

  G4double pos_MFIELD[3];
  pos_MFIELD[XCOORD] = 0.0*mm;
  pos_MFIELD[YCOORD] = 0.0*mm;
  pos_MFIELD[ZCOORD] = 0.0*mm;
  //pos_MFIELD[XCOORD] = k_move_x*mm;
  //pos_MFIELD[YCOORD] = 0.0*mm;
  //pos_MFIELD[ZCOORD] = k_pos_z*mm;
  //  pos_MFIELD[ZCOORD] = (1105.+400.)*mm;
  //  pos_MFIELD[ZCOORD] = 0.*mm;


  G4double pos_COIL1[3];
  pos_COIL1[XCOORD] = pos_MFIELD[XCOORD];
  pos_COIL1[YCOORD] = size_MFIELD[YCOORD]+size_YOKE_UD[YCOORD]*2. -(size_COIL1[YCOORD]+20.*mm);
  pos_COIL1[ZCOORD] = pos_MFIELD[ZCOORD]-size_MFIELD[ZCOORD]-(size_COIL1[ZCOORD]+20.*mm);

  G4double pos_COIL4L[3];
  pos_COIL4L[XCOORD] = pos_MFIELD[XCOORD]+size_MFIELD[XCOORD]+size_COIL4[XCOORD];
  pos_COIL4L[YCOORD] = size_MFIELD[YCOORD]/2;
  pos_COIL4L[ZCOORD] = pos_MFIELD[ZCOORD];


  G4double pos_COIL4R[3];
  pos_COIL4R[XCOORD] = pos_MFIELD[XCOORD]-size_MFIELD[XCOORD]-size_COIL4[XCOORD];
  pos_COIL4R[YCOORD] = size_MFIELD[YCOORD]/2;
  pos_COIL4R[ZCOORD] = pos_MFIELD[ZCOORD];


  G4double pos_COIL2L[3];
  pos_COIL2L[XCOORD] = pos_MFIELD[XCOORD]+size_MFIELD[XCOORD]+size_COIL4[XCOORD];
  pos_COIL2L[YCOORD] = (pos_COIL4L[YCOORD]+pos_COIL1[YCOORD])/2+(size_COIL4[YCOORD]-size_COIL1[YCOORD])/2;
  pos_COIL2L[ZCOORD] = pos_MFIELD[ZCOORD]-size_MFIELD[ZCOORD]-(size_COIL1[ZCOORD]+20.*mm);

  G4double pos_COIL2R[3];
  pos_COIL2R[XCOORD] = pos_MFIELD[XCOORD]-size_MFIELD[XCOORD]-size_COIL4[XCOORD];
  pos_COIL2R[YCOORD] = (pos_COIL4R[YCOORD]+pos_COIL1[YCOORD])/2+(size_COIL4[YCOORD]-size_COIL1[YCOORD])/2;
  pos_COIL2R[ZCOORD] = pos_MFIELD[ZCOORD]-size_MFIELD[ZCOORD]-(size_COIL1[ZCOORD]+20.*mm);


  G4double pos_COIL5[3];
  pos_COIL5[XCOORD] = pos_MFIELD[XCOORD];
  pos_COIL5[YCOORD] = size_MFIELD[YCOORD]+size_YOKE_UD[YCOORD]*2. -(size_COIL5[YCOORD]);
  pos_COIL5[ZCOORD] = pos_MFIELD[ZCOORD] + size_MFIELD[ZCOORD]+(size_COIL5[ZCOORD]+21.*mm);


  G4double pos_COIL3L[3];
  pos_COIL3L[XCOORD] = pos_MFIELD[XCOORD]+size_MFIELD[XCOORD]+size_COIL4[XCOORD];
  pos_COIL3L[YCOORD] = (pos_COIL4L[YCOORD]+pos_COIL5[YCOORD])/2+(size_COIL4[YCOORD]-size_COIL5[YCOORD])/2;
  pos_COIL3L[ZCOORD] = pos_MFIELD[ZCOORD]+size_MFIELD[ZCOORD]+(size_COIL5[ZCOORD]+21.*mm);

  G4double pos_COIL3R[3];
  pos_COIL3R[XCOORD] = pos_MFIELD[XCOORD]-size_MFIELD[XCOORD]-size_COIL4[XCOORD];
  pos_COIL3R[YCOORD] = (pos_COIL4R[YCOORD]+pos_COIL5[YCOORD])/2+(size_COIL4[YCOORD]-size_COIL5[YCOORD])/2;
  pos_COIL3R[ZCOORD] = pos_MFIELD[ZCOORD]+size_MFIELD[ZCOORD]+(size_COIL5[ZCOORD]+21.*mm);



  G4double pos_YOKE_U[3];
  pos_YOKE_U[XCOORD] = pos_MFIELD[XCOORD];
  pos_YOKE_U[YCOORD] = size_MFIELD[YCOORD]+size_YOKE_UD[YCOORD];
  pos_YOKE_U[ZCOORD] = pos_MFIELD[ZCOORD];

  G4double pos_YOKE_D[3];
  pos_YOKE_D[XCOORD] = pos_MFIELD[XCOORD];
  pos_YOKE_D[YCOORD] = -(size_MFIELD[YCOORD]+size_YOKE_UD[YCOORD]);
  pos_YOKE_D[ZCOORD] = pos_MFIELD[ZCOORD];

  G4double pos_YOKE_L[3];
  pos_YOKE_L[XCOORD] = pos_MFIELD[XCOORD] + size_MFIELD[XCOORD]+size_YOKE_LR[XCOORD]+200.*mm;
  pos_YOKE_L[YCOORD] = 0.0*mm;
  pos_YOKE_L[ZCOORD] = pos_MFIELD[ZCOORD];

  G4double pos_YOKE_R[3];
  pos_YOKE_R[XCOORD] = pos_MFIELD[XCOORD] - (size_MFIELD[XCOORD]+size_YOKE_LR[XCOORD])-200.*mm;
  pos_YOKE_R[YCOORD] = 0.0*mm;
  pos_YOKE_R[ZCOORD] = pos_MFIELD[ZCOORD];

  G4double pos_YOKE_L_GapSpace[3];
  pos_YOKE_L_GapSpace[ZCOORD] = pos_YOKE_L[ZCOORD];
  pos_YOKE_L_GapSpace[XCOORD] = pos_YOKE_L[XCOORD]-size_YOKE_LR[XCOORD]-size_YOKE_LR_GapSpace[XCOORD];
  pos_YOKE_L_GapSpace[YCOORD] = 0.0*mm;

  G4double pos_YOKE_R_GapSpace[3];
  pos_YOKE_R_GapSpace[ZCOORD] = pos_YOKE_R[ZCOORD];
  pos_YOKE_R_GapSpace[XCOORD] = pos_YOKE_R[XCOORD]+size_YOKE_LR[XCOORD]+size_YOKE_LR_GapSpace[XCOORD];
  pos_YOKE_R_GapSpace[YCOORD] = 0.0*mm;


  ///up guard
  G4double pos_UGUARD_U[3];
  pos_UGUARD_U[XCOORD] = 0.0*mm;
  pos_UGUARD_U[YCOORD] = size_UGUARD_LR[YCOORD]+size_UGUARD_UD[YCOORD];
  pos_UGUARD_U[ZCOORD] = pos_MFIELD[ZCOORD] - 820.0*mm + size_UGUARD_UD[ZCOORD];

  G4double pos_UGUARD_D[3];
  pos_UGUARD_D[XCOORD] = 0.0*mm;
  pos_UGUARD_D[YCOORD] = -(size_UGUARD_LR[YCOORD]+size_UGUARD_UD[YCOORD]);
  pos_UGUARD_D[ZCOORD] = pos_MFIELD[ZCOORD]- 820.0*mm + size_UGUARD_UD[ZCOORD];

  G4double pos_UGUARD_L[3];
  //pos_UGUARD_L[XCOORD] = 0.0*mm + (300.0*mm + size_UGUARD_LR[XCOORD]); // 150 --> gap
  pos_UGUARD_L[XCOORD] = pos_UGUARD_U[XCOORD]+size_UGUARD_UD[XCOORD] - size_UGUARD_LR[XCOORD]; //E42
  pos_UGUARD_L[YCOORD] = 0.0*mm;
  pos_UGUARD_L[ZCOORD] = pos_MFIELD[ZCOORD]- 820.0*mm + size_UGUARD_LR[ZCOORD];

  G4double pos_UGUARD_R[3];
  //pos_UGUARD_R[XCOORD] = 0.0*mm - (300.0*mm + size_UGUARD_LR[XCOORD]);
  pos_UGUARD_R[XCOORD] = pos_UGUARD_U[XCOORD] -(size_UGUARD_UD[XCOORD] - size_UGUARD_LR[XCOORD]); //E42
  pos_UGUARD_R[YCOORD] = 0.0*mm;
  pos_UGUARD_R[ZCOORD] = pos_MFIELD[ZCOORD]- 820.0*mm + size_UGUARD_LR[ZCOORD];

  G4double pos_UGUARD_R1[3];
  pos_UGUARD_R1[ZCOORD] = pos_MFIELD[ZCOORD]- 820.0*mm + size_UGUARD_LR[ZCOORD];
  pos_UGUARD_R1[XCOORD] = pos_UGUARD_U[XCOORD] -(size_UGUARD_UD[XCOORD] - size_UGUARD_R1[XCOORD]);
  pos_UGUARD_R1[YCOORD] = 0.0*mm;

  G4double pos_UGUARD_R2U[3];
  pos_UGUARD_R2U[ZCOORD] = pos_MFIELD[ZCOORD]- 820.0*mm + size_UGUARD_LR[ZCOORD];
  pos_UGUARD_R2U[XCOORD] = pos_UGUARD_U[XCOORD] -(size_UGUARD_UD[XCOORD] - size_UGUARD_R2[XCOORD]);
  pos_UGUARD_R2U[YCOORD] = size_UGUARD_R1[YCOORD]+size_UGUARD_R2[YCOORD];

  G4double pos_UGUARD_R2D[3];
  pos_UGUARD_R2D[ZCOORD] = pos_MFIELD[ZCOORD]- 820.0*mm + size_UGUARD_LR[ZCOORD];
  pos_UGUARD_R2D[XCOORD] = pos_UGUARD_U[XCOORD] -(size_UGUARD_UD[XCOORD] - size_UGUARD_R2[XCOORD]);
  pos_UGUARD_R2D[YCOORD] = -(size_UGUARD_R1[YCOORD]+size_UGUARD_R2[YCOORD]);


  G4double pos_DGUARD_U[3];
  pos_DGUARD_U[XCOORD] = pos_MFIELD[XCOORD];
  pos_DGUARD_U[YCOORD] = size_DGUARD_LR[YCOORD]+size_DGUARD_UD[YCOORD];
  //pos_DGUARD_U[YCOORD] = size_DGUARD_LR[YCOORD]+size_DGUARD_UD[YCOORD] -200.; //Hwang
  pos_DGUARD_U[ZCOORD] = pos_MFIELD[ZCOORD] + 820.0*mm - size_DGUARD_UD[ZCOORD];

  G4double pos_DGUARD_D[3];
  pos_DGUARD_D[XCOORD] = pos_MFIELD[XCOORD];
  pos_DGUARD_D[YCOORD] = -(size_DGUARD_LR[YCOORD]+size_DGUARD_UD[YCOORD]);
  //pos_DGUARD_D[YCOORD] = -(size_DGUARD_LR[YCOORD]+size_DGUARD_UD[YCOORD]) + 200.; //Hwang
  pos_DGUARD_D[ZCOORD] = pos_MFIELD[ZCOORD]+ 820.0*mm - size_DGUARD_UD[ZCOORD];

  G4double pos_DGUARD_L[3];
  //pos_DGUARD_L[XCOORD] = pos_MFIELD[XCOORD] + (size_MFIELD[XCOORD] + 50*mm + size_DGUARD_LR[XCOORD]);
  pos_DGUARD_L[XCOORD] = pos_MFIELD[XCOORD] + (550.0*mm + size_DGUARD_LR[XCOORD]); 
  pos_DGUARD_L[YCOORD] = 0.0*mm;
  pos_DGUARD_L[ZCOORD] = pos_MFIELD[ZCOORD]+ 820.0*mm - size_DGUARD_LR[ZCOORD];

  G4double pos_DGUARD_R[3];
  //pos_DGUARD_R[XCOORD] = pos_MFIELD[XCOORD] - (size_MFIELD[XCOORD] + 50*mm + size_DGUARD_LR[XCOORD]);
  pos_DGUARD_R[XCOORD] = pos_MFIELD[XCOORD] - (550*mm + size_DGUARD_LR[XCOORD]); 
  pos_DGUARD_R[YCOORD] = 0.0*mm;
  pos_DGUARD_R[ZCOORD] = pos_MFIELD[ZCOORD]+ 820.0*mm - size_DGUARD_LR[ZCOORD];


  //G4double pos_DGUARD_U[3];
  //pos_DGUARD_U[XCOORD] = pos_MFIELD[XCOORD];
  //pos_DGUARD_U[YCOORD] = size_DGUARD_LR[YCOORD]+size_DGUARD_UD[YCOORD]/2.-200.;
  //pos_DGUARD_U[ZCOORD] = pos_MFIELD[ZCOORD] + 820.0*mm - size_DGUARD_UD[ZCOORD];

  //G4double pos_DGUARD_D[3];
  //pos_DGUARD_D[XCOORD] = pos_MFIELD[XCOORD];
  //pos_DGUARD_D[YCOORD] = -(size_DGUARD_LR[YCOORD]+size_DGUARD_UD[YCOORD]/2.)+200.;
  //pos_DGUARD_D[ZCOORD] = pos_MFIELD[ZCOORD]+ 820.0*mm - size_DGUARD_UD[ZCOORD];

  //G4double pos_DGUARD_L[3];
  //pos_DGUARD_L[XCOORD] = pos_MFIELD[XCOORD] + (size_MFIELD[XCOORD]+300.0*mm + size_DGUARD_LR[XCOORD]);
  //pos_DGUARD_L[YCOORD] = 0.0*mm;
  //pos_DGUARD_L[ZCOORD] = pos_MFIELD[ZCOORD]+ 820.0*mm - size_DGUARD_LR[ZCOORD];

  //G4double pos_DGUARD_R[3];
  //pos_DGUARD_R[XCOORD] = pos_MFIELD[XCOORD] - (size_MFIELD[XCOORD]+300.0*mm + size_DGUARD_LR[XCOORD]);
  //pos_DGUARD_R[YCOORD] = 0.0*mm;
  //pos_DGUARD_R[ZCOORD] = pos_MFIELD[ZCOORD]+ 820.0*mm - size_DGUARD_LR[ZCOORD];

  char nam[128];
  int  m_model=0;


  G4double maxStep;

  //G4RotationMatrix *rotForwardSp = new G4RotationMatrix();
  G4RotationMatrix rotForwardSp;

  //G4String spectrometer_ang = getenv("Spectrometer_ang");
  //G4double fSpectrometerAngle=atof( spectrometer_ang.c_str())*deg;//out side less 100 mm. 10+5*x < 100 mm is pad_in_num
  G4double fSpectrometerAngle = 0.0*deg;
  //  G4double fSpectrometerAngle=10*deg;
  //rotForwardSp->rotateY(-fSpectrometerAngle);
  rotForwardSp.rotateY(-fSpectrometerAngle);


  /////////////
  // Construct KURAMA Magnet
  //////////////coil1U
  G4Box* Coil1_box = new G4Box("Coil1_box",
      size_COIL1[XCOORD],size_COIL1[YCOORD],size_COIL1[ZCOORD]);
  G4LogicalVolume*  Coil1_log = new G4LogicalVolume(Coil1_box, NistMan->FindOrBuildMaterial("G4_Cu"), "Coil1_log",0,0,0);
  Coil1_log->SetVisAttributes(new G4VisAttributes(G4Color::Red()));

  G4VPhysicalVolume* Coil1U_phys = new G4PVPlacement(G4Transform3D(rotForwardSp, G4ThreeVector(pos_COIL1[XCOORD],pos_COIL1[YCOORD],pos_COIL1[ZCOORD]).rotateY(fSpectrometerAngle)),
						     "Coil1U_phys",
						     Coil1_log,
						     physMother,
						     false,
						     0);
  G4VPhysicalVolume* Coil1D_phys = new G4PVPlacement(G4Transform3D(rotForwardSp, G4ThreeVector(pos_COIL1[XCOORD],-pos_COIL1[YCOORD],pos_COIL1[ZCOORD]).rotateY(fSpectrometerAngle)),
      "Coil1D_phys",
      Coil1_log,
      physMother,
      false,
      0);

  //////////////coil4RLUD
  G4Box* Coil4_box = new G4Box("Coil4_box",
      size_COIL4[XCOORD],size_COIL4[YCOORD],size_COIL4[ZCOORD]);
  G4LogicalVolume*  Coil4_log = new G4LogicalVolume(Coil4_box, NistMan->FindOrBuildMaterial("G4_Cu"), "Coil4_log",0,0,0);
  Coil4_log->SetVisAttributes(new G4VisAttributes(G4Color::Red()));


  G4VPhysicalVolume* Coil4UR_phys = new G4PVPlacement(G4Transform3D(rotForwardSp, G4ThreeVector(pos_COIL4L[XCOORD],pos_COIL4L[YCOORD],pos_COIL4L[ZCOORD]).rotateY(fSpectrometerAngle)),
      "Coil4UR_phys",
      Coil4_log,
      physMother,
      false,
      0);
  G4VPhysicalVolume* Coil4UL_phys = new G4PVPlacement(G4Transform3D(rotForwardSp, G4ThreeVector(pos_COIL4R[XCOORD],pos_COIL4R[YCOORD],pos_COIL4R[ZCOORD]).rotateY(fSpectrometerAngle)),
      "Coil4UL_phys",
      Coil4_log,
      physMother,
      false,
      0);
  G4VPhysicalVolume* Coil4DR_phys = new G4PVPlacement(G4Transform3D(rotForwardSp, G4ThreeVector(pos_COIL4L[XCOORD],-pos_COIL4L[YCOORD],pos_COIL4L[ZCOORD]).rotateY(fSpectrometerAngle)),
      "Coil4DR_phys",
      Coil4_log,
      physMother,
      false,
      0);
  G4VPhysicalVolume* Coil4DL_phys = new G4PVPlacement(G4Transform3D(rotForwardSp, G4ThreeVector(pos_COIL4R[XCOORD],-pos_COIL4R[YCOORD],pos_COIL4R[ZCOORD]).rotateY(fSpectrometerAngle)),
      "Coil4DL_phys",
      Coil4_log,
      physMother,
      false,
      0);



  //////////////coil5UD
  G4Box* Coil5_box = new G4Box("Coil5_box",
      size_COIL5[XCOORD],size_COIL5[YCOORD],size_COIL5[ZCOORD]);
  G4LogicalVolume*  Coil5_log = new G4LogicalVolume(Coil5_box, NistMan->FindOrBuildMaterial("G4_Cu"), "Coil5_log",0,0,0);
  Coil5_log->SetVisAttributes(new G4VisAttributes(G4Color::Red()));


  G4VPhysicalVolume* Coil5U_phys = new G4PVPlacement(G4Transform3D(rotForwardSp, G4ThreeVector(pos_COIL5[XCOORD],pos_COIL5[YCOORD],pos_COIL5[ZCOORD]).rotateY(fSpectrometerAngle)),
      "Coil5U_phys",
      Coil5_log,
      physMother,
      false,
      0);
  G4VPhysicalVolume* Coil5D_phys = new G4PVPlacement(G4Transform3D(rotForwardSp, G4ThreeVector(pos_COIL5[XCOORD],-pos_COIL5[YCOORD],pos_COIL5[ZCOORD]).rotateY(fSpectrometerAngle)),
      "Coil5D_phys",
      Coil5_log,
      physMother,
      false,
      0);

  //    padSolid[i] = new G4Tubs("E07 pad", pad_in[i]*mm, pad_out[i]*mm,
  //           250.*mm, 0., angle[i]*deg);

  //////////////coil6RLUD
  G4double size_COIL6[4];
  //0:in
  //1:out
  //2:z
  //3:angle
  size_COIL6[0] = 50.0*mm;
  size_COIL6[1] = size_COIL1[1]*2+size_COIL6[0];
  //  size_COIL6[1] = 330.*mm;
  size_COIL6[2] = 280./2*mm;
  size_COIL6[3] = 90.*deg;

  G4double pos_COIL6LU[3];
  G4double pos_COIL6RU[3];
  G4double pos_COIL6LD[3];
  G4double pos_COIL6RD[3];
  //LU
  pos_COIL6LU[XCOORD] = pos_MFIELD[XCOORD] +size_MFIELD[XCOORD]-size_COIL6[0];
  pos_COIL6LU[YCOORD] = pos_COIL1[YCOORD]  -(size_COIL6[0]+size_COIL1[1]);
  pos_COIL6LU[ZCOORD] = pos_MFIELD[ZCOORD] - size_MFIELD[ZCOORD]-(size_COIL6[ZCOORD]+21.*mm);
  //RU
  pos_COIL6RU[XCOORD] = pos_MFIELD[XCOORD] -size_MFIELD[XCOORD]+size_COIL6[0];
  pos_COIL6RU[YCOORD] = pos_COIL1[YCOORD]  -(size_COIL6[0]+size_COIL1[1]);
  pos_COIL6RU[ZCOORD] = pos_MFIELD[ZCOORD] - size_MFIELD[ZCOORD]-(size_COIL6[ZCOORD]+21.*mm);
  //LD
  pos_COIL6LD[XCOORD] = pos_MFIELD[XCOORD]  +size_MFIELD[XCOORD]-size_COIL6[0];
  pos_COIL6LD[YCOORD] = -pos_COIL1[YCOORD] +(size_COIL6[0]+size_COIL1[1]);
  pos_COIL6LD[ZCOORD] = pos_MFIELD[ZCOORD]  - size_MFIELD[ZCOORD]-(size_COIL6[ZCOORD]+21.*mm);
  //RD
  pos_COIL6RD[XCOORD] = pos_MFIELD[XCOORD]  -size_MFIELD[XCOORD]+size_COIL6[0];
  pos_COIL6RD[YCOORD] = -pos_COIL1[YCOORD] +(size_COIL6[0]+size_COIL1[1]);
  pos_COIL6RD[ZCOORD] = pos_MFIELD[ZCOORD]  - size_MFIELD[ZCOORD]-(size_COIL6[ZCOORD]+21.*mm);


  G4Tubs* Coil6_tub = new G4Tubs("Coil6_tubs",
      size_COIL6[0],size_COIL6[1],size_COIL6[2],0.,size_COIL6[3]);
  G4LogicalVolume*  Coil6_log = new G4LogicalVolume(Coil6_tub, NistMan->FindOrBuildMaterial("G4_Cu"), "Coil6_log",0,0,0);
  Coil6_log->SetVisAttributes(new G4VisAttributes(G4Color::Red()));
  maxStep=0.00001*mm;
  //upGuard_L_log->SetUserLimits(new G4UserLimits(maxStep));

  //G4RotationMatrix *rotcoil6lu = new G4RotationMatrix();
  //G4RotationMatrix *rotcoil6ru = new G4RotationMatrix();
  //G4RotationMatrix *rotcoil6ld = new G4RotationMatrix();
  //G4RotationMatrix *rotcoil6rd = new G4RotationMatrix();
  //rotcoil6lu->rotateY(-fSpectrometerAngle);
  //rotcoil6ru->rotateY(-fSpectrometerAngle);
  //rotcoil6ld->rotateY(-fSpectrometerAngle);
  //rotcoil6rd->rotateY(-fSpectrometerAngle);


  G4RotationMatrix rotcoil6lu;
  G4RotationMatrix rotcoil6ru;
  G4RotationMatrix rotcoil6ld;
  G4RotationMatrix rotcoil6rd;
  rotcoil6lu.rotateY(-fSpectrometerAngle);
  rotcoil6ru.rotateY(-fSpectrometerAngle);
  rotcoil6ld.rotateY(-fSpectrometerAngle);
  rotcoil6rd.rotateY(-fSpectrometerAngle);


  rotcoil6lu.rotateZ(0.*deg);
  G4VPhysicalVolume* Coil6LU_phys
    = new G4PVPlacement(G4Transform3D(rotcoil6lu,
				      G4ThreeVector(pos_COIL6LU[XCOORD],
						    pos_COIL6LU[YCOORD],
						    pos_COIL6LU[ZCOORD]).rotateY(fSpectrometerAngle)),
			"Coil6LU_phys",
			Coil6_log,
			physMother,
			false,
			0);
  //rotcoil6ru.rotateZ(-90.*deg);
  rotcoil6ru.rotateZ(90.*deg);
  G4VPhysicalVolume* Coil6RU_phys
    = new G4PVPlacement(G4Transform3D(rotcoil6ru,
				      G4ThreeVector(pos_COIL6RU[XCOORD],
						    pos_COIL6RU[YCOORD],
						    pos_COIL6RU[ZCOORD]).rotateY(fSpectrometerAngle)),
			"Coil6RU_phys",
			Coil6_log,
			physMother,
			false,
			0);

  //rotcoil6ld.rotateZ(-180.*deg);
  rotcoil6ld.rotateZ(180.*deg);
  G4VPhysicalVolume* Coil6RD_phys
    = new G4PVPlacement(G4Transform3D(rotcoil6ld,
				      G4ThreeVector(pos_COIL6RD[XCOORD],
						    pos_COIL6RD[YCOORD],
						    pos_COIL6RD[ZCOORD]).rotateY(fSpectrometerAngle)),
			"Coil6RD_phys",
			Coil6_log,
			physMother,
			false,
			0);
  //rotcoil6rd.rotateZ(-270.*deg);
  rotcoil6rd.rotateZ(270.*deg);
  G4VPhysicalVolume* Coil6LD_phys
    = new G4PVPlacement(G4Transform3D(rotcoil6rd,
				      G4ThreeVector(pos_COIL6LD[XCOORD],
						    pos_COIL6LD[YCOORD],
						    pos_COIL6LD[ZCOORD]).rotateY(fSpectrometerAngle)),
			"Coil6LD_phys",
			Coil6_log,
			physMother,
			false,
			0);
  //  rotForwardSp->rotateZ(90.*deg);
  //  rotForwardSp->rotateZ(90.*deg);



  //////////////coil8RLUD
  G4double size_COIL8[4];
  //0:in
  //1:out
  //2:z
  //3:angle
  size_COIL8[0] = 50.0*mm;
  size_COIL8[1] = size_COIL5[1]*2.+size_COIL8[0];
  size_COIL8[2] = 280./2*mm;
  size_COIL8[3] = 90.*deg;

  G4double pos_COIL8LU[3];
  G4double pos_COIL8RU[3];
  G4double pos_COIL8LD[3];
  G4double pos_COIL8RD[3];
  //LU
  pos_COIL8LU[XCOORD] = pos_MFIELD[XCOORD] +size_MFIELD[XCOORD]-size_COIL8[0];
  pos_COIL8LU[YCOORD] = pos_COIL5[YCOORD]  -(size_COIL8[0]+size_COIL5[1]);
  pos_COIL8LU[ZCOORD] = pos_MFIELD[ZCOORD] + size_MFIELD[ZCOORD]+(size_COIL8[ZCOORD]+21.*mm);
  //RU
  pos_COIL8RU[XCOORD] = pos_MFIELD[XCOORD] -size_MFIELD[XCOORD]+size_COIL8[0];
  pos_COIL8RU[YCOORD] = pos_COIL5[YCOORD]  -(size_COIL8[0]+size_COIL5[1]);
  pos_COIL8RU[ZCOORD] = pos_MFIELD[ZCOORD] + size_MFIELD[ZCOORD]+(size_COIL8[ZCOORD]+21.*mm);
  //LD
  pos_COIL8LD[XCOORD] = pos_MFIELD[XCOORD]  +size_MFIELD[XCOORD]-size_COIL8[0];
  pos_COIL8LD[YCOORD] = -pos_COIL5[YCOORD] +(size_COIL8[0]+size_COIL5[1]);
  pos_COIL8LD[ZCOORD] = pos_MFIELD[ZCOORD]  + size_MFIELD[ZCOORD]+(size_COIL8[ZCOORD]+21.*mm);
  //RD
  pos_COIL8RD[XCOORD] = pos_MFIELD[XCOORD]  -size_MFIELD[XCOORD]+size_COIL8[0];
  pos_COIL8RD[YCOORD] = -pos_COIL5[YCOORD] +(size_COIL8[0]+size_COIL5[1]);
  pos_COIL8RD[ZCOORD] = pos_MFIELD[ZCOORD]  + size_MFIELD[ZCOORD]+(size_COIL8[ZCOORD]+21.*mm);


  G4Tubs* Coil8_tub = new G4Tubs("Coil8_tubs",
      size_COIL8[0],size_COIL8[1],size_COIL8[2],0.,size_COIL8[3]);
  G4LogicalVolume*  Coil8_log = new G4LogicalVolume(Coil8_tub, NistMan->FindOrBuildMaterial("G4_Cu"), "Coil8_log",0,0,0);
  Coil8_log->SetVisAttributes(new G4VisAttributes(G4Color::Red()));
  maxStep=0.00001*mm;
  //upGuard_L_log->SetUserLimits(new G4UserLimits(maxStep));

  G4RotationMatrix rotcoil8lu;
  G4RotationMatrix rotcoil8ru;
  G4RotationMatrix rotcoil8ld;
  G4RotationMatrix rotcoil8rd;
  rotcoil8lu.rotateY(-fSpectrometerAngle);
  rotcoil8ru.rotateY(-fSpectrometerAngle);
  rotcoil8ld.rotateY(-fSpectrometerAngle);
  rotcoil8rd.rotateY(-fSpectrometerAngle);

  rotcoil8lu.rotateZ(0.*deg);
  G4VPhysicalVolume* Coil8LU_phys
    = new G4PVPlacement(G4Transform3D(rotcoil8lu,
				      G4ThreeVector(pos_COIL8LU[XCOORD],
						    pos_COIL8LU[YCOORD],
						    pos_COIL8LU[ZCOORD]).rotateY(fSpectrometerAngle)),
			"Coil8LU_phys",
			Coil8_log,
			physMother,
			false,
			0);
  //rotcoil8ru.rotateZ(-90.*deg);
  rotcoil8ru.rotateZ(90.*deg);
  G4VPhysicalVolume* Coil8RU_phys
    = new G4PVPlacement(G4Transform3D(rotcoil8ru,
				      G4ThreeVector(pos_COIL8RU[XCOORD],
						    pos_COIL8RU[YCOORD],
						    pos_COIL8RU[ZCOORD]).rotateY(fSpectrometerAngle)),
			"Coil8RU_phys",
			Coil8_log,
			physMother,
			false,
			0);

  //rotcoil8ld.rotateZ(-180.*deg);
  //rotcoil8rd.rotateZ(270.*deg);
  rotcoil8rd.rotateZ(180.*deg);
  G4VPhysicalVolume* Coil8RD_phys
    = new G4PVPlacement(G4Transform3D(rotcoil8rd,
				      G4ThreeVector(pos_COIL8RD[XCOORD],
						    pos_COIL8RD[YCOORD],
						    pos_COIL8RD[ZCOORD]).rotateY(fSpectrometerAngle)),
			"Coil8RD_phys",
			Coil8_log,
			physMother,
			false,
			0);
  //rotcoil8rd.rotateZ(-270.*deg);
  //rotcoil8ld.rotateZ(180.*deg);
  rotcoil8ld.rotateZ(270.*deg);
  G4VPhysicalVolume* Coil8LD_phys
    = new G4PVPlacement(G4Transform3D(rotcoil8ld,
				      G4ThreeVector(pos_COIL8LD[XCOORD],
						    pos_COIL8LD[YCOORD],
						    pos_COIL8LD[ZCOORD]).rotateY(fSpectrometerAngle)),
			"Coil8LD_phys",
			Coil8_log,
			
			physMother,
			false,
			0);
  //  rotForwardSp->rotateZ(90.*deg);
  //  rotForwardSp->rotateZ(90.*deg);
  


  //////////////coil7RLUD
  G4double size_COIL7[4];
  //0:in
  //1:out
  //2:z
  //3:angle
  size_COIL7[0] = 50.0*mm;
  size_COIL7[1] = size_COIL4[1]*2.+size_COIL7[0];
  size_COIL7[2] = size_COIL4[0];
  size_COIL7[3] = 90.*deg;

  G4double pos_COIL7ULU[3];
  G4double pos_COIL7URU[3];
  G4double pos_COIL7ULD[3];
  G4double pos_COIL7URD[3];
  G4double pos_COIL7DLU[3];
  G4double pos_COIL7DRU[3];
  G4double pos_COIL7DLD[3];
  G4double pos_COIL7DRD[3];
  //ULU
  //  pos_COIL7ULU[XCOORD] = pos_MFIELD[XCOORD] +size_MFIELD[XCOORD]+size_COIL7[0];
  pos_COIL7ULU[XCOORD] = pos_COIL4L[XCOORD];
  pos_COIL7ULU[YCOORD] = pos_COIL4L[YCOORD] +(size_COIL7[0]+size_COIL4[1]);
  pos_COIL7ULU[ZCOORD] = pos_MFIELD[ZCOORD] - size_COIL4[ZCOORD];
  //URU
  pos_COIL7URU[XCOORD] = pos_COIL4R[XCOORD];
  pos_COIL7URU[YCOORD] = pos_COIL4R[YCOORD] +(size_COIL7[0]+size_COIL4[1]);
  pos_COIL7URU[ZCOORD] = pos_MFIELD[ZCOORD] - size_COIL4[ZCOORD];
  //ULD
  pos_COIL7ULD[XCOORD] = pos_COIL4L[XCOORD];
  pos_COIL7ULD[YCOORD] = -pos_COIL4L[YCOORD] -(size_COIL7[0]+size_COIL4[1]);
  pos_COIL7ULD[ZCOORD] = pos_MFIELD[ZCOORD]  - size_COIL4[ZCOORD];
  //URD
  pos_COIL7URD[XCOORD] = pos_COIL4R[XCOORD];
  pos_COIL7URD[YCOORD] = -pos_COIL4R[YCOORD] -(size_COIL7[0]+size_COIL4[1]);
  pos_COIL7URD[ZCOORD] = pos_MFIELD[ZCOORD] - size_COIL4[ZCOORD];


  //DLU
  //  pos_COIL7ULU[XCOORD] = pos_MFIELD[XCOORD] +size_MFIELD[XCOORD]+size_COIL7[0];
  pos_COIL7DLU[XCOORD] = pos_COIL4L[XCOORD];
  pos_COIL7DLU[YCOORD] = pos_COIL4L[YCOORD] +(size_COIL7[0]+size_COIL4[1]);
  pos_COIL7DLU[ZCOORD] = pos_MFIELD[ZCOORD] + size_COIL4[ZCOORD];
  //DRU
  pos_COIL7DRU[XCOORD] = pos_COIL4R[XCOORD];
  pos_COIL7DRU[YCOORD] = pos_COIL4R[YCOORD] +(size_COIL7[0]+size_COIL4[1]);
  pos_COIL7DRU[ZCOORD] = pos_MFIELD[ZCOORD] + size_COIL4[ZCOORD];
  //DLD
  pos_COIL7DLD[XCOORD] = pos_COIL4L[XCOORD];
  pos_COIL7DLD[YCOORD] = -pos_COIL4L[YCOORD] -(size_COIL7[0]+size_COIL4[1]);
  pos_COIL7DLD[ZCOORD] = pos_MFIELD[ZCOORD]  + size_COIL4[ZCOORD];
  //DRD
  pos_COIL7DRD[XCOORD] = pos_COIL4R[XCOORD];
  pos_COIL7DRD[YCOORD] = -pos_COIL4R[YCOORD] -(size_COIL7[0]+size_COIL4[1]);
  pos_COIL7DRD[ZCOORD] = pos_MFIELD[ZCOORD] + size_COIL4[ZCOORD];


  G4Tubs* Coil7_tub = new G4Tubs("Coil7_tubs",
      size_COIL7[0],size_COIL7[1],size_COIL7[2],0.,size_COIL7[3]);
  G4LogicalVolume*  Coil7_log = new G4LogicalVolume(Coil7_tub, NistMan->FindOrBuildMaterial("G4_Cu"), "Coil7_log",0,0,0);
  Coil7_log->SetVisAttributes(new G4VisAttributes(G4Color::Red()));
  maxStep=0.00001*mm;
  //upGuard_L_log->SetUserLimits(new G4UserLimits(maxStep));

  G4RotationMatrix rotcoil7ulu;
  G4RotationMatrix rotcoil7uru;
  G4RotationMatrix rotcoil7uld;
  G4RotationMatrix rotcoil7urd;
  G4RotationMatrix rotcoil7dlu;
  G4RotationMatrix rotcoil7dru;
  G4RotationMatrix rotcoil7dld;
  G4RotationMatrix rotcoil7drd;
  rotcoil7ulu.rotateY(-fSpectrometerAngle);
  rotcoil7uru.rotateY(-fSpectrometerAngle);
  rotcoil7uld.rotateY(-fSpectrometerAngle);
  rotcoil7urd.rotateY(-fSpectrometerAngle);
  rotcoil7dlu.rotateY(-fSpectrometerAngle);
  rotcoil7dru.rotateY(-fSpectrometerAngle);
  rotcoil7dld.rotateY(-fSpectrometerAngle);
  rotcoil7drd.rotateY(-fSpectrometerAngle);

  
  rotcoil7ulu.rotateZ(0.*deg);
  rotcoil7ulu.rotateX(180.*deg);
  rotcoil7ulu.rotateY(90.*deg);
  G4VPhysicalVolume* Coil7ULU_phys
    = new G4PVPlacement(G4Transform3D(rotcoil7ulu,
				      G4ThreeVector(pos_COIL7ULU[XCOORD],
						    pos_COIL7ULU[YCOORD],
						    pos_COIL7ULU[ZCOORD]).rotateY(fSpectrometerAngle)),
			"Coil7ULU_phys",
			Coil7_log,
			physMother,
			false,
			0);


  rotcoil7uru.rotateZ(0.*deg);
  rotcoil7uru.rotateX(180.*deg);
  rotcoil7uru.rotateY(90.*deg);
  G4VPhysicalVolume* Coil7URU_phys
    = new G4PVPlacement(G4Transform3D(rotcoil7uru,
				      G4ThreeVector(pos_COIL7URU[XCOORD],
						    pos_COIL7URU[YCOORD],
						    pos_COIL7URU[ZCOORD]).rotateY(fSpectrometerAngle)),
			"Coil7URU_phys",
			Coil7_log,
			physMother,
			false,
			0);

  rotcoil7urd.rotateZ(0.*deg);
  //rotcoil7urd.rotateX(90.*deg);
  rotcoil7urd.rotateX(0.*deg);
  rotcoil7urd.rotateY(90.*deg);
  G4VPhysicalVolume* Coil7URD_phys
    = new G4PVPlacement(G4Transform3D(rotcoil7urd,
				      G4ThreeVector(pos_COIL7URD[XCOORD],
						    pos_COIL7URD[YCOORD],
						    pos_COIL7URD[ZCOORD]).rotateY(fSpectrometerAngle)),
			"Coil7URD_phys",
			Coil7_log,
			physMother,
			false,
			0);

  rotcoil7uld.rotateZ(0.*deg);
  //rotcoil7uld.rotateX(90.*deg);
  rotcoil7uld.rotateX(0.*deg);
  rotcoil7uld.rotateY(90.*deg);
  G4VPhysicalVolume* Coil7ULD_phys
    = new G4PVPlacement(G4Transform3D(rotcoil7uld,
				      G4ThreeVector(pos_COIL7ULD[XCOORD],
						    pos_COIL7ULD[YCOORD],
						    pos_COIL7ULD[ZCOORD]).rotateY(fSpectrometerAngle)),
			"Coil7ULD_phys",
			Coil7_log,
			physMother,
			false,
			0);


  ///down
  rotcoil7dlu.rotateX(-90.*deg);
  rotcoil7dlu.rotateY(180.*deg);
  rotcoil7dlu.rotateZ(90.*deg);
  G4VPhysicalVolume* Coil7DLU_phys
    = new G4PVPlacement(G4Transform3D(rotcoil7dlu,
				      G4ThreeVector(pos_COIL7DLU[XCOORD],
						    pos_COIL7DLU[YCOORD],
						    pos_COIL7DLU[ZCOORD]).rotateY(fSpectrometerAngle)),
			"Coil7DLU_phys",
			Coil7_log,
			physMother,
			false,
			0);



  rotcoil7dru.rotateX(-90.*deg);
  rotcoil7dru.rotateY(180.*deg);
  rotcoil7dru.rotateZ(90.*deg);
  G4VPhysicalVolume* Coil7DRU_phys
    = new G4PVPlacement(G4Transform3D(rotcoil7dru,
				      G4ThreeVector(pos_COIL7DRU[XCOORD],
						    pos_COIL7DRU[YCOORD],
						    pos_COIL7DRU[ZCOORD]).rotateY(fSpectrometerAngle)),
			"Coil7DRU_phys",
			Coil7_log,
			physMother,
			false,
			0);

  rotcoil7drd.rotateZ(0.*deg);
  rotcoil7drd.rotateX(0.*deg);
  rotcoil7drd.rotateY(-90.*deg);
  G4VPhysicalVolume* Coil7DRD_phys
    = new G4PVPlacement(G4Transform3D(rotcoil7drd,
				      G4ThreeVector(pos_COIL7DRD[XCOORD],
						    pos_COIL7DRD[YCOORD],
						    pos_COIL7DRD[ZCOORD]).rotateY(fSpectrometerAngle)),
			"Coil7DRD_phys",
			Coil7_log,
			physMother,
			false,
			0);

  rotcoil7dld.rotateZ(0.*deg);
  rotcoil7dld.rotateX(0.*deg);
  rotcoil7dld.rotateY(-90.*deg);
  G4VPhysicalVolume* Coil7DLD_phys
    = new G4PVPlacement(G4Transform3D(rotcoil7dld,
				      G4ThreeVector(pos_COIL7DLD[XCOORD],
						    pos_COIL7DLD[YCOORD],
						    pos_COIL7DLD[ZCOORD]).rotateY(fSpectrometerAngle)),
			"Coil7DLD_phys",
			Coil7_log,
			physMother,
			false,
			0);
  
  ///coil2
  //////////////coil2
  G4Box* Coil2_box = new G4Box("Coil2_box",
      size_COIL2[XCOORD],size_COIL2[YCOORD],size_COIL2[ZCOORD]);
  G4LogicalVolume*  Coil2_log = new G4LogicalVolume(Coil2_box, NistMan->FindOrBuildMaterial("G4_Cu"), "Coil2_log",0,0,0);
  Coil2_log->SetVisAttributes(new G4VisAttributes(G4Color::Red()));
  maxStep=0.00001*mm;
  //upGuard_L_log->SetUserLimits(new G4UserLimits(maxStep));

  G4VPhysicalVolume* Coil2UL_phys = new G4PVPlacement(G4Transform3D(rotForwardSp,
								    G4ThreeVector(pos_COIL2L[XCOORD],
										  pos_COIL2L[YCOORD],
										  pos_COIL2L[ZCOORD]).rotateY(fSpectrometerAngle)),
						      "Coil2UL_phys",
						      Coil2_log,
						      physMother,
						      false,
						      0);
  G4VPhysicalVolume* Coil2UR_phys = new G4PVPlacement(G4Transform3D(rotForwardSp,
								    G4ThreeVector(pos_COIL2R[XCOORD],
										  pos_COIL2R[YCOORD],
										  pos_COIL2R[ZCOORD]).rotateY(fSpectrometerAngle)),
						      "Coil2UR_phys",
						      Coil2_log,
						      physMother,
						      false,
						      0);


  G4VPhysicalVolume* Coil2DL_phys = new G4PVPlacement(G4Transform3D(rotForwardSp,
								    G4ThreeVector(pos_COIL2L[XCOORD],
										  -pos_COIL2L[YCOORD],
										  pos_COIL2L[ZCOORD]).rotateY(fSpectrometerAngle)),
						      "Coil2DL_phys",
						      Coil2_log,
						      physMother,
						      false,
						      0);
  G4VPhysicalVolume* Coil2DR_phys = new G4PVPlacement(G4Transform3D(rotForwardSp,
								    G4ThreeVector(pos_COIL2R[XCOORD],
										  -pos_COIL2R[YCOORD],
										  pos_COIL2R[ZCOORD]).rotateY(fSpectrometerAngle)),
						      "Coil2DR_phys",
						      Coil2_log,
						      physMother,
						      false,
						      0);


  ///coil3
  //////////////coil3
  G4Box* Coil3_box = new G4Box("Coil3_box",
      size_COIL3[XCOORD],size_COIL3[YCOORD],size_COIL3[ZCOORD]);
  G4LogicalVolume*  Coil3_log = new G4LogicalVolume(Coil3_box, NistMan->FindOrBuildMaterial("G4_Cu"), "Coil3_log",0,0,0);
  Coil3_log->SetVisAttributes(new G4VisAttributes(G4Color::Red()));
  maxStep=0.00001*mm;
  //upGuard_L_log->SetUserLimits(new G4UserLimits(maxStep));

  G4VPhysicalVolume* Coil3UL_phys = new G4PVPlacement(G4Transform3D(rotForwardSp,
								    G4ThreeVector(pos_COIL3L[XCOORD],
										  pos_COIL3L[YCOORD],
										  pos_COIL3L[ZCOORD]).rotateY(fSpectrometerAngle)),
						      "Coil3UL_phys",
						      Coil3_log,
						      physMother,
						      false,
						      0);
  G4VPhysicalVolume* Coil3UR_phys = new G4PVPlacement(G4Transform3D(rotForwardSp,
								    G4ThreeVector(pos_COIL3R[XCOORD],
										  pos_COIL3R[YCOORD],
										  pos_COIL3R[ZCOORD]).rotateY(fSpectrometerAngle)),
						      "Coil3UR_phys",
						      Coil3_log,
						      physMother,
						      false,
						      0);
  
  
  G4VPhysicalVolume* Coil3DL_phys = new G4PVPlacement(G4Transform3D(rotForwardSp,
								    G4ThreeVector(pos_COIL3L[XCOORD],
										  -pos_COIL3L[YCOORD],
										  pos_COIL3L[ZCOORD]).rotateY(fSpectrometerAngle)),
						      "Coil3DL_phys",
						      Coil3_log,
						      physMother,
						      false,
						      0);

  G4VPhysicalVolume* Coil3DR_phys = new G4PVPlacement(G4Transform3D(rotForwardSp,
								    G4ThreeVector(pos_COIL3R[XCOORD],
										  -pos_COIL3R[YCOORD],
										  pos_COIL3R[ZCOORD]).rotateY(fSpectrometerAngle)),
						      "Coil3DR_phys",
						      Coil3_log,
						      physMother,
						      false,
						      0);




  //-------------------- Upstream End Guard
      
    G4Box* upGuard_UD_box = new G4Box("upGuard_UD_box",
        size_UGUARD_UD[XCOORD],size_UGUARD_UD[YCOORD],size_UGUARD_UD[ZCOORD]);
    G4LogicalVolume*  upGuard_U_log = new G4LogicalVolume(upGuard_UD_box, NistMan->FindOrBuildMaterial("G4_Fe"), "upGuard_U_log",0,0,0);
    upGuard_U_log->SetVisAttributes(new G4VisAttributes(G4Color::Blue()));
    maxStep=0.00001*mm;
    //upGuard_U_log->SetUserLimits(new G4UserLimits(maxStep));
    
    G4VPhysicalVolume*  upGuard_U_phys = new G4PVPlacement(G4Transform3D(rotForwardSp,
									 G4ThreeVector(pos_UGUARD_U[XCOORD],
										       pos_UGUARD_U[YCOORD],
										       pos_UGUARD_U[ZCOORD]).rotateY(fSpectrometerAngle)),
							   "upGuard_U_phys",
							   upGuard_U_log,
							   physMother,
							   false,
							   0);
    
    G4LogicalVolume*  upGuard_D_log = new G4LogicalVolume(upGuard_UD_box, NistMan->FindOrBuildMaterial("G4_Fe"), "upGuard_D_log",0,0,0);
    upGuard_D_log->SetVisAttributes(new G4VisAttributes(G4Color::Blue()));
    maxStep=0.00001*mm;
    //upGuard_D_log->SetUserLimits(new G4UserLimits(maxStep));
    
    G4VPhysicalVolume*  upGuard_D_phys = new G4PVPlacement(G4Transform3D(rotForwardSp,
									 G4ThreeVector(pos_UGUARD_D[XCOORD],
										       pos_UGUARD_D[YCOORD],
										       pos_UGUARD_D[ZCOORD]).rotateY(fSpectrometerAngle)),
							   "upGuard_D_phys",
							   upGuard_D_log,
							   physMother,
							   false,
							   0);
    
    G4Box* upGuard_LR_box = new G4Box("upGuard_LR_box",
        size_UGUARD_LR[XCOORD],size_UGUARD_LR[YCOORD],size_UGUARD_LR[ZCOORD]);
    G4LogicalVolume*  upGuard_L_log = new G4LogicalVolume(upGuard_LR_box, NistMan->FindOrBuildMaterial("G4_Fe"), "upGuard_L_log",0,0,0);
    upGuard_L_log->SetVisAttributes(new G4VisAttributes(G4Color::Blue()));
    maxStep=0.00001*mm;
    //upGuard_L_log->SetUserLimits(new G4UserLimits(maxStep));
    G4VPhysicalVolume* upGuard_L_phys = new G4PVPlacement(G4Transform3D(rotForwardSp,
									G4ThreeVector(pos_UGUARD_L[XCOORD],
										      pos_UGUARD_L[YCOORD],
										      pos_UGUARD_L[ZCOORD]).rotateY(fSpectrometerAngle)),
							  "upGuard_L_phys",
							  upGuard_L_log,
							  physMother,
							  false,
							  0);
    
    /*
    G4LogicalVolume*  upGuard_R_log = new G4LogicalVolume(upGuard_LR_box, NistMan->FindOrBuildMaterial("G4_Fe"), "upGuard_R_log",0,0,0);
    upGuard_R_log->SetVisAttributes(new G4VisAttributes(G4Color::Blue()));
    maxStep=0.00001*mm;
    //upGuard_R_log->SetUserLimits(new G4UserLimits(maxStep));
    G4VPhysicalVolume* upGuard_R_phys = new G4PVPlacement(G4Transform3D(rotForwardSp,
									G4ThreeVector(pos_UGUARD_R[XCOORD],
										      pos_UGUARD_R[YCOORD],
										      pos_UGUARD_R[ZCOORD]).rotateY(fSpectrometerAngle)),
							  "upGuard_R_phys",							  
							  upGuard_R_log,
							  physMother,
							  false,
							  0);
    */

    G4Box* upGuard_R1_box =
      new G4Box("upGuard_R1_box",
		size_UGUARD_R1[XCOORD],size_UGUARD_R1[YCOORD],size_UGUARD_R1[ZCOORD]);
    
    G4Box* upGuard_R2_box =
      new G4Box("upGuard_R2_box",
		size_UGUARD_R2[XCOORD],size_UGUARD_R2[YCOORD],size_UGUARD_R2[ZCOORD]);
    
    
    G4LogicalVolume* upGuard_R1_log = new G4LogicalVolume(upGuard_R1_box,
							  NistMan->FindOrBuildMaterial("G4_Fe"),
							  "upGuard_R1_log",0,0,0);
    upGuard_R1_log->SetVisAttributes(new G4VisAttributes(G4Color::Blue()));
                           
    G4VPhysicalVolume* upGuard_R1_phys =
      new G4PVPlacement(G4Transform3D(rotForwardSp,
				      G4ThreeVector(pos_UGUARD_R1[XCOORD],
						    pos_UGUARD_R1[YCOORD],
						    pos_UGUARD_R1[ZCOORD]).rotateY(fSpectrometerAngle)),
			"upGuard_R1_phys",
			upGuard_R1_log,
			physMother,
			false,
			0);
    
    G4LogicalVolume* upGuard_R2U_log = new G4LogicalVolume(upGuard_R2_box,
							  NistMan->FindOrBuildMaterial("G4_Fe"),
							   "upGuard_R2U_log",0,0,0);
    upGuard_R2U_log->SetVisAttributes(new G4VisAttributes(G4Color::Blue()));
    G4VPhysicalVolume* upGuard_R2U_phys =
      new G4PVPlacement(0,
			G4ThreeVector(pos_UGUARD_R2U[XCOORD],pos_UGUARD_R2U[YCOORD],pos_UGUARD_R2U[ZCOORD]),
			"upGuard_R2U_phys",
			upGuard_R2U_log,
			physMother,
			false,
			0);
    
    G4LogicalVolume* upGuard_R2D_log = new G4LogicalVolume(upGuard_R2_box,
							   NistMan->FindOrBuildMaterial("G4_Fe"),
							   "upGuard_R2D_log",0,0,0);
    upGuard_R2D_log->SetVisAttributes(new G4VisAttributes(G4Color::Blue()));
    maxStep=0.00001*mm;
    G4VPhysicalVolume* upGuard_R2D_phys =
      new G4PVPlacement(0,
			G4ThreeVector(pos_UGUARD_R2D[XCOORD],pos_UGUARD_R2D[YCOORD],pos_UGUARD_R2D[ZCOORD]),
			"upGuard_R2D_phys",
			upGuard_R2D_log,
			physMother,
			false,
			0);


  //-------------------- Yoke
  
  G4Box* Yoke_UD_box = new G4Box("Yoke_UD_box",
      size_YOKE_UD[XCOORD],size_YOKE_UD[YCOORD],size_YOKE_UD[ZCOORD]);
  G4LogicalVolume*  Yoke_U_log = new G4LogicalVolume(Yoke_UD_box, NistMan->FindOrBuildMaterial("G4_Fe"), "Yoke_U_log",0,0,0);

  Yoke_U_log->SetVisAttributes(new G4VisAttributes(G4Color::Blue()));
  maxStep=0.00001*mm;
  //Yoke_U_log->SetUserLimits(new G4UserLimits(maxStep));
  G4VPhysicalVolume* Yoke_U_phys = new G4PVPlacement(G4Transform3D(rotForwardSp,
								   G4ThreeVector(pos_YOKE_U[XCOORD],
										 pos_YOKE_U[YCOORD],
										 pos_YOKE_U[ZCOORD]).rotateY(fSpectrometerAngle)),
						     "Yoke_U_phys",
						     Yoke_U_log,
						     physMother,
						     false,
						     0);

  G4LogicalVolume*  Yoke_D_log = new G4LogicalVolume(Yoke_UD_box, NistMan->FindOrBuildMaterial("G4_Fe"), "Yoke_D_log",0,0,0);
  Yoke_D_log->SetVisAttributes(new G4VisAttributes(G4Color::Blue()));
  maxStep=0.00001*mm;
  //Yoke_D_log->SetUserLimits(new G4UserLimits(maxStep));
  G4VPhysicalVolume* Yoke_D_phys = new G4PVPlacement(G4Transform3D(rotForwardSp,
								   G4ThreeVector(pos_YOKE_D[XCOORD],
										 pos_YOKE_D[YCOORD],
										 pos_YOKE_D[ZCOORD]).rotateY(fSpectrometerAngle)),
						     "Yoke_D_phys",
						     Yoke_D_log,
						     physMother,
						     false,
						     0);

  G4Box* Yoke_LR_box = new G4Box("Yoke_LR_box",
      size_YOKE_LR[XCOORD],size_YOKE_LR[YCOORD],size_YOKE_LR[ZCOORD]);
  G4LogicalVolume*  Yoke_L_log = new G4LogicalVolume(Yoke_LR_box, NistMan->FindOrBuildMaterial("G4_Fe"), "Yoke_L_log",0,0,0);
  Yoke_L_log->SetVisAttributes(new G4VisAttributes(G4Color::Blue()));
  maxStep=0.00001*mm;
  //Yoke_L_log->SetUserLimits(new G4UserLimits(maxStep));
  G4VPhysicalVolume* Yoke_L_phys = new G4PVPlacement(G4Transform3D(rotForwardSp,
								   G4ThreeVector(pos_YOKE_L[XCOORD],
										 pos_YOKE_L[YCOORD],
										 pos_YOKE_L[ZCOORD]).rotateY(fSpectrometerAngle)),
						     "Yoke_L_phys",
						     Yoke_L_log,
						     physMother,
						     false,
						     0);

  G4LogicalVolume*  Yoke_R_log = new G4LogicalVolume(Yoke_LR_box, NistMan->FindOrBuildMaterial("G4_Fe"), "Yoke_R_log",0,0,0);
  Yoke_R_log->SetVisAttributes(new G4VisAttributes(G4Color::Blue()));
  maxStep=0.00001*mm;
  //Yoke_R_log->SetUserLimits(new G4UserLimits(maxStep));
  G4VPhysicalVolume* Yoke_R_phys = new G4PVPlacement(G4Transform3D(rotForwardSp,
								   G4ThreeVector(pos_YOKE_R[XCOORD],
										 pos_YOKE_R[YCOORD],
										 pos_YOKE_R[ZCOORD]).rotateY(fSpectrometerAngle)),
						     "Yoke_R_phys",
						     Yoke_R_log,
						     physMother,
						     false,
						     0);
  
  //E42 Yoke addition
  G4Box* Yoke_LR_GapSpace_box = new G4Box("Yoke_LR_GapSpace_box",
                                          size_YOKE_LR_GapSpace[XCOORD],
                                          size_YOKE_LR_GapSpace[YCOORD],
                                          size_YOKE_LR_GapSpace[ZCOORD]);
  G4LogicalVolume* Yoke_L_GapSpace_log = new G4LogicalVolume(Yoke_LR_GapSpace_box,
                                                             NistMan->FindOrBuildMaterial("G4_Fe"),
                                                             "Yoke_L_GapSpace_log",0,0,0);

  Yoke_L_GapSpace_log->SetVisAttributes(new G4VisAttributes(G4Color::Blue()));
  G4VPhysicalVolume* Yoke_L_GapSpace_phys =
    new G4PVPlacement(0,
                      G4ThreeVector(pos_YOKE_L_GapSpace[XCOORD],pos_YOKE_L_GapSpace[YCOORD],pos_YOKE_L_GapSpace[ZCOORD]),
                      "Yoke_L_GapSpace_phys",
                      Yoke_L_GapSpace_log,
                      physMother,
                      false,
                      0);
  
  G4LogicalVolume* Yoke_R_GapSpace_log = new G4LogicalVolume(Yoke_LR_GapSpace_box,
                                                             NistMan->FindOrBuildMaterial("G4_Fe"),
                                                             "Yoke_R_GapSpace_log",0,0,0);

  Yoke_R_GapSpace_log->SetVisAttributes(new G4VisAttributes(G4Color::Blue()));
  G4VPhysicalVolume* Yoke_R_GapSpace_phys =
    new G4PVPlacement(0,
                      G4ThreeVector(pos_YOKE_R_GapSpace[XCOORD],pos_YOKE_R_GapSpace[YCOORD],pos_YOKE_R_GapSpace[ZCOORD]),
                      "Yoke_R_GapSpace_phys",
                      Yoke_R_GapSpace_log,
                      physMother,
                      false,
                      0);



    
  //-------------------- Downstream End Guard
  
  G4Box* downGuard_UD_box = new G4Box("downGuard_UD_box",
      size_DGUARD_UD[XCOORD],size_DGUARD_UD[YCOORD],size_DGUARD_UD[ZCOORD]);
  G4LogicalVolume*  downGuard_U_log = new G4LogicalVolume(downGuard_UD_box, NistMan->FindOrBuildMaterial("G4_Fe"), "downGuard_U_log",0,0,0);
  downGuard_U_log->SetVisAttributes(new G4VisAttributes(G4Color::Blue()));
  maxStep=0.00001*mm;
  //downGuard_U_log->SetUserLimits(new G4UserLimits(maxStep));
  G4VPhysicalVolume* downGuard_U_phys = new G4PVPlacement(G4Transform3D(rotForwardSp,
									G4ThreeVector(pos_DGUARD_U[XCOORD],
										      pos_DGUARD_U[YCOORD],
										      pos_DGUARD_U[ZCOORD]).rotateY(fSpectrometerAngle)),
							  "downGuard_U_phys",
							  downGuard_U_log,
							  physMother,
							  false,
							  0);

  G4LogicalVolume*  downGuard_D_log = new G4LogicalVolume(downGuard_UD_box, NistMan->FindOrBuildMaterial("G4_Fe"), "downGuard_D_log",0,0,0);
  downGuard_D_log->SetVisAttributes(new G4VisAttributes(G4Color::Blue()));
  maxStep=0.00001*mm;
  //downGuard_D_log->SetUserLimits(new G4UserLimits(maxStep));
  G4VPhysicalVolume* downGuard_D_phys = new G4PVPlacement(G4Transform3D(rotForwardSp,
									G4ThreeVector(pos_DGUARD_D[XCOORD],
										      pos_DGUARD_D[YCOORD],
										      pos_DGUARD_D[ZCOORD]).rotateY(fSpectrometerAngle)),
							  "downGuard_D_phys",
							  downGuard_D_log,
							  physMother,
							  false,
							  0);

  G4Box* downGuard_LR_box = new G4Box("downGuard_LR_box",
      size_DGUARD_LR[XCOORD],size_DGUARD_LR[YCOORD],size_DGUARD_LR[ZCOORD]);
  G4LogicalVolume*  downGuard_L_log = new G4LogicalVolume(downGuard_LR_box, NistMan->FindOrBuildMaterial("G4_Fe"), "downGuard_L_log",0,0,0);
  downGuard_L_log->SetVisAttributes(new G4VisAttributes(G4Color::Blue()));
  maxStep=0.00001*mm;
  //downGuard_L_log->SetUserLimits(new G4UserLimits(maxStep));
  G4VPhysicalVolume* downGuard_L_phys = new G4PVPlacement(G4Transform3D(rotForwardSp,
									G4ThreeVector(pos_DGUARD_L[XCOORD],
										      pos_DGUARD_L[YCOORD],
										      pos_DGUARD_L[ZCOORD]).rotateY(fSpectrometerAngle)),
      "downGuard_L_phys",
      downGuard_L_log,
      physMother,
      false,
      0);

  G4LogicalVolume*  downGuard_R_log = new G4LogicalVolume(downGuard_LR_box, NistMan->FindOrBuildMaterial("G4_Fe"), "downGuard_R_log",0,0,0);
  downGuard_R_log->SetVisAttributes(new G4VisAttributes(G4Color::Blue()));
  maxStep=0.00001*mm;
  //downGuard_R_log->SetUserLimits(new G4UserLimits(maxStep));
  G4VPhysicalVolume*  downGuard_R_phys = new G4PVPlacement(G4Transform3D(rotForwardSp,
									 G4ThreeVector(pos_DGUARD_R[XCOORD],
										       pos_DGUARD_R[YCOORD],
										       pos_DGUARD_R[ZCOORD]).rotateY(fSpectrometerAngle)),
      "downGuard_R_phys",
      downGuard_R_log,
      physMother,
      false,
      0);
  

  
  //Detectors inside KURAMA//
  G4ThreeVector DC1Pos ( -pos.x(), -pos.y(), - 635.4675*mm);
  G4RotationMatrix DC1Rot;

  G4ThreeVector CH1Pos ( 37.- pos.x(), 1. - pos.y(), - 387.*mm);
  G4RotationMatrix CH1Rot;

  G4ThreeVector DummyDetectorPos2 ( 0., 0., -820.0 + 100.0 + 0.01/2.0);// End HC 
  G4RotationMatrix DummyDetectorRot2;

  MakeDummyDetector2( physMother, DummyDetectorPos2, DummyDetectorRot2 );
  //MakeDC1( physMother, DC1Pos, DC1Rot);
  MakeCH1( physMother, CH1Pos, CH1Rot);

}


void DetectorConstruction::MakeDummyDetector4(G4VPhysicalVolume *pMother, G4ThreeVector &pos, G4RotationMatrix &rot)
{
  /// Dummy Detector Physics World ///
  /// Physics World ///
  // param //
  const G4double X_PW = 1000.0/2.0*mm;
  const G4double Y_PW=  800.0/2.0*mm;
  const G4double Z_PW  = 0.01/2.0*mm; 
  std::string fullNamePW = "Dummy_Detector_PhysicsWorld";
  // //
  
  G4Box* solidDetectorPW = new G4Box(G4String(fullNamePW.c_str()),
				     X_PW,
				     Y_PW,
				     Z_PW
				     );

  G4LogicalVolume* logicDetectorPW = new G4LogicalVolume(solidDetectorPW,
							 mat_air,
							 G4String(fullNamePW.c_str()));

  logicDetectorPW->SetVisAttributes(new G4VisAttributes(G4Color::Yellow()));
  
  G4RotationMatrix rot_PW = rot;
  G4ThreeVector TVPW = rot*pos;
  G4VPhysicalVolume *physMother = new G4PVPlacement(G4Transform3D(rot_PW, TVPW),
						    G4String(fullNamePW.c_str()),
						    logicDetectorPW,
						    pMother,
						    false,
						    0
						    );

  // Dummy Detector //                                                                                        
  G4Box* dummy_box = new G4Box("dummy_box", X_PW, Y_PW, Z_PW);
  G4LogicalVolume*  dummy_log = new G4LogicalVolume(dummy_box, mat_air, "dummy_log", 0, 0, 0);
  dummy_log->SetVisAttributes(G4Color::Yellow());

  G4RotationMatrix dummy_rot;
  G4VPhysicalVolume* dummy_phys = new G4PVPlacement(G4Transform3D(dummy_rot,
								  G4ThreeVector(0., 0., 0.)),
						    "dummy_phys",
						    dummy_log,
						    physMother,
						    false,
						    103);

  // sensitive detector //
  G4SDManager *SDMan = G4SDManager::GetSDMpointer();
  G4VSensitiveDetector *dcSD = SDMan -> FindSensitiveDetector("/DC");

  if(!dcSD)
    {
      DCSD* dcSD = new DCSD("/DC");
      SDMan -> AddNewDetector( dcSD );
    }

  G4VSensitiveDetector *dcSD1 = SDMan -> FindSensitiveDetector("/DC");
  dummy_log->SetSensitiveDetector( dcSD1 );

}


void DetectorConstruction::MakeDummyDetector5(G4VPhysicalVolume *pMother, G4ThreeVector &pos, G4RotationMatrix &rot)
{
  /// Dummy Detector Physics World ///
  /// Physics World ///
  // param //
  const G4double X_PW = 2200.0/2.0*mm;
  const G4double Y_PW=  1500.0/2.0*mm;
  const G4double Z_PW  = 0.01/2.0*mm; 
  std::string fullNamePW = "Dummy_Detector_PhysicsWorld";
  // //
  
  G4Box* solidDetectorPW = new G4Box(G4String(fullNamePW.c_str()),
				     X_PW,
				     Y_PW,
				     Z_PW
				     );

  G4LogicalVolume* logicDetectorPW = new G4LogicalVolume(solidDetectorPW,
							 mat_air,
							 G4String(fullNamePW.c_str()));

  logicDetectorPW->SetVisAttributes(new G4VisAttributes(G4Color::Yellow()));
  
  G4RotationMatrix rot_PW = rot;
  G4ThreeVector TVPW = rot*pos;
  G4VPhysicalVolume *physMother = new G4PVPlacement(G4Transform3D(rot_PW, TVPW),
						    G4String(fullNamePW.c_str()),
						    logicDetectorPW,
						    pMother,
						    false,
						    0
						    );

  // Dummy Detector //                                                                                        
  G4Box* dummy_box = new G4Box("dummy_box", X_PW, Y_PW, Z_PW);
  G4LogicalVolume*  dummy_log = new G4LogicalVolume(dummy_box, mat_air, "dummy_log", 0, 0, 0);
  dummy_log->SetVisAttributes(G4Color::Yellow());

  G4RotationMatrix dummy_rot;
  G4VPhysicalVolume* dummy_phys = new G4PVPlacement(G4Transform3D(dummy_rot,
								  G4ThreeVector(0., 0., 0.)),
						    "dummy_phys",
						    dummy_log,
						    physMother,
						    false,
						    104);

  // sensitive detector //
  G4SDManager *SDMan = G4SDManager::GetSDMpointer();
  G4VSensitiveDetector *dcSD = SDMan -> FindSensitiveDetector("/DC");

  if(!dcSD)
    {
      DCSD* dcSD = new DCSD("/DC");
      SDMan -> AddNewDetector( dcSD );
    }

  G4VSensitiveDetector *dcSD1 = SDMan -> FindSensitiveDetector("/DC");
  dummy_log->SetSensitiveDetector( dcSD1 );

}

void DetectorConstruction::MakeDummyDetector6(G4VPhysicalVolume *pMother, G4ThreeVector &pos, G4RotationMatrix &rot)
{
  /// Dummy Detector Physics World ///
  /// Physics World ///
  // param //
  const G4double X_PW = 3000.0/2.0*mm;
  const G4double Y_PW=  2000.0/2.0*mm;
  const G4double Z_PW  = 0.01/2.0*mm; 
  std::string fullNamePW = "Dummy_Detector_PhysicsWorld";
  // //
  
  G4Box* solidDetectorPW = new G4Box(G4String(fullNamePW.c_str()),
				     X_PW,
				     Y_PW,
				     Z_PW
				     );

  G4LogicalVolume* logicDetectorPW = new G4LogicalVolume(solidDetectorPW,
							 mat_air,
							 G4String(fullNamePW.c_str()));

  logicDetectorPW->SetVisAttributes(new G4VisAttributes(G4Color::Yellow()));
  
  G4RotationMatrix rot_PW = rot;
  G4ThreeVector TVPW = rot*pos;
  G4VPhysicalVolume *physMother = new G4PVPlacement(G4Transform3D(rot_PW, TVPW),
						    G4String(fullNamePW.c_str()),
						    logicDetectorPW,
						    pMother,
						    false,
						    0
						    );

  // Dummy Detector //                                                                                        
  G4Box* dummy_box = new G4Box("dummy_box", X_PW, Y_PW, Z_PW);
  G4LogicalVolume*  dummy_log = new G4LogicalVolume(dummy_box, mat_air, "dummy_log", 0, 0, 0);
  dummy_log->SetVisAttributes(G4Color::Yellow());

  G4RotationMatrix dummy_rot;
  G4VPhysicalVolume* dummy_phys = new G4PVPlacement(G4Transform3D(dummy_rot,
								  G4ThreeVector(0., 0., 0.)),
						    "dummy_phys",
						    dummy_log,
						    physMother,
						    false,
						    105);

  // sensitive detector //
  G4SDManager *SDMan = G4SDManager::GetSDMpointer();
  G4VSensitiveDetector *dcSD = SDMan -> FindSensitiveDetector("/DC");

  if(!dcSD)
    {
      DCSD* dcSD = new DCSD("/DC");
      SDMan -> AddNewDetector( dcSD );
    }

  G4VSensitiveDetector *dcSD1 = SDMan -> FindSensitiveDetector("/DC");
  dummy_log->SetSensitiveDetector( dcSD1 );

}
