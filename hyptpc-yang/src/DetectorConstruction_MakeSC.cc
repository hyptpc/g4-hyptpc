/*
  DetectorConstruction_MakeSC.cc
  Make SC magnet, target, tof
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
#include "BDCSD.hh"
#include "CH1SD.hh"
#include "DCSD.hh"
#include "FTOFSD.hh"
#include "ScintSD.hh"

#include "G4VisAttributes.hh"
#include "G4Colour.hh"

#include "TString.h"
#include <TMath.h>


void DetectorConstruction::MakeSCMagnet(G4VPhysicalVolume *pMother, G4ThreeVector &pos, G4RotationMatrix &rot)
{

  /// SCMagnet Physics World ///
  /// Physics World ///
  // param //
  const G4double X_PW = 1550.0/2.0*mm;
  const G4double Y_PW= 950.0/2.0*mm;
  //const G4double Z_PW  = 1200.0/2.0*mm;
  const G4double Z_PW  = 5000.0/2.0*mm; 
  std::string fullNamePW = "SCMagnet_PhysicsWorld";
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
  
  
  /// SCMagnet Core ///                                                                                      
  //  parameters  //
  G4double size_core_Width[4];
  size_core_Width[0] = 1550*mm;
  size_core_Width[1] = 1530*mm;
  size_core_Width[2] = 1480*mm;
  size_core_Width[3] = 1470*mm;

  G4double size_core_Depth[4];
  size_core_Depth[0] = 1200*mm;
  size_core_Depth[1] = 1180*mm;
  size_core_Depth[2] = 1140*mm;
  size_core_Depth[3] = 1130*mm;

  G4double size_core_Height[4];
  size_core_Height[0] = 950*mm;
  size_core_Height[1] = 890*mm;
  size_core_Height[2] = 732*mm;
  size_core_Height[3] = 712*mm;

  G4double size_core_Corner[4];
  size_core_Corner[0] = 1449.569*mm;
  size_core_Corner[1] = 1429.569*mm;
  size_core_Corner[2] = 1357.645*mm;
  size_core_Corner[3] = 1347.645*mm;

  G4double size_core_Rad[4];
  size_core_Rad[0] = 800./2.*mm;
  size_core_Rad[1] = 820./2.*mm;
  size_core_Rad[2] = 850./2.*mm;
  size_core_Rad[3] = 860./2.*mm;

  G4double size_core_Gap[4];
  size_core_Gap[0] = 300*mm;
  size_core_Gap[1] = 320*mm;
  size_core_Gap[2] = 348*mm;
  size_core_Gap[3] = 358*mm;

  //G4double size_core_SideGap[4];
  //size_core_SideGap[0] = 725*mm;
  //size_core_SideGap[1] = 745*mm;
  //size_core_SideGap[2] = 775*mm;
  //size_core_SideGap[3] = 785*mm;

  G4double size_core_SideGap[4];
  for(int i=0 ; i<4; i++)
    {
      size_core_SideGap[i]=size_core_Depth[i]/2.0*sqrt(2.)*mm;
    }

  G4double size_core_dummyBox[3];
  size_core_dummyBox[0] = 2000*mm;
  size_core_dummyBox[1] = 2000*mm;
  size_core_dummyBox[2] = 2000*mm;

  //G4double size_core_center_sub_position= 512.5*mm;
  G4double size_core_center_sub_position[4];
  for(int i=0; i<4; i++)
    {
      size_core_center_sub_position[i]=size_core_Depth[i]/2.0;
    }

  // //

  G4Box* solidBoxCF[4];
  G4Tubs* solidTubeCF[4];
  G4Box* solidBoxCF_subOuter[4];
  G4Box* solidBoxCF_subInner[4];
  G4Box* solidBoxCF_subSide[4];

  G4SubtractionSolid* solidBoxCF_Octa[4];
  G4SubtractionSolid* solidBoxCF_subCorner[4];
  G4SubtractionSolid* solidBoxCF_subSide_1[4];
  G4SubtractionSolid* solidBoxCF_subSide_2[4];
  G4SubtractionSolid* solidBoxCF_Chamber[4];
  G4RotationMatrix* rot_magnet  = new G4RotationMatrix();
  rot_magnet->rotateZ( 45*deg);

  std::string fullNameCore = "SCMagnet_CoreFrame";
  for( int i = 0; i < 4; i++)
    {
      solidBoxCF[i]           = new G4Box(G4String((fullNameCore+Form("_Main_%d",i)).c_str()),
					  size_core_Width[i]/2. ,
					  size_core_Depth[i]/2. ,
					  size_core_Height[i]/2.);
      solidBoxCF_subOuter[i]  = new G4Box(G4String((fullNameCore+Form("_SubOuter_%d",i)).c_str()),
					  size_core_dummyBox[0],
					  size_core_dummyBox[1],
					  size_core_dummyBox[2]);
      solidBoxCF_subInner[i]  = new G4Box(G4String((fullNameCore+Form("_SubInner_%d",i)).c_str()),
					  size_core_Corner[i],
					  size_core_Corner[i],
					  size_core_dummyBox[2]);
      solidBoxCF_subSide[i]   = new G4Box(G4String((fullNameCore+Form("_SubSide_%d",i)).c_str()),
					  size_core_SideGap[i]/2.,
					  size_core_SideGap[i]/2.,
					  size_core_Gap[i]/2.);
      solidTubeCF[i]          = new G4Tubs(G4String((fullNameCore+Form("_Hole_%d",i)).c_str()),
					   0, 
					   size_core_Rad[i], 
					   2000*mm, 
					   0, 
					   2.*M_PI*rad);
      
      solidBoxCF_subCorner[i] = new G4SubtractionSolid(G4String((fullNameCore+Form("_SubCorner_%d",i)).c_str()),
						       solidBoxCF_subOuter[i],
						       solidBoxCF_subInner[i],
						       0,
						       G4ThreeVector(0,0,0));
      solidBoxCF_Octa[i]      = new G4SubtractionSolid(G4String((fullNameCore+Form("_Octa_%d", i)).c_str())    ,
						       solidBoxCF[i],
						       solidBoxCF_subCorner[i],
						       rot_magnet,
						       G4ThreeVector(0,0,0));
      solidBoxCF_subSide_1[i] = new G4SubtractionSolid(G4String((fullNameCore+Form("_SubSide_1_%d",i)).c_str()),
						       solidBoxCF_Octa[i],
						       solidBoxCF_subSide[i],
						       rot_magnet,
						       G4ThreeVector(0,size_core_center_sub_position[i],0));
      solidBoxCF_subSide_2[i] = new G4SubtractionSolid(G4String((fullNameCore+Form("_SubSide_2_%d",i)).c_str()),
						       solidBoxCF_subSide_1[i],
						       solidBoxCF_subSide[i],
						       rot_magnet,
						       G4ThreeVector(0,-size_core_center_sub_position[i],0));
      solidBoxCF_Chamber[i]   = new G4SubtractionSolid(G4String((fullNameCore+Form("_Chamber_%d",i)).c_str())  ,
						       solidBoxCF_subSide_2[i], 
						       solidTubeCF[i],
						       0, 
						       G4ThreeVector(0,0,0));
    }
  
  G4SubtractionSolid* SolidCF1 = new G4SubtractionSolid(G4String((fullNameCore+"_Solid_1").c_str()), 
							solidBoxCF_Chamber[0], 
							solidBoxCF_Chamber[1], 
							0, 
							G4ThreeVector(0,0,0));
  G4SubtractionSolid* SolidCF2 = new G4SubtractionSolid(G4String((fullNameCore+"_Solid_2").c_str()), 
							solidBoxCF_Chamber[2], 
							solidBoxCF_Chamber[3], 
							0, 
							G4ThreeVector(0,0,0));
  G4UnionSolid* solidDetectorCF = new G4UnionSolid(G4String((fullNameCore).c_str()), 
						   SolidCF1, SolidCF2, 
						   0, 
						   G4ThreeVector(0,0,0));
  
  G4LogicalVolume* logicDetectorCF = new G4LogicalVolume(solidDetectorCF, 
							 NistMan->FindOrBuildMaterial("G4_Fe"), 
							 G4String(fullNameCore.c_str()));
  logicDetectorCF->SetVisAttributes(new G4VisAttributes(G4Color::Green()));
  
  G4RotationMatrix rot_CoreFrame;
  rot_CoreFrame.rotateX(90.*deg);
  G4VPhysicalVolume *phys_SDMaget_CoreFrame = new G4PVPlacement(G4Transform3D(rot_CoreFrame, 
									      G4ThreeVector(0, 0, 0)), 
								G4String(fullNameCore.c_str()), 
								logicDetectorCF, 
								physMother, 
								false, 
								0);

  ///      ///
  
  /// SCCoil Support ///
  // param //
  G4double CoilSupPos_height = 250*mm;
  G4double RadIn = 445*mm;
  G4double RadOut= 545*mm;
  G4double GapRadIn = 465*mm;
  G4double GapRadOut = 545*mm;
  G4double SupHeight = 112*mm;
  G4double GapHeight = 62*mm;
  std::string fullNameCoilSup = "SCCoilSup";
  // //
  G4Tubs* solidTube_Sup = new G4Tubs(G4String((fullNameCoilSup+"Main")).c_str(), 
				     RadIn,
				     RadOut,
				     SupHeight/2.,
				     0,
				     360*deg);
  G4Tubs* solidTube_Sub = new G4Tubs(G4String((fullNameCoilSup+"Sub")).c_str(), 
				     GapRadIn,
				     GapRadOut+20*mm,
				     GapHeight/2.,
				     0,
				     360*deg);
  G4SubtractionSolid* solidDetectorCS = new G4SubtractionSolid( G4String((fullNameCoilSup).c_str()), 
								solidTube_Sup, 
								solidTube_Sub, 
								0, 
								G4ThreeVector(0,0,0));

  G4LogicalVolume* logicDetectorCS = new G4LogicalVolume(solidDetectorCS, 
							 NistMan->FindOrBuildMaterial("G4_Fe"), 
							 G4String(fullNameCoilSup.c_str()));

  logicDetectorCS->SetVisAttributes(new G4VisAttributes(G4Color::Green()));
  //logicDetectorCS->SetVisAttributes(new G4VisAttributes(G4Color::Blue()));

  G4RotationMatrix rot_CoilSupport;
  rot_CoilSupport.rotateX(90.*deg);
  
  G4VPhysicalVolume *phys_SDCoilSupport_Up =   new G4PVPlacement(G4Transform3D(rot_CoilSupport, G4ThreeVector(0, CoilSupPos_height, 0)), 
								 G4String((fullNameCoilSup+"Up")).c_str(), 
								 logicDetectorCS, 
								 physMother, 
								 false, 
								 0);
  G4VPhysicalVolume *phys_SDCoilSupport_Down =   new G4PVPlacement(G4Transform3D(rot_CoilSupport, G4ThreeVector(0, -CoilSupPos_height, 0)), 
								   G4String((fullNameCoilSup+"Down")).c_str(), 
								   logicDetectorCS, physMother, 
								   false, 
								   0);
  /// ///


  /// SCMagnet Coil ///
  // param //
  G4double coilRad_in;
  G4double coilRad_out;
  G4double coilHeight;

  coilRad_in  = 466*mm;
  coilRad_out = 535*mm;
  coilHeight  = 60*mm;

  G4double coilPos[3];
  coilPos[0] = 0*mm;
  coilPos[1] = 250*mm;
  coilPos[2] = 0*mm;
  std::string fullNameCoil = "SCMagnet_Coil";

  G4Tubs* solidDetectorCoil = new G4Tubs(G4String(fullNameCoil.c_str()),
					 coilRad_in,
					 coilRad_out,
					 coilHeight/2.,
					 0*deg,
					 360*deg);

  G4LogicalVolume* logicDetectorCoil = new G4LogicalVolume(solidDetectorCoil,
							   NistMan->FindOrBuildMaterial("G4_Cu"),
							   G4String(fullNameCoil.c_str()));
  logicDetectorCoil->SetVisAttributes(new G4VisAttributes(G4Color::Yellow()));

  G4RotationMatrix rot_Coil;
  rot_Coil.rotateX(90.*deg);
  G4VPhysicalVolume *Phys_Coil_Up = new G4PVPlacement(G4Transform3D(rot_Coil, G4ThreeVector(coilPos[0], coilPos[1], coilPos[2])),
						      G4String((fullNameCoil+"Up").c_str()),
						      logicDetectorCoil,
						      physMother,
						      false,
						      0);
  G4VPhysicalVolume *Phys_Coil_Down = new G4PVPlacement(G4Transform3D(rot_Coil, G4ThreeVector(coilPos[0], -coilPos[1], coilPos[2])),
							G4String((fullNameCoil+"Down").c_str()),
							logicDetectorCoil,
							physMother,
							false,
							0);
  /// Detectors inside SC magnet///
  G4RotationMatrix unrotate_mat = rot.inverse();


  G4ThreeVector TPCPos ( 0.0, 0.0, 0.0*mm);
  G4RotationMatrix TPCRot;
  TPCRot.rotateX(-90.*deg);
  TPCRot.rotateY(90.*deg);
  const G4double target_offset = 143.*mm;
  //const G4double DZ_Target = 15.*mm;
  //const G4double DZ_TargetHolder  = 295.*mm; // in old code 310.*mm. it is placed out of TPC Frame.   
  const G4double DZ_TargetHolder  = 290.*mm; // 20 mm target y   
  //G4ThreeVector TargetPos ( 0., target_offset , (DZ_TargetHolder/2.));// for C target
  G4ThreeVector TargetPos ( 0., target_offset , 91.0);// for LH2 target
  //G4ThreeVector TargetPos ( 0., 0., 0.);
  G4RotationMatrix TargetRot;
  TargetRot.rotateX( -90.*deg ); // for all targets
  G4ThreeVector TOFPos ( 0., 0., 0.0*mm);
  G4RotationMatrix TOFRot;
  TOFRot.rotateX( -90.*deg );

  TargetRot = unrotate_mat*TargetRot;
  TPCRot = unrotate_mat*TPCRot;
  //TOFRot = unrotate_mat*TOFRot;

  G4ThreeVector TargetDummyPos ( 0., 0., -143.0);// for LH2 target                            
  //G4ThreeVector TargetDummyPos ( 0., 0., -650.0);// for LH2 target            
  G4RotationMatrix TargetDummyRot;
  TargetDummyRot = unrotate_mat*TargetDummyRot;

  G4ThreeVector FACPos ( 0., 0., 570.0);// for Beam Veto AC
  G4RotationMatrix FACRot;
  FACRot = unrotate_mat*FACRot;

  G4ThreeVector BDCPos ( 0., 0., -570.0);// for Beam definition counter                            
  G4RotationMatrix BDCRot;
  BDCRot = unrotate_mat*FACRot;

  MakeTOF( physMother, TOFPos, TOFRot );
  MakeHypTPC2( physMother, TPCPos, TPCRot  );
  MakeFAC( physMother, FACPos, FACRot );
  MakeBDC( physMother, BDCPos, BDCRot );


  //MakeTarget( physMother, TargetPos, TargetRot ); //E42 target
  //MakeTargetDummy( physMother, Targ etDummyPos, TargetDummyRot ); // dummy detector

  ///MakeHypTPC( physMother, TPCPos, TPCRot  );
  //MakeTargetH( physMother, TargetPos, TargetRot );
  MakeTargetHver2( physMother, TargetPos, TargetRot ); //update 2018 July

}

/*
void DetectorConstruction::MakeHypTPC(G4VPhysicalVolume *pMother, G4ThreeVector &pos, G4RotationMatrix &rot)
{
  // material //

  /// Physics world ///
  // param //
  const G4double STANG_PW     = 22.5*deg;
  const G4double DPH_PW      = 360.*deg;
  const G4double NPL_PW       = 8.;
  const G4double RIN_PW       = 0.0/2.*mm;
  const G4double ROUT_PW      = (656./2.)*mm;
  const G4double DZ_PW        =  802.*mm;
  const G4double DZ_PW_OFFSET =    0.*mm;
  std::string fullNamePW = "TPC_PhysiWorld";
  // //
  double z_PW[2] = {-DZ_PW/2.,DZ_PW/2.};
  double rmin_PW[2] = {RIN_PW, RIN_PW};
  double rmax_PW[2] = {ROUT_PW, ROUT_PW};
  G4Polyhedra* solidDetectorPW = new G4Polyhedra(G4String(fullNamePW.c_str()),
						 STANG_PW,
						 DPH_PW,
						 int(NPL_PW),
						 2,
						 z_PW,
						 rmin_PW,
						 rmax_PW);
  
  G4LogicalVolume* logicDetectorPW = new G4LogicalVolume(solidDetectorPW,
							 mat_air,
							 G4String(fullNamePW.c_str()));

  logicDetectorPW->SetVisAttributes(new G4VisAttributes(G4Color::Green()));

  // position //
  G4RotationMatrix rot_PW = rot;
  G4ThreeVector pos_PW = rot*pos;

  G4VPhysicalVolume *phyMother = new G4PVPlacement(G4Transform3D(rot_PW, pos_PW),
						   G4String(fullNamePW.c_str()),
						   logicDetectorPW,
						   pMother,
						   false,
						   0);

  /// TPC fieldcage ///
  // param //
  const G4double STANG_FC     = 22.5*deg;
  const G4double DPHI_FC      = 360.*deg;
  const G4double NPL_FC       = 8.;
  const G4double RIN_FC       = 574.5/2.*mm;
  const G4double ROUT_FC      = (654./2.)*mm;
  const G4double DZ_FC        =  800.*mm;
  const G4double DZ_FC_OFFSET =    0.*mm;
  std::string fullNameFC = "TPC_FieldCage";
  // //

  double z_FC[2] = {-DZ_FC/2.,DZ_FC/2.};
  double rmin_FC[2] = {RIN_FC, RIN_FC};
  double rmax_FC[2] = {ROUT_FC, ROUT_FC};
  G4Polyhedra* solidDetectorFC = new G4Polyhedra(G4String(fullNameFC.c_str()),
						 STANG_FC,
						 DPHI_FC,
						 int(NPL_FC),
						 2,
						 z_FC,
						 rmin_FC,
						 rmax_FC);
  
  G4LogicalVolume* logicDetectorFC = new G4LogicalVolume(solidDetectorFC,
							 mat_p10,
							 G4String(fullNameFC.c_str()));
  logicDetectorFC->SetVisAttributes(new G4VisAttributes(G4Color::Magenta()));
  G4RotationMatrix rot_FieldCage;
  //rot_FieldCage.rotateX(-90.*deg);
  G4VPhysicalVolume *Phys_TPC_Field_Cage = new G4PVPlacement(G4Transform3D(rot_FieldCage, G4ThreeVector(0, 0, 0)),
							     G4String(fullNameFC.c_str()),
							     logicDetectorFC,
							     phyMother,
							     false,
							     0);
  /// ///
  
  /// TPC frame ///
  // param //
  const G4double STANG_TPC     = 22.5*deg;
  const G4double DPHI_TPC      = 360.*deg;
  const G4double DZ_TPC        = 604.*mm;
  const G4double NPL_TPC       = 8.;
  const G4double RIN_TPC       = 0.0*mm;
  const G4double ROUT_TPC      = 574./2.*mm;
  const G4double DZ_TPC_OFFSET = 0.*mm;
  std::string fullNameTPC = "TPC_Frame";
  
  double z_TPC[2] = {-DZ_TPC/2.,DZ_TPC/2.};
  double rmin_TPC[2] = {RIN_TPC, RIN_TPC};
  double rmax_TPC[2] = {ROUT_TPC, ROUT_TPC};
  G4Polyhedra* solidDetectorTPC = new G4Polyhedra(G4String(fullNameTPC.c_str()),
						  STANG_TPC,
						  DPHI_TPC,
						  int(NPL_TPC),
						  2,
						  z_TPC,
						  rmin_TPC,
						  rmax_TPC);
  
  G4LogicalVolume* logicDetectorTPC = new G4LogicalVolume(solidDetectorTPC,
							  mat_p10,
							  G4String(fullNameTPC.c_str()));
  logicDetectorTPC->SetVisAttributes(new G4VisAttributes(G4Color::Green()));
  G4RotationMatrix rot_TPC;
  //rot_TPC.rotateX(-90.*deg);
  G4VPhysicalVolume *Phys_TPC = new G4PVPlacement(G4Transform3D(rot_TPC, G4ThreeVector(0, 0, 0)),
						  G4String(fullNameTPC.c_str()),
						  logicDetectorTPC,
						  phyMother,
						  false,
						  0);
  /// ///

  /// TPC PAD ///
  // Parameter //
  // PadParameter :                                                                                           
  // [0]: RingID                                                                                              
  // [1]: Number of pads in a ring                                                                            
  // [2]: Center radius                                                                                       
  // [3]: dTheta of a pad                                                                                     
  // [4]: Start angle of ring                                                                                 
  // [5]: length of a pad
  //Original
  /*
  const double PadParameter[32][6] = {                                                                        
    {0,         48,     14.5,   7.5     ,0.     ,9.},                                                         
    {1,         48,     24.,    7.5     ,0.     ,9.},                                                         
    {2,         72,     33.5,   5.      ,0.     ,9.},                                                         
    {3,         96,     43.,    3.75    ,0.     ,9.},                                                         
    {4,         120,    52.5,   3.      ,0.     ,9.},                                                         
    {5,         144,    62.,    2.5     ,0.     ,9.},                                                         
    {6,         168,    71.5,   2.14286 ,0.     ,9.},                                                         
    {7,         192,    81.,    1.875   ,0.     ,9.},                                                         
    {8,         216,    90.5,   1.66667 ,0.     ,9.},                                                         
    {9,         240,    100.,   1.5     ,0.     ,9.},                                                         
    {10,        208,    111.25, 1.49375 ,24.65  ,12.5},                                                       
    {11,        218,    124.25, 1.32844 ,35.2   ,12.5},                                                       
    {12,        230,    137.25, 1.2     ,42     ,12.5},                                                       
    {13,        214,    150.25, 1.09093 ,63.27  ,12.5},                                                       
    {14,        212,    163.25, 1.      ,74     ,12.5},                                                       
    {15,        214,    176.25, 0.923084,       81.23   ,12.5},                                               
    {16,        220,    189.25, 0.857182,       85.71   ,12.5},                                               
    {17,        224,    202.25, 0.801786,       90.2    ,12.5},                                               
    {18,        232,    215.25, 0.751552,       92.82   ,12.5},                                               
    {19,        238,    228.25, 0.707227,       95.84   ,12.5},                                               
    {20,        244,    241.25, 0.667869,       98.52   ,12.5},                                               
    {21,        232,    254.25, 0.632672,       106.61  ,12.5},                                               
    {22,        218,    267.25, 0.60101 ,       114.49  ,12.5},                                               
    {23,        210,    280.25, 0.573238,       119.81  ,12.5},                                               
    {24,        206,    293.25, 0.547111,       123.648 ,12.5},                                               
    {25,        202,    306.25, 0.523267,       127.15  ,12.5},                                               
    {26,        200,    319.25, 0.5014  ,       129.86  ,12.5},                                               
    {27,        196,    332.25, 0.481327,       132.83  ,12.5},                                               
    {28,        178,    345.25, 0.463371,       138.76  ,12.5},                                     
    {29,        130,    358.25, 0.446154,       151     ,12.5},                                               
    {30,        108,    371.25, 0.430185,       156.77  ,12.5},                                               
    {31,        90,     384.25, 0.415333,       161.31  ,12.5}};


  //const G4double Rad_out = 33.6*mm; //for C target
  const G4double Rad_out = 80.0/2.0*mm; //for LH2 target                                                                          
  //const G4double DZ_pad_under_target = 240*mm;                                                                
  const G4double DZ_pad_under_target = (275-118)*mm; //for LH2 target                                                                
  const G4double DZ_pad_normal_target= 550*mm;                                                                
  //const G4double PZ_pad_under_target = 155*mm; //for C target
  const G4double PZ_pad_under_target = DZ_pad_normal_target/2.0 - DZ_pad_under_target/2.0; //for LH2 target
  const G4double PZ_pad_normal_target= 0*mm;                                                                    
  const G4double tpc_centerOffset = 143.*mm;
  std::string fullNameTPCpad = "TPC_Pad";
  // //
  
  int padID = 0;                                                                                              
  G4VPhysicalVolume *tpc_ring[32];
  G4Tubs *solidTube[32];
  G4LogicalVolume *logicDetectorTPCpad[32];

  ///// Sensitive Detector /////
  G4SDManager *SDMan = G4SDManager::GetSDMpointer();
  TpcSD *tpcSD = new TpcSD( "/TPC/Tpc" );
  SDMan->AddNewDetector( tpcSD );
  //logicDetectorTPC->SetSensitiveDetector( tpcSD );      
    
  for( int i = 0; i< 32; i++)
    {                                                                                
      padID = i;                                                                                                
      G4double RMIN_pad = (PadParameter[i][2] - PadParameter[i][3]/2.)*mm;                                      
      G4double RMAX_pad = (PadParameter[i][2] + PadParameter[i][3]/2.)*mm;                                      
      G4double PZ_TPCPAD= 0;                                                                                    
      G4double DZ_TPCPAD= 0;                                                                                    
      if( RMIN_pad < Rad_out )
	{                                                                                 
	  DZ_TPCPAD = DZ_pad_under_target;                                                                        
	  PZ_TPCPAD = PZ_pad_under_target;                                                                        
	}
      else
	{                                                                                                    
	  DZ_TPCPAD = DZ_pad_normal_target;                                                                       
	  PZ_TPCPAD = PZ_pad_normal_target;                                                                       
	}                                                                                                         
      
      if( i < 10 )
	{                                                                                             
	  solidTube[i] = new G4Tubs(G4String((fullNameTPCpad+Form("_%d",i)).c_str()), 
				    RMIN_pad,
				    RMAX_pad,
				    DZ_TPCPAD/2.,
				    0,
				    360*deg);                                                          
	}
      else
	{                                                                                                    
	  solidTube[i] = new G4Tubs(G4String((fullNameTPCpad+Form("_%d",i)).c_str()), 
					     RMIN_pad,
					     RMAX_pad,
					     DZ_TPCPAD/2.,
					     (PadParameter[i][4]+90)*deg,
					     (PadParameter[i][3]*PadParameter[i][1])*deg);                           
	} 
      
      logicDetectorTPCpad[i] = new G4LogicalVolume(solidTube[i], 
						   mat_p10, 
						   G4String((fullNameTPCpad+Form("_%d",i)).c_str()));
      
      logicDetectorTPCpad[i]->SetVisAttributes(new G4VisAttributes(G4Color::Yellow()));                             
      
      G4RotationMatrix rot_TPCpad;
      rot_TPCpad.rotateZ( -90.*deg );
      G4ThreeVector TVTPCpad ( 0, tpc_centerOffset, -PZ_TPCPAD);
      TVTPCpad.rotateZ ( -90.*deg);
      tpc_ring[i] = new G4PVPlacement(G4Transform3D(rot_TPCpad, TVTPCpad),
				      G4String((fullNameTPCpad+Form("_%d",i)).c_str()),
				      logicDetectorTPCpad[i],
				      Phys_TPC,
				      false,
				      i);
      logicDetectorTPCpad[i]->SetSensitiveDetector( tpcSD );      
    }

}
*/

void DetectorConstruction::MakeTarget(G4VPhysicalVolume *pMother, G4ThreeVector &pos, G4RotationMatrix &rot)
{
  // material //

  //std::cout<<"Test: " << mat_C <<std::endl;
  /// Physics World ///
  // param //
  const G4double RIN_PW = 0.0*mm;
  const G4double ROUT_PW= 25.2*mm;
  const G4double DZ_PW  = 310.0*mm; // in old code 310.*mm. it is placed out of TPC Frame.   
  const G4double STANG_PW = 0.*deg;
  const G4double EDANG_PW = 360.*deg;
  std::string fullNamePW = "Target_PhysicsWorld";
  // //
  
  G4Tubs* solidDetectorPW = new G4Tubs(G4String(fullNamePW.c_str()),
				       RIN_PW,
				       ROUT_PW,
				       DZ_PW/2.,
				       STANG_PW,
				       EDANG_PW);
  
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

  /// Target ///
  // param //
  const G4double DX_Target = 30.*mm;
  const G4double DY_Target = 20.*mm;
  const G4double DZ_Target = 20.*mm;
  //const G4double DZ_TargetHolder  = 295.*mm; // in old code 310.*mm. it is placed out of TPC Frame.   
  const G4double DZ_TargetHolder  = 290.*mm; // in old code 310.*mm. it is placed out of TPC Frame.   
  const G4double tpc_centerOffset = 143.*mm;
  std::string fullNameTarget = "TPC_Target";
  // //
  G4Box* solidDetectorTarget = new G4Box(G4String(fullNameTarget.c_str()),
					 DX_Target/2.,
					 DY_Target/2.,
					 DZ_Target/2.);
  
  G4LogicalVolume* logicDetectorTarget = new G4LogicalVolume(solidDetectorTarget,
							     mat_C,
							     G4String(fullNameTarget.c_str()));
  logicDetectorTarget->SetVisAttributes(new G4VisAttributes(G4Color::Red()));

  G4RotationMatrix rot_Target;
  //rot_Target.rotateX( -90.*deg);
  //G4ThreeVector TVTarget ( 0, tpc_centerOffset, 0);
  G4ThreeVector TVTarget ( 0, 0, -(DZ_TargetHolder)/2.0);
  //TVTarget.rotateX( -90.*deg);
  G4VPhysicalVolume *Phys_Target = new G4PVPlacement(G4Transform3D(rot_Target, TVTarget),
						     G4String(fullNameTarget.c_str()),
						     logicDetectorTarget,
						     physMother,
						     false,
						     0
						     );

  /// Target holder ///
  // param //
  const G4double RIN_TargetHolder = 25*mm;
  const G4double ROUT_TargetHolder= 25.2*mm;

  const G4double STANG_TargetHolder = 0.*deg;
  const G4double EDANG_TargetHolder = 360.*deg;
  const G4double DZOffset_TargetHolder = DZ_TargetHolder/2. + DZ_Target/2.;
  std::string fullNameTH = "TPC_TargetHolder";
  // //
  
  G4Tubs* solidDetectorTH = new G4Tubs(G4String(fullNameTH.c_str()),
				       RIN_TargetHolder,
				       ROUT_TargetHolder,
				       DZ_TargetHolder/2.,
				       STANG_TargetHolder,
				       EDANG_TargetHolder);
  
  G4LogicalVolume* logicDetectorTH = new G4LogicalVolume(solidDetectorTH,
							 mat_p10,
							 G4String(fullNameTH.c_str()));
  logicDetectorTH->SetVisAttributes(new G4VisAttributes(G4Color::Magenta()));

  G4RotationMatrix rot_TH;
  //rot_TH.rotateX( -90.*deg);
  //G4ThreeVector TVTH ( 0, tpc_centerOffset, DZOffset_TargetHolder);
  //G4ThreeVector TVTH ( 0, 0, DZOffset_TargetHolder);
  G4ThreeVector TVTH ( 0, 0, DZ_Target/2.);
  //TVTH.rotateX ( -90.*deg);
  G4VPhysicalVolume *Phys_TargetHolder = new G4PVPlacement(G4Transform3D(rot_TH, TVTH),
							   G4String(fullNameTH.c_str()),
							   logicDetectorTH,
							   physMother,
							   false,
							   0
							   );
  /// ///
  
  ///// Sensitive Detector /////
  G4SDManager *SDMan = G4SDManager::GetSDMpointer();
  TargetSD *targetSD = new TargetSD( "/TPC/Target" );
  SDMan->AddNewDetector( targetSD );
  logicDetectorTarget->SetSensitiveDetector( targetSD );
  
}

void DetectorConstruction::MakeTOF(G4VPhysicalVolume *pMother, G4ThreeVector &pos, G4RotationMatrix &rot)
{

  /// Physics World ///
  // param //
  const G4double STANG_PW     = 22.5*deg;
  const G4double DPHI_PW      = 360.*deg;
  const G4double NPL_PW       = 8.;
  //const G4double RIN_PW       = 656.0/2.*mm;
  const G4double RIN_PW       = 688.0/2.*mm;
  //const G4double ROUT_PW      = 800./2.*mm;
  const G4double ROUT_PW      = 712./2.*mm;
  const G4double DZ_PW        =  802.*mm;
  const G4double DZ_PW_OFFSET =    0.*mm;
  std::string fullNamePW = "TOF_PhysiWorld";
  // //
  double z_PW[2] = {-DZ_PW/2.,DZ_PW/2.};
  double rmin_PW[2] = {RIN_PW, RIN_PW};
  double rmax_PW[2] = {ROUT_PW, ROUT_PW};
  G4Polyhedra* solidDetectorPW = new G4Polyhedra(G4String(fullNamePW.c_str()),
						 STANG_PW,
						 DPHI_PW,
						 int(NPL_PW),
						 2,
						 z_PW,
						 rmin_PW,
						 rmax_PW);
  
  G4LogicalVolume* logicDetectorPW = new G4LogicalVolume(solidDetectorPW,
							 mat_air,
							 G4String(fullNamePW.c_str()));
  //logicDetectorPW->SetVisAttributes(new G4VisAttributes(G4Color::Blue()));

  // position //
  G4RotationMatrix rot_PW = rot;
  G4ThreeVector pos_PW = rot*pos;

  G4VPhysicalVolume *phyMother = new G4PVPlacement(G4Transform3D(rot_PW, pos_PW),
						   G4String(fullNamePW.c_str()),
						   logicDetectorPW,
						   pMother,
						   false,
						   0);


  /// TOF ///
  // param //
  const G4double DX_TOF = 70.*mm;
  const G4double DY_TOF = 800.*mm;
  const G4double DZ_TOF = 10.*mm;
  const G4double tof_centerOffset = 350.*mm;
  std::string fullNameTOF = "TOF";
  // //
  G4Box* solidDetectorTOF = new G4Box(G4String(fullNameTOF.c_str()),
				      DX_TOF/2.,
				      DY_TOF/2.,
				      DZ_TOF/2.);
  
  G4LogicalVolume* logicDetectorTOF = new G4LogicalVolume(solidDetectorTOF,
  							  mat_Scin,
  							  G4String(fullNameTOF.c_str()));

  //beam window segments
  G4Box* solidDetectorTOF_window_all = new G4Box(G4String(fullNameTOF.c_str()),
				      DX_TOF/2.,
				      DY_TOF/2.,
				      DZ_TOF/2.);

  const G4double DX_TOF_sub = 70.0*mm;
  const G4double DY_TOF_sub = 140.0*mm;//window size +-2 cm
  const G4double DZ_TOF_sub = 10.*mm;
  G4Box* solidDetectorTOF_window_sub = new G4Box(G4String(fullNameTOF.c_str()),
						 DX_TOF_sub/2.,
						 DY_TOF_sub/2.,
						 DZ_TOF_sub/2.);

  G4VSolid* solidDetectorTOF_window = new G4SubtractionSolid(G4String(fullNameTOF.c_str()),
							     solidDetectorTOF_window_all,
							     solidDetectorTOF_window_sub);

  G4LogicalVolume* logicDetectorTOF_window = new G4LogicalVolume(solidDetectorTOF_window,
								 mat_Scin,
								 G4String(fullNameTOF.c_str()));


  logicDetectorTOF->SetVisAttributes(new G4VisAttributes(G4Color::Red()));
  logicDetectorTOF_window->SetVisAttributes(new G4VisAttributes(G4Color::Red()));

  const G4int NSeg = 4;
  const G4int NSurface = 8;

  ///// Sensitive Detector /////
  G4SDManager *SDMan = G4SDManager::GetSDMpointer();
  TofSD *tofSD = new TofSD( "/TPC/Tof" );
  SDMan->AddNewDetector( tofSD );

  G4VPhysicalVolume *Phys_TOF[NSurface*4];
  //G4LogicalVolume *logicDetectorTOF[NSurface*4];
  for(int i=0; i<NSurface; i++)
    {
      G4RotationMatrix rot_TOF;
      rot_TOF.rotateY( (45*i)*deg);
      rot_TOF.rotateX( 90.*deg);
      for(int j=0; j<NSeg; j++)
	{

	  //logicDetectorTOF[j+i*4]= new G4LogicalVolume(solidDetectorTOF,
	  //						mat_Scin,
	  //						G4String(fullNameTOF.c_str()));
	  //logicDetectorTOF[j+i*4]->SetVisAttributes(new G4VisAttributes(G4Color::Red()));
	  //logicDetectorTOF[j+i*4]->SetSensitiveDetector( tofSD );

	  G4double x_tof = -1.5*DX_TOF + DX_TOF*j;
	  G4double y_tof = 0.0;
	  G4double z_tof = tof_centerOffset;
	  G4ThreeVector TVTof ( x_tof, y_tof, z_tof );
	  TVTof.rotateY( (45*i)*deg);
	  TVTof.rotateX( 90.*deg);

	  #if 0
	  Phys_TOF[j+i*4] = new G4PVPlacement(G4Transform3D(rot_TOF, TVTof),
					      G4String(fullNameTOF.c_str()),
					      logicDetectorTOF,
					      phyMother,
					      false,
					      j+i*4);
	  #endif

	  #if 1
	  //beam window case
          if( (i == 3 && j == 3)|| (i==4 && (j == 0 || j==1 || j==2 || j==3 ))) //window              
	    {
	      Phys_TOF[j+i*4] = new G4PVPlacement(G4Transform3D(rot_TOF, TVTof),
						  G4String(fullNameTOF.c_str()),
						  logicDetectorTOF_window,
						  phyMother,
						  false,
						  j+i*4);
	    }
	  else
	    {
	      Phys_TOF[j+i*4] = new G4PVPlacement(G4Transform3D(rot_TOF, TVTof),
						  G4String(fullNameTOF.c_str()),
						  logicDetectorTOF,
						  phyMother,
						  false,
						  j+i*4);
	    }
	  #endif
	}
    }
  logicDetectorTOF->SetSensitiveDetector( tofSD );
  logicDetectorTOF_window->SetSensitiveDetector( tofSD );
}


void DetectorConstruction::MakeTargetH(G4VPhysicalVolume *pMother, G4ThreeVector &pos, G4RotationMatrix &rot)
{
  // material //

  //std::cout<<"Test: " << mat_C <<std::endl;
  /// Physics World ///
  // param //
  const G4double RIN_PW = 0.0*mm;
  const G4double ROUT_PW= 80.0/2.0*mm;
  const G4double DZ_PW  = 418.0/2.0*mm; // target tpc center (y-axis), up 350 mm, down 68 mm. 
  const G4double DZ_PW_DOWN = (118-50.0)*mm;
  const G4double STANG_PW = 0.*deg;
  const G4double EDANG_PW = 360.*deg;
  std::string fullNamePW = "Target_PhysicsWorld";
  // //
  
  G4Tubs* solidDetectorPW = new G4Tubs(G4String(fullNamePW.c_str()),
				       RIN_PW,
				       ROUT_PW,
				       DZ_PW,
				       STANG_PW,
				       EDANG_PW);
  
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

  /// ///
  
  /// Target ///
  // param //
  const G4double RIN_Target = 0.0*mm;
  const G4double ROUT_Target= 53.0/2.0*mm; // target radius 53 mm
  const G4double DZ_Target  = 100.0/2.0*mm; // target length 100 mm   
  const G4double STANG_Target = 0.*deg;
  const G4double EDANG_Target = 360.*deg;
  std::string fullNameTarget = "TPC_Target";
  // //
  G4Tubs* solidDetectorTarget = new G4Tubs(G4String(fullNameTarget.c_str()),
					   RIN_Target,
					   ROUT_Target,
					   DZ_Target,
					   STANG_Target,
					   EDANG_Target);
  
  G4LogicalVolume* logicDetectorTarget = new G4LogicalVolume(solidDetectorTarget,
							     mat_LH2,
							     G4String(fullNameTarget.c_str()));

  logicDetectorTarget->SetVisAttributes(new G4VisAttributes(G4Color::Red()));

  G4RotationMatrix rot_Target;
  G4ThreeVector TVTarget ( 0, 0, -(DZ_PW - DZ_PW_DOWN - DZ_Target));
  
  G4VPhysicalVolume *Phys_Target = new G4PVPlacement(G4Transform3D(rot_Target, TVTarget),
						     G4String(fullNameTarget.c_str()),
						     logicDetectorTarget,
						     physMother,
						     false,
						     0
						     );
  
  /// Target holder Sidewall1///

  // param //
  const G4double RIN_TargetHolder_SideWall1 = (80.0/2.0-1.5)*mm;
  const G4double ROUT_TargetHolder_SideWall1= 80.0/2.0*mm;
  const G4double DZ_TargetHolder_SideWall1 = 418.0/2.0*mm;
  const G4double STANG_TargetHolder_SideWall1 = 0.*deg;
  const G4double EDANG_TargetHolder_SideWall1 = 360.*deg;
  const G4double EDANG_TargetHolder_SideWall1_Window = 26.0*deg;
  const G4double STANG_TargetHolder_SideWall1_Window = (-270.0- EDANG_TargetHolder_SideWall1_Window/deg/2.0)*deg;
  const G4double ANG_TargetHolder_SideWall1_Window = 28.0*deg;
  //const G4double DZ_Window = DZ_TargetHolder_SideWall1+1.0;
  const G4double DZ_Window = DZ_Target;
  std::string fullNameTHSW1 = "TPC_TargetHolder_SideWall1";
  std::string fullNameTHSW1_Window = "TPC_TargetHolder_SideWall1_Window";

  G4ThreeVector TVTHSW1_Window (TVTarget);
  //G4ThreeVector TVTHSW1_Window;
  // //

  // Making window to SideWall1 //
  // Window SideWall1 //
  G4Tubs* solidDetectorTHSW1_Window1 = new G4Tubs(G4String(fullNameTHSW1_Window.c_str()),
						  RIN_TargetHolder_SideWall1-1,
						  ROUT_TargetHolder_SideWall1+1,
						  DZ_Window,
						  STANG_TargetHolder_SideWall1_Window,
						  EDANG_TargetHolder_SideWall1_Window);
  // //
  G4Tubs* solidDetectorTHSW1 = new G4Tubs(G4String(fullNameTHSW1.c_str()),
					  RIN_TargetHolder_SideWall1,
					  ROUT_TargetHolder_SideWall1,
					  DZ_TargetHolder_SideWall1,
					  STANG_TargetHolder_SideWall1,
					  EDANG_TargetHolder_SideWall1);

  G4VSolid* solidDetectorTHSW1_sub1 = new G4SubtractionSolid("Wall-Window1",
							    solidDetectorTHSW1,
							    solidDetectorTHSW1_Window1,
							    0,
							    TVTHSW1_Window
							    );
  
  G4Tubs* solidDetectorTHSW1_Window2 = new G4Tubs(G4String(fullNameTHSW1_Window.c_str()),
						  RIN_TargetHolder_SideWall1-1,
						  ROUT_TargetHolder_SideWall1+1,
						  DZ_Window,
						  STANG_TargetHolder_SideWall1_Window-2*ANG_TargetHolder_SideWall1_Window,
						  EDANG_TargetHolder_SideWall1_Window);
  
  G4VSolid* solidDetectorTHSW1_sub2 = new G4SubtractionSolid("Wall-Window2",
							     solidDetectorTHSW1_sub1,
							     solidDetectorTHSW1_Window2,
							     0,
							     TVTHSW1_Window
							     );
  
  G4Tubs* solidDetectorTHSW1_Window3 = new G4Tubs(G4String(fullNameTHSW1_Window.c_str()),
						  RIN_TargetHolder_SideWall1-1,
						  ROUT_TargetHolder_SideWall1+1,
						  DZ_Window,
						  STANG_TargetHolder_SideWall1_Window-ANG_TargetHolder_SideWall1_Window,
						  EDANG_TargetHolder_SideWall1_Window);

  G4VSolid* solidDetectorTHSW1_sub3 = new G4SubtractionSolid("Wall-Window3",
							     solidDetectorTHSW1_sub2,
							     solidDetectorTHSW1_Window3,
							     0,
							     TVTHSW1_Window
							     );

  G4Tubs* solidDetectorTHSW1_Window4 = new G4Tubs(G4String(fullNameTHSW1_Window.c_str()),
						  RIN_TargetHolder_SideWall1-1,
						  ROUT_TargetHolder_SideWall1+1,
						  DZ_Window,
						  STANG_TargetHolder_SideWall1_Window+ ANG_TargetHolder_SideWall1_Window,
						  EDANG_TargetHolder_SideWall1_Window);

  G4VSolid* solidDetectorTHSW1_sub4 = new G4SubtractionSolid("Wall-Window4",
							     solidDetectorTHSW1_sub3,
							     solidDetectorTHSW1_Window4,
							     0,
							     TVTHSW1_Window
							     );
  
  G4Tubs* solidDetectorTHSW1_Window5 = new G4Tubs(G4String(fullNameTHSW1_Window.c_str()),
						  RIN_TargetHolder_SideWall1-1,
						  ROUT_TargetHolder_SideWall1+1,
						  DZ_Window,
						  STANG_TargetHolder_SideWall1_Window + 2*ANG_TargetHolder_SideWall1_Window,
						  EDANG_TargetHolder_SideWall1_Window);
  
  G4VSolid* solidDetectorTHSW1_sub5 = new G4SubtractionSolid("Wall-Window5",
							     solidDetectorTHSW1_sub4,
							     solidDetectorTHSW1_Window5,
							     0,
							     TVTHSW1_Window
							     );
  
  G4Tubs* solidDetectorTHSW1_Window6 = new G4Tubs(G4String(fullNameTHSW1_Window.c_str()),
						  RIN_TargetHolder_SideWall1-1,
						  ROUT_TargetHolder_SideWall1+1,
						  DZ_Window,
						  STANG_TargetHolder_SideWall1_Window+180.0*degree,
						  EDANG_TargetHolder_SideWall1_Window);
  
  G4VSolid* solidDetectorTHSW1_sub6 = new G4SubtractionSolid("Wall-Window6",
							     solidDetectorTHSW1_sub5,
							     solidDetectorTHSW1_Window6,
							     0,
							     TVTHSW1_Window
							     );

  G4Tubs* solidDetectorTHSW1_Window7 = new G4Tubs(G4String(fullNameTHSW1_Window.c_str()),
						  RIN_TargetHolder_SideWall1-1,
						  ROUT_TargetHolder_SideWall1+1,
						  DZ_Window,
						  STANG_TargetHolder_SideWall1_Window+180.0*degree - 2*ANG_TargetHolder_SideWall1_Window,
						  EDANG_TargetHolder_SideWall1_Window);
  
  G4VSolid* solidDetectorTHSW1_sub7 = new G4SubtractionSolid("Wall-Window7",
							     solidDetectorTHSW1_sub6,
							     solidDetectorTHSW1_Window7,
							     0,
							     TVTHSW1_Window
							     );

  G4Tubs* solidDetectorTHSW1_Window8 = new G4Tubs(G4String(fullNameTHSW1_Window.c_str()),
						  RIN_TargetHolder_SideWall1-1,
						  ROUT_TargetHolder_SideWall1+1,
						  DZ_Window,
						  STANG_TargetHolder_SideWall1_Window+180.0*degree - ANG_TargetHolder_SideWall1_Window,
						  EDANG_TargetHolder_SideWall1_Window);
  
  G4VSolid* solidDetectorTHSW1_sub8 = new G4SubtractionSolid("Wall-Window8",
							     solidDetectorTHSW1_sub7,
							     solidDetectorTHSW1_Window8,
							     0,
							     TVTHSW1_Window
							     );

  G4Tubs* solidDetectorTHSW1_Window9 = new G4Tubs(G4String(fullNameTHSW1_Window.c_str()),
						  RIN_TargetHolder_SideWall1-1,
						  ROUT_TargetHolder_SideWall1+1,
						  DZ_Window,
						  STANG_TargetHolder_SideWall1_Window+180.0*degree + ANG_TargetHolder_SideWall1_Window,
						  EDANG_TargetHolder_SideWall1_Window);
  
  G4VSolid* solidDetectorTHSW1_sub9 = new G4SubtractionSolid("Wall-Window9",
							     solidDetectorTHSW1_sub8,
							     solidDetectorTHSW1_Window9,
							     0,
							     TVTHSW1_Window
							     );

  G4Tubs* solidDetectorTHSW1_Window10 = new G4Tubs(G4String(fullNameTHSW1_Window.c_str()),
						  RIN_TargetHolder_SideWall1-1,
						  ROUT_TargetHolder_SideWall1+1,
						  DZ_Window,
						  STANG_TargetHolder_SideWall1_Window+180.0*degree + 2*ANG_TargetHolder_SideWall1_Window,
						  EDANG_TargetHolder_SideWall1_Window);
  
  G4VSolid* solidDetectorTHSW1_sub10 = new G4SubtractionSolid("Wall-Window10",
							     solidDetectorTHSW1_sub9,
							     solidDetectorTHSW1_Window10,
							     0,
							     TVTHSW1_Window
							     );
  


  G4LogicalVolume* logicDetectorTHSW1 = new G4LogicalVolume(solidDetectorTHSW1_sub10,
							    mat_G10,
							    G4String(fullNameTHSW1.c_str()));
  logicDetectorTHSW1->SetVisAttributes(new G4VisAttributes(G4Color::Green()));

  G4RotationMatrix rot_THSW1;
  G4ThreeVector TVTHSW1 ( 0, 0, 0.);
  //TVTH.rotateX ( -90.*deg);
  G4VPhysicalVolume *Phys_TargetHolder_SideWall1 = new G4PVPlacement(G4Transform3D(rot_THSW1, TVTHSW1),
							   G4String(fullNameTHSW1.c_str()),
							   logicDetectorTHSW1,
							   physMother,
							   false,
							   0
							   );
  /// ///
 
  /// Target holder SideWall2///
  // param //
  const G4double RIN_TargetHolder_SideWall2= 53.0/2.0*mm;
  const G4double ROUT_TargetHolder_SideWall2 = (53.0/2.0+0.25)*mm;
  //const G4double DZ_TargetHolder_SideWall2 = (DZ_PW*2.0-DZ_PW_DOWN-DZ_Target*2.0)/2.0*mm;
  const G4double DZ_TargetHolder_SideWall2 = (DZ_Target*2.0 + 12.0*2.0)/2.0*mm;
  const G4double STANG_TargetHolder_SideWall2 = 0.*deg;
  const G4double EDANG_TargetHolder_SideWall2 = 360.*deg;
  std::string fullNameTHSW2 = "TPC_TargetHolder_SideWall2";
  // //
  
  G4Tubs* solidDetectorTHSW2 = new G4Tubs(G4String(fullNameTHSW2.c_str()),
					  RIN_TargetHolder_SideWall2,
					  ROUT_TargetHolder_SideWall2,
					  DZ_TargetHolder_SideWall2,
					  STANG_TargetHolder_SideWall2,
					  EDANG_TargetHolder_SideWall2);
  
  G4LogicalVolume* logicDetectorTHSW2 = new G4LogicalVolume(solidDetectorTHSW2,
							    mat_Mylar,
							    G4String(fullNameTHSW2.c_str()));
  logicDetectorTHSW2->SetVisAttributes(new G4VisAttributes(G4Color::Cyan()));

  G4RotationMatrix rot_THSW2;
  //G4ThreeVector TVTHSW2 ( 0, 0, (DZ_PW_DOWN+DZ_Target*2.0)/2.0);
  G4ThreeVector TVTHSW2 ( TVTarget);
  //TVTH.rotateX ( -90.*deg);
  G4VPhysicalVolume *Phys_TargetHolder_SideWall2 = new G4PVPlacement(G4Transform3D(rot_THSW2, TVTHSW2),
							   G4String(fullNameTHSW2.c_str()),
							   logicDetectorTHSW2,
							   physMother,
							   false,
							   0
							   );
  /// ///

  /// Target holder SideWall3///
  // param //
  const G4double RIN_TargetHolder_SideWall3= (18.0/2.0-2.0)*mm;
  const G4double ROUT_TargetHolder_SideWall3 = 18.0/2.0*mm;
  const G4double DZ_TargetHolder_SideWall3 = (DZ_PW*2.0-DZ_PW_DOWN-DZ_Target*2.0)/2.0*mm;
  const G4double STANG_TargetHolder_SideWall3 = 0.*deg;
  const G4double EDANG_TargetHolder_SideWall3 = 360.*deg;
  std::string fullNameTHSW3 = "TPC_TargetHolder_SideWall3";
  // //
  
  G4Tubs* solidDetectorTHSW3 = new G4Tubs(G4String(fullNameTHSW3.c_str()),
					  RIN_TargetHolder_SideWall3,
					  ROUT_TargetHolder_SideWall3,
					  DZ_TargetHolder_SideWall3,
					  STANG_TargetHolder_SideWall3,
					  EDANG_TargetHolder_SideWall3);
  
  G4LogicalVolume* logicDetectorTHSW3 = new G4LogicalVolume(solidDetectorTHSW3,
							    mat_G10,
							    G4String(fullNameTHSW3.c_str()));
  logicDetectorTHSW3->SetVisAttributes(new G4VisAttributes(G4Color::Blue()));

  G4RotationMatrix rot_THSW3;
  G4ThreeVector TVTHSW3 ( 0, 0, (DZ_PW_DOWN+DZ_Target*2.0)/2.0);
  //TVTH.rotateX ( -90.*deg);
  G4VPhysicalVolume *Phys_TargetHolder_SideWall3 = new G4PVPlacement(G4Transform3D(rot_THSW3, TVTHSW3),
							   G4String(fullNameTHSW3.c_str()),
							   logicDetectorTHSW3,
							   physMother,
							   false,
							   0
							   );
  /// ///

  /// Bottom plate1 ///
  // param //
  const G4double RIN_TargetHolder_BottomPlate1= 0.0*mm;
  const G4double ROUT_TargetHolder_BottomPlate1 = (RIN_TargetHolder_SideWall1/mm) *mm;
  const G4double DZ_TargetHolder_BottomPlate1 = 12.0/2.0*mm;
  const G4double STANG_TargetHolder_BottomPlate1 = 0.*deg;
  const G4double EDANG_TargetHolder_BottomPlate1 = 360.*deg;
  std::string fullNameTHBP1 = "TPC_TargetHolder_BottomPlate1";
  // //
  
  G4Tubs* solidDetectorTHBP1 = new G4Tubs(G4String(fullNameTHBP1.c_str()),
					  RIN_TargetHolder_BottomPlate1,
					  ROUT_TargetHolder_BottomPlate1,
					  DZ_TargetHolder_BottomPlate1,
					  STANG_TargetHolder_BottomPlate1,
					  EDANG_TargetHolder_BottomPlate1);
  
  G4LogicalVolume* logicDetectorTHBP1 = new G4LogicalVolume(solidDetectorTHBP1,
							    mat_G10,
							    G4String(fullNameTHBP1.c_str()));
  logicDetectorTHBP1->SetVisAttributes(new G4VisAttributes(G4Color::Black()));

  G4RotationMatrix rot_THBP1;
  G4ThreeVector TVTHBP1 ( 0, 0, -DZ_TargetHolder_SideWall1 + DZ_TargetHolder_BottomPlate1);
  //TVTH.rotateX ( -90.*deg);
  G4VPhysicalVolume *Phys_TargetHolder_BottomPlate1 = new G4PVPlacement(G4Transform3D(rot_THBP1, TVTHBP1),
							   G4String(fullNameTHBP1.c_str()),
							   logicDetectorTHBP1,
							   physMother,
							   false,
							   0
							   );
  /// ///

  /// Bottom plate2 ///
  // param //
  const G4double RIN_TargetHolder_BottomPlate2= 0.0*mm;
  const G4double ROUT_TargetHolder_BottomPlate2 = (RIN_TargetHolder_SideWall2/mm) *mm;
  const G4double DZ_TargetHolder_BottomPlate2 = DZ_TargetHolder_BottomPlate1;
  const G4double STANG_TargetHolder_BottomPlate2 = 0.*deg;
  const G4double EDANG_TargetHolder_BottomPlate2 = 360.*deg;
  std::string fullNameTHBP2 = "TPC_TargetHolder_BottomPlate2";
  // //
  
  G4Tubs* solidDetectorTHBP2 = new G4Tubs(G4String(fullNameTHBP2.c_str()),
					  RIN_TargetHolder_BottomPlate2,
					  ROUT_TargetHolder_BottomPlate2,
					  DZ_TargetHolder_BottomPlate2,
					  STANG_TargetHolder_BottomPlate2,
					  EDANG_TargetHolder_BottomPlate2);
  
  G4LogicalVolume* logicDetectorTHBP2 = new G4LogicalVolume(solidDetectorTHBP2,
							    mat_G10,
							    G4String(fullNameTHBP2.c_str()));
  logicDetectorTHBP2->SetVisAttributes(new G4VisAttributes(G4Color::Black()));

  G4RotationMatrix rot_THBP2;
  G4ThreeVector TVTHBP2 ( 0, 0, -DZ_PW + DZ_PW_DOWN + DZ_Target-DZ_Target - DZ_TargetHolder_BottomPlate2);
  //TVTH.rotateX ( -90.*deg);
  G4VPhysicalVolume *Phys_TargetHolder_BottomPlate2 = new G4PVPlacement(G4Transform3D(rot_THBP2, TVTHBP2),
							   G4String(fullNameTHBP2.c_str()),
							   logicDetectorTHBP2,
							   physMother,
							   false,
							   0
							   );
  /// ///

  /// Bottom plate3 ///
  // param //
  const G4double RIN_TargetHolder_BottomPlate3= 0.0*mm;
  const G4double ROUT_TargetHolder_BottomPlate3 = (RIN_TargetHolder_SideWall2/mm) *mm;
  const G4double DZ_TargetHolder_BottomPlate3 = DZ_TargetHolder_BottomPlate1;
  const G4double STANG_TargetHolder_BottomPlate3 = 0.*deg;
  const G4double EDANG_TargetHolder_BottomPlate3 = 360.*deg;
  std::string fullNameTHBP3 = "TPC_TargetHolder_BottomPlate3";
  // //
  G4Tubs* solidDetectorTHBP3 = new G4Tubs(G4String(fullNameTHBP3.c_str()),
					  RIN_TargetHolder_BottomPlate3,
					  ROUT_TargetHolder_BottomPlate3,
					  DZ_TargetHolder_BottomPlate3,
					  STANG_TargetHolder_BottomPlate3,
					  EDANG_TargetHolder_BottomPlate3);

  G4Tubs* solidDetectorTHBP3_hole = new G4Tubs(G4String(fullNameTHBP3.c_str()),
					       0.0*mm,
					       ROUT_TargetHolder_SideWall3,
					       DZ_TargetHolder_BottomPlate3,
					       STANG_TargetHolder_BottomPlate3,
					       EDANG_TargetHolder_BottomPlate3);


  G4VSolid* solidDetectorTHBP3_sub = new G4SubtractionSolid("Plate-hole",
							    solidDetectorTHBP3,
							    solidDetectorTHBP3_hole
							    );

  
  G4LogicalVolume* logicDetectorTHBP3 = new G4LogicalVolume(solidDetectorTHBP3_sub,
							    mat_G10,
							    G4String(fullNameTHBP3.c_str()));
  logicDetectorTHBP3->SetVisAttributes(new G4VisAttributes(G4Color::Black()));

  G4RotationMatrix rot_THBP3;
  G4ThreeVector TVTHBP3 ( 0, 0, -DZ_PW + DZ_PW_DOWN + DZ_Target + DZ_Target + DZ_TargetHolder_BottomPlate3);
  //TVTH.rotateX ( -90.*deg);
  G4VPhysicalVolume *Phys_TargetHolder_BottomPlate3 = new G4PVPlacement(G4Transform3D(rot_THBP3, TVTHBP3),
							   G4String(fullNameTHBP3.c_str()),
							   logicDetectorTHBP3,
							   physMother,
							   false,
							   0
							   );
  /// ///
    
  ///// Sensitive Detector /////
  G4SDManager *SDMan = G4SDManager::GetSDMpointer();
  TargetSD *targetSD = new TargetSD( "/TPC/Target" );
  SDMan->AddNewDetector( targetSD );
  logicDetectorTarget->SetSensitiveDetector( targetSD );
      
}


void DetectorConstruction::MakeTargetHver2(G4VPhysicalVolume *pMother, G4ThreeVector &pos, G4RotationMatrix &rot)
{
  // material //

  //std::cout<<"Test: " << mat_C <<std::endl;
  /// Physics World ///
  // param //
  const G4double RIN_PW = 0.0*mm;
  const G4double ROUT_PW= 100.0/2.0*mm;
  //const G4double DZ_PW  = 500.0/2.0*mm; // target tpc center (y-axis), up 350 mm, down 68 mm.
  const G4double DZ_PW  = 418.0/2.0*mm; // target tpc center (y-axis), up 350 mm, down 68 mm. 
  //const G4double DZ_PW_DOWN = 47*mm;
  const G4double DZ_PW_DOWN = (118-50)*mm;
  const G4double STANG_PW = 0.*deg;
  const G4double EDANG_PW = 360.*deg;
  std::string fullNamePW = "Target_PhysicsWorld";
  // //
  
  G4Tubs* solidDetectorPW = new G4Tubs(G4String(fullNamePW.c_str()),
				       RIN_PW,
				       ROUT_PW,
				       DZ_PW,
				       STANG_PW,
				       EDANG_PW);
  
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

  /// ///
  
  /// Target ///
  // param //
  const G4double RIN_Target = 0.0*mm;
  const G4double ROUT_Target= 80.0/2.0*mm; // target radius 80 mm
  const G4double DZ_Target  = 100.0/2.0*mm; // target length 100 mm   
  const G4double STANG_Target = 0.*deg;
  const G4double EDANG_Target = 360.*deg;
  std::string fullNameTarget = "TPC_Target";
  // //
  G4Tubs* solidDetectorTarget = new G4Tubs(G4String(fullNameTarget.c_str()),
					   RIN_Target,
					   ROUT_Target,
					   DZ_Target,
					   STANG_Target,
					   EDANG_Target);
  
  G4LogicalVolume* logicDetectorTarget = new G4LogicalVolume(solidDetectorTarget,
							     mat_LH2,
							     G4String(fullNameTarget.c_str()));

  logicDetectorTarget->SetVisAttributes(new G4VisAttributes(G4Color::Red()));

  G4RotationMatrix rot_Target;
  G4ThreeVector TVTarget ( 0, 0, -(DZ_PW - DZ_PW_DOWN - DZ_Target));

  G4VPhysicalVolume *Phys_Target = new G4PVPlacement(G4Transform3D(rot_Target, TVTarget),
						     G4String(fullNameTarget.c_str()),
						     logicDetectorTarget,
						     physMother,
						     false,
						     0
  						     );

  /// Target holder Sidewall1///

  // param //
  const G4double RIN_TargetHolder_SideWall1 = 94.0/2.0*mm;
  const G4double ROUT_TargetHolder_SideWall1= 100.0/2.0*mm;
  //const G4double DZ_TargetHolder_SideWall1 = 500.0/2.0*mm;
  const G4double DZ_TargetHolder_SideWall1 = 418.0/2.0*mm;
  const G4double STANG_TargetHolder_SideWall1 = 0.*deg;
  const G4double EDANG_TargetHolder_SideWall1 = 360.*deg;
  //const G4double EDANG_TargetHolder_SideWall1_Window = 26.0*deg;
  //const G4double STANG_TargetHolder_SideWall1_Window = (-270.0- EDANG_TargetHolder_SideWall1_Window/deg/2.0)*deg;
  //const G4double ANG_TargetHolder_SideWall1_Window = 28.0*deg;
  //const G4double DZ_Window = DZ_TargetHolder_SideWall1+1.0;
  //const G4double DZ_Window = DZ_Target;
  std::string fullNameTHSW1 = "TPC_TargetHolder_SideWall1";
  //std::string fullNameTHSW1_Window = "TPC_TargetHolder_SideWall1_Window";

  //G4ThreeVector TVTHSW1_Window (TVTarget);
  //G4ThreeVector TVTHSW1_Window;
  // //

  // Making window to SideWall1 //
  // Window SideWall1 //
  /*
  G4Tubs* solidDetectorTHSW1_Window1 = new G4Tubs(G4String(fullNameTHSW1_Window.c_str()),
						  RIN_TargetHolder_SideWall1-1,
						  ROUT_TargetHolder_SideWall1+1,
						  DZ_Window,
						  STANG_TargetHolder_SideWall1_Window,
						  EDANG_TargetHolder_SideWall1_Window);
  */
  // //
  G4Tubs* solidDetectorTHSW1 = new G4Tubs(G4String(fullNameTHSW1.c_str()),
					  RIN_TargetHolder_SideWall1,
					  ROUT_TargetHolder_SideWall1,
					  DZ_TargetHolder_SideWall1,
					  STANG_TargetHolder_SideWall1,
					  EDANG_TargetHolder_SideWall1);
  /*
  G4VSolid* solidDetectorTHSW1_sub1 = new G4SubtractionSolid("Wall-Window1",
							    solidDetectorTHSW1,
							    solidDetectorTHSW1_Window1,
							    0,
							    TVTHSW1_Window
							    );
  
  G4Tubs* solidDetectorTHSW1_Window2 = new G4Tubs(G4String(fullNameTHSW1_Window.c_str()),
						  RIN_TargetHolder_SideWall1-1,
						  ROUT_TargetHolder_SideWall1+1,
						  DZ_Window,
						  STANG_TargetHolder_SideWall1_Window-2*ANG_TargetHolder_SideWall1_Window,
						  EDANG_TargetHolder_SideWall1_Window);
  
  G4VSolid* solidDetectorTHSW1_sub2 = new G4SubtractionSolid("Wall-Window2",
							     solidDetectorTHSW1_sub1,
							     solidDetectorTHSW1_Window2,
							     0,
							     TVTHSW1_Window
							     );
  
  G4Tubs* solidDetectorTHSW1_Window3 = new G4Tubs(G4String(fullNameTHSW1_Window.c_str()),
						  RIN_TargetHolder_SideWall1-1,
						  ROUT_TargetHolder_SideWall1+1,
						  DZ_Window,
						  STANG_TargetHolder_SideWall1_Window-ANG_TargetHolder_SideWall1_Window,
						  EDANG_TargetHolder_SideWall1_Window);

  G4VSolid* solidDetectorTHSW1_sub3 = new G4SubtractionSolid("Wall-Window3",
							     solidDetectorTHSW1_sub2,
							     solidDetectorTHSW1_Window3,
							     0,
							     TVTHSW1_Window
							     );

  G4Tubs* solidDetectorTHSW1_Window4 = new G4Tubs(G4String(fullNameTHSW1_Window.c_str()),
						  RIN_TargetHolder_SideWall1-1,
						  ROUT_TargetHolder_SideWall1+1,
						  DZ_Window,
						  STANG_TargetHolder_SideWall1_Window+ ANG_TargetHolder_SideWall1_Window,
						  EDANG_TargetHolder_SideWall1_Window);

  G4VSolid* solidDetectorTHSW1_sub4 = new G4SubtractionSolid("Wall-Window4",
							     solidDetectorTHSW1_sub3,
							     solidDetectorTHSW1_Window4,
							     0,
							     TVTHSW1_Window
							     );
  
  G4Tubs* solidDetectorTHSW1_Window5 = new G4Tubs(G4String(fullNameTHSW1_Window.c_str()),
						  RIN_TargetHolder_SideWall1-1,
						  ROUT_TargetHolder_SideWall1+1,
						  DZ_Window,
						  STANG_TargetHolder_SideWall1_Window + 2*ANG_TargetHolder_SideWall1_Window,
						  EDANG_TargetHolder_SideWall1_Window);
  
  G4VSolid* solidDetectorTHSW1_sub5 = new G4SubtractionSolid("Wall-Window5",
							     solidDetectorTHSW1_sub4,
							     solidDetectorTHSW1_Window5,
							     0,
							     TVTHSW1_Window
							     );
  
  G4Tubs* solidDetectorTHSW1_Window6 = new G4Tubs(G4String(fullNameTHSW1_Window.c_str()),
						  RIN_TargetHolder_SideWall1-1,
						  ROUT_TargetHolder_SideWall1+1,
						  DZ_Window,
						  STANG_TargetHolder_SideWall1_Window+180.0*degree,
						  EDANG_TargetHolder_SideWall1_Window);
  
  G4VSolid* solidDetectorTHSW1_sub6 = new G4SubtractionSolid("Wall-Window6",
							     solidDetectorTHSW1_sub5,
							     solidDetectorTHSW1_Window6,
							     0,
							     TVTHSW1_Window
							     );

  G4Tubs* solidDetectorTHSW1_Window7 = new G4Tubs(G4String(fullNameTHSW1_Window.c_str()),
						  RIN_TargetHolder_SideWall1-1,
						  ROUT_TargetHolder_SideWall1+1,
						  DZ_Window,
						  STANG_TargetHolder_SideWall1_Window+180.0*degree - 2*ANG_TargetHolder_SideWall1_Window,
						  EDANG_TargetHolder_SideWall1_Window);
  
  G4VSolid* solidDetectorTHSW1_sub7 = new G4SubtractionSolid("Wall-Window7",
							     solidDetectorTHSW1_sub6,
							     solidDetectorTHSW1_Window7,
							     0,
							     TVTHSW1_Window
							     );

  G4Tubs* solidDetectorTHSW1_Window8 = new G4Tubs(G4String(fullNameTHSW1_Window.c_str()),
						  RIN_TargetHolder_SideWall1-1,
						  ROUT_TargetHolder_SideWall1+1,
						  DZ_Window,
						  STANG_TargetHolder_SideWall1_Window+180.0*degree - ANG_TargetHolder_SideWall1_Window,
						  EDANG_TargetHolder_SideWall1_Window);
  
  G4VSolid* solidDetectorTHSW1_sub8 = new G4SubtractionSolid("Wall-Window8",
							     solidDetectorTHSW1_sub7,
							     solidDetectorTHSW1_Window8,
							     0,
							     TVTHSW1_Window
							     );

  G4Tubs* solidDetectorTHSW1_Window9 = new G4Tubs(G4String(fullNameTHSW1_Window.c_str()),
						  RIN_TargetHolder_SideWall1-1,
						  ROUT_TargetHolder_SideWall1+1,
						  DZ_Window,
						  STANG_TargetHolder_SideWall1_Window+180.0*degree + ANG_TargetHolder_SideWall1_Window,
						  EDANG_TargetHolder_SideWall1_Window);
  
  G4VSolid* solidDetectorTHSW1_sub9 = new G4SubtractionSolid("Wall-Window9",
							     solidDetectorTHSW1_sub8,
							     solidDetectorTHSW1_Window9,
							     0,
							     TVTHSW1_Window
							     );

  G4Tubs* solidDetectorTHSW1_Window10 = new G4Tubs(G4String(fullNameTHSW1_Window.c_str()),
						  RIN_TargetHolder_SideWall1-1,
						  ROUT_TargetHolder_SideWall1+1,
						  DZ_Window,
						  STANG_TargetHolder_SideWall1_Window+180.0*degree + 2*ANG_TargetHolder_SideWall1_Window,
						  EDANG_TargetHolder_SideWall1_Window);
  
  G4VSolid* solidDetectorTHSW1_sub10 = new G4SubtractionSolid("Wall-Window10",
							     solidDetectorTHSW1_sub9,
							     solidDetectorTHSW1_Window10,
							     0,
							     TVTHSW1_Window
							     );
  */  


  G4LogicalVolume* logicDetectorTHSW1 = new G4LogicalVolume(solidDetectorTHSW1,
							    mat_G10,
							    G4String(fullNameTHSW1.c_str()));
  logicDetectorTHSW1->SetVisAttributes(new G4VisAttributes(G4Color::Green()));

  G4RotationMatrix rot_THSW1;
  G4ThreeVector TVTHSW1 ( 0, 0, 0.);
  //TVTH.rotateX ( -90.*deg);
  G4VPhysicalVolume *Phys_TargetHolder_SideWall1 = new G4PVPlacement(G4Transform3D(rot_THSW1, TVTHSW1),
							   G4String(fullNameTHSW1.c_str()),
							   logicDetectorTHSW1,
							   physMother,
							   false,
							   0
							   );
  /// ///
 
  /// Target holder SideWall2///
  // param //
  const G4double RIN_TargetHolder_SideWall2= 100.0/2.0*mm;
  const G4double ROUT_TargetHolder_SideWall2 = (100.0/2.0+0.25)*mm;
  //const G4double DZ_TargetHolder_SideWall2 = (DZ_PW*2.0-DZ_PW_DOWN-DZ_Target*2.0)/2.0*mm;
  const G4double DZ_TargetHolder_SideWall2 = (DZ_Target*2.0 + 12.0*2.0)/2.0*mm;
  const G4double STANG_TargetHolder_SideWall2 = 0.*deg;
  const G4double EDANG_TargetHolder_SideWall2 = 360.*deg;
  std::string fullNameTHSW2 = "TPC_TargetHolder_SideWall2";
  // //
  
  G4Tubs* solidDetectorTHSW2 = new G4Tubs(G4String(fullNameTHSW2.c_str()),
					  RIN_TargetHolder_SideWall2,
					  ROUT_TargetHolder_SideWall2,
					  DZ_TargetHolder_SideWall2,
					  STANG_TargetHolder_SideWall2,
					  EDANG_TargetHolder_SideWall2);
  
  G4LogicalVolume* logicDetectorTHSW2 = new G4LogicalVolume(solidDetectorTHSW2,
							    mat_Mylar,
							    G4String(fullNameTHSW2.c_str()));
  logicDetectorTHSW2->SetVisAttributes(new G4VisAttributes(G4Color::Cyan()));

  G4RotationMatrix rot_THSW2;
  //G4ThreeVector TVTHSW2 ( 0, 0, (DZ_PW_DOWN+DZ_Target*2.0)/2.0);
  G4ThreeVector TVTHSW2 ( TVTarget);
  //TVTH.rotateX ( -90.*deg);
  G4VPhysicalVolume *Phys_TargetHolder_SideWall2 = new G4PVPlacement(G4Transform3D(rot_THSW2, TVTHSW2),
							   G4String(fullNameTHSW2.c_str()),
							   logicDetectorTHSW2,
							   physMother,
							   false,
							   0
							   );
  /// ///

  /// Target holder SideWall3///
  // param //
  const G4double RIN_TargetHolder_SideWall3= 6.0/2.0*mm;
  const G4double ROUT_TargetHolder_SideWall3 = 9.0/2.0*mm;
  const G4double DZ_TargetHolder_SideWall3 = (DZ_PW*2.0-DZ_PW_DOWN-DZ_Target*2.0)/2.0*mm;
  const G4double STANG_TargetHolder_SideWall3 = 0.*deg;
  const G4double EDANG_TargetHolder_SideWall3 = 360.*deg;
  std::string fullNameTHSW3 = "TPC_TargetHolder_SideWall3";
  // //
  
  G4Tubs* solidDetectorTHSW3 = new G4Tubs(G4String(fullNameTHSW3.c_str()),
					  RIN_TargetHolder_SideWall3,
					  ROUT_TargetHolder_SideWall3,
					  DZ_TargetHolder_SideWall3,
					  STANG_TargetHolder_SideWall3,
					  EDANG_TargetHolder_SideWall3);
  
  G4LogicalVolume* logicDetectorTHSW3 = new G4LogicalVolume(solidDetectorTHSW3,
							    mat_G10,
							    G4String(fullNameTHSW3.c_str()));
  logicDetectorTHSW3->SetVisAttributes(new G4VisAttributes(G4Color::Blue()));

  G4RotationMatrix rot_THSW3;
  G4ThreeVector TVTHSW3 ( 0, 0, (DZ_PW_DOWN+DZ_Target*2.0)/2.0);
  //TVTH.rotateX ( -90.*deg);
  G4VPhysicalVolume *Phys_TargetHolder_SideWall3 = new G4PVPlacement(G4Transform3D(rot_THSW3, TVTHSW3),
							   G4String(fullNameTHSW3.c_str()),
							   logicDetectorTHSW3,
							   physMother,
							   false,
							   0
							   );
  /// ///

  /// Bottom plate1 ///
  // param //
  const G4double RIN_TargetHolder_BottomPlate1= 0.0*mm;
  const G4double ROUT_TargetHolder_BottomPlate1 = (RIN_TargetHolder_SideWall1/mm) *mm;
  const G4double DZ_TargetHolder_BottomPlate1 = 12.0/2.0*mm;
  const G4double STANG_TargetHolder_BottomPlate1 = 0.*deg;
  const G4double EDANG_TargetHolder_BottomPlate1 = 360.*deg;
  std::string fullNameTHBP1 = "TPC_TargetHolder_BottomPlate1";
  // //
  
  G4Tubs* solidDetectorTHBP1 = new G4Tubs(G4String(fullNameTHBP1.c_str()),
					  RIN_TargetHolder_BottomPlate1,
					  ROUT_TargetHolder_BottomPlate1,
					  DZ_TargetHolder_BottomPlate1,
					  STANG_TargetHolder_BottomPlate1,
					  EDANG_TargetHolder_BottomPlate1);
  
  G4LogicalVolume* logicDetectorTHBP1 = new G4LogicalVolume(solidDetectorTHBP1,
							    mat_G10,
							    G4String(fullNameTHBP1.c_str()));
  logicDetectorTHBP1->SetVisAttributes(new G4VisAttributes(G4Color::Black()));

  G4RotationMatrix rot_THBP1;
  G4ThreeVector TVTHBP1 ( 0, 0, -DZ_TargetHolder_SideWall1 + DZ_TargetHolder_BottomPlate1);
  //TVTH.rotateX ( -90.*deg);
  G4VPhysicalVolume *Phys_TargetHolder_BottomPlate1 = new G4PVPlacement(G4Transform3D(rot_THBP1, TVTHBP1),
							   G4String(fullNameTHBP1.c_str()),
							   logicDetectorTHBP1,
							   physMother,
							   false,
							   0
							   );
  /// ///

  /// Bottom plate2 ///
  // param //
  const G4double RIN_TargetHolder_BottomPlate2= 0.0*mm;
  const G4double ROUT_TargetHolder_BottomPlate2 = (RIN_TargetHolder_SideWall2/mm) *mm;
  const G4double DZ_TargetHolder_BottomPlate2 = DZ_TargetHolder_BottomPlate1;
  const G4double STANG_TargetHolder_BottomPlate2 = 0.*deg;
  const G4double EDANG_TargetHolder_BottomPlate2 = 360.*deg;
  std::string fullNameTHBP2 = "TPC_TargetHolder_BottomPlate2";
  // //
  
  G4Tubs* solidDetectorTHBP2 = new G4Tubs(G4String(fullNameTHBP2.c_str()),
					  RIN_TargetHolder_BottomPlate2,
					  ROUT_TargetHolder_BottomPlate2,
					  DZ_TargetHolder_BottomPlate2,
					  STANG_TargetHolder_BottomPlate2,
					  EDANG_TargetHolder_BottomPlate2);
  
  G4LogicalVolume* logicDetectorTHBP2 = new G4LogicalVolume(solidDetectorTHBP2,
							    mat_G10,
							    G4String(fullNameTHBP2.c_str()));
  logicDetectorTHBP2->SetVisAttributes(new G4VisAttributes(G4Color::Black()));

  G4RotationMatrix rot_THBP2;
  G4ThreeVector TVTHBP2 ( 0, 0, -DZ_PW + DZ_PW_DOWN + DZ_Target-DZ_Target - DZ_TargetHolder_BottomPlate2);
  //TVTH.rotateX ( -90.*deg);
  G4VPhysicalVolume *Phys_TargetHolder_BottomPlate2 = new G4PVPlacement(G4Transform3D(rot_THBP2, TVTHBP2),
							   G4String(fullNameTHBP2.c_str()),
							   logicDetectorTHBP2,
							   physMother,
							   false,
							   0
							   );
  /// ///

  /// Bottom plate3 ///
  // param //
  const G4double RIN_TargetHolder_BottomPlate3= 0.0*mm;
  const G4double ROUT_TargetHolder_BottomPlate3 = (RIN_TargetHolder_SideWall2/mm) *mm;
  const G4double DZ_TargetHolder_BottomPlate3 = DZ_TargetHolder_BottomPlate1;
  const G4double STANG_TargetHolder_BottomPlate3 = 0.*deg;
  const G4double EDANG_TargetHolder_BottomPlate3 = 360.*deg;
  std::string fullNameTHBP3 = "TPC_TargetHolder_BottomPlate3";
  // //
  G4Tubs* solidDetectorTHBP3 = new G4Tubs(G4String(fullNameTHBP3.c_str()),
					  RIN_TargetHolder_BottomPlate3,
					  ROUT_TargetHolder_BottomPlate3,
					  DZ_TargetHolder_BottomPlate3,
					  STANG_TargetHolder_BottomPlate3,
					  EDANG_TargetHolder_BottomPlate3);

  G4Tubs* solidDetectorTHBP3_hole = new G4Tubs(G4String(fullNameTHBP3.c_str()),
					       0.0*mm,
					       ROUT_TargetHolder_SideWall3,
					       DZ_TargetHolder_BottomPlate3,
					       STANG_TargetHolder_BottomPlate3,
					       EDANG_TargetHolder_BottomPlate3);


  G4VSolid* solidDetectorTHBP3_sub = new G4SubtractionSolid("Plate-hole",
							    solidDetectorTHBP3,
							    solidDetectorTHBP3_hole
							    );

  
  G4LogicalVolume* logicDetectorTHBP3 = new G4LogicalVolume(solidDetectorTHBP3_sub,
							    mat_G10,
							    G4String(fullNameTHBP3.c_str()));
  logicDetectorTHBP3->SetVisAttributes(new G4VisAttributes(G4Color::Black()));

  G4RotationMatrix rot_THBP3;
  G4ThreeVector TVTHBP3 ( 0, 0, -DZ_PW + DZ_PW_DOWN + DZ_Target + DZ_Target + DZ_TargetHolder_BottomPlate3);
  //TVTH.rotateX ( -90.*deg);
  G4VPhysicalVolume *Phys_TargetHolder_BottomPlate3 = new G4PVPlacement(G4Transform3D(rot_THBP3, TVTHBP3),
							   G4String(fullNameTHBP3.c_str()),
							   logicDetectorTHBP3,
							   physMother,
							   false,
							   0
							   );
  /// ///
    
  ///// Sensitive Detector /////
  G4SDManager *SDMan = G4SDManager::GetSDMpointer();
  TargetSD *targetSD = new TargetSD( "/TPC/Target" );
  SDMan->AddNewDetector( targetSD );
  logicDetectorTarget->SetSensitiveDetector( targetSD );
      
}



void DetectorConstruction::MakeHypTPC2(G4VPhysicalVolume *pMother, G4ThreeVector &pos, G4RotationMatrix &rot)
{
  /// Target hole ///
  // param //
  const G4double RIN_TH = 0.0*mm;
  //const G4double ROUT_TH= 80.0/2.0*mm;
  const G4double ROUT_TH= 100.0/2.0*mm;  //Change the Target Size
  const G4double DZ_TH  = (418.0)/2.0*mm; // target tpc center (y-axis), up 350 mm, down 68 mm. 
  const G4double STANG_TH = 0.*deg;
  const G4double EDANG_TH = 360.*deg;
  std::string fullNameTH = "Target_PhysicsWorld";
  // //
  
  G4Tubs* solidDetectorTH = new G4Tubs(G4String(fullNameTH.c_str()),
				       RIN_TH,
				       ROUT_TH,
				       DZ_TH,
				       STANG_TH,
				       EDANG_TH);


  /// Physics world ///
  // param //
  const G4double STANG_PW     = 22.5*deg;
  const G4double DPHI_PW      = 360.*deg;
  const G4double NPL_PW       = 8.;
  const G4double RIN_PW       = 0.0/2.*mm;
  const G4double ROUT_PW      = (656./2.)*mm;
  const G4double DZ_PW        =  802.*mm;
  const G4double DZ_PW_OFFSET =    0.*mm;
  std::string fullNamePW = "TPC_PhysiWorld";
  // //
  double z_PW[2] = {-DZ_PW/2.,DZ_PW/2.};
  double rmin_PW[2] = {RIN_PW, RIN_PW};
  double rmax_PW[2] = {ROUT_PW, ROUT_PW};
  G4Polyhedra* solidDetectorPW = new G4Polyhedra(G4String(fullNamePW.c_str()),
						 STANG_PW,
						 DPHI_PW,
						 int(NPL_PW),
						 2,
						 z_PW,
						 rmin_PW,
						 rmax_PW);

  G4ThreeVector TVTH(0.0, 143.0, 91.5);
  TVTH.rotateZ(-90.*deg);
  G4VSolid* solidDetectorPW_sub = new G4SubtractionSolid("hedra-targethall",
							 solidDetectorPW,
							 solidDetectorTH,
							 0,
							 TVTH
							 );
  
  G4LogicalVolume* logicDetectorPW = new G4LogicalVolume(solidDetectorPW_sub,
							 mat_air,
							 G4String(fullNamePW.c_str()));

  logicDetectorPW->SetVisAttributes(new G4VisAttributes(G4Color::Green()));

  // position //
  G4RotationMatrix rot_PW = rot;
  G4ThreeVector pos_PW = rot*pos;

  G4VPhysicalVolume *phyMother = new G4PVPlacement(G4Transform3D(rot_PW, pos_PW),
						   G4String(fullNamePW.c_str()),
						   logicDetectorPW,
						   pMother,
						   false,
						   0);

  /// TPC fieldcage ///
  // param //
  const G4double STANG_FC     = 22.5*deg;
  const G4double DPHI_FC      = 360.*deg;
  const G4double NPL_FC       = 8.;
  const G4double RIN_FC       = 574.5/2.*mm;
  const G4double ROUT_FC      = (654./2.)*mm;
  const G4double DZ_FC        =  800.*mm;
  const G4double DZ_FC_OFFSET =    0.*mm;
  std::string fullNameFC = "TPC_FieldCage";
  // //

  double z_FC[2] = {-DZ_FC/2.,DZ_FC/2.};
  double rmin_FC[2] = {RIN_FC, RIN_FC};
  double rmax_FC[2] = {ROUT_FC, ROUT_FC};
  G4Polyhedra* solidDetectorFC = new G4Polyhedra(G4String(fullNameFC.c_str()),
						 STANG_FC,
						 DPHI_FC,
						 int(NPL_FC),
						 2,
						 z_FC,
						 rmin_FC,
						 rmax_FC);

  G4VSolid* solidDetectorFC_sub = new G4SubtractionSolid("hedra-targethall",
							 solidDetectorFC,
							 solidDetectorTH,
							 0,
							 TVTH
							 );
  
  G4LogicalVolume* logicDetectorFC = new G4LogicalVolume(solidDetectorFC_sub,
							 mat_p10,
							 G4String(fullNameFC.c_str()));
  logicDetectorFC->SetVisAttributes(new G4VisAttributes(G4Color::Magenta()));
  G4RotationMatrix rot_FieldCage;
  //rot_FieldCage.rotateX(-90.*deg);
  G4VPhysicalVolume *Phys_TPC_Field_Cage = new G4PVPlacement(G4Transform3D(rot_FieldCage, G4ThreeVector(0, 0, 0)),
							     G4String(fullNameFC.c_str()),
							     logicDetectorFC,
							     phyMother,
							     false,
							     0);
  /// ///
  
  /// TPC frame ///
  // param //
  const G4double STANG_TPC     = 22.5*deg;
  const G4double DPHI_TPC      = 360.*deg;
  const G4double DZ_TPC        = 604.*mm;
  const G4double NPL_TPC       = 8.;
  const G4double RIN_TPC       = 0.0*mm;
  const G4double ROUT_TPC      = 574./2.*mm;
  const G4double DZ_TPC_OFFSET = 0.*mm;
  std::string fullNameTPC = "TPC_Frame";
  
  double z_TPC[2] = {-DZ_TPC/2.,DZ_TPC/2.};
  double rmin_TPC[2] = {RIN_TPC, RIN_TPC};
  double rmax_TPC[2] = {ROUT_TPC, ROUT_TPC};
  G4Polyhedra* solidDetectorTPC = new G4Polyhedra(G4String(fullNameTPC.c_str()),
						  STANG_TPC,
						  DPHI_TPC,
						  int(NPL_TPC),
						  2,
						  z_TPC,
						  rmin_TPC,
						  rmax_TPC);

  G4VSolid* solidDetectorTPC_sub = new G4SubtractionSolid("hedra-targethall",
							 solidDetectorTPC,
							 solidDetectorTH,
							 0,
							 TVTH
							 );
  
  G4LogicalVolume* logicDetectorTPC = new G4LogicalVolume(solidDetectorTPC_sub,
							  mat_p10,
							  G4String(fullNameTPC.c_str()));
  logicDetectorTPC->SetVisAttributes(new G4VisAttributes(G4Color::Green()));
  G4RotationMatrix rot_TPC;
  //rot_TPC.rotateX(-90.*deg);
  G4VPhysicalVolume *Phys_TPC = new G4PVPlacement(G4Transform3D(rot_TPC, G4ThreeVector(0, 0, 0)),
						  G4String(fullNameTPC.c_str()),
						  logicDetectorTPC,
						  phyMother,
						  false,
						  0);
  /// ///

  /// TPC PAD ///
  // Parameter //
  // PadParameter :                                                                                           
  // [0]: RingID                                                                                              
  // [1]: Number of pads in a ring                                                                            
  // [2]: Center radius                                                                                       
  // [3]: dTheta of a pad                                                                                     
  // [4]: Start angle of ring                                                                                 
  // [5]: length of a pad
  
  //Original
  /*
  const double PadParameter[32][6] = {                                                                        
    {0,         48,     14.5,   7.5     ,0.     ,9.},                                                         
    {1,         48,     24.,    7.5     ,0.     ,9.},                                                         
    {2,         72,     33.5,   5.      ,0.     ,9.},                                                         
    {3,         96,     43.,    3.75    ,0.     ,9.},                                                         
    {4,         120,    52.5,   3.      ,0.     ,9.},                                                         
    {5,         144,    62.,    2.5     ,0.     ,9.},                                                         
    {6,         168,    71.5,   2.14286 ,0.     ,9.},                                                         
    {7,         192,    81.,    1.875   ,0.     ,9.},                                                         
    {8,         216,    90.5,   1.66667 ,0.     ,9.},                                                         
    {9,         240,    100.,   1.5     ,0.     ,9.},                                                         
    {10,        208,    111.25, 1.49375 ,24.65  ,12.5},                                                       
    {11,        218,    124.25, 1.32844 ,35.2   ,12.5},                                                       
    {12,        230,    137.25, 1.2     ,42     ,12.5},                                                       
    {13,        214,    150.25, 1.09093 ,63.27  ,12.5},                                                       
    {14,        212,    163.25, 1.      ,74     ,12.5},                                                       
    {15,        214,    176.25, 0.923084,       81.23   ,12.5},                                               
    {16,        220,    189.25, 0.857182,       85.71   ,12.5},                                               
    {17,        224,    202.25, 0.801786,       90.2    ,12.5},                                               
    {18,        232,    215.25, 0.751552,       92.82   ,12.5},                                               
    {19,        238,    228.25, 0.707227,       95.84   ,12.5},                                               
    {20,        244,    241.25, 0.667869,       98.52   ,12.5},                                               
    {21,        232,    254.25, 0.632672,       106.61  ,12.5},                                               
    {22,        218,    267.25, 0.60101 ,       114.49  ,12.5},                                               
    {23,        210,    280.25, 0.573238,       119.81  ,12.5},                                               
    {24,        206,    293.25, 0.547111,       123.648 ,12.5},                                               
    {25,        202,    306.25, 0.523267,       127.15  ,12.5},                                               
    {26,        200,    319.25, 0.5014  ,       129.86  ,12.5},                                               
    {27,        196,    332.25, 0.481327,       132.83  ,12.5},                                               
    {28,        178,    345.25, 0.463371,       138.76  ,12.5},                                     
    {29,        130,    358.25, 0.446154,       151     ,12.5},                                               
    {30,        108,    371.25, 0.430185,       156.77  ,12.5},                                               
    {31,        90,     384.25, 0.415333,       161.31  ,12.5}};              
  */
  //Modified by Haein 221110
  const double PadParameter[32][6] = {                                                                        
    {0,         48,     14.75,   48     ,  0.      ,9.},                                                         
    {1,         48,     24.25,   48     ,  0.      ,9.},                                                         
    {2,         72,     33.75,   72      ,  0.      ,9.},                                                         
    {3,         96,     43.25,   96    ,  0.      ,9.},                                                         
    {4,         120,    52.75,   120      ,  0.      ,9.},                                                         
    {5,         144,    62.25,   144     ,  0.      ,9.},                                                         
    {6,         168,    71.75,   168 ,  0.      ,9.},                                                         
    {7,         192,    81.25,   192   ,  0.      ,9.},                                                         
    {8,         216,    90.75,   216 ,  0.      ,9.},                                                         
    {9,         240,    100.25,  240     ,  0.      ,9.},                                                         
    {10,        208,    111.5,   241 ,  24.65   ,12.5},                                                       
    {11,        218,    124.5,   271 ,  35.2    ,12.5},                                                       
    {12,        230,    137.5,   300     ,  42      ,12.5},                                                       
    {13,        214,    150.5,   330 ,  63.27   ,12.5},                                                       
    {14,        212,    163.5,   360     ,  74      ,12.5},                                                       
    {15,        214,    176.5,   390,  81.23   ,12.5},                                               
    {16,        220,    189.5,   420,  85.71   ,12.5},                                               
    {17,        224,    202.5,   449,  90.2    ,12.5},                                               
    {18,        232,    215.5,   479,  92.82   ,12.5},                                               
    {19,        238,    228.5,   509,  95.84   ,12.5},                                               
    {20,        244,    241.5,   539,  98.52   ,12.5},                                               
    {21,        232,    254.5,   569,  106.61  ,12.5},                                               
    {22,        218,    267.5,   599 ,  114.49  ,12.5},                                               
    {23,        210,    280.5,   628,  119.81  ,12.5},                                               
    {24,        206,    293.5,   658,  123.648 ,12.5},                                               
    {25,        202,    306.5,   688,  127.15  ,12.5},                                               
    {26,        200,    319.5,   718  ,  129.86  ,12.5},                                               
    {27,        196,    332.5,   748,  132.83  ,12.5},                                               
    {28,        178,    345.5,   777,  138.76  ,12.5},                                     
    {29,        130,    358.5,   807,  151     ,12.5},                                               
    {30,        108,    371.5,   837,  156.77  ,12.5},                                               
    {31,        90,     384.5,   867,  161.31  ,12.5}};              
  
  
  //const G4double Rad_out = 33.6*mm; //for C target
  const G4double Rad_out = 100.0/2.0*mm; //for LH2 target -> Target itself size? or contain frame of target?                                                                         
  //const G4double Rad_out = 82.0/2.0*mm; //for LH2 target                                                                          
  //const G4double DZ_pad_under_target = 240*mm;                                                                
  const G4double DZ_pad_under_target = (275-118)*mm; //for LH2 target                                                                
  const G4double DZ_pad_normal_target= 550*mm;                                                                
  //const G4double PZ_pad_under_target = 155*mm; //for C target
  const G4double PZ_pad_under_target = DZ_pad_normal_target/2.0 - DZ_pad_under_target/2.0; //for LH2 target
  const G4double PZ_pad_normal_target= 0*mm;                                                                    
  const G4double tpc_centerOffset = 143.*mm;
  std::string fullNameTPCpad = "TPC_Pad";
  // //
  
  int padID = 0;                                                                                              
  G4VPhysicalVolume *tpc_ring[32][300];
  G4Tubs *solidTube[32][300];
  G4LogicalVolume *logicDetectorTPCpad[32][300];

  ///// Sensitive Detector /////
  G4SDManager *SDMan = G4SDManager::GetSDMpointer();
  TpcSD *tpcSD = new TpcSD( "/TPC/Tpc" );
  SDMan->AddNewDetector( tpcSD );
  //logicDetectorTPC->SetSensitiveDetector( tpcSD );      


  
  for( int i = 0; i< 32; i++)
    {
      padID = i;                                                                                                
      G4double RMIN_pad = (PadParameter[i][2] - PadParameter[i][5]/2.)*mm;                                      
      G4double RMAX_pad = (PadParameter[i][2] + PadParameter[i][5]/2.)*mm;                                      
      G4double PZ_TPCPAD= 0;                                                                                    
      G4double DZ_TPCPAD= 0;                                                                                    
      if( RMIN_pad < Rad_out )
	{                                                                                 
	  DZ_TPCPAD = DZ_pad_under_target;                                                                        
	  PZ_TPCPAD = PZ_pad_under_target;                                                                        
	}
      else
	{                                                                                                    
	  DZ_TPCPAD = DZ_pad_normal_target;                                                                       
	  PZ_TPCPAD = PZ_pad_normal_target;                                                                       
	}                                                                                                         

      Double_t sTheta = 180.-(360./PadParameter[i][3])*PadParameter[i][1]/2;
      Double_t dTheta = (360./PadParameter[i][3]);
      for(int j = 0; j<(int) PadParameter[i][1] ;j++){

	
	if( i < 10 )
	{                                                                                             
	solidTube[i][j] = new G4Tubs(G4String((fullNameTPCpad+Form("_%d_%d",i,j)).c_str()), 
				     RMIN_pad,
				     RMAX_pad,
				     DZ_TPCPAD/2.,
				     (sTheta+j*dTheta)*deg,
				     dTheta*deg);                  
	 }
	
	else
	  {

	    solidTube[i][j] = new G4Tubs(G4String((fullNameTPCpad+Form("_%d_%d",i,j)).c_str()), 
				     RMIN_pad,
				     RMAX_pad,
				     DZ_TPCPAD/2.,
				     (sTheta+j*dTheta+90)*deg,
				     dTheta*deg);
	    /*
	    solidTube[i][j] = new G4Tubs(G4String((fullNameTPCpad+Form("_%d",i)).c_str()), 
				      RMIN_pad,
				      RMAX_pad,
				      DZ_TPCPAD/2.,
				      (PadParameter[i][4]+90)*deg,
				      (PadParameter[i][3]*PadParameter[i][1])*deg);                           
	    */
	  } 

      
	logicDetectorTPCpad[i][j] = new G4LogicalVolume(solidTube[i][j], 
							mat_p10, 
							G4String((fullNameTPCpad+Form("_%d_%d",i,j)).c_str()));
      
	logicDetectorTPCpad[i][j]->SetVisAttributes(new G4VisAttributes(G4Color::Yellow()));                             
      
	G4RotationMatrix rot_TPCpad;
	rot_TPCpad.rotateZ( -90.*deg );
	G4ThreeVector TVTPCpad ( 0, tpc_centerOffset, -PZ_TPCPAD);
	TVTPCpad.rotateZ ( -90.*deg);
	tpc_ring[i][j] = new G4PVPlacement(G4Transform3D(rot_TPCpad, TVTPCpad),
					   G4String((fullNameTPCpad+Form("_%d_%d",i,j)).c_str()),
					   logicDetectorTPCpad[i][j],
					   Phys_TPC,
					   false,
					   i*1000+j);
	logicDetectorTPCpad[i][j]->SetSensitiveDetector( tpcSD );      
      }
    }

  
}


void DetectorConstruction::MakeTargetDummy(G4VPhysicalVolume *pMother, G4ThreeVector &pos, G4RotationMatrix &rot)
{
  // material //

  //std::cout<<"Test: " << mat_C <<std::endl;
  /// Physics World ///
  // param //
  const G4double X_PW = 50.0/2.0*mm;
  const G4double Y_PW= 50.0/2.0*mm;
  const G4double Z_PW  = 0.01/2.0*mm; 
  std::string fullNamePW = "Target_PhysicsWorld";
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

  /// ///
  
  /// Target ///
  // param //
  const G4double X_Target = 50.0/2.0*mm;
  const G4double Y_Target= 50.0/2.0*mm; // target radius 53 mm
  const G4double Z_Target  = 0.01/2.0*mm; // target length 100 mm   
  std::string fullNameTarget = "TPC_Target";
  // //
  G4Box* solidDetectorTarget = new G4Box(G4String(fullNameTarget.c_str()),
					 X_Target,
					 Y_Target,
					 Z_Target
					 );
  
  G4LogicalVolume* logicDetectorTarget = new G4LogicalVolume(solidDetectorTarget,
							     mat_air,
							     G4String(fullNameTarget.c_str()));
  
  logicDetectorTarget->SetVisAttributes(new G4VisAttributes(G4Color::Red()));
  
  G4RotationMatrix rot_Target;
  G4ThreeVector TVTarget ( 0, 0, 0);
  
  G4VPhysicalVolume *Phys_Target = new G4PVPlacement(G4Transform3D(rot_Target, TVTarget),
						     G4String(fullNameTarget.c_str()),
						     logicDetectorTarget,
						     physMother,
						     false,
						     0
						     );
  
  /// ///

  ///// Sensitive Detector /////
  G4SDManager *SDMan = G4SDManager::GetSDMpointer();
  TargetSD *targetSD = new TargetSD( "/TPC/Target" );
  SDMan->AddNewDetector( targetSD );
  logicDetectorTarget->SetSensitiveDetector( targetSD );
      
}


void DetectorConstruction::MakeFAC(G4VPhysicalVolume *pMother, G4ThreeVector &pos, G4RotationMatrix &rot)
{
  /// FAC Physics World ///
  /// Physics World ///
  // param //
  //const G4double X_PW = 502.0/2.0*mm;
  const G4double X_PW = 1002.0/2.0*mm;
  const G4double Y_PW=  402.0/2.0*mm;
  const G4double Z_PW  = 12.0/2.0*mm; 
  std::string fullNamePW = "FAC_PhysicsWorld";
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
  
  // FAC //
  G4double maxStep = 0.00001*mm;
  //G4Box* FACSolid = new G4Box("FAC", 501./2*mm, 201./2*mm, 11./2*mm);
  G4Box* FACSolid = new G4Box("FAC", 1001./2*mm, 401./2*mm, 11./2*mm);
  G4LogicalVolume* FACLV= new G4LogicalVolume(FACSolid, mat_air, "FACLV");
  //FACLV -> SetUserLimits(new G4UserLimits(maxStep));
  G4RotationMatrix dummy_rot;
  G4VPhysicalVolume *FACPV = new G4PVPlacement(G4Transform3D(dummy_rot, G4ThreeVector(0,0,0)),
					       "FACPV",
					       FACLV, 
					       physMother, 
					       false, 
					       0);
  G4VisAttributes* FACVisAtt= new G4VisAttributes(true, G4Colour(0.,1.,0.));
  FACLV-> SetVisAttributes(FACVisAtt);

  // FAC Aerogel //
  //G4Box* FACAeroSolid = new G4Box("FACAero", 500./2*mm,200./2*mm,10./2*mm);
  G4Box* FACAeroSolid = new G4Box("FACAero", 1000./2*mm,400./2*mm,10./2*mm);
  G4LogicalVolume* FACAeroLV = new G4LogicalVolume(FACAeroSolid, mat_air, "FACAeroLV");
  //FACAeroLV -> SetUserLimits(new G4UserLimits(maxStep));
  G4VPhysicalVolume *FACAeroPV = new G4PVPlacement(G4Transform3D(dummy_rot, G4ThreeVector(0,0,0.*mm)),
						   "FACAeroPV",
						   FACAeroLV, 
						   FACPV,
						   false, 
						   1);

  G4VisAttributes* FACAeroVisAtt= new G4VisAttributes(true, G4Color::Blue());
  FACAeroLV-> SetVisAttributes(FACAeroVisAtt);

  // sensitive detector //
  G4SDManager *SDMan = G4SDManager::GetSDMpointer();
  FACSD* facSD = new FACSD("/FAC");
  SDMan -> AddNewDetector( facSD );
  FACAeroLV -> SetSensitiveDetector( facSD );


}



void DetectorConstruction::MakeBDC(G4VPhysicalVolume *pMother, G4ThreeVector &pos, G4RotationMatrix &rot)
{
  /// Dummy Detector Physics World ///
  /// Physics World ///
  // param //
  const G4double X_PW = 1000.1/2.0*mm;
  const G4double Y_PW=  400.1/2.0*mm;
  const G4double Z_PW  = 3.1/2.0*mm; 
  std::string fullNamePW = "Bdc_Detector_PhysicsWorld";
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

  // Bdc Detector //                                                                                     
  G4Box* bdc_box = new G4Box("bdc_box", X_PW, Y_PW, Z_PW);
  G4LogicalVolume*  bdc_log = new G4LogicalVolume(bdc_box, mat_air, "bdc_log", 0, 0, 0);
  bdc_log->SetVisAttributes(G4Color::Yellow());

  G4RotationMatrix bdc_rot;
  G4VPhysicalVolume* bdc_phys = new G4PVPlacement(G4Transform3D(bdc_rot,
								  G4ThreeVector(0., 0., 0.)),
						    "bdc_phys",
						    bdc_log,
						    physMother,
						    false,
						    100);

  // sensitive detector //
  G4SDManager *SDMan = G4SDManager::GetSDMpointer();
  BDCSD* bdcSD = new BDCSD("/BDC");
  SDMan -> AddNewDetector( bdcSD );
  bdc_log->SetSensitiveDetector( bdcSD );

}
