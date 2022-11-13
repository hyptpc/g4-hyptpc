/*
  DetectorConstruction.cc

  2017/8  Yang
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

DetectorConstruction::DetectorConstruction()
  : mList_(0), field_(0)
{
  NistMan = G4NistManager::Instance();
}

DetectorConstruction::~DetectorConstruction()
{
  delete mList_;
  delete field_;
}
MaterialList *DetectorConstruction::DefineMaterials()
{
  if(mList_) delete mList_;
  return new MaterialList();
}

G4MagneticField * DetectorConstruction::MakeUniformMagField( G4double fieldValue )
{
  G4MagneticField *field=0;
  G4FieldManager* fieldMan = G4TransportationManager::GetTransportationManager()->GetFieldManager();
  if( fieldValue!=0. )
    { 
      field = new G4UniformMagField( G4ThreeVector( 0., fieldValue , 0. ));
      fieldMan->CreateChordFinder( field );
    }
  fieldMan->SetDetectorField( field );
  //fieldMan->GetChordFinder()->SetDeltaChord( 1.0E-3*mm );
  G4cout<<"Field is uniform field"<<G4endl;
  return field;
}

G4MagneticField * DetectorConstruction::MakeUniformMagField2(  const std::string & filename,
							       double NormFac1, double NormFac2,
							       G4ThreeVector pos1, G4RotationMatrix rot1, 
							       G4ThreeVector pos2, G4RotationMatrix rot2 
							       )
{

  UNIField *field = new UNIField(filename, NormFac1, NormFac2 );
  field->AddPosRot(1, pos1, rot1); //SC magnet
  field->AddPosRot(2, pos2, rot2); //KURAMA magnet
  G4FieldManager* fieldMan = G4TransportationManager::GetTransportationManager()->GetFieldManager();
  fieldMan->SetDetectorField( field );
  fieldMan->CreateChordFinder( field );
  G4cout<<"Field is uniform field2"<<G4endl;

  return field;
}

G4MagneticField * DetectorConstruction::MakeMagFieldFromMap( const std::string & filename,
							     double NormFac, 
							     G4ThreeVector pos, G4RotationMatrix rot
							     )
{
  
  //G4MagneticField *field = new SCField( filename1, NormFac);
  SCField *field = new SCField( filename, NormFac);
  field->AddPosRot( pos, rot);
  G4FieldManager *fieldManager = G4TransportationManager::GetTransportationManager()->GetFieldManager();
  fieldManager->SetDetectorField( field );
  fieldManager->CreateChordFinder( field );
  fieldManager->GetChordFinder()->SetDeltaChord( 1.0E-3*mm );
  //fieldManager->GetChordFinder()->SetStepMinimum( 1.0*mm );

  return field;
}

G4MagneticField * DetectorConstruction::MakeMagFieldFromMap( const std::string & filename1,
							     const std::string & filename2,
							     double NormFac1, double NormFac2,
							     G4ThreeVector pos1, G4RotationMatrix rot1, 
							     G4ThreeVector pos2, G4RotationMatrix rot2 
							     )
{
  
  //G4MagneticField *field = new SCField( filename1, NormFac);
  SCField *field = new SCField( filename1, NormFac1);
  field->AddPosRot( pos1, rot1);
  field->AddFieldMap( filename2, NormFac2 );
  field->AddPosRot( pos2, rot2);
  G4FieldManager *fieldManager = G4TransportationManager::GetTransportationManager()->GetFieldManager();
  fieldManager->SetDetectorField( field );
  fieldManager->CreateChordFinder( field );
  fieldManager->GetChordFinder()->SetDeltaChord( 1.0E-3*mm );
  //fieldManager->GetChordFinder()->SetStepMinimum( 1.0*mm );

  return field;
}


//Main

G4VPhysicalVolume *DetectorConstruction::Construct()
{
  // Materials list //
  mList_ = DefineMaterials();
  mat_air = mList_->Air;
  mat_p10 = mList_->P10Gas;
  mat_C = mList_->C;
  mat_LH2 = mList_->LiqH2;
  mat_Scin = mList_->Scin;
  mat_G10 = mList_->G10;
  mat_Mylar = mList_->Mylar;
  mat_vacuum = mList_ -> Vacuum;
  mat_Aerogel = mList_-> Aerogel;
  mat_ArGas = mList_->ArGas;

  //field_=MakeUniformMagField( 1.0*tesla );
  
  ConfMan *confMan = ConfMan::GetConfManager();
  //int GeomFlag = confMan->GeomFlag();

  //world volume//                                                                                 
  G4LogicalVolume* logicWorld = new G4LogicalVolume( new G4Box("World", 10./2. *m, 10./2. *m, 10./2. *m), NistMan->FindOrBuildMaterial("G4_Galactic"), "World");
  logicWorld->SetVisAttributes(G4VisAttributes::Invisible);
  G4VPhysicalVolume* physiWorld = new G4PVPlacement(0,G4ThreeVector(),logicWorld,"World",0,false,0);

  //Rotation Angle
  //G4double rotation_angle = 6.5*deg;
  G4double rotation_angle = 0.*deg;

  //detector geometry//
  G4ThreeVector SCMagnetPos ( 0., 0., 0.0*mm);
  G4RotationMatrix SCMagnetRot;
  SCMagnetRot.rotateY( rotation_angle );

  G4double kurama_center_pos = 1620.0*mm;
  G4double kurama_move_x = 150.0*mm; 
  G4ThreeVector KuramaMagnetPos ( kurama_move_x, 0., kurama_center_pos);
  G4RotationMatrix KuramaMagnetRot;
  KuramaMagnetRot.rotateY( rotation_angle );

  //Make magnetic field//
  if(confMan->ExistField1() && !confMan->ExistField2())
    {
      G4cout<<"### SCMagnet Field On ###"<<G4endl;
      double fieldscale1 = confMan->GetFieldScale1();
      G4cout<<"### SCMagnet Field Scale: "<<fieldscale1 <<" ###"<<G4endl;
      field_=MakeMagFieldFromMap( confMan->FieldMapName1(), fieldscale1, SCMagnetPos, SCMagnetRot );
    }
  
  if(!confMan->ExistField1() && confMan->ExistField2())
    {
      G4cout<<"### KURAMA Field On ###"<<G4endl;
      double fieldscale2 = confMan->GetFieldScale2();
      G4cout<<"### KURAMA Field Scale: "<<fieldscale2 <<" ###"<<G4endl;
      field_=MakeMagFieldFromMap( confMan->FieldMapName2(), fieldscale2, KuramaMagnetPos+G4ThreeVector(-150.0, 0.0, 0.0), KuramaMagnetRot );
    }

  if(confMan->ExistField1() && confMan->ExistField2())
    {
      G4cout<<"### SCMagnet & KURAMA Field On ###"<<G4endl;
      double fieldscale1 = confMan->GetFieldScale1();
      G4cout<<"### SCMagnet Field Scale: "<<fieldscale1 <<" ###"<<G4endl;
      double fieldscale2 = confMan->GetFieldScale2();
      G4cout<<"### KURAMA Field Scale: "<<fieldscale2 <<" ###"<<G4endl;  
      //std::string fieldmapname2 = "/home/sbyang/g4work/lambda_search/Geant4/param/field/kurama_20160608.root";
      field_=MakeMagFieldFromMap( confMan->FieldMapName1(), confMan->FieldMapName2(), fieldscale1, fieldscale2, SCMagnetPos, SCMagnetRot, KuramaMagnetPos+G4ThreeVector(-150.0, 0.0, 0.0), KuramaMagnetRot );
    }

  if(!confMan->ExistField1() && !confMan->ExistField2())
    {
      G4cout<<"### Uniform Field On ###"<<G4endl;
      G4cout<<"### SCMagnet & KURAMA Field On ###"<<G4endl;
      double fieldscale1 = confMan->GetFieldScale1();
      G4cout<<"### SCMagnet Field Scale: "<<fieldscale1 <<" ###"<<G4endl;
      double fieldscale2 = confMan->GetFieldScale2();
      G4cout<<"### KURAMA Field Scale: "<<fieldscale2 <<" ###"<<G4endl;  
      //field_=MakeUniformMagField( 1.0*tesla );
      field_=MakeUniformMagField2( confMan->FieldMapName1(), fieldscale1, fieldscale2, SCMagnetPos, SCMagnetRot, KuramaMagnetPos, KuramaMagnetRot );
    }

  // Dummy Detectors //
  G4ThreeVector DummyDetectorPos1 ( 0., 0., 600.0+0.01/2.0);// End HC 
  G4RotationMatrix DummyDetectorRot1;
  DummyDetectorRot1.rotateY(rotation_angle);
  MakeDummyDetector1( physiWorld, DummyDetectorPos1, DummyDetectorRot1 );

  //SC Spectrometer
  MakeSCMagnet( physiWorld, SCMagnetPos, SCMagnetRot );

  return physiWorld;
}

G4bool DetectorConstruction::IsVolumeStopper( G4VPhysicalVolume *physVol ) const
{
  // Check SksMagnet                                                                                                                     
  G4String name = physVol->GetName();
  if( name=="SCMagnet_CoreFrame" || name=="SCCoilSupUp" ||
      name=="SCCoilSupDown" || name=="SCMagnet_CoilUp" || name=="SCMagnet_CoilDown")
    {
      //G4cout<<"particle hit: "<<name<<" Stepping action stop!"<<G4endl;                                                                
      return true;
    }
  else
    return false;
}

G4bool DetectorConstruction::IsVolumeStopper_KURAMA( G4VPhysicalVolume *physVol ) const
{
  // Check SksMagnet                                                                                                                     
  G4String name = physVol->GetName();
  if( name=="Coil1U_phys" || name=="Coil1D_phys" || name=="Coil4UR_phys" || name== "Coil4UL_phys" ||
      name=="Coil4DR_phys" || name=="Coil4DL_phys" || name=="Coil5U_phys" ||name=="Coil5D_phys" ||
      name=="Coil6LU_phys" || name=="Coil6RU_phys" || name=="Coil6RD_phys" || name=="Coil6LD_phys" ||
      name=="Coil8LU_phys" || name=="Coil8RU_phys" || name=="Coil8RD_phys" || name=="Coil8LD_phys" ||
      name=="Coil7ULU_phys" || name=="Coil7URU_phys" || name=="Coil7URD_phys" || name=="Coil7ULD_phys"||
      name=="Coil7DLU_phys" || name=="Coil7DRU_phys" || name=="Coil7DRD_phys" || name=="Coil7DLD_phys"||
      name=="Coil2UL_phys" || name=="Coil2UR_phys" || name=="Coil2DL_phys" || name=="Coil2DR_phys"||
      name=="Coil3UL_phys" || name=="Coil3UR_phys" || name=="Coil3DL_phys" || name=="Coil3DR_phys"||
      name=="upGuard_U_phys" || name=="upGuard_D_phys" || name=="upGuard_L_phys" || name=="upGuard_R_phys"||
      name=="Yoke_U_phys" || name=="Yoke_D_phys" || name=="Yoke_L_phys" || name=="Yoke_R_phys"||
      name=="downGuard_U_phys" || name=="downGuard_D_phys" || name=="downGuard_L_phys" || name=="downGuard_R_phys"||
      name=="upGuard_R1_phys" || name=="upGuard_R2D_phys" || name=="upGuard_R2U_phys" ||
      name=="Yoke_L_GapSpace_phys" || name=="Yoke_R_GapSpace_phys"
      )
    {
      //G4cout<<"particle hit: "<<name<<" Stepping action stop!"<<G4endl;                                                                
      return true;
    }
  else
    return false;
}
