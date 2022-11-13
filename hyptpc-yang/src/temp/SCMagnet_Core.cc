/**
 *  @file
 *  @brief  GsimSCMagnet_Core
 *  $Id: GsimSCMagnet_Core.cc,v 1.3 2006/12/06 19:14:30 nanjo Exp $ 
 *  $Log: GsimSCMagnet_Core.cc,v $
 *  Revision 1.3  2006/12/06 19:14:30  nanjo
 *  CVS variables.
 * 
 */
#include "GsimDetector/GsimSCMagnet_Core.h"
#include "GsimKernel/GsimUtil.h"
#include "GsimPersistency/GsimMessage.h"
#include "TROOT.h"

#include "G4Box.hh"
#include "G4Tubs.hh"
#include "G4IntersectionSolid.hh"
#include "G4SubtractionSolid.hh"
#include "G4UnionSolid.hh"
#include "G4LogicalVolume.hh"
#include "G4VisAttributes.hh"
#include "G4PVPlacement.hh"

GsimSCMagnet_Core::GsimSCMagnet_Core(std::string name,
		       GsimDetector* motherDetector,
		       G4ThreeVector transV,G4ThreeVector rotV,
		       int userFlag)
  : GsimDetector(name,motherDetector,transV,rotV,userFlag)
{
#ifdef GSIMDEBUG
  GsimMessage::getInstance()->debugEnter(__PRETTY_FUNCTION__);
#endif
  
  m_className = "GsimSCMagnet_Core";

  
  
  m_numberOfParameters=6;
  m_parameterArray= new double[m_numberOfParameters];
  m_parameterArray[0] =10.*mm; // x of box
  m_parameterArray[1] =10.*mm; // y of box
  m_parameterArray[2] =10.*mm; // z of box
  m_parameterArray[3] =0.*mm;  // x pos of tube
  m_parameterArray[4] =0.*mm;  // y pos of tube
  m_parameterArray[5] =3.*mm;  // r of tube
  
  m_outerMaterialName = "G4_Al";
  m_outerColorName    ="yellow";
  
  
  
  update();
  

#ifdef GSIMDEBUG
  GsimMessage::getInstance()->debugExit(__PRETTY_FUNCTION__);
#endif
}

GsimSCMagnet_Core::~GsimSCMagnet_Core()
{
#ifdef GSIMDEBUG
  GsimMessage::getInstance()->debugEnter(__PRETTY_FUNCTION__);
#endif
  ;
#ifdef GSIMDEBUG
  GsimMessage::getInstance()->debugExit(__PRETTY_FUNCTION__);
#endif
}

void GsimSCMagnet_Core::constructPhysiDetector()
{
#ifdef GSIMDEBUG
  GsimMessage::getInstance()->debugEnter(__PRETTY_FUNCTION__);
#endif
  G4VPhysicalVolume* physiMother=getMotherPhysiDetector();
  std::string fullName = getFullName();

  //W     : 1550, 1530, 1480, 1470,
  //D     : 1200, 1180, 1140, 1130,
  //H     : 950, 890, 732, 712,
  //Corner: 1449.569, 1429.569, 1357.645, 1347.645
  //R     : 800, 820, 850, 860
  //GAP   : 300, 320, 348, 358
  
  G4double size_Width[4];
  size_Width[0] = 1550*mm;
  size_Width[1] = 1530*mm;
  size_Width[2] = 1480*mm;
  size_Width[3] = 1470*mm;

  G4double size_Depth[4];
  size_Depth[0] = 1200*mm;
  size_Depth[1] = 1180*mm;
  size_Depth[2] = 1140*mm;
  size_Depth[3] = 1470*mm;

  G4double size_Height[4];
  size_Height[0] = 950*mm;
  size_Height[1] = 890*mm;
  size_Height[2] = 732*mm;
  size_Height[3] = 712*mm;

  G4double size_Corner[4];
  size_Corner[0] = 1449.569*mm;
  size_Corner[1] = 1429.569*mm;
  size_Corner[2] = 1357.645*mm;
  size_Corner[3] = 1347.645*mm;

  G4double size_Rad[4];
  size_Rad[0] = 800./2.*mm;
  size_Rad[1] = 820./2.*mm;
  size_Rad[2] = 850./2.*mm;
  size_Rad[3] = 860./2.*mm;

  G4double size_Gap[4];
  size_Gap[0] = 300*mm;
  size_Gap[1] = 320*mm;
  size_Gap[2] = 348*mm;
  size_Gap[3] = 358*mm;

  G4double size_SideGap[4];
  size_SideGap[0] = 725*mm;
  size_SideGap[1] = 745*mm;
  size_SideGap[2] = 775*mm;
  size_SideGap[3] = 785*mm;

  
  G4double size_dummyBox[3];
  size_dummyBox[0] = 2000*mm;
  size_dummyBox[1] = 2000*mm;
  size_dummyBox[2] = 2000*mm;

  G4double size_center_sub_position= 512.5*mm;
  

  if(physiMother) {

    G4Box* solidBox[4];
    G4Tubs* solidTube[4];
    G4Box* solidBox_subOuter[4];
    G4Box* solidBox_subInner[4];
    G4Box* solidBox_subSide[4];

    G4SubtractionSolid* solidBox_Octa[4];
    G4SubtractionSolid* solidBox_subCorner[4];
    G4SubtractionSolid* solidBox_subSide_1[4];
    G4SubtractionSolid* solidBox_subSide_2[4];
    G4SubtractionSolid* solidBox_Chamber[4];
    G4RotationMatrix* rot  = new G4RotationMatrix();
    rot->rotateZ( 45*deg);
    for( int i = 0; i < 4; i++){
      solidBox[i]           = new G4Box(G4String((fullName+Form("_Main_%d",i)).c_str())     ,size_Width[i]/2. ,size_Depth[i]/2. ,size_Height[i]/2.);
      solidBox_subOuter[i]  = new G4Box(G4String((fullName+Form("_SubOuter_%d",i)).c_str()) ,size_dummyBox[0] ,size_dummyBox[1] ,size_dummyBox[2]);
      solidBox_subInner[i]  = new G4Box(G4String((fullName+Form("_SubInner_%d",i)).c_str()) ,size_Corner[i]   ,size_Corner[i]   ,size_dummyBox[2]);
      solidBox_subSide[i]   = new G4Box(G4String((fullName+Form("_SubSide_%d",i)).c_str())  ,size_SideGap[i]/2.,size_SideGap[i]/2.,size_Gap[i]/2.);
      solidTube[i]          = new G4Tubs(G4String((fullName+Form("_Hole_%d",i)).c_str()),0, size_Rad[i], 2000*mm, 0, 2.*M_PI*rad);

      solidBox_subCorner[i] = new G4SubtractionSolid(G4String((fullName+Form("_SubCorner_%d",i)).c_str()),solidBox_subOuter[i]  ,solidBox_subInner[i] ,0   ,G4ThreeVector(0,0,0)); 
      solidBox_Octa[i]      = new G4SubtractionSolid(G4String((fullName+Form("_Octa_%d", i)).c_str())    ,solidBox[i]           ,solidBox_subCorner[i],rot,G4ThreeVector(0,0,0));
      solidBox_subSide_1[i] = new G4SubtractionSolid(G4String((fullName+Form("_SubSide_1_%d",i)).c_str()),solidBox_Octa[i]      ,solidBox_subSide[i]  ,rot,G4ThreeVector(0,size_center_sub_position,0));
      solidBox_subSide_2[i] = new G4SubtractionSolid(G4String((fullName+Form("_SubSide_2_%d",i)).c_str()),solidBox_subSide_1[i] ,solidBox_subSide[i]  ,rot,G4ThreeVector(0,-size_center_sub_position,0));
      solidBox_Chamber[i]   = new G4SubtractionSolid(G4String((fullName+Form("_Chamber_%d",i)).c_str())  ,solidBox_subSide_2[i] ,solidTube[i]         ,0, G4ThreeVector(0,0,0));
    }
    
    G4SubtractionSolid* Solid1 = new G4SubtractionSolid(G4String((fullName+"_Solid_1").c_str()), solidBox_Chamber[1], solidBox_Chamber[2], 0, G4ThreeVector(0,0,0));  
    G4UnionSolid* solidDetector = new G4UnionSolid( G4String((fullName).c_str()), solidBox_Chamber[1], Solid1, 0, G4ThreeVector(0,0,0));
    
    
    G4LogicalVolume* logicDetector
      = new G4LogicalVolume(solidDetector,
			    GsimUtil::getUtil()->findMaterial(m_outerMaterialName),
			    G4String(fullName.c_str()));
    G4VisAttributes* visAtt
      = new G4VisAttributes(GsimUtil::getUtil()->findColor(m_outerColorName));
    logicDetector->SetVisAttributes(visAtt);
    
    m_physiDetector = 
      new G4PVPlacement(m_rotationMatrix,
			G4ThreeVector(),
			G4String(fullName.c_str()),
			logicDetector,
			physiMother,
			false,0
			);
  }
#ifdef GSIMDEBUG
  GsimMessage::getInstance()->debugExit(__PRETTY_FUNCTION__);
#endif
}

void GsimSCMagnet_Core::updateParametersEtc()
{
#ifdef GSIMDEBUG
  GsimMessage::getInstance()->debugEnter(__PRETTY_FUNCTION__);
#endif
  if(!m_physiDetector) return;
  
  G4LogicalVolume* logicDetector = m_physiDetector->GetLogicalVolume();
  if(!logicDetector) {
    GsimMessage::getInstance()
      ->report("warning","No logicDetector is constructed.");
    return; 
  }
  G4VSolid* solidDetector = logicDetector->GetSolid();
  if(!solidDetector) {
    GsimMessage::getInstance()
      ->report("warning","No solidDetector is constructed.");
    return;
  }
  
  
  delete solidDetector;
  solidDetector=0;
  
  std::string fullName = getFullName();


  //W     : 1550, 1530, 1480, 1470,
  //D     : 1200, 1180, 1140, 1130,
  //H     : 950, 890, 732, 712,
  //Corner: 1449.569, 1429.569, 1357.645, 1347.645
  //R     : 800, 820, 850, 860
  //GAP   : 300, 320, 348, 358
  
  G4double size_Width[4];
  size_Width[0] = 1550*mm;
  size_Width[1] = 1530*mm;
  size_Width[2] = 1480*mm;
  size_Width[3] = 1470*mm;

  G4double size_Depth[4];
  size_Depth[0] = 1200*mm;
  size_Depth[1] = 1180*mm;
  size_Depth[2] = 1140*mm;
  size_Depth[3] = 1470*mm;

  G4double size_Height[4];
  size_Height[0] = 950*mm;
  size_Height[1] = 890*mm;
  size_Height[2] = 732*mm;
  size_Height[3] = 712*mm;

  G4double size_Corner[4];
  size_Corner[0] = 1449.569*mm;
  size_Corner[1] = 1429.569*mm;
  size_Corner[2] = 1357.645*mm;
  size_Corner[3] = 1347.645*mm;

  G4double size_Rad[4];
  size_Rad[0] = 800./2.*mm;
  size_Rad[1] = 820./2.*mm;
  size_Rad[2] = 850./2.*mm;
  size_Rad[3] = 860./2.*mm;

  G4double size_Gap[4];
  size_Gap[0] = 300*mm;
  size_Gap[1] = 320*mm;
  size_Gap[2] = 348*mm;
  size_Gap[3] = 358*mm;

  G4double size_SideGap[4];
  size_SideGap[0] = 725*mm;
  size_SideGap[1] = 745*mm;
  size_SideGap[2] = 775*mm;
  size_SideGap[3] = 785*mm;

  
  G4double size_dummyBox[3];
  size_dummyBox[0] = 2000*mm;
  size_dummyBox[1] = 2000*mm;
  size_dummyBox[2] = 2000*mm;

  G4double size_center_sub_position= 512.5*mm;
  
  G4Box* solidBox[4];
  G4Tubs* solidTube[4];
  G4Box* solidBox_subOuter[4];
  G4Box* solidBox_subInner[4];
  G4Box* solidBox_subSide[4];
  
  G4SubtractionSolid* solidBox_Octa[4];
  G4SubtractionSolid* solidBox_subCorner[4];
  G4SubtractionSolid* solidBox_subSide_1[4];
  G4SubtractionSolid* solidBox_subSide_2[4];
  G4SubtractionSolid* solidBox_Chamber[4];
  G4RotationMatrix* rot  = new G4RotationMatrix();
  rot->rotateZ( 45*deg);
  for( int i = 0; i < 4; i++){
    solidBox[i]           = new G4Box(G4String((fullName+Form("_Main_%d",i)).c_str())     ,size_Width[i]/2. ,size_Depth[i]/2. ,size_Height[i]/2.);
    solidBox_subOuter[i]  = new G4Box(G4String((fullName+Form("_SubOuter_%d",i)).c_str()) ,size_dummyBox[0] ,size_dummyBox[1] ,size_dummyBox[2]);
    solidBox_subInner[i]  = new G4Box(G4String((fullName+Form("_SubInner_%d",i)).c_str()) ,size_Corner[i]   ,size_Corner[i]   ,size_dummyBox[2]);
    solidBox_subSide[i]   = new G4Box(G4String((fullName+Form("_SubSide_%d",i)).c_str())  ,size_SideGap[i]/2.,size_SideGap[i]/2.,size_Gap[i]/2.);
    solidTube[i]          = new G4Tubs(G4String((fullName+Form("_Hole_%d",i)).c_str()),0, size_Rad[i], 2000*mm, 0, 2.*M_PI*rad);
    
    solidBox_subCorner[i] = new G4SubtractionSolid(G4String((fullName+Form("_SubCorner_%d",i)).c_str()),solidBox_subOuter[i]  ,solidBox_subInner[i] ,0   ,G4ThreeVector(0,0,0)); 
    solidBox_Octa[i]      = new G4SubtractionSolid(G4String((fullName+Form("_Octa_%d", i)).c_str())    ,solidBox[i]           ,solidBox_subCorner[i],rot,G4ThreeVector(0,0,0));
    solidBox_subSide_1[i] = new G4SubtractionSolid(G4String((fullName+Form("_SubSide_1_%d",i)).c_str()),solidBox_Octa[i]      ,solidBox_subSide[i]  ,rot,G4ThreeVector(0,size_center_sub_position,0));
    solidBox_subSide_2[i] = new G4SubtractionSolid(G4String((fullName+Form("_SubSide_2_%d",i)).c_str()),solidBox_subSide_1[i] ,solidBox_subSide[i]  ,rot,G4ThreeVector(0,-size_center_sub_position,0));
    solidBox_Chamber[i]   = new G4SubtractionSolid(G4String((fullName+Form("_Chamber_%d",i)).c_str())  ,solidBox_subSide_2[i] ,solidTube[i]         ,0, G4ThreeVector(0,0,0));
  }
  
  G4SubtractionSolid* Solid1 = new G4SubtractionSolid(G4String((fullName+"_Solid_1").c_str()), solidBox_Chamber[1], solidBox_Chamber[2], 0, G4ThreeVector(0,0,0));
  
  solidDetector = new G4UnionSolid( G4String((fullName).c_str()), solidBox_Chamber[1], Solid1);
    
  logicDetector->SetSolid(solidDetector);
#ifdef GSIMDEBUG
  GsimMessage::getInstance()->debugExit(__PRETTY_FUNCTION__);
#endif
}


void GsimSCMagnet_Core::extractTwoPlanesDXF(std::list<TVector3>& v0List,
				     std::list<TVector3>& v1List,
				     std::list<TVector3>& v2List,
				     std::list<TVector3>& v3List)
{
#ifdef GSIMDEBUG
  GsimMessage::getInstance()->debugEnter(__PRETTY_FUNCTION__);
#endif

  double dx=m_parameterArray[0]/2.;
  double dy=m_parameterArray[1]/2.;
  double dz=m_parameterArray[2]/2.;

  TVector3 v[2][4];
  v[0][0]=TVector3(+dx,+dy,+dz);
  v[0][1]=TVector3(-dx,+dy,+dz);
  v[0][2]=TVector3(-dx,-dy,+dz);
  v[0][3]=TVector3(+dx,-dy,+dz);
  v[1][0]=TVector3(+dx,+dy,-dz);
  v[1][1]=TVector3(-dx,+dy,-dz);
  v[1][2]=TVector3(-dx,-dy,-dz);
  v[1][3]=TVector3(+dx,-dy,-dz);
  
  for(int i=0;i<4;i++) {
    v0List.push_back(v[0][i]);
    v1List.push_back(v[1][i]);
  }


  double r1=m_parameterArray[5];
  double phi0=0;
  double dphi=2*M_PI;
  double sx=m_parameterArray[3];
  double sy=m_parameterArray[3];

  for(int i=0;i<64;i++) {
    double phi=dphi/64.*i + phi0;

    double x=r1*cos(phi)+sx;
    double y=r1*sin(phi)+sy;
    double z=+dz;
    v2List.push_back( TVector3(x,y,z) );
    z=-dz;
    v3List.push_back( TVector3(x,y,z) );
    
  }
  
#ifdef GSIMDEBUG
  GsimMessage::getInstance()->debugExit(__PRETTY_FUNCTION__);
#endif
}
