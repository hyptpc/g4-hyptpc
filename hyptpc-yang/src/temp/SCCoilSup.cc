/**
 *  @file
 *  @brief  GsimSCCoilSup
 *  $Id: GsimSCCoilSup.cc,v 1.3 2006/12/06 19:14:30 nanjo Exp $ 
 *  $Log: GsimSCCoilSup.cc,v $
 *  Revision 1.3  2006/12/06 19:14:30  nanjo
 *  CVS variables.
 * 
 */
#include "GsimDetector/GsimSCCoilSup.h"
#include "GsimKernel/GsimUtil.h"
#include "GsimPersistency/GsimMessage.h"

#include "G4Box.hh"
#include "G4Tubs.hh"
#include "G4IntersectionSolid.hh"
#include "G4SubtractionSolid.hh"
#include "G4LogicalVolume.hh"
#include "G4VisAttributes.hh"
#include "G4PVPlacement.hh"

GsimSCCoilSup::GsimSCCoilSup(std::string name,
		       GsimDetector* motherDetector,
		       G4ThreeVector transV,G4ThreeVector rotV,
		       int userFlag)
  : GsimDetector(name,motherDetector,transV,rotV,userFlag)
{
#ifdef GSIMDEBUG
  GsimMessage::getInstance()->debugEnter(__PRETTY_FUNCTION__);
#endif
  
  m_className = "GsimSCCoilSup";
  
  m_numberOfParameters=6;
  m_parameterArray= new double[m_numberOfParameters];
  m_parameterArray[0] =10.*mm; // x of box
  m_parameterArray[1] =10.*mm; // y of box
  m_parameterArray[2] =10.*mm; // z of box
  m_parameterArray[3] =0.*mm;  // x pos of tube
  m_parameterArray[4] =0.*mm;  // y pos of tube
  m_parameterArray[5] =3.*mm;  // r of tube
  
  m_outerMaterialName = "G4_Fe";
  m_outerColorName    ="yellow";
  
  
  
  update();
  

#ifdef GSIMDEBUG
  GsimMessage::getInstance()->debugExit(__PRETTY_FUNCTION__);
#endif
}

GsimSCCoilSup::~GsimSCCoilSup()
{
#ifdef GSIMDEBUG
  GsimMessage::getInstance()->debugEnter(__PRETTY_FUNCTION__);
#endif
  ;
#ifdef GSIMDEBUG
  GsimMessage::getInstance()->debugExit(__PRETTY_FUNCTION__);
#endif
}

void GsimSCCoilSup::constructPhysiDetector()
{
#ifdef GSIMDEBUG
  GsimMessage::getInstance()->debugEnter(__PRETTY_FUNCTION__);
#endif
  G4VPhysicalVolume* physiMother=getMotherPhysiDetector();
  std::string fullName = getFullName();
  if(physiMother) {

    G4double RadIn = 445*mm;
    G4double RadOut= 545*mm;
    G4double GapRadIn = 465*mm;
    G4double GapRadOut = 545*mm;

    G4double SupHeight = 112*mm;
    G4double GapHeight = 62*mm;


    G4Tubs* solidTube_Sup = new G4Tubs(G4String((fullName+"Main")).c_str(),
				       RadIn,RadOut,SupHeight/2.,0,360*deg);
    G4Tubs* solidTube_Sub = new G4Tubs(G4String((fullName+"Sub")).c_str(),
				       GapRadIn,GapRadOut+20*mm,GapHeight/2.,0,360*deg);
    G4SubtractionSolid* solidDetector = new G4SubtractionSolid( G4String((fullName).c_str()),
							  solidTube_Sup,
							  solidTube_Sub,
							  0, G4ThreeVector(0,0,0));    
    
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

void GsimSCCoilSup::updateParametersEtc()
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
  
  G4double RadIn = 445*mm;
  G4double RadOut= 545*mm;
  G4double GapRadIn = 465*mm;
  G4double GapRadOut = 545*mm;
  
  G4double SupHeight = 112*mm;
  G4double GapHeight = 62*mm;
  
  
  G4Tubs* solidTube_Sup = new G4Tubs(G4String((fullName+"Main")).c_str(),
				     RadIn,RadOut,SupHeight/2.,0,360*deg);
  G4Tubs* solidTube_Sub = new G4Tubs(G4String((fullName+"Sub")).c_str(),
				     GapRadIn,GapRadOut+20*mm,GapHeight/2.,0,360*deg);
  solidDetector = new G4SubtractionSolid( G4String((fullName).c_str()),
					  solidTube_Sup,
					  solidTube_Sub,
					  0, G4ThreeVector(0,0,0));    
    
  logicDetector->SetSolid(solidDetector);
#ifdef GSIMDEBUG
  GsimMessage::getInstance()->debugExit(__PRETTY_FUNCTION__);
#endif
}


void GsimSCCoilSup::extractTwoPlanesDXF(std::list<TVector3>& v0List,
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
