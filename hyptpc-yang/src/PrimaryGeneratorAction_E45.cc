/*
  PrimaryGeneratorAction.cc

  2017/10  Yang
*/

#include "PrimaryGeneratorAction.hh"
#include "DetectorConstruction.hh"
#include "AnalysisManager.hh"
#include "ConfMan.hh"

#include "globals.hh"
#include "G4Event.hh"
#include "G4ParticleGun.hh"
#include "G4ParticleTypes.hh"
#include "G4ParticleDefinition.hh"
#include "G4ParticleTable.hh"
#include "G4LorentzVector.hh"
#include "G4RandomDirection.hh"
#include "G4SystemOfUnits.hh"

#include "TFile.h"
#include "TTree.h"

#include "Randomize.hh"

#include "EvtGen/EvtGen.hh"
#include "EvtGenBase/EvtParticle.hh"
#include "EvtGenBase/EvtParticleFactory.hh"
#include "EvtGenBase/EvtPatches.hh"
#include "EvtGenBase/EvtPDL.hh"
#include "EvtGenBase/EvtRandom.hh"
#include "EvtGenBase/EvtReport.hh"
#include "EvtGenBase/EvtHepMCEvent.hh"
//#include "EvtGenBase/EvtStdlibRandomEngine.hh"
#include "EvtGenBase/EvtSimpleRandomEngine.hh"
#include "EvtGenBase/EvtAbsRadCorr.hh"
#include "EvtGenBase/EvtDecayBase.hh"
#include "EvtGenExternal/EvtExternalGenList.hh"

#include "EvtGenBase/EvtVector4R.hh"
#include "EvtGenBase/EvtVectorParticle.hh"
#include "EvtGenBase/EvtParticle.hh"
#include "EvtGenBase/EvtScalarParticle.hh"
#include "EvtGenBase/EvtDecayTable.hh"
#include "EvtGenBase/EvtPDL.hh"
#include "EvtGenBase/EvtStdHep.hh"
#include "EvtGenBase/EvtSecondary.hh"
#include "EvtGenBase/EvtReport.hh"
#include "EvtGenBase/EvtId.hh"
#include "EvtGenBase/EvtRandom.hh"
#include "EvtGenBase/EvtRandomEngine.hh"
#include "EvtGenBase/EvtParticleFactory.hh"
#include "CLHEP/Vector/LorentzVector.h"
#include "EvtGenBase/EvtStatus.hh"
#include "EvtGenBase/EvtAbsRadCorr.hh"
#include "EvtGenBase/EvtRadCorr.hh"

#include <cmath>
#include <iostream>
#include <fstream>
#include <string>
#include <list>
#include <zlib.h>
//#include <boost/iostreams/filtering_streambuf.hpp>
//#include <boost/iostreams/copy.hpp>
//#include <boost/iostreams/filter/gzip.hpp>

//beam pi-
void PrimaryGeneratorAction::GenerateBeamPiM(G4Event* anEvent, G4ThreeVector D, G4ThreeVector P)
{
  anaMan_->SetBeam(1, D, P);
  double mass_pim = 0.13957061;
  //particleGun -> SetParticleDefinition (particleTable -> FindParticle("kaonT1"));
  //particleGun -> SetParticleDefinition (particleTable -> FindParticle("kaon-"));
  particleGun -> SetParticleDefinition (particleTable -> FindParticle("pi-"));
  //std::cout<<"Beam energy GenerateTest 2" << std::endl;
  //G4ThreeVector dir = G4RandomDirection();
  G4ThreeVector beamx = D + GaussPosition(7.6, 4.2, 0.0);
  G4ThreeVector beamp ( P.x(), P.y(), P.z());
  G4ThreeVector beampu =  beamp/beamp.mag();
  G4double energy = (sqrt(mass_pim*mass_pim+beamp.mag2()) - mass_pim )*GeV;
  //std::cout<<"Beam energy: "<< energy /GeV << std::endl;

  particleGun->SetParticleMomentumDirection ( beampu );
  //particleGun->SetParticleMomentum ( beamp );
  particleGun->SetParticleTime ( 0.0 );
  particleGun->SetParticlePosition( beamx );
  particleGun->SetParticleEnergy( energy );
  particleGun->GeneratePrimaryVertex( anEvent);

  anaMan_->SetBeam(1, beamx, beamp);
}

//beam pi+
void PrimaryGeneratorAction::GenerateBeamPiP(G4Event* anEvent, G4ThreeVector D, G4ThreeVector P)
{
  anaMan_->SetBeam(1, D, P);
  double mass_pip = 0.13957061;
  //particleGun -> SetParticleDefinition (particleTable -> FindParticle("kaonT1"));
  //particleGun -> SetParticleDefinition (particleTable -> FindParticle("kaon-"));
  particleGun -> SetParticleDefinition (particleTable -> FindParticle("pi+"));
  //std::cout<<"Beam energy GenerateTest 2" << std::endl;
  //G4ThreeVector dir = G4RandomDirection();
  G4ThreeVector beamx = D + GaussPosition(7.6, 4.2, 0.0);
  G4ThreeVector beamp ( P.x(), P.y(), P.z());
  G4ThreeVector beampu =  beamp/beamp.mag();
  G4double energy = (sqrt(mass_pip*mass_pip+beamp.mag2()) - mass_pip )*GeV;
  //std::cout<<"Beam energy: "<< energy /GeV << std::endl;

  particleGun->SetParticleMomentumDirection ( beampu );
  //particleGun->SetParticleMomentum ( beamp );
  particleGun->SetParticleTime ( 0.0 );
  particleGun->SetParticlePosition( beamx );
  particleGun->SetParticleEnergy( energy );
  particleGun->GeneratePrimaryVertex( anEvent);
  anaMan_->SetBeam(1, beamx, beamp);

}





//elastic scattering pi- p --> pi- p
void PrimaryGeneratorAction::GenerateElasticPiM(G4Event* anEvent, G4ThreeVector D, G4ThreeVector P)
{

  //G4cout<< "[PrimaryGeneratorAction]:Elastic Generator Start!!" << G4endl;
  G4double Mi1 = 0.13957061; //pion mass
  G4double Mi2 = 0.938272; //proton mass
  G4double Mf1 = 0.13957061; //pion mass
  G4double Mf2 = 0.938272; //proton mass
 

  //beam profile generation
  G4ThreeVector LPos = GaussPosition(7.6, 4.2, 54.0) + D; 
  bool gaussposition_status = false;
  while( gaussposition_status )
    {
      double hz = LPos.z() + 143.0;
      double r = 27.0;
      double theta = asin( fabs(hz)/r );
      double dx = r*cos(theta);
      if( LPos.y() < 50.0 && LPos.y() > -50.0 &&  hz > -1.*r && hz < r && LPos.x() < dx && LPos.x() > -1.*dx )
	{
	  gaussposition_status = true;
	}
      else
	{
	  LPos = GaussPosition(7.6, 4.2, 54.0) + D;
	}
    }
  
  //beam momentum
  G4ThreeVector Pb = P;

  double cosx = 0.0; //set differential cross-section
  G4ThreeVector LPf1, LPf2;
  bool status = Scattering2Body(Mi1, Mi2, Mf1, Mf2, Pb, LPf1, LPf2, cosx);
  if( status == false) 
    {
      G4cout<<"[PrimaryGeneratorAction]: Elastic scattering generator problem"<<G4endl;
    }

  G4ThreeVector gPos(LPos), gPf1(LPf1), gPf2(LPf2);

  //G4cout<<"Beam momentum: "<< Pb << " cosx: "<< cosx <<" gPos: " << gPos <<G4endl;


  EvtVector4R evx4_pion, evp4_pion;
  evx4_pion.set(1, gPos.x());
  evx4_pion.set(2, gPos.y());
  evx4_pion.set(3, gPos.z());
  evx4_pion.set(0, 0.0);
  
  evp4_pion.set(1, gPf1.x());
  evp4_pion.set(2, gPf1.y());
  evp4_pion.set(3, gPf1.z());
  evp4_pion.set(0, sqrt(Mf1*Mf1 + gPf1.mag2()));

  EvtVector4R evx4_proton, evp4_proton;
  evx4_proton.set(1, gPos.x());
  evx4_proton.set(2, gPos.y());
  evx4_proton.set(3, gPos.z());
  evx4_proton.set(0, 0.0);
  
  evp4_proton.set(1, gPf2.x());
  evp4_proton.set(2, gPf2.y());
  evp4_proton.set(3, gPf2.z());
  evp4_proton.set(0, sqrt(Mf2*Mf2 + gPf2.mag2()));

  makeGun(anEvent, -211, evx4_pion, evp4_pion);  //pi-
  makeGun(anEvent, 2212, evx4_proton, evp4_proton);  //proton
  anaMan_->SetEvtGen(1, -211, 0, 0, 0, 0, 0, evx4_pion, evp4_pion);
  anaMan_->SetEvtGen(2, 2212, 0, 0, 0, 0, 1, evx4_proton, evp4_proton);
}

//elastic scattering pi+ p --> pi+ p
void PrimaryGeneratorAction::GenerateElasticPiP(G4Event* anEvent, G4ThreeVector D, G4ThreeVector P)
{

  //G4cout<< "[PrimaryGeneratorAction]:Elastic Generator Start!!" << G4endl;
  G4double Mi1 = 0.13957061; //pion mass
  G4double Mi2 = 0.938272; //proton mass
  G4double Mf1 = 0.13957061; //pion mass
  G4double Mf2 = 0.938272; //proton mass
 

  //beam profile generation
  G4ThreeVector LPos = GaussPosition(7.6, 4.2, 54.0) + D; 
  bool gaussposition_status = true;
  while( gaussposition_status )
    {
      double hz = LPos.z() + 143.0;
      double r = 27.0;
      double theta = asin( fabs(hz)/r );
      double dx = r*cos(theta);
      if( LPos.y() < 50.0 && LPos.y() > -50.0 &&  hz > -1.*r && hz < r && LPos.x() < dx && LPos.x() > -1.*dx )
	{
	  gaussposition_status = false;
	}
      else
	{
	  LPos = GaussPosition(7.6, 4.2, 54.0) + D;
	}
    }
  
  //beam momentum
  G4ThreeVector Pb = GaussMomentum(0.0, 0.0, 0.013) + P;
  anaMan_->SetBeam(1, LPos, Pb);  


  double cosx = 0.0; //set differential cross-section
  double comp = 0.0;
  //////////////////////pi+pi differential cross section
  G4double coeff[10]={1., 1.518, 2.332, 2.800, 2.272, 1.317, 1.084, 0.267, 0.105, 0.075}; //p = 1.6 GeV/c

  while(1)
    {
      comp = -100.;
      cosx=-100.;
      cosx = CLHEP::RandFlat::shoot(-1.,1.);
      //cosx = -1.;
      comp = CLHEP::RandFlat::shoot(0.,1.);
      //G4cout<<"test of legendre function"<<G4endl;                                                   
      G4double cross=0.;
      for(G4int jj=0; jj<10; jj++)
	{
	  //G4cout<<Legendre(jj,-1.)<<G4endl;                                                          
	  cross += coeff[jj]*Legendre(jj,cosx);
	}
      //G4cout<<"cross:comp --> "<<cross/5.9 <<" : "<< comp <<G4endl;                               
      
      if(comp <= (cross/15.05) )
	{
	  break;
	}
      else if(cross/15.05 > 1.)
	{
	  G4cout<<"cross is larger than normalization"<<G4endl;
	  exit(-1);
	}
    }
  
  //cosx = 0.9; // test
  G4ThreeVector LPf1, LPf2;
  bool status = Scattering2Body(Mi1, Mi2, Mf1, Mf2, Pb, LPf1, LPf2, cosx);
  if( status == false) 
    {
      G4cout<<"[PrimaryGeneratorAction]: Elastic scattering generator problem"<<G4endl;
    }
  
  G4ThreeVector gPos(LPos), gPf1(LPf1), gPf2(LPf2);
  
  //G4cout<<"Beam momentum: "<< Pb << " cosx: "<< cosx <<" gPos: " << gPos <<G4endl;
  
  
  EvtVector4R evx4_pion, evp4_pion;
  evx4_pion.set(1, gPos.x());
  evx4_pion.set(2, gPos.y());
  evx4_pion.set(3, gPos.z());
  evx4_pion.set(0, 0.0);
  
  evp4_pion.set(1, gPf1.x());
  evp4_pion.set(2, gPf1.y());
  evp4_pion.set(3, gPf1.z());
  evp4_pion.set(0, sqrt(Mf1*Mf1 + gPf1.mag2()));

  EvtVector4R evx4_proton, evp4_proton;
  evx4_proton.set(1, gPos.x());
  evx4_proton.set(2, gPos.y());
  evx4_proton.set(3, gPos.z());
  evx4_proton.set(0, 0.0);
  
  evp4_proton.set(1, gPf2.x());
  evp4_proton.set(2, gPf2.y());
  evp4_proton.set(3, gPf2.z());
  evp4_proton.set(0, sqrt(Mf2*Mf2 + gPf2.mag2()));

  anaMan_->SetEvtGen(1, 211, 0, 0, 0, 0, 1, evx4_pion, evp4_pion);
  anaMan_->SetEvtGen(2, 2212, 0, 0, 0, 0, 2, evx4_proton, evp4_proton);
  anaMan_->SetCosTheta(cosx);

  makeGun(anEvent,  211, evx4_pion, evp4_pion);  //pi+
  makeGun(anEvent, 2212, evx4_proton, evp4_proton);  //proton



}


void PrimaryGeneratorAction::GeneratePiPSchannel(G4Event* anEvent, EvtGen *evtGenerator, G4ThreeVector D, G4ThreeVector P)
{
  //beam profile generator
  //beam profile generation
  G4ThreeVector LPos = GaussPosition(7.6, 4.2, 54.0) + D; 
  bool gaussposition_status = true;
  while( gaussposition_status )
    {
      double hz = LPos.z() + 143.0;
      double r = 27.0;
      double theta = asin( fabs(hz)/r );
      double dx = r*cos(theta);
      if( LPos.y() < 50.0 && LPos.y() > -50.0 &&  hz > -1.*r && hz < r && LPos.x() < dx && LPos.x() > -1.*dx )
	{
	  gaussposition_status = false;
	}
      else
	{
	  LPos = GaussPosition(7.6, 4.2, 54.0) + D;
	}
    }

  G4ThreeVector Pb = GaussMomentum(0.0, 0.0, 0.013) + P;
  anaMan_->SetBeam(1, LPos, Pb);  
 
  //G4cout<<"Start Generate Test2"<<G4endl;
  /// mother particle momentum //
  double mass_pi = 0.13957061;
  double mass_proton = 0.938272081;
  G4LorentzVector lv_beam;
  G4LorentzVector lv_target;
  G4LorentzVector lv_particle;
  double pbeam = P.mag();

  lv_beam.setX(Pb.x());
  lv_beam.setY(Pb.y());
  lv_beam.setZ(Pb.z());
  lv_beam.setE(sqrt(mass_pi*mass_pi + pbeam*pbeam)); 

  lv_target.setX(0.0);
  lv_target.setY(0.0);
  lv_target.setZ(0.0);
  lv_target.setE(sqrt(mass_proton*mass_proton)); 
  
  lv_particle = lv_beam + lv_target;

  // make mother particle //
  EvtParticle* res(0);
  static EvtId RES = EvtPDL::getId(std::string("ResPiBeam++"));
  G4LorentzVector LvRes;
  G4ThreeVector TVp (lv_particle.x(), lv_particle.y(), lv_particle.z());
  G4ThreeVector TVx (LPos.x(), LPos.y(), LPos.z());
  //double mass_res = EvtPDL::getMass(RES);
  double mass_res = sqrt((lv_beam.e()+lv_target.e())*(lv_beam.e()+lv_target.e()) - pbeam*pbeam);

  // check total energy //
  if(mass_res < 0.14*2 + 0.94 )
    {
      G4cout<<"### Beam momentum is not enough to generate Res ###"<<G4endl;
      return;
    }
  /* 
  G4cout<<"########################### Test  ##############################"<<G4endl;
  G4cout<<"Momentum of K-: "<<pbeam << " GeV/c" <<G4endl;
  G4cout<<"Invariant mass of K + p: "<<lv_particle.m() << " GeV/c2" <<G4endl;
  G4cout<<"Momentum of Res: "<< TVp.mag() << " GeV/c"<<G4endl;
  G4cout<<"Mass of Res: "<< mass_res << " GeV/c2" <<G4endl;
  G4cout<<"################################################################" << G4endl;
  */
  LvRes.setVect(TVp);
  LvRes.setE(sqrt(mass_res*mass_res+TVp.mag2()));
  //LvRes.setE(sqrt(1.664*1.664+TVp.mag2()));

  EvtVector4R pInit_res( LvRes.e(), LvRes.vect().x(), LvRes.vect().y(), LvRes.vect().z() );
  res = EvtParticleFactory::particleFactory(RES, pInit_res);



  GenerateDecay(anEvent, evtGenerator, res, LPos);
}
