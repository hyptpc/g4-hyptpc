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


//beam K-
void PrimaryGeneratorAction::GenerateBeamKaonM(G4Event* anEvent, G4ThreeVector D, G4ThreeVector P)
{
  //anaMan_->SetBeam(1, D, P);
  double mass_kaonm = 0.493677;
  //particleGun -> SetParticleDefinition (particleTable -> FindParticle("kaonT1"));
  //particleGun -> SetParticleDefinition (particleTable -> FindParticle("kaon-"));
  particleGun -> SetParticleDefinition (particleTable -> FindParticle("kaon-"));
  //std::cout<<"Beam energy GenerateTest 2" << std::endl;

  //G4ThreeVector dir = G4RandomDirection();
  //G4ThreeVector beamx = D + GaussPosition(7.6, 4.2, 0.0);
  //G4ThreeVector beamp ( P.x(), P.y(), P.z());

  G4ThreeVector beamx = D;
  G4ThreeVector beamp = P;

  G4ThreeVector beampu =  beamp/beamp.mag();
  G4double energy = (sqrt(mass_kaonm*mass_kaonm+beamp.mag2()) - mass_kaonm )*GeV;
  //std::cout<<"Beam energy: "<< energy /GeV << std::endl;

  particleGun->SetParticleMomentumDirection ( beampu );
  //particleGun->SetParticleMomentum ( beamp );
  particleGun->SetParticleTime ( 0.0 );
  particleGun->SetParticlePosition( beamx );
  particleGun->SetParticleEnergy( energy );
  particleGun->GeneratePrimaryVertex( anEvent);
  anaMan_->SetBeam(1, beamx, beamp);

}

void PrimaryGeneratorAction::GenerateScatteredKPlus(G4Event* anEvent, G4ThreeVector D, G4ThreeVector P)
{
  //anaMan_->SetBeam(1, D, P);
  double mass_kaonm = 0.493677;
  particleGun -> SetParticleDefinition (particleTable -> FindParticle("kaon+"));

  //beam rotate angle
  //G4double rotate_angle = -3.13*degree; //1.8 GeV/c && field = 1.0 case
  //G4double rotate_angle = -2.664*degree; //1.8 GeV/c && field = 0.85 case
  //G4double rotate_angle = -2.507*degree; //1.8 GeV/c && field = 0.8 case
  //G4double rotate_angle = -2.35*degree; //1.8 GeV/c && field = 0.75 case
  //G4double rotate_angle = 0.0*degree;                                             
  G4double rotate_angle = 3.135*degree; //1.8 GeV/c && field = -1.0 case

  G4ThreeVector TVbeam_rotation = BeamMomRotate( P, rotate_angle);
  //G4ThreeVector TVbeam_rotation = P; //test
  G4ThreeVector pos_shift = GaussPosition(3.0, 3.0, 20.0); //C target

  bool flag = true;
  while ( flag )
    {
      if( fabs(pos_shift.x()) < 15.0 && fabs(pos_shift.y()) < 10.0 && fabs(pos_shift.z()) < 10.0 )
	{
	  flag = false;
	}
      else
	{
	  pos_shift = GaussPosition(3.0, 3.0, 20.0); //C target
	}
    }

  G4ThreeVector beam_box_x = pos_shift + D;
  //G4ThreeVector beam_box_x = D; //test

  anaMan_->SetBeam(1, beam_box_x, TVbeam_rotation);

  double phi = (G4UniformRand()) * 2.0 * acos(-1.) ;
  double theta = (G4UniformRand()) * 30.0 * 3.1415 / 180.0 ;
  double pmag = G4UniformRand() * 0.2 + 1.2;

  G4ThreeVector beam_rotation_unit = TVbeam_rotation/TVbeam_rotation.mag();
  G4ThreeVector beam_rotation_orthogonal_unit = beam_rotation_unit.orthogonal();

  G4ThreeVector p_unit = beam_rotation_unit;
  p_unit.rotate( theta, beam_rotation_orthogonal_unit); //theta rotation
  p_unit.rotate( phi, beam_rotation_unit ); //phi rotation


  G4ThreeVector p = p_unit*pmag;

  
  G4ThreeVector beamx = beam_box_x;
  G4ThreeVector beamp ( p.x(), p.y(), p.z() );
  G4ThreeVector beampu =  beamp/beamp.mag();

  /*
  G4double energy = (sqrt(mass_kaonm*mass_kaonm+beamp.mag2()) - mass_kaonm )*GeV;
  //std::cout<<"Beam energy: "<< energy /GeV << std::endl;
  particleGun->SetParticleMomentumDirection ( beampu );
  //particleGun->SetParticleMomentum ( beamp );
  particleGun->SetParticleTime ( 0.0 );
  particleGun->SetParticlePosition( beamx );
  particleGun->SetParticleEnergy( energy );
  particleGun->GeneratePrimaryVertex( anEvent);
  */


  G4ThreeVector gPos(beamx), gPf1(beamp);
  //G4ThreeVector gPos(beamx), gPf1(TVbeam_rotation); // test
  
  EvtVector4R evx4_kaon, evp4_kaon;
  evx4_kaon.set(1, gPos.x());
  evx4_kaon.set(2, gPos.y());
  evx4_kaon.set(3, gPos.z());
  evx4_kaon.set(0, 0.0);
  
  evp4_kaon.set(1, gPf1.x());
  evp4_kaon.set(2, gPf1.y());
  evp4_kaon.set(3, gPf1.z());
  evp4_kaon.set(0, sqrt(mass_kaonm*mass_kaonm + gPf1.mag2()));


  anaMan_->SetEvtGen(1, 321, 0, 0, 0, 0, 1, evx4_kaon, evp4_kaon);
  anaMan_->SetCosTheta(cos(theta));
  makeGun(anEvent, 321, evx4_kaon, evp4_kaon); //K+
  //makeGun(anEvent, -321, evx4_kaon, evp4_kaon); //K-

}



//beam K- 1.8 GeV/c Profile from Root
void PrimaryGeneratorAction::GenerateBeamKaonMBeamRoot(G4Event* anEvent)
{
  //target V0+1175 E13 target V0+1657
  //target size x: 30 z: 10 y: 15

  double mass_kaonm = 0.493677;
  particleGun -> SetParticleDefinition (particleTable -> FindParticle("kaon-"));
  G4ThreeVector X, P;
  double inu, inv;
  double z_target, z_diff, z_diff_w;
  double profile_x, profile_y, profile_z;
  bool flag = true;

					       
  while ( flag )
    {
      ReadBeamProfile(X, P);
      
      
      if(P.z() == 0.0)
	{
	  inu  = 0.0;
	  inv  = 0.0;
	}
      else
	{
	  inu = P.x()/P.z();
	  inv = P.y()/P.z();
	}

      z_target = (G4UniformRand() - 0.5) * 20.0 ; //target z size 20 mm
      z_diff = 1175.0 - 1657.0;
      z_diff_w = z_target + z_diff;
      
      profile_x = X.x() + z_diff_w*inu;
      profile_y = X.y() + z_diff_w*inv;
      profile_z = z_target -143.0;

      if(profile_x < 15.0 && profile_x > -15.0 && profile_y > -10.0 && profile_y < 10.0)
	{
	  flag = false;
	}

    }


  //G4cout<<"[PrimaryGeneratorAction]beam u and v: "<<inu<<" "<<inv<<G4endl;


  G4ThreeVector beamx (profile_x, profile_y, profile_z);
  G4ThreeVector beamp ( P.x(), P.y(), P.z());

  G4ThreeVector beampu =  beamp/beamp.mag();
  G4double energy = (sqrt(mass_kaonm*mass_kaonm+beamp.mag2()) - mass_kaonm )*GeV;
  //std::cout<<"Beam energy: "<< energy /GeV << std::endl;

  particleGun->SetParticleMomentumDirection ( beampu );
  //particleGun->SetParticleMomentum ( beamp );
  particleGun->SetParticleTime ( 0.0 );
  particleGun->SetParticlePosition( beamx );
  particleGun->SetParticleEnergy( energy );
  particleGun->GeneratePrimaryVertex( anEvent);
  anaMan_->SetBeam(1, beamx, beamp);

}

//E42 Test beam profile from Root
void PrimaryGeneratorAction::GenerateHDiBaryonBeamRoot_test(G4Event* anEvent)
{
  //target V0+1175 E13 target V0+1657
  //target size x: 30 z: 10 y: 15

  double mass_kaonm = 0.493677;
  particleGun -> SetParticleDefinition (particleTable -> FindParticle("kaon-"));
  G4ThreeVector X, P;
  double inu, inv;
  double z_target, z_diff, z_diff_w;
  double profile_x, profile_y, profile_z;
  bool flag = true;

					       
  while ( flag )
    {
      ReadBeamProfile(X, P);
      
      
      if(P.z() == 0.0)
	{
	  inu  = 0.0;
	  inv  = 0.0;
	}
      else
	{
	  inu = P.x()/P.z();
	  inv = P.y()/P.z();
	}

      z_target = (G4UniformRand() - 0.5) * 20.0 ; //target z size 20 mm
      z_diff = 1175.0 - 1657.0;
      z_diff_w = z_target + z_diff;
      
      profile_x = X.x() + z_diff_w*inu;
      profile_y = X.y() + z_diff_w*inv;
      profile_z = z_target -143.0;

      if(profile_x < 15.0 && profile_x > -15.0 && profile_y > -10.0 && profile_y < 10.0)
	{
	  flag = false;
	}

    }


  //G4cout<<"[PrimaryGeneratorAction]beam u and v: "<<inu<<" "<<inv<<G4endl;


  G4ThreeVector beamx (profile_x, profile_y, profile_z);
  //G4ThreeVector beamp ( P.x(), P.y(), P.z());
  G4ThreeVector beamp ( 0.0, 0.0, 1.8);
  anaMan_->SetBeam(1, beamx, beamp);
  
  //
  double cosx = 0.0; //set differential cross-section
  //cosx = (G4UniformRand()-0.5)*2.0; //phase space
  cosx = 1.0-(G4UniformRand())*0.5; //test
  G4double Mi1 = mass_kaonm; //Kaon-
  G4double Mi2 = 0.938272; //proton
  G4double Mf1 = mass_kaonm;  //Kaon+
  G4double Mf2 = 1.3217;   //Xi-
  G4ThreeVector Pb = beamp;

  //cosx = 0.9; // test
  G4ThreeVector LPf1, LPf2;
  bool status = Scattering2Body(Mi1, Mi2, Mf1, Mf2, Pb, LPf1, LPf2, cosx);
  if( status == false) 
    {
      G4cout<<"[PrimaryGeneratorAction]: Elastic scattering generator problem"<<G4endl;
    }
  
  G4ThreeVector gPos(beamx), gPf1(LPf1), gPf2(LPf2);
  
  //G4cout<<"Beam momentum: "<< Pb << " cosx: "<< cosx <<" gPos: " << gPos <<G4endl;
  
  
  EvtVector4R evx4_kaon, evp4_kaon;
  evx4_kaon.set(1, gPos.x());
  evx4_kaon.set(2, gPos.y());
  evx4_kaon.set(3, gPos.z());
  evx4_kaon.set(0, 0.0);
  
  evp4_kaon.set(1, gPf1.x());
  evp4_kaon.set(2, gPf1.y());
  evp4_kaon.set(3, gPf1.z());
  evp4_kaon.set(0, sqrt(Mf1*Mf1 + gPf1.mag2()));

  EvtVector4R evx4_xi, evp4_xi;
  evx4_xi.set(1, gPos.x());
  evx4_xi.set(2, gPos.y());
  evx4_xi.set(3, gPos.z());
  evx4_xi.set(0, 0.0);
  
  evp4_xi.set(1, gPf2.x());
  evp4_xi.set(2, gPf2.y());
  evp4_xi.set(3, gPf2.z());
  evp4_xi.set(0, sqrt(Mf2*Mf2 + gPf2.mag2()));

  anaMan_->SetEvtGen(1, 321, 0, 0, 0, 0, 1, evx4_kaon, evp4_kaon);
  anaMan_->SetEvtGen(2, 3312, 0, 0, 0, 0, 2, evx4_xi, evp4_xi);
  anaMan_->SetCosTheta(cosx);

  makeGun(anEvent,  321, evx4_kaon, evp4_kaon);  //K+
  makeGun(anEvent, 3312, evx4_xi, evp4_xi);  //Xi-

  /*
  G4ThreeVector beampu =  beamp/beamp.mag();
  G4double energy = (sqrt(mass_kaonm*mass_kaonm+beamp.mag2()) - mass_kaonm )*GeV;
  //std::cout<<"Beam energy: "<< energy /GeV << std::endl;

  particleGun->SetParticleMomentumDirection ( beampu );
  //particleGun->SetParticleMomentum ( beamp );
  particleGun->SetParticleTime ( 0.0 );
  particleGun->SetParticlePosition( beamx );
  particleGun->SetParticleEnergy( energy );
  particleGun->GeneratePrimaryVertex( anEvent);
  */
}


//Test H-dibaryon by EvtGen
void PrimaryGeneratorAction::GenerateHDiBaryon_test(G4Event* anEvent, EvtGen *evtGenerator, G4ThreeVector D, G4ThreeVector P)
{
  G4ThreeVector beamx = D;
  G4ThreeVector beamp = P;

  double sigmax_beam = 3.0;
  double sigmay_beam = 3.0;
  double zlength = 20.0;

  //beam rotate angle
  G4double rotate_angle = -3.13*degree; //1.8 GeV/c && field = 1.0 case
  //G4double rotate_angle = -2.664*degree; //1.8 GeV/c && field = 0.85 case
  //G4double rotate_angle = -2.507*degree; //1.8 GeV/c && field = 0.8 case
  //G4double rotate_angle = -2.35*degree; //1.8 GeV/c && field = 0.75 case
  //G4double rotate_angle = 0.0*degree;                                             
  //G4double rotate_angle = 3.135*degree; //1.8 GeV/c && field = -1.0 case

  G4ThreeVector beamp_rotate = BeamMomRotate( P, rotate_angle);
  G4ThreeVector pos_shift = GaussPosition(sigmax_beam, sigmay_beam, zlength); //C target

  bool flag = true;
  while ( flag )
    {
      if( fabs(pos_shift.x()) < 15.0 && fabs(pos_shift.y()) < 10.0 && fabs(pos_shift.z()) < 10.0 )
        {
          flag = false;
        }
      else
        {
          pos_shift = GaussPosition(sigmax_beam, sigmay_beam, zlength); //C target                                     
        }
    }

  G4ThreeVector beamx_shift = pos_shift + beamx;


  double mass_km = 0.493677;
  double mass_target = 0.938272081*2.0;
  G4LorentzVector lv_beam;
  G4LorentzVector lv_target;
  G4LorentzVector lv_particle;
  double pbeam = beamp_rotate.mag();

  lv_beam.setX(beamp_rotate.x());
  lv_beam.setY(beamp_rotate.y());
  lv_beam.setZ(beamp_rotate.z());
  lv_beam.setE(sqrt(mass_km*mass_km + pbeam*pbeam)); 

  lv_target.setX(0.0);
  lv_target.setY(0.0);
  lv_target.setZ(0.0);
  lv_target.setE(sqrt(mass_target*mass_target)); 
  
  lv_particle = lv_beam + lv_target;

  // make mother particle //
  EvtParticle* lam1663(0);
  //G4cout<<"Test Started: "<<G4endl;
  static EvtId LAM1663 = EvtPDL::getId(std::string("MP+"));
  G4LorentzVector LvLam1663;
  G4ThreeVector TVp (lv_particle.x(), lv_particle.y(), lv_particle.z());

  double mass_lam1663 = sqrt((lv_beam.e()+lv_target.e())*(lv_beam.e()+lv_target.e()) - pbeam*pbeam);

  // check total energy //
  if(mass_lam1663 < 2.0*1.115683 + 0.49 )
    {
      G4cout<<"### Beam momentum is not enough to generate Lam1663 ###"<<G4endl;
      return;
    }
  /*
  G4cout<<"########################### Test  ##############################"<<G4endl;
  G4cout<<"Momentum of K-: "<<pbeam << " GeV/c" <<G4endl;
  G4cout<<"Invariant mass of K + p: "<<lv_particle.m() << " GeV/c2" <<G4endl;
  G4cout<<"Momentum of Lam1663: "<< TVp.mag() << " GeV/c"<<G4endl;
  G4cout<<"Mass of Lam1663: "<< mass_lam1663 << " GeV/c2" <<G4endl;
  G4cout<<"################################################################" << G4endl;
  */
  LvLam1663.setVect(TVp);
  LvLam1663.setE(sqrt(mass_lam1663*mass_lam1663+TVp.mag2()));
  //LvLam1663.setE(sqrt(1.664*1.664+TVp.mag2()));

  EvtVector4R pInit_lam1663( LvLam1663.e(), LvLam1663.vect().x(), LvLam1663.vect().y(), LvLam1663.vect().z() );
  lam1663 = EvtParticleFactory::particleFactory(LAM1663, pInit_lam1663);
  GenerateDecay(anEvent, evtGenerator, lam1663, beamx_shift);
  anaMan_->SetBeam(1, beamx_shift, beamp_rotate);
}
