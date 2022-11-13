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



void PrimaryGeneratorAction::GenerateTest(G4Event* anEvent, G4ThreeVector D, G4ThreeVector P)
{
  double mass_km = 0.493677;
  //particleGun -> SetParticleDefinition (particleTable -> FindParticle("kaonT1"));
  //particleGun -> SetParticleDefinition (particleTable -> FindParticle("kaon-"));
  particleGun -> SetParticleDefinition (particleTable -> FindParticle("pi-"));
  //std::cout<<"Beam energy GenerateTest 2" << std::endl;
  //G4ThreeVector dir = G4RandomDirection();
  G4ThreeVector beamx ( D.x(), D.y(), D.z());
  G4ThreeVector beamp ( P.x(), P.y(), P.z());
  G4ThreeVector beampu =  beamp/beamp.mag();
  G4double energy = (sqrt(mass_km*mass_km+beamp.mag2()) - mass_km )*GeV;
  //std::cout<<"Beam energy: "<< energy /GeV << std::endl;

  particleGun->SetParticleMomentumDirection ( beampu );
  //particleGun->SetParticleMomentum ( beamp );
  particleGun->SetParticleTime ( 0.0 );
  particleGun->SetParticlePosition( beamx );
  particleGun->SetParticleEnergy( energy );
  particleGun->GeneratePrimaryVertex( anEvent);
}


void PrimaryGeneratorAction::GenerateTest2(G4Event* anEvent, EvtGen *evtGenerator, G4ThreeVector D, G4ThreeVector P)
{
  //G4cout<<"Start Generate Test2"<<G4endl;
  /// mother particle momentum //
  double mass_km = 0.493677;
  double mass_proton = 0.938272081;
  G4LorentzVector lv_beam;
  G4LorentzVector lv_target;
  G4LorentzVector lv_particle;
  double pbeam = P.mag();

  lv_beam.setX(P.x());
  lv_beam.setY(P.y());
  lv_beam.setZ(P.z());
  lv_beam.setE(sqrt(mass_km*mass_km + pbeam*pbeam)); 

  lv_target.setX(0.0);
  lv_target.setY(0.0);
  lv_target.setZ(0.0);
  lv_target.setE(sqrt(mass_proton*mass_proton)); 
  
  lv_particle = lv_beam + lv_target;

  // make mother particle //
  EvtParticle* lam1663(0);
  static EvtId LAM1663 = EvtPDL::getId(std::string("Lambda(1663)0"));
  G4LorentzVector LvLam1663;
  G4ThreeVector TVp (lv_particle.x(), lv_particle.y(), lv_particle.z());
  G4ThreeVector TVx (D.x(), D.y(), D.z());
  //double mass_lam1663 = EvtPDL::getMass(LAM1663);
  double mass_lam1663 = sqrt((lv_beam.e()+lv_target.e())*(lv_beam.e()+lv_target.e()) - pbeam*pbeam);

  // check total energy //
  if(mass_lam1663 < 1.115683 + 0.547862 )
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
  GenerateDecay(anEvent, evtGenerator, lam1663, D);
}

//beam K- 2018.10.17
void PrimaryGeneratorAction::GenerateBeamKaonMBr(G4Event* anEvent, G4ThreeVector D, G4ThreeVector P)
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

  //G4ThreeVector pos_shift = GaussPosition(7.6, 2.5, 0.0); //C target
  //G4ThreeVector beamx = D + pos_shift;
  //G4ThreeVector beamx = D;

  G4ThreeVector X_file, P_file;
  ReadBeamProfile(X_file, P_file);

  G4ThreeVector beamx(X_file.x() + D.x(), X_file.y() + D.y(), D.z());

  double beam_file_mag = P_file.mag();
  double beam_mag_set = P.mag()/0.9*beam_file_mag; //for 0.9 GeV/c beam profile file

  //G4ThreeVector p_dir = GaussMomentum(-0.009, 0.02, 0.0, 0.0005); //u and v, rad unit 
  //G4ThreeVector p_dir = GaussMomentum(0.00, 0.02, 0.0, 0.0005); //u and v, rad unit 
  G4ThreeVector p_dir(P_file.x()/beam_file_mag, P_file.y()/beam_file_mag, P_file.z()/beam_file_mag) ; //from beam profile file 

  //G4ThreeVector beamp (0., 0., beam_file_mag);
  G4ThreeVector beamp (beam_mag_set*p_dir.x(), beam_mag_set*p_dir.y(), beam_mag_set*p_dir.z());
  //G4ThreeVector beamp = P_file;

  //beam rotate angle
  //G4double rotate_angle = -3.13*degree; //1.8 GeV/c && field = 1.0 case
  G4double rotate_angle = 0.0*degree; //1.8 GeV/c && field = 0.0 case
  //G4double rotate_angle = -2.664*degree; //1.8 GeV/c && field = 0.85 case                   
  //G4double rotate_angle = 7.7*degree; //0.73 GeV/c && field = 0.8 case                     
 

  G4ThreeVector beamp_rotate = BeamMomRotate( beamp, rotate_angle);

  //G4ThreeVector beampu =  beamp/beamp.mag();
  G4ThreeVector beampu =  beamp_rotate/beamp_rotate.mag();
  //G4double energy = (sqrt(mass_kaonm*mass_kaonm+beamp.mag2()) - mass_kaonm )*GeV;
  G4double energy = (sqrt(mass_kaonm*mass_kaonm+beamp_rotate.mag2()) - mass_kaonm )*GeV;
  //std::cout<<"Beam energy: "<< energy /GeV << std::endl;

  particleGun->SetParticleMomentumDirection ( beampu );
  //particleGun->SetParticleMomentum ( beamp );
  particleGun->SetParticleTime ( 0.0 );
  particleGun->SetParticlePosition( beamx );
  particleGun->SetParticleEnergy( energy );
  particleGun->GeneratePrimaryVertex( anEvent);
  anaMan_->SetBeam(1, beamx, beamp_rotate);

}




void PrimaryGeneratorAction::GenerateLambdaStarBr(G4Event* anEvent, EvtGen *evtGenerator, G4ThreeVector D, G4ThreeVector P)
{
  G4ThreeVector beamx = D;

  G4ThreeVector X_file, P_file;
  ReadBeamProfile(X_file, P_file);

  double z = GaussZ(54.0); //at target center;
  //G4cout<<"read file1 x: "<<X_file.x()<<" y: "<<X_file.y()<<" z: " << z <<G4endl; 

  bool flag = false;
  //double target_posz = 1300.0 - 143.0 + z; //in experimental hall
  double target_posz = 1300.0 + z; //in experimental hall

  double beam_vx = 0.0;
  double beam_vy = 0.0;
  double beam_vz = 0.0;
  double beam_vz_target = 0.0;

  while ( !flag )
    {
      double u = P_file.x()/P_file.z();
      double v = P_file.y()/P_file.z();

      double beam_vx_temp = X_file.x() + target_posz*u + D.x();
      double beam_vy_temp = X_file.y() + target_posz*v + D.y();
      double beam_vz_target_temp = z; //center is SC magnet center
      //G4cout<<"read file2 x: "<<beam_vx_temp<<" y: "<<beam_vy_temp<<" z: " << beam_vz_target_temp<<G4endl;
      
      if( fabs(beam_vy_temp) < 50.0 && sqrt(beam_vx_temp*beam_vx_temp + beam_vz_target_temp*beam_vz_target_temp) < 54.0/2.0) 
	{
	  beam_vx = beam_vx_temp;
	  beam_vy = beam_vy_temp;
	  beam_vz_target = beam_vz_target_temp;
	  flag = true;
	}
      else
	{
	  ReadBeamProfile(X_file, P_file);
	  z = GaussZ(54.0);
	  //target_posz = 1300.0 - 143.0 +z;
	  target_posz = 1300.0 +z;
	}
    }

  beam_vz = beam_vz_target + D.z();
  //G4cout<<"read file3 x: "<<beam_vx<<" y: "<<beam_vy<<" z: " << beam_vz<<G4endl;
  G4ThreeVector Xb(beam_vx, beam_vy, beam_vz);
  
  double beam_file_mag = P_file.mag();
  double beam_mag_set = P.mag()/0.9; //for beam 0.9 GeV/c
  G4ThreeVector beamp(P_file.x()*beam_mag_set, P_file.y()*beam_mag_set, P_file.z()*beam_mag_set);

  //beam rotate angle                                                                                
  G4double rotate_angle = -7.8*degree; //1.8 GeV/c && field = 1.0 case                             
  //G4double rotate_angle = -2.664*degree; //1.8 GeV/c && field = 0.85 case                          
  //G4double rotate_angle = -2.507*degree; //1.8 GeV/c && field = 0.8 case                           
  //G4double rotate_angle = -2.35*degree; //1.8 GeV/c && field = 0.75 case                           
  //G4double rotate_angle = 0.0*degree;                                                              
  //G4double rotate_angle = 3.135*degree; //1.8 GeV/c && field = -1.0 case                           

  G4ThreeVector beamp_rotate = BeamMomRotate( beamp, rotate_angle);
  G4ThreeVector Pb = beamp_rotate;

  anaMan_->SetBeam(1, Xb, beamp_rotate);

  //G4cout<<"Start Generate Test2"<<G4endl;
  /// mother particle momentum //
  double mass_km = 0.493677;
  double mass_proton = 0.938272081;
  G4LorentzVector lv_beam;
  G4LorentzVector lv_target;
  G4LorentzVector lv_particle;
  double pbeam = Pb.mag();

  lv_beam.setX(Pb.x());
  lv_beam.setY(Pb.y());
  lv_beam.setZ(Pb.z());
  lv_beam.setE(sqrt(mass_km*mass_km + pbeam*pbeam)); 

  lv_target.setX(0.0);
  lv_target.setY(0.0);
  lv_target.setZ(0.0);
  lv_target.setE(sqrt(mass_proton*mass_proton)); 
  
  lv_particle = lv_beam + lv_target;

  // make mother particle //
  EvtParticle* lam1663(0);
  static EvtId LAM1663 = EvtPDL::getId(std::string("Lambda(1663)0"));
  G4LorentzVector LvLam1663;
  G4ThreeVector TVp (lv_particle.x(), lv_particle.y(), lv_particle.z());
  G4ThreeVector TVx (Xb.x(), Xb.y(), Xb.z());
  //double mass_lam1663 = EvtPDL::getMass(LAM1663);
  double mass_lam1663 = sqrt((lv_beam.e()+lv_target.e())*(lv_beam.e()+lv_target.e()) - pbeam*pbeam);

  // check total energy //
  if(mass_lam1663 < 1.115683 + 0.547862 )
    {
      //G4cout<<"### Beam momentum is not enough to generate Lam1663 ###"<<G4endl;
      //G4cout<<"Beam Momentum: "<< pbeam << " GeV/c" << G4endl;
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
  GenerateDecay(anEvent, evtGenerator, lam1663, Xb);
}

void PrimaryGeneratorAction::GenerateLambdaStarBrCS(G4Event* anEvent, EvtGen *evtGenerator, G4ThreeVector D, G4ThreeVector P)
{
  G4ThreeVector beamx = D;

  G4ThreeVector X_file, P_file;
  ReadBeamProfile(X_file, P_file);

  double z = GaussZ(54.0); //at target center;
  //G4cout<<"read file1 x: "<<X_file.x()<<" y: "<<X_file.y()<<" z: " << z <<G4endl; 

  bool flag = false;
  //double target_posz = 1300.0 - 143.0 + z; //in experimental hall
  double target_posz = 1300.0 + z; //in experimental hall

  double beam_vx = 0.0;
  double beam_vy = 0.0;
  double beam_vz = 0.0;
  double beam_vz_target = 0.0;

  while ( !flag )
    {
      double u = P_file.x()/P_file.z();
      double v = P_file.y()/P_file.z();

      double beam_vx_temp = X_file.x() + target_posz*u + D.x();
      double beam_vy_temp = X_file.y() + target_posz*v + D.y();
      double beam_vz_target_temp = z; //center is SC magnet center
      //G4cout<<"read file2 x: "<<beam_vx_temp<<" y: "<<beam_vy_temp<<" z: " << beam_vz_target_temp<<G4endl;
      
      if( fabs(beam_vy_temp) < 50.0 && sqrt(beam_vx_temp*beam_vx_temp + beam_vz_target_temp*beam_vz_target_temp) < 54.0/2.0) 
	{
	  beam_vx = beam_vx_temp;
	  beam_vy = beam_vy_temp;
	  beam_vz_target = beam_vz_target_temp;
	  flag = true;
	}
      else
	{
	  ReadBeamProfile(X_file, P_file);
	  z = GaussZ(54.0);
	  //target_posz = 1300.0 - 143.0 +z;
	  target_posz = 1300.0 +z;
	}
    }

  double beam_file_mag = P_file.mag();
  double beam_mag_set = P.mag()/0.9; //for beam 0.9 GeV/c
  G4ThreeVector beamp(P_file.x()*beam_mag_set, P_file.y()*beam_mag_set, P_file.z()*beam_mag_set);

  beam_vz = beam_vz_target + D.z();
  //G4cout<<"read file3 x: "<<beam_vx<<" y: "<<beam_vy<<" z: " << beam_vz<<G4endl;
  G4ThreeVector Xb(beam_vx, beam_vy, beam_vz);

  //beam rotate angle                                                                                
  G4double rotate_angle = -7.8*degree; //1.8 GeV/c && field = 1.0 case                             
  //G4double rotate_angle = -2.664*degree; //1.8 GeV/c && field = 0.85 case                          
  //G4double rotate_angle = -2.507*degree; //1.8 GeV/c && field = 0.8 case                           
  //G4double rotate_angle = -2.35*degree; //1.8 GeV/c && field = 0.75 case                           
  //G4double rotate_angle = 0.0*degree;                                                              
  //G4double rotate_angle = 3.135*degree; //1.8 GeV/c && field = -1.0 case                           

  G4ThreeVector beamp_rotate = BeamMomRotate( beamp, rotate_angle);
  G4ThreeVector Pb = beamp_rotate;

  //cross-section
  double cross_section = bw_swave(beamp.mag());
  double max_cross_section = 1.3923;
  double comp = CLHEP::RandFlat::shoot(0.,1.0);
  //double comp = CLHEP::RandFlat::shoot(0., max_cross_section);
  if(comp > cross_section/max_cross_section)
    {
      anaMan_->SetBeam(1, Xb, beamp_rotate);
      return;
    }
  anaMan_->SetBeam(1, Xb, beamp_rotate);
  //G4cout<<"Start Generate Test2"<<G4endl;
  /// mother particle momentum //
  double mass_km = 0.493677;
  double mass_proton = 0.938272081;
  G4LorentzVector lv_beam;
  G4LorentzVector lv_target;
  G4LorentzVector lv_particle;
  double pbeam = Pb.mag();

  lv_beam.setX(Pb.x());
  lv_beam.setY(Pb.y());
  lv_beam.setZ(Pb.z());
  lv_beam.setE(sqrt(mass_km*mass_km + pbeam*pbeam)); 

  lv_target.setX(0.0);
  lv_target.setY(0.0);
  lv_target.setZ(0.0);
  lv_target.setE(sqrt(mass_proton*mass_proton)); 
  
  lv_particle = lv_beam + lv_target;

  // make mother particle //
  EvtParticle* lam1663(0);
  static EvtId LAM1663 = EvtPDL::getId(std::string("Lambda(1663)0"));
  G4LorentzVector LvLam1663;
  G4ThreeVector TVp (lv_particle.x(), lv_particle.y(), lv_particle.z());
  G4ThreeVector TVx (Xb.x(), Xb.y(), Xb.z());
  //double mass_lam1663 = EvtPDL::getMass(LAM1663);
  double mass_lam1663 = sqrt((lv_beam.e()+lv_target.e())*(lv_beam.e()+lv_target.e()) - pbeam*pbeam);

  // check total energy //
  if(mass_lam1663 < 1.115683 + 0.547862 )
    {
      //G4cout<<"### Beam momentum is not enough to generate Lam1663 ###"<<G4endl;
      //G4cout<<"Beam Momentum: "<< pbeam << " GeV/c" << G4endl;
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
  GenerateDecay(anEvent, evtGenerator, lam1663, Xb);
}






double PrimaryGeneratorAction::bw_swave(double beamp)
{
  double m_lambda = 1.11568;
  double m_proton = 0.938272;
  double m_kaon = 0.493677;
  double m_eta = 0.547853;
  double e_cm = sqrt(m_kaon*m_kaon + m_proton*m_proton + 2.0*m_proton*sqrt(m_kaon*m_kaon+beamp*beamp));

  if(e_cm <= m_lambda + m_eta) return 0.;

  double pl2  = 0.00715;
  double width_a = 0.004;
  double width_b = 0.005;
  double width_c = 0.024;
  double p_eta_cm = sqrt((e_cm*e_cm-m_eta*m_eta-m_lambda*m_lambda)*(e_cm*e_cm-m_eta*m_eta-m_lambda*m_lambda) - 4.0*m_eta*m_eta*m_lambda*m_lambda)/2.0/e_cm ;
  double p_beam_cm = sqrt((e_cm*e_cm-m_kaon*m_kaon-m_proton*m_proton)*(e_cm*e_cm-m_kaon*m_kaon-m_proton*m_proton)- 4.0*m_kaon*m_kaon*m_proton*m_proton)/2.0/e_cm;
  double par0 = 1.05836e-03;
  double m_lam1670 = 1.67;

  double Gb = width_b*sqrt((e_cm - m_eta - m_lambda)/(m_lam1670 - m_eta - m_lambda));
  double Gt = Gb + width_a + width_c;
  //double Gt = Gb + width_a + par[1];                                                                             

  //double cs = par[0]*pl2*width_a*Gb/((e_cm-m_lam1670)*(e_cm-m_lam1670) + Gt*Gt/4.0);                             
  //double cs = par[0]*pl2*width_a*Gb/((e_cm-m_lam1670)*(e_cm-m_lam1670) + width_c*width_c/4.0);                   
  double cs = par0*p_eta_cm/p_beam_cm/((e_cm-m_lam1670)*(e_cm-m_lam1670) + width_c*width_c/4.0);
  return cs;
}
