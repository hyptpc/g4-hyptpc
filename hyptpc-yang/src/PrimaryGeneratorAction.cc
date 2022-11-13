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

PrimaryGeneratorAction::PrimaryGeneratorAction( DetectorConstruction *det, 
						AnalysisManager *analysisManager,
						EvtGen *evtgen)
  : G4VUserPrimaryGeneratorAction(), det_(det), 
    anaMan_(analysisManager),
    evtgen_(evtgen)
{

  G4int n_particle = 1;
  bp_nAccess = 0; 
  particleGun = new G4ParticleGun(n_particle);
  particleTable = G4ParticleTable::GetParticleTable();

}

PrimaryGeneratorAction::~PrimaryGeneratorAction()
{
  //if(BP_) delete BP_;
  DeleteGuns();
}

void PrimaryGeneratorAction::GeneratePrimaries( G4Event *anEvent )
{
  //std::cout<<"Beam generated "<< std::endl;
  G4int event_num = anEvent->GetEventID()+1;
  //G4cout<<"Events: "<< event_num <<G4endl;
  if(event_num% 1000 == 0)
    {
      G4cout<<"Event# "<<event_num<<G4endl;
    }

  ConfMan *confMan = ConfMan::GetConfManager();
  int reactionMode = confMan->ReactionMode();

  //vertex point and beam momentum
  bpx_ = confMan->GetBeamPX();
  bpy_ = confMan->GetBeamPY();
  bpz_ = confMan->GetBeamPZ();
  bvx_ = confMan->GetBeamVX();
  bvy_ = confMan->GetBeamVY();
  bvz_ = confMan->GetBeamVZ();
  G4ThreeVector D(bvx_, bvy_, bvz_);
  G4ThreeVector P(bpx_, bpy_, bpz_);
  //anaMan_->SetBeam(1, D, P);

  //G4cout<<"Reaction Mode: "<< reactionMode << G4endl;
  switch(reactionMode)
    {
    case 1: GenerateTest(anEvent, D, P); break;
    case 2: GenerateTest2(anEvent, evtgen_, D, P); break;
    case 3: GenerateElasticPiM(anEvent, D, P); break;
    case 4: GenerateElasticPiP(anEvent, D, P); break;
    case 5: GenerateBeamPiM(anEvent, D, P); break;
    case 6: GenerateBeamPiP(anEvent, D, P); break;
    case 7: GeneratePiPSchannel(anEvent, evtgen_, D, P); break;
    case 8: GenerateBeamKaonM(anEvent, D, P); break;
    case 9: GenerateBeamKaonMBeamRoot(anEvent); break; //read beam profiel from Root file
    case 10: GenerateHDiBaryonBeamRoot_test(anEvent); break; //read beam profiel from Root file
    case 11: GenerateHDiBaryon_test(anEvent, evtgen_, D, P); break; //K-pp -> HK+ reaction
    case 12: GenerateScatteredKPlus(anEvent, D, P); break; //Random K+ generation
      //case 13: GenerateXi_test(anEvent, evtgen_, D, P); break; //K-p -> Xi-K+ reaction
    case 14: GenerateBeamKaonMBr(anEvent, D, P); break; // K- beam beam file
    case 15: GenerateLambdaStarBr(anEvent, evtgen_, D, P); break; //K-p -> Lambda(1663) beam file
    case 16: GenerateLambdaStarBrCS(anEvent, evtgen_, D, P); break; //K-p -> Lambda(1663) beam file cross-section
    }
}

void PrimaryGeneratorAction::DeleteGuns()
{
}

G4ParticleGun * PrimaryGeneratorAction::chooseGun( G4int Pid )
{
  //G4String name=PIDParticleName(Pid);
  //return chooseGun( name );
}

G4ParticleGun * PrimaryGeneratorAction::chooseGun( const G4String & name  )
{
}


void PrimaryGeneratorAction::GenerateDecay(G4Event* anEvent, EvtGen *evtGenerator, EvtParticle* particle, G4ThreeVector D)
{

                                                                     
  static EvtStdHep evtstdhep;
  static EvtSecondary evtsecondary;

  EvtId        list_of_stable[10];
  EvtParticle* stable_parent[10];

  list_of_stable[0]=EvtId(-1,-1);
  stable_parent[0]=0;

  evtsecondary.init();
  evtstdhep.init();
  evtGenerator -> generateDecay(particle);
  particle->makeStdHep(evtstdhep,evtsecondary,list_of_stable);

  int npart = evtstdhep.getNPart();

  bool generate_flag = false;
  if(npart<100) generate_flag = true;
  while(!generate_flag)
    {
      G4cout<<"!!!Particles in EvtGen is more than 99."<<G4endl;
      evtsecondary.init();
      evtstdhep.init();
      evtGenerator -> generateDecay(particle);
      particle->makeStdHep(evtstdhep,evtsecondary,list_of_stable);
      if(evtstdhep.getNPart() < 100) generate_flag = true;
    }

  int j;
  int istat;
  int partnum;
  double px,py,pz,e,m;
  double x,y,z,t;

  EvtVector4R p4,x4;
  int n_beam = 0; // number of beams

  
  
  for(int i=0;i<evtstdhep.getNPart();i++)
    {
      j=i+1;
      int jmotherfirst=evtstdhep.getFirstMother(i)+1;
      int jmotherlast=evtstdhep.getLastMother(i)+1;
      int jdaugfirst=evtstdhep.getFirstDaughter(i)+1;
      int jdauglast=evtstdhep.getLastDaughter(i)+1;
      
      partnum=evtstdhep.getStdHepID(i);
  
      istat=evtstdhep.getIStat(i);
  
      p4=evtstdhep.getP4(i);
      x4=evtstdhep.getX4(i);
   
      px=p4.get(1);
      py=p4.get(2);
      pz=p4.get(3);
      e=p4.get(0);
   
      x=x4.get(1)+D.x();
      y=x4.get(2)+D.y();
      z=x4.get(3)+D.z();
      //t=x4.get(0)+D.t();
      t=x4.get(0)*0.001/2.999792458 *1.e1; //mm/c --> ns
      m=p4.mass();

      EvtVector4R evx4, evp4;
      evx4.set(1, x);
      evx4.set(2, y);
      evx4.set(3, z);
      evx4.set(0, t);

      evp4.set(1, px);
      evp4.set(2, py);
      evp4.set(3, pz);
      evp4.set(0, e);

      bool beam_flag = false;
      int tr_id = -1;
      //if(jdaugfirst ==0 && jdauglast ==0 && (partnum == 2212 || partnum == -211) && jmotherfirst == 2)
      if(jdaugfirst ==0 && jdauglast ==0)
	{
	  makeGun(anEvent, partnum, evx4, evp4);
	  beam_flag = true;
	  n_beam++;
	}
      if(beam_flag == true)
	{
	  tr_id = n_beam;
	}
      anaMan_->SetEvtGen(j, partnum, jmotherfirst, jmotherlast, jdaugfirst, jdauglast, tr_id, evx4, evp4);
    }

#if 0
  G4cout<<"############# Particle decay table  ##############"<<G4endl;
  G4cout<<"Npart: "<<npart<<G4endl;
  for(int i=0;i<evtstdhep.getNPart();i++)
    {
      j=i+1;
      int jmotherfirst=evtstdhep.getFirstMother(i)+1;
      int jmotherlast=evtstdhep.getLastMother(i)+1;
      int jdaugfirst=evtstdhep.getFirstDaughter(i)+1;
      int jdauglast=evtstdhep.getLastDaughter(i)+1;

      p4=evtstdhep.getP4(i);
      x4=evtstdhep.getX4(i);
   
      px=p4.get(1);
      py=p4.get(2);
      pz=p4.get(3);
      e=p4.get(0);
   
      x=x4.get(1)+D.x();
      y=x4.get(2)+D.y();
      z=x4.get(3)+D.z();
      t=x4.get(0);
      m=p4.mass();
      
      partnum=evtstdhep.getStdHepID(i);
      G4cout<<"ID: " << j<< "  Particle Num: "<<partnum<<"  mf: "<<jmotherfirst<< "  ml: "<<jmotherlast << "  df: "<<jdaugfirst << "  dl: "<<jdauglast<<G4endl;
      G4cout<< "   p: "<<(float)sqrt(px*px+py*py+pz*pz) << " e: " << (float)e << " t: "<< (float)t<< " m: "<< (float)m <<G4endl;
    }

  G4cout<<"##################################################"<<G4endl;
#endif
  particle->deleteTree();

}

void PrimaryGeneratorAction::makeGun(G4Event* anEvent, int partnum, EvtVector4R x4, EvtVector4R p4)
{
  particleGun -> SetParticleDefinition (particleTable -> FindParticle(partnum));
  G4ThreeVector beamx ( x4.get(1), x4.get(2), x4.get(3));
  G4ThreeVector beamp ( p4.get(1), p4.get(2), p4.get(3));
  G4ThreeVector beampu = beamp/beamp.mag();
  G4double energy = (p4.get(0) - p4.mass())*GeV;
  //G4double energy = beamp.mag()*GeV;
  double time = x4.get(0); 
  G4double g4time =  time * ns;

  //G4cout<<"[Test]pid: "<<partnum << " mass: "<<p4.mass() <<G4endl;
  /*
  std::cout<<"Beam energy: "<< energy << std::endl;
  G4cout<<"###############  beam generated  #################" << G4endl;
  G4cout<<"particle: "<<partnum<< " vertex(x,y,z) ("<< beamx.x() << ", " << beamx.y() << ", " << beamx.z() << ") "
  	<< "energy: " <<energy/GeV<< " momentum(x,y,z) ("<< p4.get(1) << ", " << p4.get(2) << ", " <<p4.get(3) << ") "<< beamp.mag()<< " time: "<<time
  	<< G4endl;
  G4cout<<"##################################################" << G4endl;
  */
  particleGun->SetParticleMomentumDirection ( beampu );
  //particleGun->SetParticleMomentum ( beamp );
  particleGun->SetParticleTime ( g4time );
  particleGun->SetParticlePosition( beamx );
  particleGun->SetParticleEnergy( energy );
  //if(partnum == 321 || partnum == 2212 || partnum == -211 )particleGun->GeneratePrimaryVertex( anEvent);
  //if(partnum == 321 )particleGun->GeneratePrimaryVertex( anEvent);
  particleGun->GeneratePrimaryVertex( anEvent);
}


//Gauss Position
G4ThreeVector PrimaryGeneratorAction::GaussPosition( double sigx, double sigy, double hz )
{
  double x = 0. , y = 0., z = 0.0;//target center
  x += G4RandGauss::shoot(0.0, sigx);; 
  y += G4RandGauss::shoot(0.0, sigy);; 
  z += (G4UniformRand() - 0.5) * hz;; 
  
  return G4ThreeVector (x, y, z);

}

//Gauss Momentum
G4ThreeVector PrimaryGeneratorAction::GaussMomentum( double sigx, double sigy, double sigz )
{
  double x = 0. , y = 0., z = 0.0;//target center
  x += G4RandGauss::shoot(0.0, sigx);; 
  y += G4RandGauss::shoot(0.0, sigy);; 
  z += G4RandGauss::shoot(0.0, sigz);; 
  
  return G4ThreeVector (x, y, z);
}


//Read Beam Profile
void PrimaryGeneratorAction::ReadBeamProfile( G4ThreeVector & X, G4ThreeVector & P )
{

  ConfMan *conf = ConfMan::GetConfManager();
  G4String beamfilename = conf->GetBeamProfileFile();
  TFile *beam_file = new TFile(beamfilename, "R");
  TTree *beam_tree = (TTree*)beam_file->Get("tr");

  int ntK18;
  double pointInx[1];
  double pointIny[1];
  double pointInz[1];
  double pInx[1];
  double pIny[1];
  double pInz[1];


  beam_tree->SetBranchAddress("ntK18",&ntK18);
  beam_tree->SetBranchAddress("pointInx",  &pointInx);
  beam_tree->SetBranchAddress("pointIny",  &pointIny);
  beam_tree->SetBranchAddress("pointInz",  &pointInz);

  beam_tree->SetBranchAddress("pInx",  &pInx);
  beam_tree->SetBranchAddress("pIny",  &pIny);
  beam_tree->SetBranchAddress("pInz",  &pInz);

  bp_file_ndata = beam_tree->GetEntries();
  if(bp_file_ndata == bp_nAccess) bp_nAccess = 0;
  
  beam_tree->GetEntry(bp_nAccess);
 
  G4ThreeVector TVp(pInx[0], pIny[0], pInz[0]);
  G4ThreeVector TVx(pointInx[0], pointIny[0], pointInz[0]);


  X=TVx;
  P=TVp;
  bp_nAccess++;

  //G4cout<<"[PrimaryGeneratorAction]nAccess: "<< bp_nAccess <<G4endl;
  //G4cout<<"[PrimaryGeneratorAction]x and p: "<< TVx <<"   "<< TVp <<G4endl;

  beam_file->Close();

}

//scattering 2body
bool PrimaryGeneratorAction::Scattering2Body( double Mi1, double Mi2, double Mf1, double Mf2, const G4ThreeVector & Pini, G4ThreeVector & Pf1,  G4ThreeVector & Pf2, double cosx )
{
  //  std::cout << "Pini=" << Pini/GeV << "GeV/c" << std::endl;                                   
  
  double Mc=sqrt(Mi1*Mi1 + Mi2*Mi2 + 2.*Mi2*sqrt(Mi1*Mi1+Pini.mag2()) );
  return Decay2Body( Mc, Mf1, Mf2, Pini, Pf1, Pf2, cosx );
}

bool PrimaryGeneratorAction::Decay2Body( double Mini, double Mf1, double Mf2, const G4ThreeVector & Pini, G4ThreeVector & Pf1,  G4ThreeVector & Pf2, double cosx )
{
  if(Mini<Mf1+Mf2)
    {
      std::cerr << "Mini < Mf1+Mf2 Mini=" << Mini/GeV << "GeV/c2 "
		<< "Mf1=" <<  Mf1/GeV << "GeV/c2 "
		<< "Mf2=" <<  Mf2/GeV << "GeV/c2 " << std::endl;
      return false;
    }
  
  G4ThreeVector beta( Pini/sqrt(Mini*Mini+Pini.mag2()) );
  
  double Ecmf1=(Mini*Mini-Mf2*Mf2+Mf1*Mf1)/Mini*0.5;
  double Ecmf2=(Mini*Mini-Mf1*Mf1+Mf2*Mf2)/Mini*0.5;
  double Pcm=sqrt((Mini*Mini-(Mf1+Mf2)*(Mf1+Mf2))*
                  (Mini*Mini-(Mf1-Mf2)*(Mf1-Mf2)))/Mini*0.5;
  
  //G4ThreeVector UnitDir=generator.GenerateDirection();
  
  double cost=cosx;
  double sint=sqrt(1.-cost*cost);
  double phi=G4UniformRand()*acos(-1.)*2.;
  double cosp=cos(phi), sinp=sin(phi);

  G4ThreeVector UnitDir( sint*cosp, sint*sinp, cost );

  G4ThreeVector Pcmf1 =  Pcm*UnitDir;
  G4ThreeVector Pcmf2 = -Pcm*UnitDir;
  
  G4LorentzVector LVf1( Pcmf1, Ecmf1), LVf2( Pcmf2, Ecmf2 );
  LVf1.boost(beta); Pf1=LVf1.vect();
  LVf2.boost(beta); Pf2=LVf2.vect();
  
  //  std::cout << "CosTCM=" << cost << " PhiCM=" << phi/degree                                   
  //        << " degree" << std::endl;                                                            
  //  std::cout << "PCM=" << Pcm/GeV << " " << Pcmf1/GeV << " -->"                                
  //                << Pf1/GeV << " " << Pf1.mag()/GeV << std::endl;                              
  
  return true;
}

//Legendre function
G4double PrimaryGeneratorAction::Legendre(G4int order, G4double x)
{
  G4double output=0;
  switch( order  )
    {
    case 0:
      output=1.;
      break;
    case 1:
      output= x;
      break;
    case 2:
      output= 1./2.*( 3.*x*x - 1. );
      break;
    case 3:
      output= 1./2.*( 5.*x*x*x - 3.*x );
      break;
    case 4:
      output= 1./8.*( 35.*pow(x,4) - 30.*pow(x,2) + 3. );
      break;
    case 5:
      output= 1./8.*( 63.*pow(x,5) - 70.*pow(x,3) + 15.*x );
      break;
    case 6:
      output= 1./16.*( 231.*pow(x,6) - 315.*pow(x,4) + 105.*pow(x,2) - 5. );
      break;
    case 7:
      output= 1./16.*( 429.*pow(x,7) - 693.*pow(x,5) + 315.*pow(x,3) - 35.*x );
      break;
    case 8:
      output= 1./128.*( 6435.*pow(x,8) - 12012.*pow(x,6) + 6930.*pow(x,4) - 1260.*pow(x,2) + 35. );
      break;
    case 9:
      output= 1./128.*( 12155.*pow(x,9) - 25740.*pow(x,7) + 18018.*pow(x,5) - 4620.*pow(x,3) + 315.*x );
      break;
    case 10:
      output= 1./256.*( 46189.*pow(x,10) - 109395.*pow(x,8) + 90090.*pow(x,6) - 30030.*pow(x,4) + 3465.*pow(x,2) - 63. );
      break;

    default:
      G4cout<<"error Legendre order:"<<order<<G4endl;
    }
  
  return output;
}

G4ThreeVector PrimaryGeneratorAction::BeamMomRotate(G4ThreeVector beam, G4double x)
{
  G4ThreeVector beama;
  G4ThreeVector yunit(0., 1., 0.); //y direction rotation, y direction B field
  beam.rotate(x, yunit);

  beama = beam;
  
  return beama;
}
 
G4ThreeVector PrimaryGeneratorAction::GaussMomentum( double umean, double usig, double vmean, double vsig)
{
  double u_mean = tan(umean);
  double v_mean = tan(vmean);
  double u_sig = tan(usig);
  double v_sig = tan(vsig);
  double u = G4RandGauss::shoot(u_mean, u_sig) ;;
  double v = G4RandGauss::shoot(v_mean, v_sig) ;;
  
  double pmag = sqrt(u*u + v*v + 1.0);
  G4ThreeVector pdir (u/pmag, v/pmag, 1.0/pmag);
  
  //std::cout<<"U: "<<u<<" V: "<<v<<std::endl;
  return pdir;
  
}

double PrimaryGeneratorAction::GaussZ( double hz)
{
  double z = (G4UniformRand() - 0.5) * hz;; 
  return z;
}
