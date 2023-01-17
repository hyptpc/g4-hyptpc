#include "SACDetectorConstruction.hh"
#include "SACPrimaryGeneratorAction.hh"
#include "G4LogicalVolumeStore.hh"
#include "G4LogicalVolume.hh"
#include "G4Box.hh"
#include "G4RunManager.hh"
#include "G4ParticleGun.hh"
#include "G4ParticleTable.hh"
#include "G4ParticleDefinition.hh"
#include "G4SystemOfUnits.hh"


#include "TFile.h"
#include "TTree.h"
#include "Randomize.hh"

SACPrimaryGeneratorAction::SACPrimaryGeneratorAction(const G4String &par1, const G4String &par2)
  : G4VUserPrimaryGeneratorAction(), parameter1(par1), parameter2(par2)
{
  G4int n_particle = 1;
  
  fParticleGun = new G4ParticleGun(n_particle);

}

SACPrimaryGeneratorAction::~SACPrimaryGeneratorAction()
{
  
  delete fParticleGun;
}

void SACPrimaryGeneratorAction::GeneratePrimaries(G4Event* anEvent)
{

  G4int event_num = anEvent->GetEventID()+1;
  if(event_num% 1000 == 0)
    {
      G4cout<<"Event# "<<event_num<<G4endl;
    }



  pa1 = std::stod(parameter1);
  pa2 = std::stod(parameter2);
  G4double x_moving = pa1*mm;
  G4double y_moving = pa2*mm;
  
  //homogeneous test------------------------------------
  
  G4double momentum = 1.2; //GeV/c
  //G4double momentum = 4;
  G4double beam_size = 125;

  
  G4double tight_size_x = 10;
  G4double tight_size_y = 10;

  
  if(particle=="kaon"){
    energy = (sqrt(mass_kaonm*mass_kaonm+momentum*momentum) - mass_kaonm )*GeV;
    fParticleGun -> SetParticleDefinition (particleTable -> FindParticle("kaon-"));
  }
  if(particle=="pion"){
    energy = (sqrt(mass_pim*mass_pim+momentum*momentum) - mass_pim )*GeV;
    fParticleGun->SetParticleDefinition (particleTable -> FindParticle("pi-"));
  }
  if(particle=="muon"){
    energy = 4*GeV;
    fParticleGun->SetParticleDefinition (particleTable -> FindParticle("mu-"));
  }

  if(particle=="positron"){
    energy = 1081*MeV;
    fParticleGun->SetParticleDefinition (particleTable -> FindParticle("e+"));
  } 
    
  
  fParticleGun->SetParticleMomentumDirection ( G4ThreeVector(0,0,1) );

  G4double x = tight_size_x*0.5-G4UniformRand()*tight_size_x*mm;
  G4double y = tight_size_y*0.5-G4UniformRand()*tight_size_y*mm;
  fParticleGun->SetParticlePosition(G4ThreeVector(x+x_moving,y+y_moving,-5*cm) );

  fParticleGun->SetParticleEnergy(energy);
  fParticleGun->GeneratePrimaryVertex(anEvent);


  


}
