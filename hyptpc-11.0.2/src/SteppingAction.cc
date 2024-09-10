// -*- C++ -*-

#include "SteppingAction.hh"

#include <unordered_map>

#include <G4Material.hh>
#include <G4ParticleDefinition.hh>
#include <G4ParticleTypes.hh>
#include <G4SteppingManager.hh>
#include <G4Step.hh>
#include <G4StepPoint.hh>
#include <G4Track.hh>
#include <G4TrackStatus.hh>
#include <G4VPhysicalVolume.hh>

#include "ConfMan.hh"
#include "PrintHelper.hh"
#include "AnaManager.hh"

namespace
{
auto& gAnaMan = AnaManager::GetInstance();
const auto& gConf = ConfMan::GetInstance();
}

//_____________________________________________________________________________
SteppingAction::SteppingAction()
  : G4UserSteppingAction()
{
}

//_____________________________________________________________________________
SteppingAction::~SteppingAction()
{
}

//_____________________________________________________________________________
void
SteppingAction::UserSteppingAction(const G4Step* theStep)
{
  static const G4bool KillStepInIron = gConf.Get<G4bool>("KillStepInIron");

  auto theTrack = theStep->GetTrack();
  auto theParticle = theTrack->GetParticleDefinition();
  auto particleName = theParticle->GetParticleName();
  auto particlePdgCode = theParticle->GetPDGEncoding();
  auto prePoint = theStep->GetPreStepPoint();
  auto prePV = prePoint->GetPhysicalVolume();
  auto prePVName = prePV->GetName();
  auto postPoint = theStep->GetPostStepPoint();
  auto theProcess = postPoint->GetProcessDefinedStep()->GetProcessName();
  // check if it is alive
  //  if(theTrack->GetTrackStatus() != fAlive) { return; }

  // check if it is primary
  //  if(theTrack->GetParentID() != 0) { return; }

  // check if it is NOT muon
  //  auto definition = theTrack->GetDefinition();
  //  if((definition == G4MuonPlus::MuonPlusDefinition()) ||
  //      (definition == G4MuonMinus::MuonMinusDefinition()))
  //  { return; }

  //  G4cout<<"start stepping action:"<<prePVName<<G4endl;

  // -- check decay particle -----
  std::pair<G4String, G4String> previous_particle = gAnaMan.GetPreviousParticle();
  G4int generator = gAnaMan.GetNextGenerator();
  if (previous_particle.second == "Decay" && previous_particle.first == gAnaMan.GetFocusParticle(generator) ){
    gAnaMan.SetDecayParticleCode( particlePdgCode );
  }
  gAnaMan.SetPreviousParticle(particleName, theProcess);
  
  
#ifdef DEBUG
  PrintHelper helper(3, std::ios::fixed, G4cout);
  auto time = prePoint->GetGlobalTime();
  auto track_id = theTrack->GetTrackID();
  std::stringstream particle_ss;
  particle_ss << particleName << "(" << track_id << ") ";
  if(theProcess != "eIoni" &&
     theProcess != "hIoni" &&
     theProcess != "msc" &&
     theProcess != "eBeam" &&
     theProcess != "Transportation"){
    G4cout << "   " << time/CLHEP::ns << " ns : "
	   << particle_ss.str() << theProcess << G4endl;
  }

  auto secondary = theStep->GetSecondaryInCurrentStep();
  for(const auto& s : *secondary){
    auto particle = s->GetDefinition();
    auto name = particle->GetParticleName();
    // auto type = particle->GetParticleType();
    if(true
       || (particleName == "lambda" && name == "proton")){
      G4cout << "   " << particle_ss.str() << "\tP" << prePoint->GetMomentum()
	     << " X" << prePoint->GetPosition()
	     << "\t-> " << name << " P" << s->GetMomentum() << G4endl;
    }
  }

  if(false
     && particleName == "proton"){
    auto preMaterial = prePoint->GetMaterial();
    G4double edep = theStep->GetTotalEnergyDeposit();
    G4cout << "   " << particleName << " " << theProcess
	   << " " << theTrack->GetTrackStatus()
	   << " " << preMaterial->GetName()
	   << " x=" << prePoint->GetPosition()
	   << " p=" << prePoint->GetMomentum()
	   << " edep=" << edep
	   << G4endl;
  }
#endif

  if(KillStepInIron){
    auto preMaterial = prePoint->GetMaterial();
    if(preMaterial->GetName() == "Iron"){
      theTrack->SetTrackStatus(fStopAndKill);
      return;
    }
  }

  if(particleName == "e-" || particleName == "e+"){
    theTrack->SetTrackStatus(fStopAndKill);
    return;
  }

  // if(prePVName.contains("Coil") || prePVName.contains("Guard")){
  //   theTrack->SetTrackStatus(fStopAndKill);
  //   return;
  // }
}
