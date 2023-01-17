#include "SACStackingAction.hh"

#include "G4VProcess.hh"

#include "G4ParticleDefinition.hh"
#include "G4VDiscreteProcess.hh"
#include "G4ParticleChange.hh"
#include "G4ParticleTypes.hh"
#include "G4Track.hh"
#include "G4ios.hh"
#include "G4ClassificationOfNewTrack.hh"

int gCerenkovCounter;
double wavecheck[10000];


SACStackingAction::SACStackingAction()
  : G4UserStackingAction(),
    fScintillationCounter(0), fCerenkovCounter(0)
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

SACStackingAction::~SACStackingAction()
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4ClassificationOfNewTrack
SACStackingAction::ClassifyNewTrack(const G4Track * aTrack)
{

  G4ClassificationOfNewTrack classification = fWaiting;
  const G4double h = 6.628e-34;
  const G4double c = 3.0e+8;
  
  if(aTrack->GetDefinition() == G4OpticalPhoton::OpticalPhotonDefinition())
    { // particle is optical photon
      if(aTrack->GetTouchable()->GetVolume()->GetCopyNo() ==9999){ //Get Cherenkov photon only from aerogel
	if(aTrack->GetParentID()>0)
	  { // particle is secondary
	    if(aTrack->GetCreatorProcess()->GetProcessName() == "Scintillation")
	      fScintillationCounter++;
	    if(aTrack->GetCreatorProcess()->GetProcessName() == "Cerenkov"){
	      if(((h*c)/(aTrack->GetTotalEnergy()*1.6e-13))*(1e+9)>190&&((h*c)/(aTrack->GetTotalEnergy()*1.6e-13))*(1e+9)<900){ //Only consider photons with wavelength that PMT can detect
		fCerenkovCounter++;
		gCerenkovCounter++;
		wavecheck[gCerenkovCounter] = ((h*c)/(aTrack->GetTotalEnergy()*1.6e-13))*(1e+9);

	      }
	      else return fKill;
	    }
	  }
      }
    }
		   
  return fUrgent;
}


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void SACStackingAction::NewStage()
{
  std::cout << "Number of Cerenkov photons produced in this event : "
	    << fCerenkovCounter << std::endl;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void SACStackingAction::PrepareNewEvent()
{
  fCerenkovCounter = 0;
  gCerenkovCounter = 0;
  for(int i=0;i<10000;i++)wavecheck[i] = 0;

}
