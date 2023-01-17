#ifndef SACStackingAction_H
#define SACStackingAction_H 1

#include "globals.hh"
#include "G4UserStackingAction.hh"



class G4HCofThisEvent;


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

class SACStackingAction : public G4UserStackingAction
{
  public:
    SACStackingAction();
    virtual ~SACStackingAction();

  public:
    virtual G4ClassificationOfNewTrack ClassifyNewTrack(const G4Track* aTrack);
    virtual void NewStage();
    virtual void PrepareNewEvent();


  private:
    G4int fScintillationCounter;
    G4int fCerenkovCounter;

  
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif
