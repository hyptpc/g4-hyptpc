// ====================================================================
//   CH1SD.hh
//
// ====================================================================
#ifndef CH1_SD_H
#define CH1_SD_H

#include "G4VSensitiveDetector.hh"
#include "CH1Hit.hh"

class G4HCofThisEvent;
class G4Step;
class G4TouchableHistory;

class CH1SD : public G4VSensitiveDetector {
  private:
  //G4THitsCollection<CH1Hit>* hitsCollection;
  CH1HitsCollection *Ch1Collection;

  public:
    CH1SD( G4String name );
    virtual ~CH1SD();

    // virtual methods
    virtual G4bool ProcessHits(G4Step* aStep, G4TouchableHistory* ROhist);
    virtual void Initialize(G4HCofThisEvent* HCTE);
    virtual void EndOfEvent(G4HCofThisEvent* HCTE);

    virtual void DrawAll();
    virtual void PrintAll();

};

#endif
