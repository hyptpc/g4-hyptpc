// ====================================================================
//   ACSD.hh
//
// ====================================================================
#ifndef AC_SD_H
#define AC_SD_H

#include "G4VSensitiveDetector.hh"
#include "ACHit.hh"

class G4HCofThisEvent;
class G4Step;
class G4TouchableHistory;

class ACSD : public G4VSensitiveDetector {
  private:
    G4THitsCollection<ACHit>* hitsCollection;

  public:
    ACSD(const G4String& name);
    virtual ~ACSD();

    // virtual methods
    virtual G4bool ProcessHits(G4Step* aStep, G4TouchableHistory* ROhist);
    virtual void Initialize(G4HCofThisEvent* HCTE);
    virtual void EndOfEvent(G4HCofThisEvent* HCTE);

    virtual void DrawAll();
    virtual void PrintAll();

};

#endif
