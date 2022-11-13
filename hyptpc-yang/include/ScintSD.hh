// ====================================================================
//   ScintSD.hh
//
// ====================================================================
#ifndef SCINT_SD_H
#define SCINT_SD_H

#include "G4VSensitiveDetector.hh"
#include "ScintHit.hh"

class G4HCofThisEvent;
class G4Step;
class G4TouchableHistory;

class ScintSD : public G4VSensitiveDetector {
  private:
    G4THitsCollection<ScintHit>* hitsCollection;

  public:
    ScintSD(const G4String& name);
    virtual ~ScintSD();

    // virtual methods
    virtual G4bool ProcessHits(G4Step* aStep, G4TouchableHistory* ROhist);
    virtual void Initialize(G4HCofThisEvent* HCTE);
    virtual void EndOfEvent(G4HCofThisEvent* HCTE);

    virtual void DrawAll();
    virtual void PrintAll();

};

#endif
