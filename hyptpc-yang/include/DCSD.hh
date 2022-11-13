// ====================================================================
//   DCSD.hh
//
// ====================================================================
#ifndef DC_SD_H
#define DC_SD_H

#include "G4VSensitiveDetector.hh"
#include "DCHit.hh"

class G4HCofThisEvent;
class G4Step;
class G4TouchableHistory;

class DCSD : public G4VSensitiveDetector {
  private:
  DCHitsCollection *DcCollection;
  //    G4THitsCollection<DCHit>* hitsCollection;

  public:
    DCSD( G4String name );
    virtual ~DCSD();

    // virtual methods
    virtual G4bool ProcessHits(G4Step* aStep, G4TouchableHistory* ROhist);
    virtual void Initialize(G4HCofThisEvent* HCTE);
    virtual void EndOfEvent(G4HCofThisEvent* HCTE);

    virtual void DrawAll();
    virtual void PrintAll();

};

#endif
