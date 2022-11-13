// ====================================================================
//   FTOFSD.hh
//
// ====================================================================
#ifndef FTOF_SD_H
#define FTOF_SD_H

#include "G4VSensitiveDetector.hh"
#include "FTOFHit.hh"

class G4HCofThisEvent;
class G4Step;
class G4TouchableHistory;

class FTOFSD : public G4VSensitiveDetector {
  private:
  //G4THitsCollection<FTOFHit>* hitsCollection;
  FTOFHitsCollection *FtofCollection;

  public:
    FTOFSD( G4String name );
    virtual ~FTOFSD();

    // virtual methods
    virtual G4bool ProcessHits(G4Step* aStep, G4TouchableHistory* ROhist);
    virtual void Initialize(G4HCofThisEvent* HCTE);
    virtual void EndOfEvent(G4HCofThisEvent* HCTE);

    virtual void DrawAll();
    virtual void PrintAll();

};

#endif
