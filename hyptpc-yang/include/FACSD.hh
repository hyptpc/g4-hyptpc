// ====================================================================
//   FACSD.hh
//
// ====================================================================
#ifndef FAC_SD_H
#define FAC_SD_H 1

#include "G4VSensitiveDetector.hh"
#include "FACHit.hh"

class G4HCofThisEvent;
class G4Step;
class G4TouchableHistory;

class FACSD : public G4VSensitiveDetector 
{
private:
  FACHitsCollection *FacCollection;

public:
  FACSD( G4String name );
  ~FACSD();

  // virtual methods
  G4bool ProcessHits(G4Step* aStep, G4TouchableHistory* ROhist);
  void Initialize(G4HCofThisEvent* HCTE);
  void EndOfEvent(G4HCofThisEvent* HCTE);
  
  void DrawAll();
  void PrintAll();
  
};

#endif
