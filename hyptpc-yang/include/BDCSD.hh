// ====================================================================
//   BDCSD.hh
//
// ====================================================================
#ifndef BDC_SD_H
#define BDC_SD_H 1

#include "G4VSensitiveDetector.hh"
#include "BDCHit.hh"

class G4HCofThisEvent;
class G4Step;
class G4TouchableHistory;

class BDCSD : public G4VSensitiveDetector 
{
private:
  BDCHitsCollection *BdcCollection;

public:
  BDCSD( G4String name );
  ~BDCSD();

  // virtual methods
  G4bool ProcessHits(G4Step* aStep, G4TouchableHistory* ROhist);
  void Initialize(G4HCofThisEvent* HCTE);
  void EndOfEvent(G4HCofThisEvent* HCTE);
  
  void DrawAll();
  void PrintAll();
  
};

#endif
