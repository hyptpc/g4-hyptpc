// -*- C++ -*-

#ifndef BAC_SD_HH
#define BAC_SD_HH

#include <G4VSensitiveDetector.hh>

#include "BACHit.hh"

class G4HCofThisEvent;
class G4Step;
class G4TouchableHistory;

//_____________________________________________________________________________
class BACSD : public G4VSensitiveDetector
{
public:
  BACSD(const G4String& name);
  virtual ~BACSD();

private:
  G4THitsCollection<BACHit>* hitsCollection;

public:
  virtual G4bool ProcessHits(G4Step* aStep, G4TouchableHistory* ROhist);
  virtual void   Initialize(G4HCofThisEvent* HCTE);
  virtual void   EndOfEvent(G4HCofThisEvent* HCTE);
  virtual void   DrawAll();
  virtual void   PrintAll();
};

#endif
