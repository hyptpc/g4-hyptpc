// -*- C++ -*-

#ifndef KVC_SD_HH
#define KVC_SD_HH

#include <G4VSensitiveDetector.hh>

#include "KVCHit.hh"

class G4HCofThisEvent;
class G4Step;
class G4TouchableHistory;

//_____________________________________________________________________________
class KVCSD : public G4VSensitiveDetector
{
public:
  KVCSD(const G4String& name);
  virtual ~KVCSD();

private:
  G4THitsCollection<KVCHit>* hitsCollection;

public:
  virtual G4bool ProcessHits(G4Step* aStep, G4TouchableHistory* ROhist);
  virtual void   Initialize(G4HCofThisEvent* HCTE);
  virtual void   EndOfEvent(G4HCofThisEvent* HCTE);
  virtual void   DrawAll();
  virtual void   PrintAll();
};

#endif
