// -*- C++ -*-

#ifndef CVC_SD_HH
#define CVC_SD_HH

#include <G4VSensitiveDetector.hh>

#include "CVCHit.hh"

class G4HCofThisEvent;
class G4Step;
class G4TouchableHistory;

//_____________________________________________________________________________
class CVCSD : public G4VSensitiveDetector
{
public:
  static G4String ClassName();
  CVCSD(const G4String& name);
  virtual ~CVCSD();

private:
  G4THitsCollection<CVCHit>* m_hits_collection;

public:
  virtual G4bool ProcessHits(G4Step* aStep, G4TouchableHistory* ROhist);
  virtual void   Initialize(G4HCofThisEvent* HCTE);
  virtual void   EndOfEvent(G4HCofThisEvent* HCTE);
  virtual void   DrawAll();
  virtual void   PrintAll();
};

//_____________________________________________________________________________
inline G4String
CVCSD::ClassName()
{
  static G4String s_name("CVCSD");
  return s_name;
}

#endif