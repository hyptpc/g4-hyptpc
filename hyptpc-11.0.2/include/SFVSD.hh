// -*- C++ -*-

#ifndef SFV_SD_HH
#define SFV_SD_HH

#include <G4VSensitiveDetector.hh>

#include "SFVHit.hh"

class G4HCofThisEvent;
class G4Step;
class G4TouchableHistory;

//_____________________________________________________________________________
class SFVSD : public G4VSensitiveDetector
{
public:
  static G4String ClassName();
  SFVSD(const G4String& name);
  virtual ~SFVSD();

private:
  G4THitsCollection<SFVHit>* m_hits_collection;

public:
  virtual G4bool ProcessHits(G4Step* aStep, G4TouchableHistory* ROhist);
  virtual void   Initialize(G4HCofThisEvent* HCTE);
  virtual void   EndOfEvent(G4HCofThisEvent* HCTE);
  virtual void   DrawAll();
  virtual void   PrintAll();
};

//_____________________________________________________________________________
inline G4String
SFVSD::ClassName()
{
  static G4String s_name("SFVSD");
  return s_name;
}

#endif