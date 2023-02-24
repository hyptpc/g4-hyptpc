// -*- C++ -*-

#ifndef TPC_PAD_SD_HH
#define TPC_PAD_SD_HH

#include <G4VSensitiveDetector.hh>

#include "TPCHit.hh"

class G4HCofThisEvent;
class G4Step;
class G4TouchableHistory;

//_____________________________________________________________________________
class TPCSD : public G4VSensitiveDetector
{
public:
  static G4String ClassName();
  TPCSD(const G4String& name);
  virtual ~TPCSD();

private:
  G4THitsCollection<TPCHit>* hitsCollection;
  G4double select_plane;
  G4int num_plane;
  G4int num_deadarea;
  G4int select_dead;

  G4int m_gem_discharge;
  G4int m_gem_fix_dead;
  G4int m_gem_dead_plane;
  G4int m_gem_dead_plane_division;
  G4double m_dead_area;

public:
  G4int ntrk;
  virtual G4bool ProcessHits(G4Step* aStep, G4TouchableHistory* ROhist);
  virtual void   Initialize(G4HCofThisEvent* HCTE);
  virtual void   EndOfEvent(G4HCofThisEvent* HCTE);
  virtual void   DrawAll();
  virtual void   PrintAll();
};

//_____________________________________________________________________________
inline G4String
TPCSD::ClassName()
{
  static G4String s_name("TPCSD");
  return s_name;
}
#endif
