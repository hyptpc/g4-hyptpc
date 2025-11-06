// -*- C++ -*-

#ifndef SAC3_SD_HH
#define SAC3_SD_HH

#include <G4VSensitiveDetector.hh>

#include "SAC3Hit.hh"

class G4HCofThisEvent;
class G4Step;
class G4TouchableHistory;

//_____________________________________________________________________________
class SAC3SD : public G4VSensitiveDetector
{
public:
  static G4String ClassName();
  SAC3SD(const G4String& name);
  virtual ~SAC3SD();

private:
  G4THitsCollection<SAC3Hit>* m_hits_collection;
  G4double                   m_refractive_index;

public:
  G4double GetRefractiveIndex() const { return m_refractive_index; }
  void     SetRefractiveIndex(G4double index){ m_refractive_index = index; }

public:
  virtual G4bool ProcessHits(G4Step* aStep, G4TouchableHistory* ROhist);
  virtual void   Initialize(G4HCofThisEvent* HCTE);
  virtual void   EndOfEvent(G4HCofThisEvent* HCTE);
  virtual void   DrawAll();
  virtual void   PrintAll();
};

//_____________________________________________________________________________
inline G4String
SAC3SD::ClassName()
{
  static G4String s_name("SAC3SD");
  return s_name;
}

#endif
