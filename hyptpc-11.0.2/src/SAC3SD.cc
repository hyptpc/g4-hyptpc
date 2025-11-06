// -*- C++ -*-

#include "SAC3SD.hh"

#include <G4Step.hh>
#include <G4TouchableHistory.hh>
#include <G4Track.hh>
#include <G4VPhysicalVolume.hh>
#include <G4VTouchable.hh>

#include "FuncName.hh"
#include "SAC3Hit.hh"

//_____________________________________________________________________________
SAC3SD::SAC3SD(const G4String& name)
  : G4VSensitiveDetector(name),
    m_hits_collection(),
    m_refractive_index()
{
  collectionName.insert("hit");
}

//_____________________________________________________________________________
SAC3SD::~SAC3SD()
{
}

//_____________________________________________________________________________
void
SAC3SD::Initialize(G4HCofThisEvent* HCTE)
{
  m_hits_collection = new G4THitsCollection<SAC3Hit>(SensitiveDetectorName,
                                                    collectionName[0]);
  HCTE->AddHitsCollection(GetCollectionID(0), m_hits_collection);
}

//_____________________________________________________________________________
G4bool
SAC3SD::ProcessHits(G4Step* aStep, G4TouchableHistory* /* ROhist */)
{
  const auto preStepPoint = aStep->GetPreStepPoint();
  const auto aTrack = aStep->GetTrack();
  const auto Definition = aTrack->GetDefinition();
  const G4String particleName = Definition->GetParticleName();
  const G4String particleType = Definition->GetParticleType();

  if(preStepPoint->GetStepStatus() != fGeomBoundary)
    return false;
  if(Definition->GetPDGCharge() == 0.)
    return false;

  // if(particleName == "e-")
  //   return false;
  // if(particleName == "e+")
  //   return false;
  // if(particleName != "kaon+")
  //   return false;
  // if(particleName != "pi-" && particleName != "pi+")
  //   return false;
  // if(particleName != "pi+" && particleName != "pi-" &&
  //     particleName != "proton")
  //   return false;
  // if(particleType == "lepton")
  //   return false;

  // const G4double momentum_threshold =
  //   hit->GetMass() / std::sqrt(m_refractive_index*m_refractive_index - 1.);
  // if(hit->GetMomentum().mag() < momentum_threshold)
  //   return false;

  // if(hit->GetParticleName().contains("proton")){
  //   hit->Print();
  //   G4cout << momentum_threshold << G4endl;
  // }

  m_hits_collection->insert(new SAC3Hit(SensitiveDetectorName, aStep));

  return true;
}

//_____________________________________________________________________________
void
SAC3SD::EndOfEvent(G4HCofThisEvent* /* HCTE */)
{
}

//_____________________________________________________________________________
void
SAC3SD::DrawAll()
{
}

//_____________________________________________________________________________
void
SAC3SD::PrintAll()
{
  m_hits_collection->PrintAllHits();
}
