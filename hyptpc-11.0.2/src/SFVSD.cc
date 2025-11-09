// -*- C++ -*-

#include "SFVSD.hh"

#include <G4Step.hh>
#include <G4TouchableHistory.hh>
#include <G4Track.hh>
#include <G4VPhysicalVolume.hh>
#include <G4VTouchable.hh>

#include "SFVHit.hh"
#include "FuncName.hh"

//_____________________________________________________________________________
SFVSD::SFVSD(const G4String& name)
    : G4VSensitiveDetector(name),
      m_hits_collection()
{
    collectionName.insert("hit");
}

//_____________________________________________________________________________
SFVSD::~SFVSD()
{
}

//_____________________________________________________________________________
void SFVSD::Initialize(G4HCofThisEvent* HCTE)
{
    m_hits_collection = new G4THitsCollection<SFVHit>(SensitiveDetectorName,
                                                      collectionName[0]);
    HCTE->AddHitsCollection(GetCollectionID(0), m_hits_collection);
}

//_____________________________________________________________________________
G4bool
SFVSD::ProcessHits(G4Step* aStep, G4TouchableHistory* /* ROhist */)
{
    const auto preStepPoint = aStep->GetPreStepPoint();

    if (preStepPoint->GetStepStatus() != fGeomBoundary)
        return false;

    auto aHit = new SFVHit(SensitiveDetectorName, aStep);
    // aHit->Print();

    m_hits_collection->insert(aHit);

    return true;
}

//_____________________________________________________________________________
void SFVSD::EndOfEvent(G4HCofThisEvent* /* HCTE */)
{
}

//_____________________________________________________________________________
void SFVSD::DrawAll()
{
}

//_____________________________________________________________________________
void SFVSD::PrintAll()
{
    m_hits_collection->PrintAllHits();
}