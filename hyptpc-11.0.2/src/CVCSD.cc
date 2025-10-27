// -*- C++ -*-

#include "CVCSD.hh"

#include <G4Step.hh>
#include <G4TouchableHistory.hh>
#include <G4Track.hh>
#include <G4VPhysicalVolume.hh>
#include <G4VTouchable.hh>

#include "CVCHit.hh"
#include "FuncName.hh"

//_____________________________________________________________________________
CVCSD::CVCSD(const G4String& name)
    : G4VSensitiveDetector(name),
      m_hits_collection()
{
    collectionName.insert("hit");
}

//_____________________________________________________________________________
CVCSD::~CVCSD()
{
}

//_____________________________________________________________________________
void CVCSD::Initialize(G4HCofThisEvent* HCTE)
{
    m_hits_collection = new G4THitsCollection<CVCHit>(SensitiveDetectorName,
                                                      collectionName[0]);
    HCTE->AddHitsCollection(GetCollectionID(0), m_hits_collection);
}

//_____________________________________________________________________________
G4bool
CVCSD::ProcessHits(G4Step* aStep, G4TouchableHistory* /* ROhist */)
{
    const auto preStepPoint = aStep->GetPreStepPoint();

    if (preStepPoint->GetStepStatus() != fGeomBoundary)
        return false;

    auto aHit = new CVCHit(SensitiveDetectorName, aStep);
    // aHit->Print();

    m_hits_collection->insert(aHit);

    return true;
}

//_____________________________________________________________________________
void CVCSD::EndOfEvent(G4HCofThisEvent* /* HCTE */)
{
}

//_____________________________________________________________________________
void CVCSD::DrawAll()
{
}

//_____________________________________________________________________________
void CVCSD::PrintAll()
{
    m_hits_collection->PrintAllHits();
}