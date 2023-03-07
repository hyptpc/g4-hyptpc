// -*- C++ -*-

#include "KVCSD.hh"

#include <G4Step.hh>
#include <G4TouchableHistory.hh>
#include <G4Track.hh>
#include <G4VPhysicalVolume.hh>
#include <G4VTouchable.hh>

#include "KVCHit.hh"

//_____________________________________________________________________________
KVCSD::KVCSD(const G4String& name)
  : G4VSensitiveDetector(name)
{
  collectionName.insert("hit");
}

//_____________________________________________________________________________
KVCSD::~KVCSD()
{
}

//_____________________________________________________________________________
void
KVCSD::Initialize(G4HCofThisEvent* HCTE)
{
  hitsCollection = new G4THitsCollection<KVCHit>(SensitiveDetectorName,
						    collectionName[0]);
  G4int hcid = GetCollectionID(0);
  HCTE->AddHitsCollection(hcid, hitsCollection);
}

//_____________________________________________________________________________
G4bool
KVCSD::ProcessHits(G4Step* aStep, G4TouchableHistory* /* ROhist */)
{
  const G4StepPoint* preStepPoint= aStep-> GetPreStepPoint();
  G4String particleName;

  if(aStep-> GetTrack()-> GetDefinition()-> GetPDGCharge() == 0.)
    return false;

  particleName = aStep-> GetTrack()-> GetDefinition()-> GetParticleName();
  /*
  // e+/e- rejection
  if(particleName == "e-")
    return false;
  if(particleName == "e+")
    return false;
  */
  const G4Track* aTrack = aStep->GetTrack();
  G4String particleType;
  particleType = aTrack->GetDefinition()->GetParticleType();

  //    if(particleType == "lepton")
  //      return false;
  //
  //  if((particleName != "kaon+")){
  //    return false;
  //  }

  //  if((particleName != "pi-") && (particleName != "pi+")){
  //    return false;
  //  }
  //  if((particleName != "pi+")&& (particleName != "pi-")
  //     && (particleName != "proton")){
  //    return false;
  //  }
  if(preStepPoint-> GetStepStatus() != fGeomBoundary) return false;

  G4TouchableHistory* theTouchable
    = (G4TouchableHistory*)(aStep->GetPreStepPoint()->GetTouchable());

  G4ThreeVector VertexPosition = aTrack->GetVertexPosition();
  G4ThreeVector VertexMomentum = aTrack->GetVertexMomentumDirection();
  G4double VertexEnergy = aTrack -> GetVertexKineticEnergy(); // Ek = sqrt(p^2+m^2)-m

  G4ThreeVector pos= preStepPoint-> GetPosition();
  G4ThreeVector mom= preStepPoint-> GetMomentum();
  G4double tof= preStepPoint-> GetGlobalTime();
  G4int tid =  aStep-> GetTrack()-> GetTrackID();
  G4int pid =  aStep-> GetTrack()-> GetDefinition() -> GetPDGEncoding();
  G4double mass =  aStep-> GetTrack()-> GetDynamicParticle() -> GetMass();
  G4int qq =  aStep-> GetTrack()-> GetDynamicParticle() -> GetCharge();
  G4double tlength = aStep->GetTrack()-> GetTrackLength();
  G4int parentID =  aStep-> GetTrack()-> GetParentID();
  G4int copy_number = theTouchable->GetCopyNumber();

  auto ahit = new KVCHit(pos, mom, tof, tid, pid, copy_number, mass, qq,
                         parentID,
                         VertexPosition, VertexMomentum, VertexEnergy,tlength);
  hitsCollection-> insert(ahit);

  return true;
}

//_____________________________________________________________________________
void
KVCSD::EndOfEvent(G4HCofThisEvent* /* HCTE */)
{
}

//_____________________________________________________________________________
void
KVCSD::DrawAll()
{
}

//_____________________________________________________________________________
void
KVCSD::PrintAll()
{
  hitsCollection-> PrintAllHits();
}
