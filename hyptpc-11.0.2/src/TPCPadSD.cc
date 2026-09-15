// -*- C++ -*-

#include "TPCPadSD.hh"

#include <cmath>
#include <cstdio>

#include <G4DynamicParticle.hh>
#include <G4Step.hh>
#include <G4SystemOfUnits.hh>
#include <G4TouchableHistory.hh>
#include <G4Track.hh>
#include <G4VPhysicalVolume.hh>
#include <Randomize.hh>

#include "TPCdEdx.hh"
#include "TPCPadHit.hh"
#include "TPCPadHelper.hh"

//_____________________________________________________________________________
TPCPadSD::TPCPadSD(const G4String& name)
  : G4VSensitiveDetector(name),
    m_hits_collection()
{
  collectionName.insert("hit");
}

//_____________________________________________________________________________
TPCPadSD::~TPCPadSD(void)
{
}

//_____________________________________________________________________________
void
TPCPadSD::Initialize(G4HCofThisEvent* HCTE)
{
  ntrk = 0;
  m_hits_collection = new G4THitsCollection<TPCPadHit>(SensitiveDetectorName,
                                                       collectionName[0]);
  G4int hcid = GetCollectionID(0);
  HCTE->AddHitsCollection(hcid, m_hits_collection);
}

//_____________________________________________________________________________
G4bool
TPCPadSD::ProcessHits(G4Step* aStep, G4TouchableHistory* /* ROhist */)
{
  const G4StepPoint* preStepPoint = aStep->GetPreStepPoint();
  const G4Track* aTrack = aStep->GetTrack();
  G4int copyNo = preStepPoint->GetPhysicalVolume()->GetCopyNo();

  if (preStepPoint->GetStepStatus() != fGeomBoundary)
    return false;

  if (aTrack->GetDefinition()->GetPDGCharge() == 0.)
    return false;

  // Optional particle filters (uncomment refs + conditions together):
  // const G4String& particleName = aTrack->GetDefinition()->GetParticleName();
  // const G4String& particleType = aTrack->GetDefinition()->GetParticleType();

  G4TouchableHistory* theTouchable =
    (G4TouchableHistory*)(preStepPoint->GetTouchable());
  G4VPhysicalVolume* physVol = theTouchable->GetVolume();

  G4ThreeVector pos  = preStepPoint->GetPosition();
  G4double      hitx = pos.getX();
  G4double      hitz = pos.getZ();

  G4ThreeVector VertexPosition = aTrack->GetVertexPosition();
  G4ThreeVector VertexMomentum = aTrack->GetVertexMomentumDirection();
  G4double      VertexEnergy   = aTrack->GetVertexKineticEnergy(); // Ek = sqrt(p^2+m^2)-m

  G4ThreeVector mom          = preStepPoint->GetMomentum();
  G4double      tof          = preStepPoint->GetGlobalTime();
  G4double      beta         = preStepPoint->GetBeta();
  G4int         tid          = aTrack->GetTrackID();
  G4int         pid          = aTrack->GetDefinition()->GetPDGEncoding();
  G4double      mass         = aTrack->GetDynamicParticle()->GetMass();
  G4int         charge       = aTrack->GetDefinition()->GetPDGCharge();
  G4double      tlength      = aTrack->GetTrackLength();
  G4double      slength      = aTrack->GetStepLength();
  G4int         parentID     = aTrack->GetParentID();
  G4int         parentID_pid =
    aTrack->GetDynamicParticle()->GetDefinition()->GetPDGEncoding();

  G4int iLay_copyNo = copyNo;
  G4int iLay        = 0;
  G4int iRow        = 0;
  G4int iPad        = TPCPadHelper::FindPad(hitz, hitx, iLay, iRow);
  // valid layers 0..31; reject gap/OOB sentinels and invalid layer
  if (iPad < 0 || iLay > 31)
    return false;

  G4double PadLen = TPCPadHelper::GetLength(iLay);
  G4ThreeVector PadPos(hitx, 0, hitz - TPCPadHelper::GetZTarget());
  G4ThreeVector MomT(mom.x(), 0, mom.z());
  G4double alpha = PadPos.theta() - MomT.theta();
  G4double PathT = PadLen * 1. / std::cos(alpha);
  G4double Pitch = mom.y() / MomT.mag();
  G4double Path  = PathT * std::sqrt(1 + Pitch * Pitch);
  slength        = Path;

  G4double edepMean = TPCdEdx::dEdx(mass, beta) * Path;
  G4double edepSig  = TPCdEdx::dEdxSigma(mass, mom.mag()) * Path;
  if (edepSig / edepMean < 0.01 || edepSig / edepMean > 0.5)
    edepSig = 0.2 * edepMean;
  G4double edep = G4RandGauss::shoot(edepMean, edepSig);
  if (edep < 0.1 * edepMean)
    edep = 0.1 * edepMean;

#ifdef DEBUG
  G4double      radius  = std::hypot(hitx, hitz - TPCPadHelper::GetZTarget());
  G4ThreeVector Point   = TPCPadHelper::GetPosition(iPad);
  G4int         iPad_re = TPCPadHelper::FindPadID(Point.z(), Point.x());
  G4cout << "hitx = " << hitx
         << ", pointx " << Point.x()
         << ", hitz = " << hitz
         << ", pointz " << Point.z()
         << ", radius = " << radius
         << ", iPad = " << iPad
         << ", iPad_re = " << iPad_re
         << ", iLay_copyNo = " << iLay_copyNo
         << ", iLay = " << iLay << G4endl;
  G4cout << "dx = " << hitx - Point.x()
         << ", dz = " << hitz - Point.z() << G4endl;
#endif

  G4String name = physVol->GetName();
  if (name == "TPC_PV") {
    name = "PadPV-1";
  }
  sscanf(name, "PadPV%d", &iLay_copyNo);

  TPCPadHit* ahit = new TPCPadHit(pos, mom, tof, tid, pid, iLay, iRow, beta,
                                  edep, parentID, tlength, mass, charge,
                                  VertexPosition, VertexMomentum, VertexEnergy,
                                  slength, parentID_pid);
  m_hits_collection->insert(ahit);
  return true;
}

//_____________________________________________________________________________
void
TPCPadSD::EndOfEvent(G4HCofThisEvent* /* HCTE */)
{
}

//_____________________________________________________________________________
void
TPCPadSD::DrawAll(void)
{
}

//_____________________________________________________________________________
void
TPCPadSD::PrintAll(void)
{
  m_hits_collection->PrintAllHits();
}
