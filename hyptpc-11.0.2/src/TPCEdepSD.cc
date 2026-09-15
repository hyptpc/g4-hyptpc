// -*- C++ -*-

#include "TPCEdepSD.hh"

#include <cmath>
#include <cstdio>

#include <G4DynamicParticle.hh>
#include <G4Step.hh>
#include <G4SystemOfUnits.hh>
#include <G4TouchableHistory.hh>
#include <G4Track.hh>
#include <G4VPhysicalVolume.hh>

#include "ConfMan.hh"
#include "TPCEdepHit.hh"
#include "TPCPadHelper.hh"

namespace
{
const auto& gConf = ConfMan::GetInstance();
}

//_____________________________________________________________________________
TPCEdepSD::TPCEdepSD(const G4String& name)
  : G4VSensitiveDetector(name),
    m_hits_collection()
{
  collectionName.insert("hit");
}

//_____________________________________________________________________________
TPCEdepSD::~TPCEdepSD(void)
{
}

//_____________________________________________________________________________
void
TPCEdepSD::Initialize(G4HCofThisEvent* HCTE)
{
  ntrk = 0;
  m_hits_collection = new G4THitsCollection<TPCEdepHit>(SensitiveDetectorName,
                                                        collectionName[0]);
  G4int hcid = GetCollectionID(0);
  HCTE->AddHitsCollection(hcid, m_hits_collection);
}

//_____________________________________________________________________________
G4bool
TPCEdepSD::ProcessHits(G4Step* aStep, G4TouchableHistory* /* ROhist */)
{
  const G4StepPoint* preStepPoint  = aStep->GetPreStepPoint();
  const G4StepPoint* postStepPoint = aStep->GetPostStepPoint();
  const G4Track*     aTrack        = aStep->GetTrack();
  G4int              copyNo        = preStepPoint->GetPhysicalVolume()->GetCopyNo();
  G4int              copyNo_post   = postStepPoint->GetPhysicalVolume()->GetCopyNo();

  if (aTrack->GetDefinition()->GetPDGCharge() == 0.)
    return false;

  G4TouchableHistory* theTouchable =
    (G4TouchableHistory*)(preStepPoint->GetTouchable());
  G4VPhysicalVolume* physVol = theTouchable->GetVolume();

  G4ThreeVector pos  = preStepPoint->GetPosition();
  G4double      hitx = pos.getX();
  G4double      hitz = pos.getZ();

  G4ThreeVector VertexPosition = aTrack->GetVertexPosition();
  G4ThreeVector VertexMomentum = aTrack->GetVertexMomentumDirection();
  G4double VertexEnergy = aTrack->GetVertexKineticEnergy(); // Ek = sqrt(p^2+m^2)-m

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

  G4int iLay_copyNo = copyNo - 2000;
  G4int iLay_pad    = 0;
  G4int iRow        = 0;
  G4int iPad        = TPCPadHelper::FindPad(hitz, hitx, iLay_pad, iRow);
  // valid layers 0..31; reject gap/OOB sentinels and invalid layer
  if (iPad < 0 || iLay_pad > 31)
    return false;
  G4int iLay = iLay_copyNo;

  G4double PadLen = TPCPadHelper::GetLength(iLay);
  G4ThreeVector PadPos(hitx, 0, hitz - TPCPadHelper::GetZTarget());
  G4ThreeVector MomT(mom.x(), 0, mom.z());
  G4double      alpha  = PadPos.theta() - MomT.theta();
  G4double      PathT  = PadLen * 1. / std::cos(alpha);
  G4double      Pitch  = mom.y() / MomT.mag();
  G4double      Path   = PathT * std::sqrt(1 + Pitch * Pitch);
  slength = Path;

  static const G4double conversion_factor =
    gConf.GetOrDefault<G4double>("TpcConversionFactor",
                                 TPCPadHelper::kDefaultConversionFactor);
  G4double cmTomm = 10.;
  G4double edep   = aStep->GetTotalEnergyDeposit() * conversion_factor * cmTomm;

#ifdef DEBUG
  // Compare geometry layer vs FindPadID-derived layer (not self-compare of iLay).
  const G4int    iLay_pad = TPCPadHelper::GetLayerID(iPad);
  const G4double radius   = std::hypot(hitx, hitz - TPCPadHelper::GetZTarget());

  G4ThreeVector  pos_post       = postStepPoint->GetPosition();
  const G4double hitx_post      = pos_post.getX();
  const G4double hitz_post      = pos_post.getZ();
  const G4int    iPad_post      = TPCPadHelper::FindPadID(hitz_post, hitx_post);
  const G4int    iLay_post      = TPCPadHelper::GetLayerID(iPad_post);
  const G4int    iLay_post_copy = copyNo_post - 2000;

  if (iLay_pad < 0 || iLay_pad > 31 || iLay_pad != iLay_copyNo) {
    G4cout << "pre strange layer"
           << " radius=" << radius
           << " iPad=" << iPad
           << " iLay_copyNo=" << iLay_copyNo
           << " iLay_pad=" << iLay_pad
           << " hitx=" << hitx
           << " hitz=" << hitz
           << G4endl;
  }
  // Post: only when still in TPCEdep copy range, or FindPadID returned a sentinel.
  const G4bool post_in_edep =
    (copyNo_post >= 2000 && copyNo_post < 2000 + 32);
  if (iPad_post < 0
      || (post_in_edep
          && (iLay_post < 0 || iLay_post > 31 || iLay_post != iLay_post_copy))) {
    G4cout << "post strange layer"
           << " radius=" << std::hypot(hitx_post, hitz_post - TPCPadHelper::GetZTarget())
           << " iPad=" << iPad_post
           << " iLay_copyNo=" << iLay_post_copy
           << " iLay_pad=" << iLay_post
           << " hitx=" << hitx_post
           << " hitz=" << hitz_post
           << G4endl;
  }
#endif

  G4String name = physVol->GetName();
  if (name == "TPC_PV") {
    name = "EdepPV-1";
  }
  sscanf(name, "EdepPV%d", &iLay_copyNo);

  TPCEdepHit* ahit = new TPCEdepHit(pos, mom, tof, tid, pid, iLay, iRow, beta,
                                    edep, parentID, tlength, mass, charge,
                                    VertexPosition, VertexMomentum, VertexEnergy,
                                    slength, parentID_pid);
  m_hits_collection->insert(ahit);
  return true;
}

//_____________________________________________________________________________
void
TPCEdepSD::EndOfEvent(G4HCofThisEvent* /* HCTE */)
{
}

//_____________________________________________________________________________
void
TPCEdepSD::DrawAll(void)
{
}

//_____________________________________________________________________________
void
TPCEdepSD::PrintAll(void)
{
  m_hits_collection->PrintAllHits();
}
