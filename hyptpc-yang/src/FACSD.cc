// ====================================================================
//   FACSD.cc
//
// ====================================================================
#include "G4VPhysicalVolume.hh"
#include "G4Step.hh"
#include "G4Track.hh"
#include "G4VTouchable.hh"
#include "G4TouchableHistory.hh"
#include "FACSD.hh"
#include "FACHit.hh"

////////////////////////////////////////////////
  FACSD::FACSD( G4String name)
: G4VSensitiveDetector(name)
  ////////////////////////////////////////////////
{
  collectionName.insert( G4String( "FacCollection" ));
}

/////////////////////////////////
FACSD::~FACSD()
  /////////////////////////////////
{
}


////////////////////////////////////////////////
void FACSD::Initialize(G4HCofThisEvent* HCTE)
  ////////////////////////////////////////////////
{
  static int HCID = -1;
  // create hit collection(s)
  FacCollection = new FACHitsCollection( SensitiveDetectorName,
      collectionName[0]);

  // push H.C. to "Hit Collection of This Event"
  if( HCID<0 )
    {
      HCID = GetCollectionID(0);
    }
  HCTE-> AddHitsCollection(HCID, FacCollection);
}


///////////////////////////////////////////////////////////
G4bool FACSD::ProcessHits(G4Step* aStep,
    G4TouchableHistory* ROhist)
  ///////////////////////////////////////////////////////////
{

  // get step information from "PreStepPoint"
  const G4StepPoint* preStepPoint= aStep-> GetPreStepPoint();
  G4String particleName;

  //if(aStep-> GetTrack()-> GetDefinition()-> GetPDGCharge() == 0.)
  //  return false;

  particleName = aStep-> GetTrack()-> GetDefinition()-> GetParticleName();
  /*
  // e+/e- rejection
  if( particleName == "e-")
  return false;
  if( particleName == "e+")
  return false;
  */
  const G4Track* aTrack = aStep->GetTrack();
  G4String particleType;
  particleType = aTrack->GetDefinition()->GetParticleType();

  //if(particleType == "lepton")
  //  return false;
  //
  //  if( (particleName != "kaon+")){}
  //    return false;
  //  {}

  //  if( (particleName != "pi-") && (particleName != "pi+")){}
  //    return false;
  //  {}
  //  if((particleName != "pi+")&& (particleName != "pi-")
  //     && (particleName != "proton")){}
  //    return false;
  //  {}
  if(preStepPoint-> GetStepStatus() != fGeomBoundary) return false;
  G4TouchableHistory* theTouchable
    = (G4TouchableHistory*)(aStep->GetPreStepPoint()->GetTouchable());
  G4VPhysicalVolume* physVol = theTouchable->GetVolume();


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
  //  G4double slength = aStep->GetTrack()-> GetStepLength();
  //  G4cout<<mass<<":"<<tlength<<":"<<slength<<G4endl;
  G4int parentID =  aStep-> GetTrack()-> GetParentID();
  //    G4cout <<"test: "<<parentID << G4endl;
  // Get Pad number
  G4int iDet;
  G4String name = physVol->GetName();
  //  G4cout << name << G4endl;
  iDet = preStepPoint -> GetPhysicalVolume()->GetCopyNo();
  //  sscanf(name,"ACPV%d",&iDet);

  //G4cout<<"[FACSD]"<< name <<" trid: "<< tid <<" collection_size: " <<G4endl;

  // create a new hit and push them to "Hit Coleltion"
  FACHit* ahit= new FACHit(pos, mom, tof, tid, pid, iDet,mass,qq,parentID, VertexPosition, VertexMomentum, VertexEnergy,tlength);
  FacCollection-> insert(ahit);
  return true;
}

////////////////////////////////////////////////////
void FACSD::EndOfEvent(G4HCofThisEvent* HCTE)
  ////////////////////////////////////////////////////
{
}

////////////////////////////
void FACSD::DrawAll()
  ////////////////////////////
{
}

/////////////////////////////
void FACSD::PrintAll()
  /////////////////////////////
{
  FacCollection-> PrintAllHits();
}
