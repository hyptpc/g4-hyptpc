#include "AeroHit.hh"
//#include "BACDetectorConstruction_167_2.hh"

#include "G4VVisManager.hh"
#include "G4Circle.hh"
#include "G4VisAttributes.hh"
#include "G4AttDef.hh"
#include "G4AttValue.hh"
#include "G4UIcommand.hh"
#include "G4UnitsTable.hh"
#include "G4SystemOfUnits.hh"
#include "G4ios.hh"
#include "G4AttDefStore.hh"

#include <iomanip>

G4Allocator<AeroHit> AeroHitAllocator;

AeroHit::AeroHit()
  :xyz(0.,0.,0.),tof(0.)
{}

AeroHit::AeroHit(G4ThreeVector& axyz, G4double t)
  :xyz(axyz),tof(t)
{}

AeroHit::~AeroHit()
{}

AeroHit::AeroHit(G4ThreeVector &axyz, G4double t, G4int pid, G4double angle)
  :xyz(axyz), tof(t),particleID(pid),angleget(angle)
{}
