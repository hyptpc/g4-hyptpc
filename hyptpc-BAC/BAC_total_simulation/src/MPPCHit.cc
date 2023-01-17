#include "MPPCHit.hh"
//#include "BACDetectorConstruction.hh"

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

G4Allocator<MPPCHit> MPPCHitAllocator;

MPPCHit::MPPCHit()
  :xyz(0.,0.,0.),tof(0.)
{}

MPPCHit::MPPCHit(G4ThreeVector& axyz, G4double t)
  :xyz(axyz),tof(t)
{}

MPPCHit::~MPPCHit()
{}

MPPCHit::MPPCHit(G4ThreeVector &axyz, G4ThreeVector &wxyz, G4double t, G4int pid, G4double Wavelength, G4int cn)
  :xyz(axyz), worldxyz(wxyz), tof(t),particleID(pid), wavelengthMP(Wavelength), copynum(cn)
{}

