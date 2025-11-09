// -*- C++ -*-

#include "SFVHit.hh"

#include <CLHEP/Units/PhysicalConstants.h>
#include <G4Step.hh>

G4Allocator<SFVHit> SFVHitAllocator;

//_____________________________________________________________________________
SFVHit::SFVHit(const G4String& name, G4Step* step)
  : G4VHit(),
    VHitInfo(name, step)
{
}

//_____________________________________________________________________________
SFVHit::~SFVHit()
{
}

//_____________________________________________________________________________
void
SFVHit::Draw()
{
}

//_____________________________________________________________________________
void
SFVHit::Print()
{
  VHitInfo::Print();
}